from pathlib import Path
import gzip
import re
import pandas as pd

EXPR = Path("results/GSE67462/03_expression_for_EDA.csv")
META = Path("results/GSE67462/02_sample_metadata.csv")
SOFT = Path("Data/GPL19972_family.soft.gz")
GTF = Path("Data/GSE67520/mm9.refGene.gtf.gz")
OUT = Path("results/Dynamics/z6_gse67462_mechanistic_interpretation")
OUT.mkdir(parents=True, exist_ok=True)

def clean_id(x):
    x = str(x).strip().strip('"').strip("'").strip()
    return x

def accession_key(x):
    x = clean_id(x)
    x = x.split("///")[0].split("//")[0].strip()
    return re.sub(r"\.\d+$", "", x)

def attrs_dict(s):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', s))

# GTF: transkrypt -> symbol genu
tx_to_gene = {}
gene_symbols = set()

with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9 or parts[2] != "transcript":
            continue
        a = attrs_dict(parts[8])
        gene = a.get("gene_name") or a.get("gene_id")
        tx = a.get("transcript_id")
        if gene:
            gene_symbols.add(gene)
            if tx:
                tx_to_gene[accession_key(tx)] = gene

symbol_lookup = {g.casefold(): g for g in gene_symbols}

# Adnotacja GPL
platform = {}
inside = False
header = None

with gzip.open(SOFT, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.rstrip("\r\n")
        if line == "!platform_table_begin":
            inside, header = True, None
            continue
        if line == "!platform_table_end":
            break
        if not inside or not line or line.startswith("!"):
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        fields += [""] * max(0, len(header) - len(fields))
        row = dict(zip(header, fields))
        platform[clean_id(row.get("ID", ""))] = row

probe_to_gene = {}
source = {}

for probe, row in platform.items():
    acc = clean_id(row.get("GB_ACC", ""))
    gene = tx_to_gene.get(accession_key(acc))
    if gene:
        probe_to_gene[probe] = gene
        source[probe] = "transcript_accession"
        continue

    desc = row.get("Description", "")
    candidates = re.findall(r"\(([^()]*)\)", desc)
    matches = {
        symbol_lookup[c.strip().casefold()]
        for c in candidates
        if c.strip().casefold() in symbol_lookup
    }
    if len(matches) == 1:
        probe_to_gene[probe] = next(iter(matches))
        source[probe] = "description_symbol"
    elif len(matches) > 1:
        source[probe] = "ambiguous_description"
    else:
        source[probe] = "unmapped"

# Ekspresja i normalizacja ID_REF
expr = pd.read_csv(EXPR, dtype=str)
if "ID_REF" not in expr.columns:
    raise SystemExit("STOP: brak kolumny ID_REF")

expr["probe_clean"] = expr["ID_REF"].map(clean_id)
meta = pd.read_csv(META, dtype=str)
samples = meta["sample"].tolist()

missing_samples = [s for s in samples if s not in expr.columns]
if missing_samples:
    raise SystemExit(f"STOP: brak próbek w ekspresji: {missing_samples}")

expr["gene"] = expr["probe_clean"].map(probe_to_gene)
mapped = expr[expr["gene"].notna()].copy()

for s in samples:
    mapped[s] = pd.to_numeric(mapped[s], errors="coerce")

gene_expr = mapped.groupby("gene", sort=False)[samples].median()
gene_expr.index.name = "gene"
gene_expr.to_csv(OUT / "01_gene_expression_by_sample.csv")

# Raport dla wszystkich sond występujących w macierzy ekspresji
rows = []
for probe in expr["probe_clean"].dropna().unique():
    row = platform.get(probe, {})
    rows.append({
        "probe_id": probe,
        "GB_ACC": row.get("GB_ACC", ""),
        "Description": row.get("Description", ""),
        "gene_symbol": probe_to_gene.get(probe, ""),
        "mapping_source": source.get(probe, "probe_not_in_platform_table"),
    })
pd.DataFrame(rows).to_csv(OUT / "02_probe_mapping_report.csv", index=False)

# Średnia ekspresja w grupie etapowej — dwa powtórzenia na etap
long = gene_expr.T.copy()
long["sample"] = long.index
long = long.merge(meta, on="sample", how="left")
stage_expr = long.groupby("group")[gene_expr.index.tolist()].mean().T
stage_expr.index.name = "gene"
stage_expr.to_csv(OUT / "03_gene_expression_by_stage.csv")

report = pd.DataFrame(rows)
print("\n=== WYNIK MAPOWANIA ===")
print("Sond w macierzy:", expr["probe_clean"].nunique())
print("Sond znalezionych w GPL:",
      expr["probe_clean"].isin(platform).sum())
print("Sond przypisanych do genu:",
      expr["probe_clean"].isin(probe_to_gene).sum())
print("Unikalne geny ekspresji:", len(gene_expr))
print("\nŹródła mapowania:")
print(report["mapping_source"].value_counts().to_string())
print("\nWymiary macierzy ekspresji:", gene_expr.shape)
print("Wymiary ekspresji etapowej:", stage_expr.shape)
print("Brakujące wartości ekspresji:", int(gene_expr.isna().sum().sum()))
print("\nZapisano w:", OUT)
