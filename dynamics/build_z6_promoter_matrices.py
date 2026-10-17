from pathlib import Path
import csv
import gzip
import re
import tarfile
import heapq
from collections import defaultdict

import pandas as pd

ROOT = Path(".")
EXPR = ROOT / "results/GSE67462/03_expression_for_EDA.csv"
META = ROOT / "results/GSE67462/02_sample_metadata.csv"
SOFT = ROOT / "Data/GPL19972_family.soft.gz"
GTF = ROOT / "Data/GSE67520/mm9.refGene.gtf.gz"
TAR = ROOT / "Data/GSE67520_RAW.tar"

OUT = ROOT / "results/Dynamics/z6_gse67462_mechanistic_interpretation"
OUT.mkdir(parents=True, exist_ok=True)

def stop(message):
    raise SystemExit("STOP: " + message)

for p in [EXPR, META, SOFT, GTF, TAR]:
    if not p.exists():
        stop(f"Brak pliku: {p}")

def accession_key(value):
    value = str(value).strip()
    value = value.split("///")[0].split("//")[0].strip()
    return re.sub(r"\.\d+$", "", value)

def parse_attrs(s):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', s))

# 1. GTF: symbole genów, transkrypty i promotory strand-aware
print("1/5 Czytam GTF...")
tx_to_gene = {}
gene_symbols = set()
promoters = defaultdict(list)

with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "transcript":
            continue

        chrom, _, _, start, end, _, strand, _, attr = fields
        a = parse_attrs(attr)
        gene = a.get("gene_name") or a.get("gene_id")
        tx = a.get("transcript_id")
        if not gene:
            continue

        gene_symbols.add(gene)
        if tx:
            tx_to_gene[accession_key(tx)] = gene

        start, end = int(start), int(end)
        # GTF: 1-based inclusive; przedziały wynikowe: BED 0-based, half-open.
        if strand == "+":
            left, right = max(0, start - 1 - 2000), start - 1 + 500
        elif strand == "-":
            left, right = max(0, end - 500), end + 2000
        else:
            continue

        if right > left:
            promoters[gene].append((chrom, left, right))

# Scal nakładające się promotory tego samego genu, aby nie liczyć
# tego samego piku wielokrotnie przez alternatywne transkrypty.
merged_promoters = []
for gene, intervals in promoters.items():
    by_chrom = defaultdict(list)
    for chrom, start, end in intervals:
        by_chrom[chrom].append((start, end))

    for chrom, intervals2 in by_chrom.items():
        intervals2.sort()
        left, right = intervals2[0]
        for start, end in intervals2[1:]:
            if start <= right:
                right = max(right, end)
            else:
                merged_promoters.append((chrom, left, right, gene))
                left, right = start, end
        merged_promoters.append((chrom, left, right, gene))

print("  Geny z promotorami:", len(promoters))
print("  Zmergowane przedziały:", len(merged_promoters))

# 2. Adnotacja platformy: accession -> gene; opis jest fallbackiem.
print("2/5 Czytam adnotację platformy...")
platform = {}
in_table = False
header = None

with gzip.open(SOFT, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.rstrip("\r\n")
        if line == "!platform_table_begin":
            in_table, header = True, None
            continue
        if line == "!platform_table_end":
            break
        if not in_table or not line or line.startswith("!"):
            continue

        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        if len(fields) < len(header):
            fields += [""] * (len(header) - len(fields))
        row = dict(zip(header, fields))
        platform[row.get("ID", "")] = row

# Normalizujemy symbole z GTF tylko na potrzeby dopasowania bez względu
# na wielkość liter; końcowy symbol zachowujemy dokładnie jak w GTF.
symbol_lookup = {g.casefold(): g for g in gene_symbols}
probe_to_gene = {}
mapping_source = {}

for probe, row in platform.items():
    accession = row.get("GB_ACC", "").strip()
    description = row.get("Description", "")
    gene = tx_to_gene.get(accession_key(accession))

    if gene:
        probe_to_gene[probe] = gene
        mapping_source[probe] = "transcript_accession"
        continue

    candidates = re.findall(r"\(([^()]*)\)", description)
    matches = {symbol_lookup[c.strip().casefold()]
               for c in candidates if c.strip().casefold() in symbol_lookup}
    if len(matches) == 1:
        gene = next(iter(matches))
        probe_to_gene[probe] = gene
        mapping_source[probe] = "description_symbol"
    elif len(matches) > 1:
        mapping_source[probe] = "ambiguous_description"
    else:
        mapping_source[probe] = "unmapped"

# 3. Ekspresja: sondy -> geny; wiele sond na gen agregujemy medianą.
print("3/5 Czytam ekspresję i metadane...")
expr = pd.read_csv(EXPR, dtype={"ID_REF": str})
if "ID_REF" not in expr.columns:
    stop("Brak kolumny ID_REF w macierzy ekspresji.")

meta = pd.read_csv(META, dtype=str)
if not {"sample", "group"}.issubset(meta.columns):
    stop("Metadane muszą zawierać kolumny sample i group.")

sample_cols = [c for c in expr.columns if c in set(meta["sample"])]
if len(sample_cols) != len(meta):
    missing = sorted(set(meta["sample"]) - set(expr.columns))
    stop(f"Próbki z metadanych nie występują w ekspresji: {missing}")

expr["gene"] = expr["ID_REF"].map(probe_to_gene)
mapped = expr[expr["gene"].notna()].copy()
for c in sample_cols:
    mapped[c] = pd.to_numeric(mapped[c], errors="coerce")

gene_expr = mapped.groupby("gene", sort=False)[sample_cols].median()
gene_expr.index.name = "gene"
gene_expr.to_csv(OUT / "01_gene_expression_by_sample.csv")

mapping_rows = []
for probe in expr["ID_REF"].dropna().unique():
    row = platform.get(probe, {})
    mapping_rows.append({
        "probe_id": probe,
        "GB_ACC": row.get("GB_ACC", ""),
        "Description": row.get("Description", ""),
        "gene_symbol": probe_to_gene.get(probe, ""),
        "mapping_source": mapping_source.get(probe, "probe_not_in_platform_table"),
    })
pd.DataFrame(mapping_rows).to_csv(OUT / "02_probe_mapping_report.csv", index=False)

expr_by_stage = gene_expr.T.copy()
expr_by_stage["sample"] = expr_by_stage.index
expr_by_stage = expr_by_stage.merge(meta, on="sample", how="left")
stage_expr = expr_by_stage.groupby("group")[gene_expr.index.tolist()].mean().T
stage_expr.index.name = "gene"
stage_expr.to_csv(OUT / "03_gene_expression_by_stage.csv")

print("  Sond w ekspresji:", expr["ID_REF"].nunique())
print("  Sond z mapowaniem:", expr["ID_REF"].isin(probe_to_gene).sum())
print("  Geny po agregacji:", len(gene_expr))
print("  Grupy:", ", ".join(meta["group"].drop_duplicates()))

# 4. Wczytaj piki z TAR i policz nakładanie na promotory metodą sweep-line.
def count_overlaps(promoter_intervals, peak_intervals):
    # Zwraca liczbę unikalnych pików nakładających się na promotory każdego genu.
    events = []
    for i, (start, end, gene) in enumerate(promoter_intervals):
        events.append((start, 0, end, i, gene))
    for i, (start, end) in enumerate(peak_intervals):
        events.append((start, 1, end, i, None))
    events.sort(key=lambda x: (x[0], x[1]))

    active_prom = {}
    active_peaks = {}
    heap_prom, heap_peaks = [], []
    counts = defaultdict(int)

    for start, kind, end, idx, gene in events:
        while heap_prom and heap_prom[0][0] <= start:
            _, old = heapq.heappop(heap_prom)
            active_prom.pop(old, None)
        while heap_peaks and heap_peaks[0][0] <= start:
            _, old = heapq.heappop(heap_peaks)
            active_peaks.pop(old, None)

        if kind == 0:  # nowy promotor
            for peak_id in active_peaks:
                counts[gene] += 1
            active_prom[idx] = (end, gene)
            heapq.heappush(heap_prom, (end, idx))
        else:  # nowy pik
            for _, active_gene in active_prom.values():
                counts[active_gene] += 1
            active_peaks[idx] = end
            heapq.heappush(heap_peaks, (end, idx))

    return counts

prom_by_chrom = defaultdict(list)
for chrom, start, end, gene in merged_promoters:
    prom_by_chrom[chrom].append((start, end, gene))

def modality_stage(filename):
    base = Path(filename).name
    m = re.search(
        r"_(Flag|Oct4|K4me1|K27ac|K4me3|K27me3|RNAPII)_(d\d+|iPSC)_peaks\.",
        base, re.I
    )
    if not m:
        return None
    mark, stage = m.group(1).lower(), m.group(2).lower()
    mark_map = {
        "oct4": "OCT4",
        "k4me1": "H3K4me1",
        "k27ac": "H3K27ac",
        "k4me3": "H3K4me3",
        "k27me3": "H3K27me3",
        "rnapii": "RNAPII",
        "flag": "FLAG_tagged_exogenous_OCT4",
    }
    stage = "iPSC" if stage == "ipsc" else "day" + stage[1:]
    return mark_map[mark], stage

print("4/5 Czytam pliki pików i liczę nakładanie...")
records = []
with tarfile.open(TAR, "r") as tar:
    members = [m for m in tar.getmembers() if m.isfile()]
    peak_members = [m for m in members if modality_stage(m.name)]

    for mi, member in enumerate(peak_members, 1):
        modality, stage = modality_stage(member.name)
        raw = tar.extractfile(member)
        if raw is None:
            continue

        peaks_by_chrom = defaultdict(list)
        with gzip.GzipFile(fileobj=raw) as gz:
            for line in gz:
                parts = line.decode("utf-8", errors="replace").strip().split()
                if len(parts) < 3:
                    continue
                try:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                except (ValueError, IndexError):
                    continue
                if end > start:
                    peaks_by_chrom[chrom].append((start, end))

        counts = defaultdict(int)
        n_peaks = sum(map(len, peaks_by_chrom.values()))
        for chrom, peaks in peaks_by_chrom.items():
            if chrom in prom_by_chrom:
                part = count_overlaps(prom_by_chrom[chrom], peaks)
                for gene, count in part.items():
                    counts[gene] += count

        for gene in gene_symbols:
            records.append({
                "gene": gene,
                "stage": stage,
                "modality": modality,
                "peak_count": counts.get(gene, 0),
                "peak_present": int(counts.get(gene, 0) > 0),
                "n_peaks_in_file": n_peaks,
                "source_file": Path(member.name).name,
            })
        print(f"  [{mi}/{len(peak_members)}] {modality} {stage}: {n_peaks} pików")

# 5. Zapis wyników oraz raportów.
print("5/5 Zapisuję wyniki...")
peaks_df = pd.DataFrame(records)
peaks_df.to_csv(OUT / "04_promoter_peak_counts_long.csv", index=False)

presence = peaks_df.pivot_table(
    index="gene", columns=["modality", "stage"],
    values="peak_present", aggfunc="max", fill_value=0
)
presence.to_csv(OUT / "05_promoter_peak_presence_matrix.csv")

count_matrix = peaks_df.pivot_table(
    index="gene", columns=["modality", "stage"],
    values="peak_count", aggfunc="sum", fill_value=0
)
count_matrix.to_csv(OUT / "06_promoter_peak_count_matrix.csv")

summary = (
    peaks_df.groupby(["modality", "stage"], as_index=False)
    .agg(
        genes_with_peak=("peak_present", "sum"),
        genes_tested=("gene", "nunique"),
        total_peak_overlaps=("peak_count", "sum"),
        peaks_in_source_file=("n_peaks_in_file", "first"),
    )
)
summary.to_csv(OUT / "07_modality_coverage_summary.csv", index=False)

print("\nGOTOWE. Wyniki:", OUT)
print("Pliki:")
for p in sorted(OUT.glob("*.csv")):
    print(" ", p.name, f"({p.stat().st_size:,} B)")
print("\nUWAGA: peak_count to liczba nakładających się przedziałów pików,")
print("a peak_present to wskaźnik 0/1. Żadna z tych miar nie jest")
print("bezpośrednim pomiarem ilościowym siły wiązania.")


