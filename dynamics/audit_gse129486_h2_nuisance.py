"""Audit predeclared H2 nuisance axes in GSE129486 without H3 comparison."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"
TPM = DATA / "GSE129486_rnaseq-data-1_gene-tpm.tsv.gz"
META = DATA / "GSE129486_rnaseq-data-1_metadata.tsv.gz"
EXPECTED = {TPM.name: "6e3d7860f4f38d95830a15b8dd570d58f9226170df6002343e432a05c96fcbf0", META.name: "78a60e353461ab819672e478a9f38b20822fdfc493a1773677ddef3eb167ff04"}

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def load_mygene_map(ids):
    import mygene
    info = mygene.MyGeneInfo()
    clean = [str(x).split(".")[0].upper() for x in ids]
    out = {}
    for i in range(0, len(clean), 1000):
        result = info.querymany(clean[i:i + 1000], scopes="ensembl.gene", fields="symbol", species="human", as_dataframe=False, returnall=False, step=1000, verbose=False)
        for row in result:
            if row.get("notfound"):
                continue
            q, symbol = row.get("query"), row.get("symbol")
            if q and symbol:
                out[str(q).upper()] = str(symbol).upper()
    return out

def load_inputs():
    meta = pd.read_csv(META, sep="\t")
    raw = pd.read_csv(TPM, sep="\t")
    gene_col = raw.columns[0]
    sample_cols = [c for c in raw.columns if c != gene_col]
    raw[gene_col] = raw[gene_col].astype(str).str.split(".").str[0].str.upper()
    mapping = load_mygene_map(raw[gene_col].tolist())
    raw["HGNC"] = raw[gene_col].map(mapping)
    raw = raw[raw["HGNC"].notna()].copy()
    expr = raw.groupby("HGNC")[sample_cols].mean().T
    expr.index = expr.index.astype(str)
    meta["sample"] = meta["sample"].astype(str)
    meta = meta.set_index("sample")
    common = expr.index.intersection(meta.index)
    return expr.loc[common], meta.loc[common]

def pathway_scores(expr_hgnc):
    """Score PROGENy with the current decoupler API."""
    import decoupler as dc
    net = dc.op.progeny(organism="human", top=100)
    # Current decoupler Method API takes data/net; source/target/weight are
    # network columns and must not be passed as kwargs to the low-level ULM.
    result = dc.mt.ulm(data=expr_hgnc, net=net, tmin=1)
    if isinstance(result, tuple):
        result = result[0]
    if isinstance(result, dict):
        result = result.get("estimate", result.get("score", result))
    if not isinstance(result, pd.DataFrame):
        result = pd.DataFrame(result)
    if "sample" in result.columns:
        result = result.set_index("sample")
    result.index = result.index.astype(str)
    wanted = [c for c in ["JAK-STAT", "NFkB", "NF-kB", "NFKB"] if c in result.columns]
    return result.loc[:, wanted].copy() if wanted else pd.DataFrame(index=expr_hgnc.index)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"status": "ok", "scope": "GSE129486 nuisance-axis audit only; no Yamanaka comparison, H3 similarity, trajectory-agreement statistic, or H3 decision", "input_hashes": {}}
    for path in (TPM, META):
        actual = sha256(path)
        report["input_hashes"][path.name] = {"sha256": actual, "expected": EXPECTED[path.name], "matches": actual == EXPECTED[path.name]}
    expr, meta = load_inputs()
    scores = pathway_scores(expr)
    joined = meta.join(scores, how="inner")
    report["data"] = {"n_samples": int(len(joined)), "n_mapped_hgnc": int(expr.shape[1]), "n_cell_lines": int(meta["cell_line"].nunique()), "n_stimulations": int(meta["stimulation"].nunique()), "n_timepoints": int(meta["time"].nunique())}
    report["nuisance_axes"] = {}
    for axis in scores.columns:
        s = joined[axis].astype(float)
        report["nuisance_axes"][axis] = {"mean": float(s.mean()), "std": float(s.std(ddof=1)), "min": float(s.min()), "max": float(s.max()), "mean_by_time": {str(k): float(v) for k, v in joined.groupby("time")[axis].mean().items()}, "mean_by_stimulation": {str(k): float(v) for k, v in joined.groupby("stimulation")[axis].mean().items()}, "mean_by_cell_line": {str(k): float(v) for k, v in joined.groupby("cell_line")[axis].mean().items()}, "time_slope": float(np.polyfit(joined["time"].astype(float), s, 1)[0]) if joined["time"].nunique() > 1 else None}
    report["interpretation"] = "These scores quantify predeclared nuisance-axis activity in GSE129486. They do not by themselves establish confounding of the Yamanaka representation."
    report["guardrail"] = "No H3 similarity or biological decision is computed."
    (OUT / "H3_GSE129486_H2_NUISANCE_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
