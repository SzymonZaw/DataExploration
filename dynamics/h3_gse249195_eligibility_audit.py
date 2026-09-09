from __future__ import annotations

import gzip, hashlib, json, re, urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"
GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE249195&file=GSE249195_raw_counts.txt.gz&format=file"
DATA_FILE = DATA / "GSE249195_raw_counts.txt.gz"
EXPECTED_SAMPLES = 24
EXPECTED_TIMEPOINTS = np.array([0.0, 1.0, 2.0, 4.0, 7.0, 14.0])
EXPECTED_REPLICATES = 4
MIN_PROGENY_OVERLAP = 0.90
MIN_DOROTHEA_OVERLAP = 0.90
MIN_PC1_ABS_SPEARMAN = 0.80
MIN_DIRECTIONAL_GENE_FRACTION = 0.60
MIN_MEDIAN_ABS_GENE_SPEARMAN = 0.60
MIN_TRANSITION_MEDIAN_ABS_LOG2FC = 0.50
H2_NUISANCE_FAMILIES = {"JAK-STAT", "NFkB", "TNFa", "p53", "Hypoxia"}
SAMPLE_RE = re.compile(r"^S1_D(\d+(?:\.\d+)?)_R(\d+)$", re.I)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def acquire() -> dict:
    DATA.mkdir(parents=True, exist_ok=True)
    if not DATA_FILE.exists():
        urllib.request.urlretrieve(GEO_URL, DATA_FILE)
    return {"url": GEO_URL, "local_path": str(DATA_FILE), "sha256": sha256(DATA_FILE), "bytes": DATA_FILE.stat().st_size}


def read_counts() -> pd.DataFrame:
    # GEO supplementary file is expected to be a tab-delimited gene-by-sample matrix.
    d = pd.read_csv(DATA_FILE, sep="\t", compression="gzip", low_memory=False)
    sample_cols = [c for c in d.columns if SAMPLE_RE.fullmatch(str(c))]
    if len(sample_cols) != EXPECTED_SAMPLES:
        raise ValueError(f"Expected {EXPECTED_SAMPLES} sample columns, found {len(sample_cols)}")
    gene_col = d.columns[0]
    x = d[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    x.index = d[gene_col].astype(str).str.strip()
    x = x.groupby(level=0).sum()
    return x


def design_audit(columns: list[str]) -> dict:
    rows = []
    for c in columns:
        m = SAMPLE_RE.fullmatch(str(c))
        if not m:
            continue
        rows.append((float(m.group(1)), int(m.group(2))))
    d = pd.DataFrame(rows, columns=["day", "replicate"])
    counts = d.groupby("day").replicate.nunique().sort_index()
    complete = np.array_equal(counts.index.to_numpy(float), EXPECTED_TIMEPOINTS) and np.all(counts.to_numpy() == EXPECTED_REPLICATES)
    return {
        "n_samples": int(len(d)),
        "timepoints": counts.index.astype(float).tolist(),
        "replicates_per_timepoint": {str(k): int(v) for k, v in counts.items()},
        "expected_timepoints": EXPECTED_TIMEPOINTS.tolist(),
        "complete_design": bool(complete),
        "independent_trajectories": EXPECTED_REPLICATES,
        "technical_replicates_treated_as_independent": False,
        "time_span_days": float(EXPECTED_TIMEPOINTS[-1] - EXPECTED_TIMEPOINTS[0]),
    }


def map_hgnc(x: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    ids = pd.Index(x.index.astype(str))
    ensembl_fraction = sum(bool(re.fullmatch(r"ENSG\d+(?:\.\d+)?", i.upper())) for i in ids) / max(len(ids), 1)
    if ensembl_fraction < 0.80:
        y = x.copy()
        y.index = ids.str.upper()
        before = len(y)
        y = y.groupby(level=0).mean()
        return y, {"input_rows": before, "mapped_rows": len(y), "mapping_fraction": 1.0, "mapping_mode": "already_symbol_like", "duplicate_symbol_collisions": int(before - len(y))}

    try:
        import mygene
    except ImportError as exc:
        raise RuntimeError("mygene is required for Ensembl-to-HGNC mapping") from exc
    mg = mygene.MyGeneInfo()
    clean = sorted({i.split(".", 1)[0].upper() for i in ids})
    mapping: dict[str, str] = {}
    for start in range(0, len(clean), 1000):
        res = mg.querymany(clean[start:start + 1000], scopes="ensembl.gene", fields="symbol", species="human", as_dataframe=False, verbose=False)
        for r in res:
            if not r.get("notfound") and r.get("query") and r.get("symbol"):
                mapping[str(r["query"]).upper()] = str(r["symbol"]).upper()
    symbols = [mapping.get(i.split(".", 1)[0].upper(), "") for i in ids]
    keep = np.array([bool(s) for s in symbols])
    y = x.loc[keep].copy()
    y.index = np.asarray(symbols)[keep]
    before = len(y)
    y = y.groupby(level=0).mean()
    return y, {
        "input_rows": int(len(ids)),
        "mapped_rows": int(len(y)),
        "mapping_fraction": float(keep.mean()),
        "mapping_mode": "ensembl_to_hgnc_mygene",
        "duplicate_symbol_collisions": int(before - len(y)),
    }


def network_overlap(expr: pd.DataFrame) -> dict:
    import decoupler as dc
    out = {}
    progeny = dc.op.progeny(organism="human", top=100)
    dorothea = dc.op.dorothea(organism="human", levels=["A", "B", "C"])
    for family, net in (("PROGENy", progeny), ("DoRothEA", dorothea)):
        targets = pd.Index(net["target"].astype(str).str.upper().unique())
        overlap = len(expr.index.intersection(targets)) / max(len(targets), 1)
        out[family] = {"network_targets": int(len(targets)), "common_genes": int(len(expr.index.intersection(targets))), "overlap_fraction": float(overlap)}
    return out


def log_cpm(expr: pd.DataFrame) -> pd.DataFrame:
    lib = expr.sum(axis=0).replace(0, np.nan)
    return np.log2(expr.div(lib, axis=1) * 1e6 + 1.0)


def temporal_geometry(expr: pd.DataFrame, columns: list[str]) -> dict:
    from scipy.stats import spearmanr
    lcpm = log_cpm(expr[columns])
    meta = []
    for c in columns:
        m = SAMPLE_RE.fullmatch(str(c))
        if m:
            meta.append((c, float(m.group(1))))
    md = pd.DataFrame(meta, columns=["sample", "day"]).set_index("sample")
    day_mean = lcpm.T.join(md).groupby("day").mean(numeric_only=True).T
    x = day_mean.columns.to_numpy(float)
    rhos = []
    for _, row in day_mean.iterrows():
        if np.std(row.to_numpy(float)) == 0:
            continue
        rhos.append(float(spearmanr(x, row.to_numpy(float)).statistic))
    rhos = np.asarray(rhos)
    directional = np.abs(rhos) >= 0.80
    transition = day_mean.iloc[:, -1] - day_mean.iloc[:, 0]
    return {
        "pc1_abs_spearman": None,
        "gene_trajectory_count": int(len(rhos)),
        "directional_gene_fraction_abs_spearman_ge_0_80": float(directional.mean()) if len(rhos) else 0.0,
        "median_abs_gene_spearman": float(np.median(np.abs(rhos))) if len(rhos) else 0.0,
        "median_abs_day14_day0_log2fc": float(np.median(np.abs(transition.to_numpy(float)))),
        "note": "Gene-level temporal geometry is evaluated on log2-CPM day means. PC1 criterion is intentionally reported as pending because it requires a locked PCA implementation shared with the project-wide representation code.",
    }


def h2_nuisance_audit(expr: pd.DataFrame, columns: list[str]) -> dict:
    # H2 is a separation audit, not a post-hoc exclusion based on the candidate result.
    import decoupler as dc
    progeny = dc.op.progeny(organism="human", top=100)
    n = progeny.copy()
    n["target"] = n["target"].astype(str).str.upper()
    common = expr.index.intersection(n["target"])
    r = dc.mt.ulm(data=expr.loc[common, columns].T, net=n[n.target.isin(common)], tmin=5)
    if isinstance(r, tuple):
        r = r[0]
    if not isinstance(r, pd.DataFrame) or not {"source", "score"}.issubset(r.columns):
        return {"status": "AUDIT_UNAVAILABLE", "reason": "Unexpected decoupler ULM output format"}
    activity = r.pivot(index=r.index, columns="source", values="score")
    time = np.array([float(SAMPLE_RE.fullmatch(str(c)).group(1)) for c in columns])
    from scipy.stats import spearmanr
    rows = []
    for family in H2_NUISANCE_FAMILIES:
        if family not in activity.columns:
            continue
        rho = float(spearmanr(time, activity[family].to_numpy(float)).statistic)
        rows.append({"family": family, "spearman_time": rho, "abs_spearman_time": abs(rho)})
    return {"status": "COMPUTED", "nuisance_families_present": rows, "interpretation": "Strong nuisance-axis dynamics do not by themselves reject the control; they trigger documented H2 review before eligibility is frozen."}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    acquisition = acquire()
    raw = read_counts()
    design = design_audit([str(c) for c in raw.columns])
    expr, mapping = map_hgnc(raw)
    overlap = network_overlap(expr)
    geometry = temporal_geometry(expr, [str(c) for c in raw.columns])
    h2 = h2_nuisance_audit(expr, [str(c) for c in raw.columns])

    representation_pass = overlap["PROGENy"]["overlap_fraction"] >= MIN_PROGENY_OVERLAP and overlap["DoRothEA"]["overlap_fraction"] >= MIN_DOROTHEA_OVERLAP
    geometry_pass = (
        geometry["directional_gene_fraction_abs_spearman_ge_0_80"] >= MIN_DIRECTIONAL_GENE_FRACTION
        and geometry["median_abs_gene_spearman"] >= MIN_MEDIAN_ABS_GENE_SPEARMAN
        and geometry["median_abs_day14_day0_log2fc"] >= MIN_TRANSITION_MEDIAN_ABS_LOG2FC
    )
    design_pass = design["complete_design"] and design["time_span_days"] >= 7.0
    pc1_pending = True
    eligibility = "ELIGIBILITY_PENDING_PC1_AUDIT" if design_pass and representation_pass and geometry_pass else "INELIGIBLE"

    result = {
        "status": "ok",
        "accession": "GSE249195",
        "source": "NCBI GEO",
        "acquisition": acquisition,
        "design": design,
        "mapping": mapping,
        "representation_overlap": overlap,
        "temporal_geometry": geometry,
        "h2_nuisance_audit": h2,
        "locked_thresholds": {
            "min_progeny_overlap": MIN_PROGENY_OVERLAP,
            "min_dorothea_overlap": MIN_DOROTHEA_OVERLAP,
            "min_directional_gene_fraction_abs_spearman_ge_0_80": MIN_DIRECTIONAL_GENE_FRACTION,
            "min_median_abs_gene_spearman": MIN_MEDIAN_ABS_GENE_SPEARMAN,
            "min_median_abs_day14_day0_log2fc": MIN_TRANSITION_MEDIAN_ABS_LOG2FC,
            "min_time_span_days": 7.0,
            "pc1_abs_spearman_threshold": MIN_PC1_ABS_SPEARMAN,
        },
        "decision": {
            "status": eligibility,
            "design_pass": bool(design_pass),
            "representation_pass": bool(representation_pass),
            "directional_geometry_pass": bool(geometry_pass),
            "pc1_audit_complete": bool(not pc1_pending),
            "reason": "GSE249195 is a promising directional days-scale differentiation control, but it must pass the locked PC1 temporal-geometry check and documented H2 review before similarity is allowed." if eligibility != "INELIGIBLE" else "Candidate failed one or more locked design/representation/temporal-geometry gates; similarity is not permitted.",
        },
        "guardrail": "No similarity/null analysis is run by this eligibility script. Eligibility must be frozen before any outcome-inspecting similarity analysis.",
    }
    (OUT / "H3_GSE249195_ELIGIBILITY_AUDIT.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
