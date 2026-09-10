"""Robustness audit for the GSE242421 frozen-module orthogonal Z4 test.

Reuses the already derived peak-by-sample CPM matrix and re-derives gene
activity from the published peak-gene links at progressively stricter
absolute-correlation thresholds. No target-derived fitting or feature
selection is performed.
"""
from __future__ import annotations
import argparse, json, re, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DEFAULT_PEAK = Path("results/Dynamics/z4_gse242421_gene_activity/02_peak_sample_cpm.tsv.gz")
DEFAULT_INTEGRATION = Path("Data/GSE242421/scATAC_scRNA_integration.zip")
DEFAULT_MODULES = Path("results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv")
DEFAULT_REFERENCE = Path("results/Dynamics/z6_gse67462_temporal_modules/02_temporal_modules.csv")
DEFAULT_MAPPING = Path("results/Dynamics/z4_frozen_orthology_mapping/01_frozen_human_mouse_mapping.csv")
DEFAULT_OUT = Path("results/Dynamics/z4_gse242421_peak_link_robustness")
TIME_ORDER = ["D0","D2","D4","D6","D8","D10","D12","D14"]
ALL_ORDER = TIME_ORDER + ["iPSC"]
THRESHOLDS = [0.45, 0.50, 0.60, 0.70]


def norm(x): return str(x).strip().upper()

def day(x):
    m = re.search(r"D(\d+(?:\.\d+)?)", str(x), re.I)
    return float(m.group(1)) if m else np.nan

def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan
    return float(spearmanr(x[m], y[m]).statistic)

def ref_direction(ref, module):
    g = ref[pd.to_numeric(ref.module, errors="coerce") == module]
    if g.empty: return None
    label = str(g.trajectory_label.iloc[0]).lower()
    if label == "late_rising": return 1.0
    if label == "early_declining": return -1.0
    return None

def read_links(path):
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        target = next((n for n in names if n.endswith("peak_gene_links_fdr1e-4.tsv")), None)
        if target is None:
            raise FileNotFoundError("peak_gene_links_fdr1e-4.tsv not found in integration zip")
        with z.open(target) as f:
            links = pd.read_csv(f, sep="\t")
    cols = {c.lower(): c for c in links.columns}
    peak_col = next((cols[k] for k in ["peak","peak_id","peakid"] if k in cols), None)
    gene_col = next((cols[k] for k in ["gene","gene_name","gene_symbol"] if k in cols), None)
    corr_col = next((cols[k] for k in ["correlation","corr","score","rho","r"] if k in cols), None)
    if not all([peak_col, gene_col, corr_col]):
        raise ValueError(f"Cannot identify peak/gene/correlation columns: {list(links.columns)}")
    out = links[[peak_col, gene_col, corr_col]].copy()
    out.columns = ["peak", "gene", "correlation"]
    out["peak"] = out.peak.astype(str)
    out["gene"] = out.gene.map(norm)
    out["correlation"] = pd.to_numeric(out.correlation, errors="coerce")
    return out.dropna(subset=["peak","gene","correlation"])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peak-cpm", type=Path, default=DEFAULT_PEAK)
    ap.add_argument("--integration-zip", type=Path, default=DEFAULT_INTEGRATION)
    ap.add_argument("--modules", type=Path, default=DEFAULT_MODULES)
    ap.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    ap.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--permutations", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=424242)
    ap.add_argument("--thresholds", type=float, nargs="+", default=THRESHOLDS)
    a = ap.parse_args(); a.output.mkdir(parents=True, exist_ok=True)
    for p in [a.peak_cpm, a.integration_zip, a.modules, a.reference, a.mapping]:
        if not p.exists(): raise FileNotFoundError(p)

    peak = pd.read_csv(a.peak_cpm, sep="\t", compression="gzip", index_col=0)
    peak.index = peak.index.astype(str)
    peak = peak[~peak.index.duplicated(keep="first")]
    if any(c not in peak.columns for c in ALL_ORDER):
        raise ValueError("Peak CPM matrix missing required samples")
    links = read_links(a.integration_zip)
    modules = pd.read_csv(a.modules)
    mapping = pd.read_csv(a.mapping)
    mapping["mouse"] = mapping.mouse_gene.map(norm); mapping["human"] = mapping.human_gene.map(norm)
    mouse_to_human = dict(zip(mapping.mouse, mapping.human))
    modules["gene_norm"] = modules.gene.map(norm)
    ref = pd.read_csv(a.reference)
    time = np.array([day(x) for x in TIME_ORDER], float)
    rng = np.random.default_rng(a.seed)
    module_ids = sorted(pd.to_numeric(modules.module, errors="coerce").dropna().astype(int).unique())
    rows=[]
    for threshold in a.thresholds:
        l = links[links.correlation.abs() >= threshold].copy()
        l = l[l.peak.isin(peak.index)]
        # Same binary peak->gene aggregation rule as the primary derivation.
        l = l.drop_duplicates(["peak","gene"])
        for module in module_ids:
            mg = modules.loc[pd.to_numeric(modules.module, errors="coerce") == module, "gene_norm"].unique()
            human = sorted({mouse_to_human[g] for g in mg if g in mouse_to_human})
            gene_peaks = l[l.gene.isin(human)]
            validated = sorted(gene_peaks.gene.unique())
            if not validated:
                continue
            # Gene activity = mean accessibility of linked peaks for each gene,
            # followed by mean across genes, matching the frozen derivation logic.
            gene_scores = []
            for gene in validated:
                peaks = gene_peaks.loc[gene_peaks.gene == gene, "peak"].unique()
                gene_scores.append(peak.loc[peaks, ALL_ORDER].mean(axis=0).to_numpy(float))
            score = np.vstack(gene_scores).mean(axis=0)
            obs = spearman(time, score[:len(TIME_ORDER)])
            direction = ref_direction(ref, module)
            p = np.nan; concordant = np.nan
            if direction is not None and np.isfinite(obs):
                null = np.array([spearman(time, rng.permutation(score[:len(TIME_ORDER)])) for _ in range(a.permutations)])
                signed_obs = direction * obs
                p = float((np.sum(direction * null >= signed_obs) + 1) / (len(null) + 1))
                concordant = bool(np.sign(obs) == np.sign(direction) and obs != 0)
            rows.append({"threshold":threshold,"module":module,"frozen_mouse_genes":len(mg),"validated_human_genes":len(validated),"coverage_fraction":len(validated)/max(len(mg),1),"target_time_spearman":obs,"frozen_reference_direction":direction,"direction_concordant":concordant,"permutation_p_directional":p,"n_links":len(l)})
    summary = pd.DataFrame(rows).sort_values(["threshold","module"])
    summary.to_csv(a.output/"01_threshold_module_summary.csv", index=False)
    key = summary[summary.module.isin([4,6])].copy()
    key.to_csv(a.output/"02_key_modules_4_6.csv", index=False)
    passes = key.direction_concordant.fillna(False) & (key.permutation_p_directional < 0.05)
    by_thr = []
    for threshold, g in key.groupby("threshold"):
        by_thr.append({"threshold":float(threshold),"modules_tested":int(g.direction_concordant.notna().sum()),"modules_passed":int(passes.loc[g.index].sum()),"both_4_6_pass":bool(passes.loc[g.index].all()) if len(g)==2 else False})
    robustness = pd.DataFrame(by_thr)
    robustness.to_csv(a.output/"03_robustness_by_threshold.csv", index=False)
    all_pass = bool(not robustness.empty and robustness.both_4_6_pass.all())
    decision = "Z4_ORTHO_ROBUST" if all_pass else "Z4_ORTHO_THRESHOLD_SENSITIVE"
    result={"candidate":"GSE242421","thresholds":[float(x) for x in a.thresholds],"primary_threshold":0.45,"key_modules":[4,6],"permutations":a.permutations,"gene_activity_derivation":"reuse_peak_sample_cpm_with_binary_peak_gene_links","module_fitting_on_target":False,"feature_selection_on_target":False,"decision":decision,"interpretation":"ROBUSTNESS_AUDIT_ONLY"}
    (a.output/"04_summary.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    print("GSE242421 Z4 PEAK-LINK ROBUSTNESS AUDIT")
    print(f"thresholds: {', '.join(map(str,a.thresholds))}")
    print(f"decision: {decision}")
    print(f"output: {a.output}")

if __name__ == "__main__": main()
