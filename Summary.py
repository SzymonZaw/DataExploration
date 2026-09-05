from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "SUMMARY"
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "GSE148158": {
        "title": "A PIANO (Proper, Insufficient, Aberrant, and NO reprogramming) response to the Yamanaka factors in the initial stages of human iPSC reprogramming",
        "type": "RNA-seq / expression profiling by high-throughput sequencing",
        "platform": "GPL16791 - Illumina HiSeq 2500",
        "samples": 13,
        "description": "Human fibroblasts (BJ), hESC controls, GFP controls and cells undergoing OSKM reprogramming. The dataset describes early transcriptional responses during induced pluripotency.",
        "interpretation": "This dataset is useful for examining how global gene-expression profiles change from fibroblast and control states toward the early OSKM reprogramming state. The input file already contains GEO normalized counts; log2(x+1) is used only for exploratory analyses and PCA.",
    },
    "GSE28688": {
        "title": "Molecular insights into induced pluripotency mediated by the OCT4, SOX2, KLF and c-MYC gene regulatory network",
        "type": "Expression profiling by array",
        "platform": "GPL6883 - Illumina HumanRef-8 v3.0 expression beadchip",
        "samples": 14,
        "description": "Human HFF1 fibroblasts measured before and 24, 48 and 72 hours after transduction, together with H1/H9 hESC controls and iPS2/iPS4 samples.",
        "interpretation": "This dataset allows exploration of the temporal transcriptional response after reprogramming-factor transduction and comparison with established pluripotent states. Replicate pairs make within-group consistency particularly useful to inspect.",
    },
    "GSE52052": {
        "title": "Comparison of gene expression at day 11 of iPSC reprogramming with different reprogramming cocktails",
        "type": "Expression profiling by array",
        "platform": "GPL14550 - Agilent-028004 SurePrint G3 Human GE 8x60K Microarray",
        "samples": 6,
        "description": "Day-11 reprogramming samples generated with OSK, OSK+let-7 inhibitor, OSKM or OSK+LIN-41, plus a GFP control and H1 hESC reference.",
        "interpretation": "This dataset compares global expression at the same reprogramming time point under different molecular cocktails. It is especially useful for visualizing how strongly the cocktails separate transcriptional states, while remembering that most conditions have only one sample.",
    },
}

GROUPS = {
    "GSE148158": {
        "BJ_2 [re-analysis]": "BJ_fibroblast", "BJ_1 [re-analysis]": "BJ_fibroblast",
        "BJ_3 [re-analysis]": "BJ_fibroblast", "BJ_4 [re-analysis]": "BJ_fibroblast",
        "H1_2 [re-analysis]": "hESC", "H9 [re-analysis]": "hESC", "H1 [re-analysis]": "hESC",
        "BJ_GFP48": "GFP_48h", "BJ_GFP48b": "GFP_48h",
        "BJ_GFP72": "GFP_72h", "BJ_GFP72b": "GFP_72h",
        "OSKM48": "OSKM_48h", "OSKM72": "OSKM_72h",
    },
    "GSE28688": {
        "HFF1-a": "HFF1", "HFF1-b": "HFF1",
        "HFF1-24 h post-transduction-a": "24h", "HFF1-24 h post-transduction-b": "24h",
        "HFF1-48 h post-transduction-a": "48h", "HFF1-48 h post-transduction-b": "48h",
        "HFF1-72 h post-transduction-a": "72h", "HFF1-72 h post-transduction-b": "72h",
        "H1": "H1_hESC", "H9": "H9_hESC",
        "iPS2 from HFF1-a": "iPS2", "iPS2 from HFF1-b": "iPS2",
        "iPS4 from HFF1-a": "iPS4", "iPS4 from HFF1-b": "iPS4",
    },
    "GSE52052": {
        "HDF_GFP(+)_day11": "GFP_control",
        "HDF_OSK+control_inh_TRA(+)_day11": "OSK_control_inhibitor",
        "HDF_OSK+let-7_inh_TRA(+)_day11": "OSK_let7_inhibitor",
        "HDF_OSKM_TRA(+)_day11": "OSKM",
        "HDF_OSK+LIN-41_TRA(+)_day11": "OSK_LIN41",
        "H1_hESC": "H1_hESC",
    },
}


def read_csv(path, index_col=0):
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, index_col=index_col)
    except Exception:
        return None


def fmt(value, digits=3):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "n/a"
    return f"{value:.{digits}f}"


def summarize_dataset(name):
    d = DATASETS[name]
    folder = RESULTS / name
    lines = []
    lines.append(f"{name}")
    lines.append("=" * len(name))
    lines.append(f"Title: {d['title']}")
    lines.append(f"Experiment: {d['type']}")
    lines.append(f"Platform: {d['platform']}")
    lines.append(f"Expected samples from GEO design: {d['samples']}")
    lines.append(f"What the data represent: {d['description']}")
    lines.append(f"Why this dataset matters: {d['interpretation']}")
    lines.append("")

    if not folder.exists():
        lines.append("RESULTS NOT FOUND: run the corresponding exploration script first.")
        return lines

    meta = read_csv(folder / "01_sample_metadata.csv", index_col=None)
    qc = read_csv(folder / "03_sample_QC.csv")
    corr = read_csv(folder / "04_sample_correlation.csv")
    var = read_csv(folder / "06_feature_variance.csv")
    pca = read_csv(folder / "07_PCA_coordinates.csv")

    lines.append("Exploration performed")
    lines.append("----------------------")
    lines.append("- Distribution check: boxplot of expression values per sample.")
    lines.append("- Sample QC: mean, median, standard deviation, minimum, maximum and missing values.")
    lines.append("- Sample correlation: Pearson correlation between complete sample expression profiles.")
    lines.append("- Variability: ranking features by variance after log2(x+1).")
    lines.append("- PCA: dimensionality reduction using the most variable features.")
    lines.append("- Top-variable feature heatmap: standardized patterns across samples.")
    lines.append("")

    if meta is not None and "group" in meta.columns:
        counts = meta["group"].value_counts().sort_index()
        lines.append("Sample groups observed")
        lines.append("----------------------")
        for group, count in counts.items():
            lines.append(f"- {group}: {count} sample(s)")
        lines.append("")

    if qc is not None and not qc.empty:
        lines.append("QC summary")
        lines.append("----------")
        if "missing" in qc.columns:
            lines.append(f"- Maximum missing values in one sample: {int(qc['missing'].max())}")
        if "mean" in qc.columns:
            lines.append(f"- Sample mean range: {fmt(qc['mean'].min())} to {fmt(qc['mean'].max())}")
        if "sd" in qc.columns:
            lines.append(f"- Sample SD range: {fmt(qc['sd'].min())} to {fmt(qc['sd'].max())}")
        lines.append("")

    if corr is not None and not corr.empty:
        numeric = corr.apply(pd.to_numeric, errors="coerce")
        vals = numeric.to_numpy(dtype=float)
        if vals.size:
            off_diag = vals[~np.eye(vals.shape[0], dtype=bool)]
            off_diag = off_diag[np.isfinite(off_diag)]
            if len(off_diag):
                lines.append("Correlation summary")
                lines.append("-------------------")
                lines.append(f"- Pearson r between different samples: min {fmt(off_diag.min())}, median {fmt(np.median(off_diag))}, max {fmt(off_diag.max())}")
                lines.append("- High correlation suggests similar global expression profiles; lower correlation suggests stronger biological or technical separation.")
                lines.append("")

    if var is not None and not var.empty:
        variance_col = var.columns[0]
        vv = pd.to_numeric(var[variance_col], errors="coerce").dropna()
        if len(vv):
            lines.append("Feature variability")
            lines.append("-------------------")
            lines.append(f"- Features ranked by variance: {len(vv):,}")
            lines.append(f"- Highest observed variance: {fmt(vv.iloc[0])}")
            if len(vv) >= 10:
                lines.append("- The highest-variance features are the strongest contributors to visible sample differences and PCA structure.")
            lines.append("")

    if pca is not None and not pca.empty and "PC1" in pca.columns:
        lines.append("PCA interpretation")
        lines.append("------------------")
        if "PC2" in pca.columns:
            lines.append("- PC1 and PC2 summarize major axes of variation among samples.")
        if "group" in pca.columns:
            group_counts = pca["group"].value_counts()
            singleton = int((group_counts == 1).sum())
            if singleton:
                lines.append(f"- {singleton} group(s) contain only one sample; separation of those groups is descriptive rather than evidence of reproducibility.")
        lines.append("- Samples that cluster together have more similar global expression patterns; separated samples have larger overall expression differences.")
        lines.append("")

    lines.append("Important interpretation limits")
    lines.append("-------------------------------")
    lines.append("- These scripts perform exploratory data analysis, not differential-expression testing.")
    lines.append("- A visible PCA separation does not by itself prove a biological mechanism or statistical significance.")
    lines.append("- Correlation and clustering reflect the complete expression profile and can be influenced by technical effects as well as biology.")
    if name == "GSE52052":
        lines.append("- GSE52052 has one sample for each main reprogramming condition, so group-level biological conclusions should be treated as descriptive.")
    if name == "GSE148158":
        lines.append("- GSE148158 contains singleton OSKM time points; these are useful for visualization but do not provide replicate-based estimates at those exact time points.")
    if name == "GSE28688":
        lines.append("- GSE28688 has replicate pairs for the HFF1/time-course and iPS groups, which supports inspection of within-group consistency.")

    return lines


all_lines = [
    "SUMMARY OF THE THREE GENE-EXPRESSION EXPLORATIONS",
    "=" * 48,
    "",
    "Purpose",
    "-------",
    "This report explains what the three GEO datasets contain and what was obtained from the exploratory scripts. The analysis focuses on data quality, sample similarity, variability and global structure of the expression data.",
    "",
    "How to read the results",
    "------------------------",
    "Boxplots show the distribution of expression values in each sample. QC tables quantify basic sample characteristics and missing values. Correlation matrices show how similar whole-sample expression profiles are. Variance tables identify features that change the most across samples. PCA compresses thousands of expression measurements into a few axes so that major sample relationships can be visualized. Heatmaps show the expression patterns of the most variable features.",
    "",
    "Important distinction",
    "----------------------",
    "The three datasets are not identical types of measurements. GSE148158 is RNA-seq-derived normalized expression data, whereas GSE28688 and GSE52052 are microarray expression datasets. Their expression values therefore should not be interpreted as directly comparable numerical measurements across datasets. Each dataset is primarily interpreted within its own experiment.",
    "",
]

for name in DATASETS:
    section = summarize_dataset(name)
    all_lines.extend(section)
    all_lines.append("")

all_lines.extend([
    "Overall conclusion",
    "==================",
    "Together, the three explorations provide complementary views of early human iPSC reprogramming. GSE148158 emphasizes early OSKM-associated transcriptional responses, GSE28688 provides a time-course from fibroblast cells through early post-transduction stages and established pluripotent controls, and GSE52052 compares different reprogramming cocktails at day 11.",
    "",
    "The main output of the current project is therefore a quality-controlled exploratory picture of how samples relate to one another globally: which samples are similar, which are separated, which features are most variable, and whether the biological grouping is visible in PCA/correlation structure. The results should be used as a foundation for later statistical analyses rather than as proof of differential expression or causality.",
])

text = "\n".join(all_lines) + "\n"
(OUT / "Summary.txt").write_text(text, encoding="utf-8")

print(text)
print(f"Summary saved to: {OUT / 'Summary.txt'}")
