from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "SUMMARY"
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "GSE148158": ("A PIANO (Proper, Insufficient, Aberrant, and NO reprogramming) response to the Yamanaka factors in the initial stages of human iPSC reprogramming", "Human; RNA-seq / high-throughput sequencing", "GPL16791 - Illumina HiSeq 2500", 13, "Very early human OSKM transcriptional response; BJ fibroblasts, hESC and GFP controls."),
    "GSE28688": ("Molecular insights into induced pluripotency mediated by the OCT4, SOX2, KLF and c-MYC gene regulatory network", "Human; expression profiling by array", "GPL6883 - Illumina HumanRef-8 v3.0", 14, "Early human time course: fibroblasts, 24/48/72 h after OSKM transduction, hESC and iPS controls."),
    "GSE52052": ("Comparison of gene expression at day 11 of iPSC reprogramming with different reprogramming cocktails", "Human; expression profiling by array", "GPL14550 - Agilent-028004 SurePrint G3 Human GE 8x60K", 6, "Day-11 comparison of GFP control, OSK-based cocktails, OSKM and H1 hESC."),
    "GSE67462": ("Expression data from OSKM-mediated 2nd reprogramming cells and the corresponding iPS cell line", "Mouse; expression profiling by array; time series", "GPL19972 - Affymetrix Mouse Gene 1.0 ST Array", 18, "Second reprogramming sampled at day 0, 1, 3, 5, 7, 11, 15, 18 and iPSC, with two replicates per stage."),
    "GSE67520": ("Oct4 binding and Histone modification profiling during OSKM-mediated 2nd reprogramming", "Mouse; ChIP-seq / genome binding and occupancy profiling", "GPL13112 - Illumina HiSeq 2000; GPL17021 - Illumina HiSeq 2500", 71, "Regulatory dataset profiling Oct4, histone modifications, RNAPII and input controls during second reprogramming."),
}

GSE52052_GROUPS = {
    "HDF_GFP(+)_day11": "GFP_control",
    "HDF_OSK+control_inh_TRA(+)_day11": "OSK_control_inhibitor",
    "HDF_OSK+let-7_inh_TRA(+)_day11": "OSK_let7_inhibitor",
    "HDF_OSKM_TRA(+)_day11": "OSKM",
    "HDF_OSK+LIN-41_TRA(+)_day11": "OSK_LIN41",
    "H1_hESC": "H1_hESC",
}

GSE52052_GSM_GROUPS = {
    "GSM1258008": "GFP_control",
    "GSM1258009": "OSK_control_inhibitor",
    "GSM1258010": "OSK_let7_inhibitor",
    "GSM1258011": "OSKM",
    "GSM1258012": "OSK_LIN41",
    "GSM1258013": "H1_hESC",
}


def group_for_gse52052(sample):
    s = str(sample).strip()
    for gsm, group in GSE52052_GSM_GROUPS.items():
        if gsm in s:
            return group
    for key, group in GSE52052_GROUPS.items():
        if s == key or s.startswith(key) or key in s:
            return group
    return "unclassified"


def read_csv(path, index_col=0):
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, index_col=index_col)
    except Exception:
        return None


def fmt(x):
    try:
        return f"{float(x):.3f}"
    except (TypeError, ValueError):
        return "n/a"


def expression_section(name, lines):
    folder = RESULTS / name
    meta = read_csv(folder / "01_sample_metadata.csv", index_col=None)
    qc = read_csv(folder / "03_sample_QC.csv")
    corr = read_csv(folder / "04_sample_correlation.csv")
    var = read_csv(folder / "06_feature_variance.csv")
    if var is None and name == "GSE52052":
        var = read_csv(folder / "07_probe_variance.csv")
    pca = read_csv(folder / "07_PCA_coordinates.csv")
    if pca is None and name == "GSE52052":
        pca = read_csv(folder / "08_PCA_coordinates.csv")
    if name == "GSE52052" and meta is not None and "sample" in meta.columns:
        meta = meta.copy()
        meta["group"] = meta["sample"].map(group_for_gse52052)

    lines += ["Exploration performed", "----------------------",
              "- Expression distributions / boxplot.",
              "- Per-sample QC: mean, median, SD and missing values.",
              "- Pearson sample correlation.",
              "- Feature variance ranking.",
              "- PCA using the most variable features.",
              "- Heatmap of top variable features."]
    if name == "GSE148158":
        lines.append("- Supplied matrix is GEO-normalized; no raw-count normalization was added.")
    elif name == "GSE28688":
        lines += ["- Non-normalized matrix was transformed with log2(x+1).",
                  "- Exploratory quantile normalization was used for visualization.",
                  "- RAW archive was identified as GPL6883 BGX probe annotation rather than expression data."]
    elif name == "GSE52052":
        lines += ["- Agilent RAW files were parsed using gProcessedSignal and non-control probes.",
                  "- Exploratory quantile normalization was generated; this does not replace GEO/GeneSpring processing."]
    elif name == "GSE67462":
        lines += ["- Processed expression values were read from the GEO Series Matrix.",
                  "- RAW archive contains 18 Affymetrix CEL files and was inventoried separately."]
    lines.append("")

    if meta is not None and "group" in meta.columns:
        counts = meta["group"].fillna("unclassified").astype(str).value_counts().sort_index()
        lines += ["Sample groups observed", "----------------------"]
        lines += [f"- {g}: {n} sample(s)" for g, n in counts.items()]
        if "unclassified" in counts.index:
            lines.append("- Unclassified samples should be checked before interpreting group-level plots.")
        lines.append("")

    if qc is not None and not qc.empty:
        lines += ["QC summary", "----------"]
        if "missing" in qc.columns:
            lines.append(f"- Maximum missing values in one sample: {int(qc['missing'].max())}")
        if "mean" in qc.columns:
            lines.append(f"- Sample mean range: {fmt(qc['mean'].min())} to {fmt(qc['mean'].max())}")
        if "sd" in qc.columns:
            lines.append(f"- Sample SD range: {fmt(qc['sd'].min())} to {fmt(qc['sd'].max())}")
        lines.append("")

    if corr is not None and not corr.empty:
        a = corr.apply(pd.to_numeric, errors="coerce").to_numpy(float)
        if a.ndim == 2 and a.shape[0] == a.shape[1]:
            x = a[~np.eye(a.shape[0], dtype=bool)]
            x = x[np.isfinite(x)]
            if len(x):
                lines += ["Correlation summary", "-------------------",
                          f"- Pearson r between different samples: min {fmt(x.min())}, median {fmt(np.median(x))}, max {fmt(x.max())}",
                          "- Higher correlation indicates greater similarity of whole-sample expression profiles.", ""]

    if var is not None and not var.empty:
        col = "variance" if "variance" in var.columns else var.columns[0]
        x = pd.to_numeric(var[col], errors="coerce").dropna()
        if len(x):
            lines += ["Feature variability", "-------------------",
                      f"- Features ranked by variance: {len(x):,}",
                      f"- Highest observed variance: {fmt(x.iloc[0])}", ""]

    if pca is not None and not pca.empty and "PC1" in pca.columns:
        lines += ["PCA interpretation", "------------------",
                  "- PC1 and PC2 summarize major axes of variation among samples."]
        if "group" in pca.columns:
            gc = pca["group"].fillna("unclassified").astype(str).value_counts()
            singletons = int((gc == 1).sum())
            if singletons:
                lines.append(f"- {singletons} group(s) contain only one sample; separation is descriptive rather than evidence of reproducibility.")
        lines += ["- PCA separation can reflect biology, technical effects or both.", ""]

    lines += ["Important interpretation limits", "-------------------------------",
              "- These scripts perform exploratory data analysis, not differential-expression testing.",
              "- PCA, correlation and clustering do not by themselves prove mechanism or statistical significance."]
    if name == "GSE148158":
        lines.append("- OSKM 48 h and 72 h are singleton conditions, so replicate-based inference at those exact time points is not possible.")
    if name == "GSE28688":
        lines.append("- Replicate pairs support inspection of within-group consistency.")
    if name == "GSE52052":
        lines.append("- Most conditions have one sample, so condition-level conclusions are descriptive.")
    if name == "GSE67462":
        lines.append("- Two replicates are available at each reprogramming stage, making this dataset suitable for later replicate-aware time-course analysis.")
    lines.append("")


def chipseq_section(lines):
    folder = RESULTS / "GSE67520"
    inv = read_csv(folder / "01_RAW_file_inventory.csv", index_col=None)
    peaks = read_csv(folder / "06_peak_summary_sorted.csv", index_col=None)
    if peaks is None:
        peaks = read_csv(folder / "02_peak_summary.csv", index_col=None)

    lines += ["Exploration performed", "----------------------",
              "- RAW archive inventory with experimental mark and reprogramming stage.",
              "- BED/BROADPEAK genomic interval parsing.",
              "- Peak counts, lengths, chromosome counts and covered base pairs.",
              "- Peak counts by stage and chromatin mark, including stage dynamics.",
              "- No differential peak calling was performed.", ""]
    if inv is not None and not inv.empty:
        lines += ["RAW archive summary", "-------------------", f"- Files inventoried: {len(inv)}"]
        if "mark" in inv.columns:
            c = inv["mark"].fillna("unknown").value_counts().sort_index()
            lines.append("- Categories: " + ", ".join(f"{k} ({v})" for k, v in c.items()))
        lines.append("")
    if peaks is not None and not peaks.empty:
        lines += ["Peak summary", "------------", f"- Parsed peak files: {len(peaks)}"]
        if "peaks" in peaks.columns:
            total = pd.to_numeric(peaks["peaks"], errors="coerce").fillna(0).sum()
            lines.append(f"- Total parsed peaks across files: {int(total):,}")
        if "mark" in peaks.columns:
            lines.append("- Marks observed: " + ", ".join(sorted(peaks["mark"].dropna().astype(str).unique())))
        if "stage" in peaks.columns:
            order = [f"day{x}" for x in [0, 1, 3, 5, 7, 11, 15, 18]] + ["iPSC"]
            observed = set(peaks["stage"].astype(str))
            lines.append("- Stages observed: " + ", ".join(s for s in order if s in observed))
        lines.append("")
    lines += ["Important interpretation limits", "-------------------------------",
              "- GSE67520 is not a gene-expression dataset; ChIP-seq peaks represent genomic regions enriched for the assayed factor or mark.",
              "- More peaks does not directly mean higher gene expression.",
              "- Interpretation should consider peak location, signal, controls and genomic annotation.",
              "- Differential binding and causal regulatory analysis are outside the current EDA.", ""]


lines = [
    "SUMMARY OF THE FIVE GEO DATASET EXPLORATIONS",
    "=============================================", "",
    "Purpose", "-------",
    "This report summarizes the five GEO datasets currently included in the project. Four datasets provide gene-expression measurements and one provides ChIP-seq regulatory information.", "",
    "How to read the results", "------------------------",
    "Expression boxplots show distributions; QC tables summarize sample characteristics; correlation matrices show whole-sample similarity; variance tables identify variable features; PCA summarizes major axes of variation; heatmaps show patterns among the most variable features.", "",
    "GSE67520 is different: it contains ChIP-seq peak calls rather than expression values, so its exploration focuses on genomic enrichment and changes across stages and chromatin marks.", "",
    "Important distinction between datasets", "--------------------------------------------",
    "Expression values should be interpreted primarily within each experiment. GSE148158 is RNA-seq-derived; GSE28688 and GSE52052 are human microarrays; GSE67462 is a mouse Affymetrix microarray. Different technologies, species, preprocessing and value scales make raw numeric expression values from different datasets non-equivalent.", "",
    "GSE67520 is complementary rather than directly comparable: it provides regulatory information about Oct4 binding, histone modifications and RNAPII. Integration should therefore use genomic annotation and biological hypotheses rather than compare peak counts with expression values.", "",
]

for name, (title, experiment, platform, samples, description) in DATASETS.items():
    lines += [name, "=" * len(name), f"Title: {title}", f"Experiment: {experiment}", f"Platform: {platform}", f"Expected samples from GEO design: {samples}", f"What the data represent: {description}", ""]
    folder = RESULTS / name
    if not folder.exists():
        lines += ["RESULTS NOT FOUND: run the corresponding exploration script first.", ""]
    elif name == "GSE67520":
        chipseq_section(lines)
    else:
        expression_section(name, lines)

lines += ["Overall conclusion", "==================",
          "The five datasets now provide complementary views of cellular reprogramming at transcriptional and regulatory levels: very early human responses (GSE148158), an early human time course (GSE28688), day-11 cocktail comparison (GSE52052), a longer mouse second-reprogramming trajectory (GSE67462), and ChIP-seq regulatory profiling (GSE67520).", "",
          "The current project provides an exploratory foundation for quality control, hypothesis generation and identification of candidate biological transitions. These outputs are not by themselves evidence of differential expression, differential binding, mechanism or causality.", "",
          "The next stage should use replicate-aware statistical testing and multiple-testing correction. GSE67462 is particularly suitable for time-course modeling. GSE67520 can later be extended with differential binding and genomic annotation so regulatory findings can be connected to expression changes. Cross-dataset integration must account for species, platform, preprocessing and experimental design.", "",
          "Raw archives that are too large for GitHub remain local; the analysis scripts and lightweight result summaries are the reproducible project components."]

text = "\n".join(lines) + "\n"
(OUT / "Summary.txt").write_text(text, encoding="utf-8")
print(text)
print(f"Summary saved to: {OUT / 'Summary.txt'}")
