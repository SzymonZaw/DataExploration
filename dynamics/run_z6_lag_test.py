from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("results/Dynamics/z6_gse67462_mechanistic_interpretation")
EXPR_FILE = BASE / "03_gene_expression_by_stage.csv"
SIGNAL_FILE = BASE / "09_counts_clean.csv"

# Wyłącznie uporządkowane etapy reprogramowania.
stages = ["day0", "day1", "day3", "day5",
          "day7", "day11", "day15", "day18"]

modalities = ["H3K27ac", "H3K4me3", "RNAPII",
              "OCT4", "H3K27me3"]

expr = pd.read_csv(EXPR_FILE, index_col=0)
signal = pd.read_csv(SIGNAL_FILE, index_col=0)

# Upewnij się, że wszystkie potrzebne etapy są dostępne.
missing_expr = [s for s in stages if s not in expr.columns]
if missing_expr:
    raise SystemExit(f"Brak etapów ekspresji: {missing_expr}")

needed = [f"{m}__{s}" for m in modalities for s in stages]
missing_signal = [c for c in needed if c not in signal.columns]
if missing_signal:
    raise SystemExit(f"Brak kolumn sygnału: {missing_signal}")

common = expr.index.intersection(signal.index)
if len(common) == 0:
    raise SystemExit("STOP: brak wspólnych symboli genów.")

expr = expr.loc[common, stages].apply(pd.to_numeric, errors="coerce")
signal = signal.loc[common, needed].apply(
    pd.to_numeric, errors="coerce"
).fillna(0)

# Kolejne etapy: sygnał w t, ekspresja w t oraz w t+1.
source_stages = stages[:-1]
target_stages = stages[1:]

rows = []
for modality in modalities:
    cols = [f"{modality}__{s}" for s in source_stages]
    sig = np.log1p(signal[cols].to_numpy(dtype=float))

    for i, gene in enumerate(common):
        x = pd.Series(sig[i, :])
        concurrent_y = pd.Series(
            expr.loc[gene, source_stages].to_numpy(dtype=float)
        )
        lagged_y = pd.Series(
            expr.loc[gene, target_stages].to_numpy(dtype=float)
        )

        valid_c = x.notna() & concurrent_y.notna()
        valid_l = x.notna() & lagged_y.notna()

        # Korelacja nie jest informatywna przy stałym sygnale.
        unique_signal = int(x.nunique())
        if unique_signal < 2:
            concurrent_rho = np.nan
            lagged_rho = np.nan
        else:
            concurrent_rho = x[valid_c].corr(
                concurrent_y[valid_c], method="spearman"
            )
            lagged_rho = x[valid_l].corr(
                lagged_y[valid_l], method="spearman"
            )

        gain = (
            lagged_rho - concurrent_rho
            if pd.notna(lagged_rho) and pd.notna(concurrent_rho)
            else np.nan
        )

        rows.append({
            "gene": gene,
            "modality": modality,
            "concurrent_rho": concurrent_rho,
            "lagged_rho": lagged_rho,
            "lead_gain": gain,
            "n_timepoints": int(valid_l.sum()),
            "signal_unique_values": unique_signal,
            "candidate_gain_ge_010_rho_ge_030": bool(
                pd.notna(gain) and gain >= 0.10
                and lagged_rho >= 0.30
                and valid_l.sum() >= 4
            ),
        })

scores = pd.DataFrame(rows)
scores.to_csv(BASE / "10_z6_gene_lag_scores.csv", index=False)

summary = (
    scores.groupby("modality")
    .agg(
        genes_tested=("gene", "size"),
        genes_with_valid_gain=("lead_gain", "count"),
        median_concurrent_rho=("concurrent_rho", "median"),
        median_lagged_rho=("lagged_rho", "median"),
        median_lead_gain=("lead_gain", "median"),
        candidate_genes=(
            "candidate_gain_ge_010_rho_ge_030", "sum"
        ),
    )
    .reset_index()
)
summary.to_csv(BASE / "11_z6_modality_summary.csv", index=False)

print("\n=== Z6: PODSUMOWANIE REGULATORÓW ===")
print(summary.to_string(index=False))
print("\nWspólne geny:", len(common))
print("Etapy w teście:", " -> ".join(stages))
print("Liczba porównań czasowych:", len(source_stages))
print("Zapisano:", BASE / "10_z6_gene_lag_scores.csv")
print("Zapisano:", BASE / "11_z6_modality_summary.csv")

print("\n=== 20 NAJWIĘKSZYCH WZROSTÓW KORELACJI ===")
valid = scores.dropna(subset=["lead_gain", "lagged_rho"])
top = valid.sort_values("lead_gain", ascending=False).head(20)
print(top[[
    "gene", "modality", "concurrent_rho",
    "lagged_rho", "lead_gain"
]].to_string(index=False))
