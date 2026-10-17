from pathlib import Path
import pandas as pd
import numpy as np

base = Path("results/Dynamics/z6_gse67462_mechanistic_interpretation")

# Kolejność kolumn odpowiada kolejności zapisu macierzy.
layout = [
    ("FLAG_tagged_exogenous_OCT4", 8),
    ("H3K27ac", 9),
    ("H3K27me3", 9),
    ("H3K4me1", 9),
    ("H3K4me3", 9),
    ("OCT4", 9),
    ("RNAPII", 9),
]

def load_matrix(filename):
    path = base / filename
    raw = pd.read_csv(path, index_col=0, low_memory=False)

    if "stage" not in raw.index:
        raise ValueError(
            f"Nie znaleziono wiersza 'stage' w {filename}. "
            f"Pierwsze indeksy: {list(raw.index[:5])}"
        )

    stages = raw.loc["stage"].astype(str).tolist()
    modalities = []
    for modality, n in layout:
        modalities.extend([modality] * n)

    if len(modalities) != len(raw.columns):
        raise ValueError(
            f"Niezgodna liczba kolumn: {len(raw.columns)}; "
            f"oczekiwano {len(modalities)}"
        )

    # Usuń techniczne wiersze zapisane z nagłówkiem macierzy.
    data = raw.loc[~raw.index.isin(["modality", "stage", "gene"])].copy()
    data.index = data.index.astype(str)
    data = data.apply(pd.to_numeric, errors="coerce").fillna(0)

    data.columns = pd.MultiIndex.from_arrays(
        [modalities, stages], names=["modality", "stage"]
    )
    return data

presence = load_matrix("05_promoter_peak_presence_matrix.csv")
counts = load_matrix("06_promoter_peak_count_matrix.csv")

print("\n=== KONTROLA MACIERZY ===")
print("Obecność:", presence.shape)
print("Zliczenia:", counts.shape)
print("Te same geny:", presence.index.equals(counts.index))
print("Te same kolumny:", presence.columns.equals(counts.columns))

print("\n=== POKRYCIE GENÓW WEDŁUG REGULATORA I ETAPU ===")
for modality in presence.columns.get_level_values("modality").unique():
    for stage in presence[modality].columns:
        x = presence[(modality, stage)]
        y = counts[(modality, stage)]
        print(
            f"{modality:14s} {stage:5s} "
            f"geny_z_sygnalem={(x > 0).sum():5d} "
            f"suma_zliczen={int(y.sum()):7d}"
        )

# Sprawdź, czy macierz obecności jest zgodna z macierzą zliczeń.
expected_presence = counts.gt(0).astype(int)
observed_presence = presence.gt(0).astype(int)
print("\nZgodność obecności ze zliczeniami:",
      bool(expected_presence.equals(observed_presence)))

# Zapisz macierze z jednoznacznymi nazwami kolumn.
def flatten_columns(df):
    out = df.copy()
    out.columns = [f"{m}__{s}" for m, s in out.columns]
    return out

flatten_columns(presence).to_csv(base / "08_presence_clean.csv")
flatten_columns(counts).to_csv(base / "09_counts_clean.csv")
print("\nZapisano 08_presence_clean.csv i 09_counts_clean.csv")

