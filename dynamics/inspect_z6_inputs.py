from pathlib import Path
import pandas as pd

out = Path("results/Dynamics/z6_gse67462_mechanistic_interpretation")

files = [
    "03_gene_expression_by_stage.csv",
    "05_promoter_peak_presence_matrix.csv",
    "06_promoter_peak_count_matrix.csv",
    "07_modality_coverage_summary.csv",
]

for name in files:
    path = out / name
    print(f"\n{'=' * 70}\n{name}")
    if not path.exists():
        print("BRAK PLIKU")
        continue

    df = pd.read_csv(path, index_col=0 if name !=
                     "07_modality_coverage_summary.csv" else None)
    print("Wymiary:", df.shape)
    print("Kolumny:", list(df.columns[:20]))
    print(df.head(5).to_string())

    if name == "05_promoter_peak_presence_matrix.csv":
        print("\nSuma obecności sygnału w kolumnach:")
        print(df.sum(numeric_only=True).to_string())

    if name == "06_promoter_peak_count_matrix.csv":
        print("\nSumy zliczeń w kolumnach:")
        print(df.sum(numeric_only=True).to_string())

print("\nKontrola zakończona.")
