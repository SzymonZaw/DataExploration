from pathlib import Path
import gzip
import csv

root = Path("Data/GSE67520_extracted")
out = Path("results/Dynamics/z6_gse67520_peak_file_audit.csv")
out.parent.mkdir(parents=True, exist_ok=True)

rows = []

for path in sorted(root.glob("*.gz")):
    with gzip.open(path, "rt", errors="replace") as f:
        first_lines = []
        n = 0
        chroms = set()
        malformed = 0

        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            n += 1

            if len(first_lines) < 2:
                first_lines.append(line.rstrip("\n")[:250])

            if len(fields) >= 3:
                chroms.add(fields[0])
            else:
                malformed += 1

    name = path.name
    if "_Flag_" in name:
        marker = "FLAG_tagged_exogenous_OCT4"
    elif "_Oct4_" in name:
        marker = "OCT4_total"
    elif "_K4me1_" in name:
        marker = "H3K4me1"
    elif "_K27ac_" in name:
        marker = "H3K27ac"
    elif "_K4me3_" in name:
        marker = "H3K4me3"
    elif "_K27me3_" in name:
        marker = "H3K27me3"
    elif "_RNAPII_" in name:
        marker = "RNAPII"
    else:
        marker = "UNKNOWN"

    rows.append({
        "file": name,
        "marker": marker,
        "peak_rows": n,
        "chromosome_count": len(chroms),
        "rows_with_fewer_than_3_columns": malformed,
        "example_row_1": first_lines[0] if first_lines else "",
        "example_row_2": first_lines[1] if len(first_lines) > 1 else "",
    })

with out.open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print("Liczba plikow:", len(rows))
print("Zapisano:", out)
print("\nLiczba wierszy pikow wedlug markera:")
totals = {}
for r in rows:
    totals.setdefault(r["marker"], []).append(r["peak_rows"])

for marker, values in sorted(totals.items()):
    print(f"{marker:30s} pliki={len(values):2d}, min={min(values):7d}, "
          f"mediana={sorted(values)[len(values)//2]:7d}, max={max(values):7d}")
