from pathlib import Path
import gzip
import statistics

root = Path("Data/GSE67520_extracted")

files = [
    "GSM1648624_Oct4_d0_peaks.broadPeak.gz",
    "GSM1648625_Oct4_d1_peaks.broadPeak.gz",
    "GSM1648616_Flag_d0_peaks.broadPeak.gz",
    "GSM1648617_Flag_d1_peaks.broadPeak.gz",
    "GSM1648623_Flag_d18_peaks.broadPeak.gz",
]

for filename in files:
    path = root / filename
    widths, scores, signals, pvalues, qvalues = [], [], [], [], []
    bad = 0
    examples = []

    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                start, end = int(fields[1]), int(fields[2])
                if end < start:
                    raise ValueError("end < start")
                widths.append(end - start)

                if len(fields) >= 5:
                    scores.append(float(fields[4]))
                if len(fields) >= 7:
                    signals.append(float(fields[6]))
                if len(fields) >= 8:
                    pvalues.append(float(fields[7]))
                if len(fields) >= 9:
                    qvalues.append(float(fields[8]))

                if len(examples) < 2:
                    examples.append(line.rstrip()[:220])
            except (ValueError, IndexError):
                bad += 1

    def summary(values):
        if not values:
            return "brak"
        vals = sorted(values)
        return (
            f"mediana={statistics.median(vals):.2f}, "
            f"min={vals[0]:.2f}, max={vals[-1]:.2f}"
        )

    print(f"\n=== {filename} ===")
    print(f"Poprawne szerokosci: {len(widths)}; bledne wiersze: {bad}")
    print("Szerokosc piku (bp):", summary(widths))
    print("Score:", summary(scores))
    print("Signal value:", summary(signals))
    print("p-value pole:", summary(pvalues))
    print("q-value pole:", summary(qvalues))
    for row in examples:
        print("Przyklad:", row)
