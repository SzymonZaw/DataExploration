"""Experimental seed for mathematical dynamics of OSKM reprogramming."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "Dynamics"
OUT.mkdir(parents=True, exist_ok=True)

# Existing EDA outputs. Different PCA spaces are NOT assumed comparable yet.
DATASETS = {
    "GSE28688": RESULTS / "GSE28688" / "non_normalized" / "07_PCA_coordinates.csv",
    "GSE148158": RESULTS / "GSE148158" / "07_PCA_coordinates.csv",
    "GSE52052": RESULTS / "GSE52052" / "07_PCA_coordinates.csv",
    "GSE67462": RESULTS / "GSE67462" / "07_PCA_coordinates.csv",
    "GSE297234": RESULTS / "GSE297234" / "07_PCA_coordinates.csv",
}


def load_pca(path):
    if not path.exists():
        return None
    df = pd.read_csv(path, index_col=0)
    cols = [c for c in df.columns if str(c).startswith("PC")]
    return df[cols].apply(pd.to_numeric, errors="coerce") if cols else None


def orient_state(pc1):
    """Deterministic PC1 sign convention; sign has no biological meaning."""
    x = pd.Series(pc1, dtype=float).copy()
    if x.notna().sum() == 0:
        return x
    i = x.abs().idxmax()
    return -x if x.loc[i] < 0 else x


def estimate_velocity(values, times=None):
    """Finite-difference estimate of local state velocity."""
    x = np.asarray(values, dtype=float)
    t = np.arange(len(x), dtype=float) if times is None else np.asarray(times, dtype=float)
    return np.full(len(x), np.nan) if len(x) < 2 else np.gradient(x, t)


def stability_indicators(values, window=5):
    """Exploratory rolling variance and lag-1 autocorrelation."""
    s = pd.Series(values, dtype=float)
    min_periods = max(3, window // 2)
    variance = s.rolling(window, min_periods=min_periods).var()
    autocorr = s.rolling(window, min_periods=min_periods).apply(
        lambda x: pd.Series(x).autocorr(lag=1) if pd.Series(x).std() else np.nan,
        raw=False,
    )
    return variance, autocorr


def candidate_library(x):
    """Small interpretable basis for a future symbolic-regression/GP stage."""
    x = np.asarray(x, dtype=float)
    return pd.DataFrame({
        "1": np.ones_like(x),
        "x": x,
        "x2": x**2,
        "x3": x**3,
        "sin_x": np.sin(x),
        "cos_x": np.cos(x),
        "log_abs_x": np.log1p(np.abs(x)),
    })


def process(name, path):
    coords = load_pca(path)
    if coords is None:
        return None
    state = orient_state(coords["PC1"])
    velocity = estimate_velocity(state)
    variance, autocorr = stability_indicators(state)
    out = pd.DataFrame({
        "sample": state.index.astype(str),
        "dataset": name,
        "state_PC1": state.to_numpy(),
        "velocity": velocity,
        "rolling_variance": variance.to_numpy(),
        "rolling_autocorrelation": autocorr.to_numpy(),
    })
    return out


def main():
    availability = []
    results = []
    for name, path in DATASETS.items():
        result = process(name, path)
        availability.append({"dataset": name, "PCA_file_found": result is not None, "path": str(path)})
        if result is not None:
            results.append(result)

    pd.DataFrame(availability).to_csv(OUT / "01_dataset_availability.csv", index=False)
    if not results:
        print("No existing PCA trajectories found. Run the EDA scripts first.")
        return

    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(OUT / "02_state_dynamics_exploratory.csv", index=False)

    first = results[0]
    library = candidate_library(first["state_PC1"])
    library.insert(0, "sample", first["sample"].to_numpy())
    library.to_csv(OUT / "03_symbolic_candidate_library.csv", index=False)

    report = [
        "Experimental mathematical-dynamics prototype",
        "",
        "Current layer:",
        "  1. Reuse existing PCA coordinates from GEO EDA.",
        "  2. Use PC1 only as an exploratory state coordinate.",
        "  3. Estimate local velocity.",
        "  4. Estimate rolling variance and lag-1 autocorrelation.",
        "  5. Prepare an interpretable candidate library for symbolic regression.",
        "",
        "Target hypothesis:",
        "  Independent reprogramming experiments may share a generalisable latent",
        "  dynamical structure after appropriate cross-study alignment.",
        "",
        "Next scientific steps:",
        "  * build a shared latent space rather than comparing raw PCA axes,",
        "  * encode biological time explicitly,",
        "  * use replicate-aware modelling,",
        "  * discover equations with symbolic regression / Genetic Programming,",
        "  * hold out an entire dataset for external validation,",
        "  * quantify uncertainty and robustness to preprocessing choices.",
        "",
        "This prototype makes no claim of critical transition, causality, universality,",
        "or quantum/relativistic mechanism.",
    ]
    (OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")
    print(f"Dynamics prototype results written to: {OUT}")


if __name__ == "__main__":
    main()
