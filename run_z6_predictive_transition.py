"""Run the official Z6 predictive cellular-state transition benchmark.

The scientific benchmark implementation lives in dynamics.model_benchmark.
This wrapper only redirects its outputs to the dedicated Z6 result directory.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from dynamics import run_model_benchmark as benchmark

ROOT = Path(__file__).resolve().parent
Z6_OUT = ROOT / "results" / "Dynamics" / "z6_predictive_transition"


def main() -> None:
    p = argparse.ArgumentParser(description="Run the frozen Z6 predictive transition benchmark")
    p.add_argument("--max-genes", type=int, default=2000)
    p.add_argument("--state-dim", type=int, default=8)
    p.add_argument("--hidden-dim", type=int, default=128)
    p.add_argument("--epochs", type=int, default=250)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--prefix-fraction", type=float, default=.6)
    p.add_argument("--permutation-n", type=int, default=1000)
    p.add_argument("--seeds", type=int, nargs="+", default=list(benchmark.SEEDS))
    a = p.parse_args()

    Z6_OUT.mkdir(parents=True, exist_ok=True)
    benchmark.OUT = Z6_OUT
    benchmark.run(
        max_genes=a.max_genes,
        state_dim=a.state_dim,
        hidden_dim=a.hidden_dim,
        epochs=a.epochs,
        lr=a.lr,
        prefix_fraction=a.prefix_fraction,
        seeds=tuple(a.seeds),
        permutation_n=a.permutation_n,
    )


if __name__ == "__main__":
    main()
