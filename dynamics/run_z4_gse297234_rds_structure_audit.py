"""Python launcher for the GSE297234 Z4 RDS structure audit.

Normal usage requires no arguments. The R implementation owns the audit logic;
this module only provides a short, consistent Python entry point.
"""
from __future__ import annotations

import subprocess
from pathlib import Path


SCRIPT = Path(__file__).with_name("run_z4_gse297234_rds_structure_audit.R")


def main() -> None:
    result = subprocess.run(
        ["Rscript", str(SCRIPT)],
        check=False,
    )
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
