"""Canonical local paths for Z4 audits.

All Z4 command-line tools should use these defaults unless an explicit path is
provided. Keeping paths here avoids long, error-prone PowerShell commands.
"""
from pathlib import Path

DATA = Path("Data")
RESULTS = Path("results") / "Dynamics"

GSE297234_RDS = [
    DATA / "GSE297234_GM00731_SEVOSKM.rds",
    DATA / "GSE297234_HFIB_COMBINED_SEVOSKM.rds",
]

GSE297234_METADATA = RESULTS / "z4_external_metadata_audit" / "GSE297234" / "metadata.csv"
GSE28688_METADATA = RESULTS / "z4_external_metadata_audit" / "GSE28688" / "metadata.csv"
GEO_METADATA_DIR = DATA / "GEO_Z4"
