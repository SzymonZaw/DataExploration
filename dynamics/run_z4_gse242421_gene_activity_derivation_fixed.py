"""Compatibility wrapper for the GSE242421 gene-activity derivation.

The published scATAC analysis product stores Matrix Market as a gzip-compressed
member inside scATAC.zip. The original derivation script passed the compressed
stream directly to scipy.io.mmread, which correctly rejected it because the
Matrix Market banner is inside the gzip layer. This wrapper patches that loader
without changing the derivation logic or frozen modules.
"""

from __future__ import annotations

import gzip

from scipy.io import mmread
from scipy import sparse

from . import run_z4_gse242421_gene_activity_derivation as derivation


def load_mtx_from_zip(zf, member):
    with zf.open(member) as raw_fh:
        if member.lower().endswith(".gz"):
            with gzip.GzipFile(fileobj=raw_fh, mode="rb") as fh:
                matrix = mmread(fh)
        else:
            matrix = mmread(raw_fh)
    return sparse.csc_matrix(matrix, dtype="float64")


def main():
    derivation.load_mtx_from_zip = load_mtx_from_zip
    derivation.main()


if __name__ == "__main__":
    main()
