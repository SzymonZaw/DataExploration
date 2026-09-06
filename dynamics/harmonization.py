"""Gene-level biological harmonization (Stage 2.6).

This module is the intended home for Ensembl/RefSeq normalization, BioMart and
MyGene mapping, platform annotation handling, and construction of the common
human gene space. Existing implementations remain in Dynamics.py for this
non-breaking first refactor step.
"""

from Dynamics import (
    biomart_ensembl_to_human,
    biomart_refseq_to_human,
    clean_refseq,
    direct_gene_matrix,
    extract_ensembl,
    extract_refseq,
    map_rows,
    merge_direct_and_mapping,
    mygene_human_ensembl,
    mygene_mouse_refseq,
    normalize_gene_symbol,
    read_feature_matrix,
    stage2_6,
)

__all__ = [
    "biomart_ensembl_to_human",
    "biomart_refseq_to_human",
    "clean_refseq",
    "direct_gene_matrix",
    "extract_ensembl",
    "extract_refseq",
    "map_rows",
    "merge_direct_and_mapping",
    "mygene_human_ensembl",
    "mygene_mouse_refseq",
    "normalize_gene_symbol",
    "read_feature_matrix",
    "stage2_6",
]
