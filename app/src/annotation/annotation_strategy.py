from typing import Optional

from Bio.SeqRecord import SeqRecord
from app.src.io.genbank_parser import parse_cds_features, parse_mat_peptides


"""
Module: annotation_strategy.py

Purpose:
    Detect which feature type a GenBank record uses for gene-level annotation.

Returns:
    "mat_peptide"  — record has mat_peptide features (e.g. FMDV)
    "CDS"          — record has CDS features (e.g. PRRSV)
    None           — record has no annotation → needs tblastn lifting

Notes:
    - mat_peptide takes priority over CDS because some viruses (e.g. FMDV) have a
      single whole-genome polyprotein CDS that carries no gene-level names; the real
      names live in mat_peptide sub-features.
    - This module only detects the feature type. It does NOT perform extraction.
"""


def get_feature_type(record: SeqRecord) -> Optional[str]:
    """
    Detect the annotation feature type used by a GenBank record.

    Priority:
        1. mat_peptide present → "mat_peptide"
        2. CDS present         → "CDS"
        3. Neither             → None  (record needs tblastn lifting)

    Args:
        record: A Biopython SeqRecord

    Returns:
        "mat_peptide", "CDS", or None
    """
    if parse_mat_peptides(record):
        return "mat_peptide"

    if parse_cds_features(record):
        return "CDS"

    return None
