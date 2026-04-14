from typing import Tuple

from Bio.SeqRecord import SeqRecord
from app.src.io.genbank_parser import parse_cds_features, parse_mat_peptides


"""
Module: annotation_strategy.py

Purpose:
    Decide which extraction strategy and feature type to use for a query genome:
    - ("direct", "mat_peptide"): query has mat_peptide features (e.g. FMDV polyprotein)
    - ("direct", "CDS"):         query has CDS features (e.g. PRRSV)
    - ("minimap", None):         query has no annotation -> coordinate transfer from reference

Notes:
    - mat_peptide takes priority over CDS because some viruses (e.g. FMDV) have a single
      polyprotein CDS that is not useful for naming; gene names live in mat_peptide instead.
    - This module only decides the strategy, it does NOT perform extraction.
    - Alias-based naming and feature normalization should be handled separately.
"""


def choose_strategy(query_record: SeqRecord) -> Tuple[str, str]:
    """
    Select extraction strategy and feature type for a query genome.

    Priority:
        1. mat_peptide present -> ("direct", "mat_peptide")
        2. CDS present         -> ("direct", "CDS")
        3. Neither             -> ("minimap", None)

    Args:
        query_record: Query genome (GenBank record)

    Returns:
        Tuple of (strategy, feature_type):
            strategy     - "direct" or "minimap"
            feature_type - "mat_peptide", "CDS", or None
    """
    if parse_mat_peptides(query_record):
        return "direct", "mat_peptide"

    if parse_cds_features(query_record):
        return "direct", "CDS"

    return "minimap", None