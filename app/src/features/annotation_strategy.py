from typing import Literal, Optional

from Bio.SeqRecord import SeqRecord
from app.src.io.genbank_parser import parse_cds_features, parse_mat_peptides


"""
Module: annotation_strategy.py

Purpose:
    Determine the annotation strategy for a query GenBank record given a reference.

Public API:
    get_strategy(query_record, ref_feature_type) → "direct" | "tblastn"

    "direct"  — query already has meaningful gene-level annotation; extract and
                normalize names without alignment.
    "tblastn" — query lacks annotation (or has only a polyprotein shell);
                must lift coordinates from the reference via tblastn.

Internal helper:
    get_feature_type(record) — detect which feature type a record uses.
    Still exposed for reference-side detection in main.py / streamlit_app.py.

Notes:
    - mat_peptide takes priority over CDS because some viruses (e.g. FMDV) encode
      a single whole-genome polyprotein CDS with no gene-level names; the real
      names live in mat_peptide sub-features.
    - A query with exactly one CDS whose name contains "polyprotein" (or is unknown)
      is treated as an unannotated shell and routed to tblastn, regardless of ref type.
"""


def get_feature_type(record: SeqRecord) -> Optional[str]:
    """
    Detect the annotation feature type used by a GenBank record.

    Priority:
        1. mat_peptide present → "mat_peptide"
        2. CDS present         → "CDS"
        3. Neither             → None

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


def get_strategy(
    query_record: SeqRecord,
    ref_feature_type: str,
) -> Literal["direct", "tblastn"]:
    """
    Decide whether to use direct extraction or tblastn lifting for a query record.

    "direct"  — query has gene-level annotation compatible with direct extraction.
    "tblastn" — query lacks annotation or has only a polyprotein shell;
                coordinates must be lifted from the reference.

    Args:
        query_record:     Query GenBank record to evaluate.
        ref_feature_type: Feature type of the reference ("CDS" or "mat_peptide").

    Returns:
        "direct" or "tblastn"
    """
    query_type = get_feature_type(query_record)

    if query_type is None:
        return "tblastn"

    # A single CDS whose name looks like a polyprotein shell (no gene-level names)
    # should be lifted via tblastn regardless of whether the ref uses CDS or mat_peptide.
    cds_list = parse_cds_features(query_record)
    if query_type == "CDS" and len(cds_list) == 1:
        name = (cds_list[0].get("name") or "").lower()
        if "polyprotein" in name or name in ("", "unknown"):
            return "tblastn"

    return "direct"
