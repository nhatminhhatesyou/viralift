from typing import Dict, Literal, Optional

from Bio.SeqRecord import SeqRecord
from app.src.io.genbank_parser import (
    _LOOKUP_QUALIFIER_KEYS,
    parse_cds_features,
    parse_cds_features_basic,
    parse_mat_peptides,
)


"""
Module: annotation_strategy.py

Purpose:
    Determine the annotation strategy for a query GenBank record given a reference.

Public API:
    get_strategy(query_record, ref_feature_type, alias_lookup=None) → "direct" | "tblastn"

    "direct"  — query already has meaningful gene-level annotation; extract and
                normalize names without alignment.
    "tblastn" — query lacks annotation (or has only a polyprotein shell);
                must lift coordinates from the reference via tblastn.

Internal helpers:
    get_feature_type(record)              — detect which feature type a record uses.
                                            Still exposed for reference-side detection
                                            in main.py / streamlit_app.py.
    _all_names_ignored(record, alias_lookup) — True when every qualifier value across
                                            all CDS features maps to IGNORED_SENTINEL
                                            or is empty (i.e. no informative gene names).

Notes:
    - mat_peptide takes priority over CDS because some viruses (e.g. FMDV) encode
      a single whole-genome polyprotein CDS with no gene-level names; the real
      names live in mat_peptide sub-features.
    - When alias_lookup is provided, a CDS record is treated as an unannotated shell
      (routed to tblastn) if ALL qualifier values on ALL its CDS features map to
      IGNORED_SENTINEL — this replaces the old hardcoded "polyprotein" / len==1 check
      and covers any virus whose uninformative names are declared in its alias config.
    - When alias_lookup is None (no config available), a legacy fallback applies:
      exactly one CDS whose name contains "polyprotein", is empty, or is "unknown".
"""


def _all_names_ignored(record: SeqRecord, alias_lookup: Dict) -> bool:
    """
    Return True if every qualifier value across all CDS features maps to
    IGNORED_SENTINEL or is empty — indicating the record has no informative
    gene-level names and should be treated as unannotated.

    Args:
        record:       A Biopython SeqRecord (CDS-type).
        alias_lookup: Alias lookup dict loaded from the virus alias config.

    Returns:
        True  → all CDS features carry only ignored/empty names → route to tblastn.
        False → at least one feature has an informative name → route to direct.
    """
    from app.src.alias.gene_alias import lookup_field_value, IGNORED_SENTINEL

    cds_list = parse_cds_features_basic(record)
    if not cds_list:
        return False

    for feat in cds_list:
        for field in _LOOKUP_QUALIFIER_KEYS:
            value = feat.get(field)
            if not value:
                continue
            canonical = lookup_field_value(value, alias_lookup)
            if canonical != IGNORED_SENTINEL:
                return False  # found at least one informative name

    return True  # every non-empty qualifier value resolved to IGNORED_SENTINEL


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
    alias_lookup: Optional[Dict] = None,
) -> Literal["direct", "tblastn"]:
    """
    Decide whether to use direct extraction or tblastn lifting for a query record.

    "direct"  — query has gene-level annotation compatible with direct extraction.
    "tblastn" — query lacks annotation or has only a polyprotein shell;
                coordinates must be lifted from the reference.

    Args:
        query_record:     Query GenBank record to evaluate.
        ref_feature_type: Feature type of the reference ("CDS" or "mat_peptide").
        alias_lookup:     Optional alias lookup dict. When provided, a CDS record is
                          treated as a shell if ALL its qualifier values map to
                          IGNORED_SENTINEL (config-driven check). When None, a
                          legacy hardcoded fallback is used instead.

    Returns:
        "direct" or "tblastn"
    """
    query_type = get_feature_type(query_record)

    if query_type is None:
        return "tblastn"

    if query_type == "CDS":
        if alias_lookup:
            # Config-driven: treat as unannotated shell if every qualifier value on
            # every CDS feature is either empty or maps to IGNORED_SENTINEL.
            if _all_names_ignored(query_record, alias_lookup):
                return "tblastn"
        else:
            # Legacy fallback (no alias config available): single CDS with a
            # hardcoded uninformative name → unannotated shell.
            cds_list = parse_cds_features(query_record)
            if len(cds_list) == 1:
                name = (cds_list[0].get("name") or "").lower()
                if "polyprotein" in name or name in ("", "unknown"):
                    return "tblastn"

    return "direct"
