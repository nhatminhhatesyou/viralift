from typing import Dict, Literal, Optional, Tuple

from Bio.SeqRecord import SeqRecord
from app.src.io.genbank_parser import (
    parse_cds_features,
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
    get_feature_type(record)              — detect which feature types exist.
                                            Still exposed for reference-side detection
                                            in main.py / streamlit_app.py.
    get_best_feature_type(record, alias_lookup=None)
                                          — choose the most informative feature
                                            level for extraction.
Notes:
    - get_feature_type keeps the legacy existence-based priority
      (mat_peptide before CDS) for compatibility.
    - get_best_feature_type scores both mat_peptide and CDS, then chooses the
      level with the most useful gene-level names.
    - When alias_lookup is None (no config available), a legacy fallback applies:
      exactly one CDS whose name contains "polyprotein", is empty, or is "unknown".
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


def _parse_features_for_type(record: SeqRecord, feature_type: str):
    if feature_type == "mat_peptide":
        return parse_mat_peptides(record)
    if feature_type == "CDS":
        return parse_cds_features(record)
    return []


def _legacy_feature_type(record: SeqRecord) -> Optional[str]:
    """
    Fallback selector used when no alias lookup is available.

    This preserves the old behavior for unregistered viruses, including the
    single-polyprotein CDS shell check.
    """
    if parse_mat_peptides(record):
        return "mat_peptide"

    cds_list = parse_cds_features(record)
    if not cds_list:
        return None

    if len(cds_list) == 1:
        name = (cds_list[0].get("name") or "").lower()
        if "polyprotein" in name or name in ("", "unknown"):
            return None

    return "CDS"


def score_feature_type(record: SeqRecord, feature_type: str, alias_lookup: Optional[Dict] = None) -> int:
    """
    Score whether a feature level carries useful gene-level annotation.

    Higher is better. Alias-resolved names dominate raw names so a smaller set
    of recognized gene-level annotations is preferred over many unknown labels.
    Raw names are still weakly useful because the resolver can ask the user to
    map them. Ignored or ambiguous names count against the feature level.
    """
    features = _parse_features_for_type(record, feature_type)
    if not features:
        return 0

    if not alias_lookup:
        return len(features)

    from app.src.alias.gene_alias import apply_alias_to_feature

    resolved_count = 0
    raw_count = 0
    ignored_or_ambiguous_count = 0
    for feature in features:
        resolved = apply_alias_to_feature(feature, alias_lookup)
        name_source = resolved.get("name_source")

        if name_source in ("alias", "alias_conflict_resolved"):
            resolved_count += 1
        elif name_source == "raw":
            raw_count += 1
        elif name_source in ("ignored", "ambiguous"):
            ignored_or_ambiguous_count += 1

    return (resolved_count * 100) + raw_count - ignored_or_ambiguous_count


def get_best_feature_type(record: SeqRecord, alias_lookup: Optional[Dict] = None) -> Optional[str]:
    """
    Choose the most informative annotation feature level for a record.

    When an alias lookup exists, both mat_peptide and CDS are scored by how many
    useful names they contain. If neither level is informative, return None so
    the caller can route to tblastn. Without an alias lookup, preserve the legacy
    feature-existence behavior.
    """
    if not alias_lookup:
        return _legacy_feature_type(record)

    scores = {
        "mat_peptide": score_feature_type(record, "mat_peptide", alias_lookup),
        "CDS": score_feature_type(record, "CDS", alias_lookup),
    }

    best_type, best_score = max(scores.items(), key=lambda item: item[1])
    if best_score <= 0:
        return None

    return best_type


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
                          Retained for API compatibility; selection is based on
                          query feature usefulness.
        alias_lookup:     Optional alias lookup dict used to score CDS and
                          mat_peptide usefulness. When None, a legacy fallback is
                          used instead.

    Returns:
        "direct" or "tblastn"
    """
    query_type = get_best_feature_type(query_record, alias_lookup)

    if query_type is None:
        return "tblastn"

    return "direct"


def choose_strategy(
    record: SeqRecord,
    alias_lookup: Optional[Dict] = None,
) -> Tuple[Literal["direct", "tblastn"], Optional[str]]:
    """
    Compatibility wrapper for older validation scripts.

    Returns both the strategy and the selected feature type.
    """
    feature_type = get_best_feature_type(record, alias_lookup)
    strategy: Literal["direct", "tblastn"] = "direct" if feature_type else "tblastn"
    return strategy, feature_type
