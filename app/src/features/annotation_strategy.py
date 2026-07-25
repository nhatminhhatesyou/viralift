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
    select_feature_type(record, alias_lookup=None) → "CDS" | "mat_peptide" | None
        Single entry point for feature type selection.
        With alias_lookup: scores both levels, returns the most informative.
        Without alias_lookup: existence check with polyprotein-shell guard.
        Returns None when the record has no usable gene-level annotation.

    get_strategy(query_record, alias_lookup=None) → ("direct" | "tblastn", feature_type | None)
        Combines type selection and routing into one call.
        Returns ("direct", feature_type) when the query has usable annotation,
        or ("tblastn", None) when coordinates must be lifted from the reference.
"""


def _parse_features_for_type(record: SeqRecord, feature_type: str):
    if feature_type == "mat_peptide":
        return parse_mat_peptides(record)
    if feature_type == "CDS":
        return parse_cds_features(record)
    return []


def _score_feature_type(
    record: SeqRecord,
    feature_type: str,
    alias_lookup: Dict,
) -> int:
    """
    Score how informative a feature level is given an alias lookup.

    Higher is better. Alias-resolved names dominate raw names so a smaller set
    of recognised gene-level annotations is preferred over many unknown labels.
    Excluded names count against the feature level.
    """
    features = _parse_features_for_type(record, feature_type)
    if not features:
        return 0

    from app.src.alias.gene_alias import apply_alias_to_feature

    resolved_count = raw_count = excluded_count = 0
    for feature in features:
        resolved = apply_alias_to_feature(feature, alias_lookup)
        name_source = resolved.get("name_source")
        if name_source in ("alias", "alias_conflict_resolved"):
            resolved_count += 1
        elif name_source == "raw":
            raw_count += 1
        elif name_source in ("excluded", "ignored", "ambiguous"):
            excluded_count += 1

    return (resolved_count * 100) + raw_count - excluded_count


def select_feature_type(
    record: SeqRecord,
    alias_lookup: Optional[Dict] = None,
    allowed_types: Optional[Tuple[str, ...]] = None,
) -> Optional[str]:
    """
    Choose which feature level (CDS or mat_peptide) to use for a record.

    With alias_lookup: scores both levels; returns the most informative, or None
    if neither level contains any useful gene-level names (scores ≤ 0).

    Without alias_lookup: simple existence check with a polyprotein-shell guard —
    returns None when the only CDS is a whole-genome polyprotein placeholder.

    Args:
        record:       A Biopython SeqRecord.
        alias_lookup: Optional {normalised_name: canonical} dict.
        allowed_types: Restrict the choice to these levels. Callers pass the
            reference's own level so a query cannot switch to a level the
            reference does not describe.

    Returns:
        "mat_peptide", "CDS", or None.
    """
    candidates = allowed_types or ("mat_peptide", "CDS")

    if alias_lookup:
        scores = {
            name: _score_feature_type(record, name, alias_lookup)
            for name in candidates
        }
        if not scores:
            return None
        best_type, best_score = max(scores.items(), key=lambda x: x[1])
        return best_type if best_score > 0 else None

    # No alias lookup — existence check with polyprotein-shell guard.
    if "mat_peptide" in candidates and parse_mat_peptides(record):
        return "mat_peptide"
    if "CDS" not in candidates:
        return None
    cds_list = parse_cds_features(record)
    if not cds_list:
        return None
    if len(cds_list) == 1:
        name = (cds_list[0].get("name") or "").lower()
        if "polyprotein" in name or name in ("", "unknown"):
            return None
    return "CDS"


def get_strategy(
    query_record: SeqRecord,
    alias_lookup: Optional[Dict] = None,
    allowed_types: Optional[Tuple[str, ...]] = None,
) -> Tuple[Literal["direct", "tblastn"], Optional[str]]:
    """
    Decide whether to use direct extraction or tblastn lifting for a query record.

    Calls select_feature_type once and returns both the routing decision and the
    selected feature type so the caller does not need to call select_feature_type
    again.

    Args:
        query_record: Query GenBank record to evaluate.
        alias_lookup: Optional alias lookup dict used to score feature levels.

    Returns:
        ("direct", feature_type) — query has usable gene-level annotation.
        ("tblastn", None)        — query lacks annotation; lift from reference.
    """
    feature_type = select_feature_type(query_record, alias_lookup, allowed_types)
    if feature_type is None:
        return "tblastn", None
    return "direct", feature_type
