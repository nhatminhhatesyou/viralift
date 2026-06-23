from typing import Dict, List, Optional

from Bio.SeqRecord import SeqRecord

from app.src.alias.gene_alias import apply_alias_to_features
from app.src.io.genbank_parser import parse_cds_features, parse_mat_peptides
from app.src.lifting.base import LiftedFeature


"""
Module: direct_extractor.py

Purpose:
    Extract feature sequences directly from query annotation coordinates
    and normalize names using the alias lookup — no alignment needed.

Use case:
    Query records that already carry gene-level annotation (CDS or mat_peptide).
    Much faster than tblastn lifting; just parse, rename, and slice sequence.
"""


def direct_extract_with_alias(
    query_record: SeqRecord,
    query_feature_type: str,
    ref_features: List[Dict],
    alias_lookup: Dict[str, str],
) -> List[LiftedFeature]:
    """
    Extract annotated features from a query record and normalize their names.

    Args:
        query_record:       Query genome with existing annotation.
        query_feature_type: "CDS" or "mat_peptide".
        ref_features:       Reference features (alias-normalized).
                            Used to populate ref_start / ref_end by name match.
        alias_lookup:       Alias normalization lookup {raw_name: canonical}.

    Returns:
        List of LiftedFeature objects with method="direct". Status values:
        - ok: name resolved and exists in the reference feature set
        - unresolved_name: name was not resolved by the alias lookup
        - ambiguous_name: name is known ambiguous and needs user resolution
        - not_in_reference: name resolved, but the selected reference lacks it
    """
    if query_feature_type == "mat_peptide":
        query_features = parse_mat_peptides(query_record)
    else:
        query_features = parse_cds_features(query_record)

    if alias_lookup:
        query_features = apply_alias_to_features(query_features, alias_lookup)

    # canonical name → ref feature (for ref_start / ref_end)
    ref_by_name: Dict[str, Dict] = {f["name"]: f for f in ref_features}
    query_features = _collapse_duplicate_canonicals(query_features, ref_by_name)

    ignored: List[str] = []
    results: List[LiftedFeature] = []
    for qf in query_features:
        if qf.get("name_source") == "ignored":
            ignored.append(qf.get("raw_name") or qf["name"])
            continue

        name = qf["name"]
        name_source = qf.get("name_source")
        ref_match = ref_by_name.get(name)

        if name_source == "ambiguous":
            status = "ambiguous_name"
        elif alias_lookup and name_source == "raw":
            status = "unresolved_name"
        elif ref_match is None:
            status = "not_in_reference"
        else:
            status = "ok"

        start  = qf["start"]
        end    = qf["end"]
        strand = qf.get("strand", "+")

        seq = query_record.seq[start - 1: end]
        if strand == "-":
            seq = seq.reverse_complement()

        results.append(LiftedFeature(
            name=name,
            source_name=(
                qf.get("raw_name")
                if name_source in ("alias", "alias_conflict_resolved")
                else None
            ),
            ref_start=ref_match["start"] if ref_match else None,
            ref_end=ref_match["end"]   if ref_match else None,
            strand=strand,
            method="direct",
            query_start=start,
            query_end=end,
            sequence=str(seq),
            coverage=1.0,
            status=status,
        ))

    if ignored:
        print(
            f"  [WARN] {query_record.id}: {len(ignored)} feature(s) skipped "
            f"(ignored by alias config): {ignored}"
        )

    return results


def _collapse_duplicate_canonicals(
    features: List[Dict],
    ref_by_name: Dict[str, Dict],
) -> List[Dict]:
    """
    Keep one direct feature per canonical name when annotations contain
    overlapping duplicate entries.

    Some viral records annotate a polyprotein both as a full-length feature and
    as a shorter nested/sub-feature while both resolve to the same canonical
    name. For direct extraction, emitting both creates duplicated output rows.
    If the canonical exists in the reference, keep the feature whose length is
    closest to the reference feature. If it does not exist in the reference,
    keep the longest feature, which preserves the broadest query annotation.
    """
    best_by_name: Dict[str, Dict] = {}
    order_by_name: Dict[str, int] = {}

    for order, feature in enumerate(features):
        name = feature.get("name")
        if not name:
            continue

        existing = best_by_name.get(name)
        if existing is None:
            best_by_name[name] = feature
            order_by_name[name] = order
            continue

        chosen = _choose_better_duplicate(feature, existing, ref_by_name.get(name))
        if chosen is feature:
            best_by_name[name] = feature

    return [
        best_by_name[name]
        for name, _ in sorted(order_by_name.items(), key=lambda item: item[1])
    ]


def _choose_better_duplicate(
    candidate: Dict,
    current: Dict,
    ref_match: Optional[Dict],
) -> Dict:
    candidate_len = int(candidate.get("length") or _feature_len(candidate))
    current_len = int(current.get("length") or _feature_len(current))

    if ref_match:
        ref_len = int(ref_match.get("length") or _feature_len(ref_match))
        candidate_delta = abs(candidate_len - ref_len)
        current_delta = abs(current_len - ref_len)
        if candidate_delta != current_delta:
            return candidate if candidate_delta < current_delta else current

    if candidate_len != current_len:
        return candidate if candidate_len > current_len else current

    candidate_start = int(candidate.get("start") or 0)
    current_start = int(current.get("start") or 0)
    if candidate_start != current_start:
        return candidate if candidate_start < current_start else current

    return current


def _feature_len(feature: Dict) -> int:
    start = feature.get("start")
    end = feature.get("end")
    if start is None or end is None:
        return 0
    return int(end) - int(start) + 1
