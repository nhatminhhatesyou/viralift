import hashlib
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.SeqRecord import SeqRecord

from app.src.alias.gene_alias import normalize_text


"""
Module: alias_payload.py

Purpose:
    Extract structured feature information from GenBank records
    and build JSON payloads intended for LLM-assisted alias building.

Two payload types:
    - "map_aliases":     virus is known, some feature names are unresolved
    - "build_alias_map": virus is new, need to build alias map from scratch

Main entry point:
    run_alias_pipeline() — processes ref + query records, returns per-record
    results and an aggregated LLM payload if anything is unresolved.

Notes:
    - Payloads are compact and intended for optional review-only LLM calls.
    - Nucleotide sequences are never included.
    - Only naming-relevant fields are extracted per feature.
"""


# ---------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------

def _extract_feature_info(feature: Dict) -> Dict:
    return {
        "name": feature.get("name"),
        "gene": feature.get("gene"),
        "product": feature.get("product"),
        "length": feature.get("length"),
        "order": feature.get("order"),
    }


def _get_organism(record: SeqRecord) -> str:
    return record.annotations.get("organism", "") or record.description or record.id


def _unique_ordered(values: List[str]) -> List[str]:
    seen = set()
    result = []
    for v in values:
        if v not in seen:
            seen.add(v)
            result.append(v)
    return result


# ---------------------------------------------------------------------
# Payload builders
# ---------------------------------------------------------------------

def build_map_aliases_payload(
    record: SeqRecord,
    unresolved_features: List[Dict],
    existing_canonical_names: List[str],
    feature_type: str,
) -> Dict:
    """
    Build payload for virus known, some feature names not yet in alias map.
    """
    return {
        "task": "map_aliases",
        "record_id": record.id,
        "organism": _get_organism(record),
        "feature_type": feature_type,
        "existing_canonical_names": existing_canonical_names,
        "unresolved_features": [_extract_feature_info(f) for f in unresolved_features],
    }


def build_new_virus_payload(
    record: SeqRecord,
    all_features: List[Dict],
    feature_type: str,
) -> Dict:
    """
    Build payload for virus not yet in registry. LLM builds alias map from scratch.
    """
    return {
        "task": "build_alias_map",
        "record_id": record.id,
        "organism": _get_organism(record),
        "feature_type": feature_type,
        "all_features": [_extract_feature_info(f) for f in all_features],
    }


def build_uncertain_suggestion_review_payload(
    virus_name: str,
    canonical_names: List[str],
    suggestions: List[Dict],
    ignored_names: Optional[List[str]] = None,
    ambiguous_names: Optional[List[str]] = None,
) -> Dict:
    """
    Build a compact payload for LLM review of uncertain alias suggestions.

    This payload is intentionally row-level and small: it contains only the
    deterministic suggestion context needed to advise the user. It never
    includes nucleotide sequence or full GenBank records.
    """
    return {
        "task": "review_uncertain_alias_suggestions",
        "virus": virus_name,
        "available_canonicals": canonical_names,
        "ignored_names": ignored_names or [],
        "ambiguous_names": ambiguous_names or [],
        "instructions": (
            "Review only the supplied uncertain alias rows. Recommend save_alias "
            "only when the raw alias clearly belongs to one available canonical. "
            "Treat available_canonicals as authoritative. If matching_available_canonical "
            "is present, prefer that canonical unless there is strong evidence of a true "
            "shared/ambiguous name. query_name is weak context and may describe a parent "
            "or combined ORF; do not reject a specific raw alias solely because query_name differs. "
            "Generic means a broad description without a specific gene/ORF identifier, "
            "for example 'polyprotein' or 'replicase polyprotein' alone. Names such as "
            "'polyprotein 1a', 'polyprotein 1b', and 'polyprotein 1ab' are specific "
            "and should usually map to ORF1a, ORF1b, and ORF1ab when coordinate evidence "
            "or matching_available_canonical supports that target. "
            "Use move_to_ambiguous for names that appear shared across genes, "
            "ignore for generic descriptions, and skip when evidence is weak."
        ),
        "suggestions": [
            _extract_suggestion_review_info(row, canonical_names)
            for row in suggestions
        ],
    }


def build_unresolved_name_review_payload(
    virus_name: str,
    canonical_names: List[str],
    unknown_items: Dict[str, Dict],
    ambiguous_items: Optional[Dict[str, Dict]] = None,
    ignored_names: Optional[List[str]] = None,
) -> Dict:
    """
    Build a payload for LLM review of names that did not get coordinate-backed
    alias suggestions.

    This is used on the manual Name review page. It is intentionally advisory:
    the model can suggest a mapping, but the UI still requires user approval.
    """
    unresolved_rows = []
    for representative, info in (unknown_items or {}).items():
        unresolved_rows.append(
            _extract_unresolved_review_info(
                representative,
                info,
                canonical_names,
                is_ambiguous=False,
            )
        )
    for representative, info in (ambiguous_items or {}).items():
        unresolved_rows.append(
            _extract_unresolved_review_info(
                representative,
                info,
                canonical_names,
                is_ambiguous=True,
            )
        )

    return {
        "task": "review_unresolved_names",
        "virus": virus_name,
        "available_canonicals": canonical_names,
        "ignored_names": ignored_names or [],
        "instructions": (
            "Review unresolved query annotation names that were not resolved by the "
            "alias config and may not have coordinate-backed tblastn suggestions. "
            "Recommend save_alias only when the representative name or one of its "
            "candidate qualifiers clearly belongs to exactly one available canonical. "
            "Recommend ignore for broad descriptions, comments, locus-like values, "
            "or parent annotations that should not become a reusable alias. Recommend "
            "skip when the evidence is insufficient. Use move_to_ambiguous when the "
            "same raw name is genuinely shared across multiple genes. Treat "
            "available_canonicals as authoritative and do not invent new canonicals. "
            "If a combined ORF name such as 'ORF1a/1b', 'ORF1a/b', or 'contains ORF1a "
            "and ORF1b' appears and ORF1ab is available, prefer ORF1ab; otherwise skip "
            "and let the user add a new canonical first."
        ),
        "suggestions": unresolved_rows,
    }


# ---------------------------------------------------------------------
# Feature helpers
# ---------------------------------------------------------------------

def get_unresolved_features(normalized_features: List[Dict]) -> List[Dict]:
    """Return features that were not resolved by alias lookup (name_source == 'raw')."""
    return [f for f in normalized_features if f.get("name_source") == "raw"]


def _extract_suggestion_review_info(row: Dict, canonical_names: Optional[List[str]] = None) -> Dict:
    return {
        "review_id": row.get("llm_review_id"),
        "raw_value": row.get("raw_value"),
        "field": row.get("field"),
        "matching_available_canonical": _matching_available_canonical(
            row.get("raw_value"),
            canonical_names or [],
        ),
        "deterministic_action": row.get("suggested_action"),
        "deterministic_confidence": row.get("confidence"),
        "deterministic_reason": row.get("reason"),
        "deterministic_score": row.get("score"),
        "canonical_candidate": row.get("canonical_name"),
        "query_feature_type": row.get("query_feature_type"),
        "query_name": row.get("query_name"),
        "support_count": row.get("support_count"),
        "support_records": row.get("support_records"),
        "iou": row.get("iou"),
        "coverage": row.get("coverage"),
        "identity": row.get("identity"),
    }


def _extract_unresolved_review_info(
    representative: str,
    info: Dict,
    canonical_names: List[str],
    is_ambiguous: bool,
) -> Dict:
    candidates = [str(value) for value in info.get("candidates", []) if value]
    matching = (
        _matching_available_canonical(representative, canonical_names)
        or _first_matching_candidate(candidates, canonical_names)
    )
    return {
        "review_id": _unresolved_review_id(representative),
        "raw_value": representative,
        "candidate_values": candidates,
        "matching_available_canonical": matching,
        "deterministic_action": "ambiguous" if is_ambiguous else "unresolved",
        "deterministic_confidence": "low",
        "deterministic_reason": (
            "raw name maps to multiple canonicals"
            if is_ambiguous
            else "raw name was not resolved by alias lookup"
        ),
        "canonical_candidate": matching,
        "support_count": len(info.get("records", [])),
        "support_records": list(info.get("records", []))[:12],
        "is_ambiguous": is_ambiguous,
    }


def _unresolved_review_id(representative: str) -> str:
    digest = hashlib.sha1(normalize_text(representative).encode("utf-8")).hexdigest()[:12]
    return f"unresolved_{digest}"


def _first_matching_candidate(candidates: List[str], canonical_names: List[str]) -> Optional[str]:
    for candidate in candidates:
        matched = _matching_available_canonical(candidate, canonical_names)
        if matched:
            return matched
    return None


def _matching_available_canonical(raw_value: str, canonical_names: List[str]) -> Optional[str]:
    raw_norm = normalize_text(raw_value or "")
    if not raw_norm:
        return None
    matches = []
    for canonical in canonical_names:
        canonical_norm = normalize_text(canonical)
        if not canonical_norm:
            continue
        if raw_norm == canonical_norm:
            matches.append((0, canonical))
        elif raw_norm in {
            f"{canonical_norm}protein",
            f"{canonical_norm}polyprotein",
            f"{canonical_norm}gene",
            f"{canonical_norm}cds",
        }:
            matches.append((1, canonical))
        elif _descriptive_orf_polyprotein_matches(raw_norm, canonical_norm):
            matches.append((2, canonical))
        elif _combined_orf_matches(raw_value or "", canonical_norm):
            matches.append((3, canonical))
    if not matches:
        return None
    return sorted(matches, key=lambda item: (item[0], len(item[1])))[0][1]


def _descriptive_orf_polyprotein_matches(raw_norm: str, canonical_norm: str) -> bool:
    if not raw_norm.startswith("polyprotein") and "polyprotein" not in raw_norm:
        return False

    canonical_match = re.fullmatch(r"orf(\d+[a-z]*)", canonical_norm)
    if not canonical_match:
        return False

    raw_match = re.search(r"polyprotein(?:orf)?(\d+[a-z]*)$", raw_norm)
    return bool(raw_match and raw_match.group(1) == canonical_match.group(1))


def _combined_orf_matches(raw_value: str, canonical_norm: str) -> bool:
    if canonical_norm != "orf1ab":
        return False
    raw = (raw_value or "").lower()
    raw_compact = normalize_text(raw)
    if raw_compact in {"orf1a/1b", "orf1a/b"}:
        return True
    if re.search(r"orf\s*1a\s*/\s*(?:orf\s*)?1?b", raw):
        return True
    return bool(
        re.search(r"contains\s+orf\s*1a\s+and\s+orf\s*1b", raw)
        or re.search(r"orf\s*1a\s+and\s+orf\s*1b", raw)
    )


# ---------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------

def run_alias_pipeline(
    ref_record: SeqRecord,
    query_records: List[SeqRecord],
    registry_path: Path,
    use_ref_naming: bool = True,
) -> Tuple[List[Dict], Optional[Dict]]:
    """
    Run alias normalization pipeline over all query records.

    Uses the reference record to detect the virus and load the alias config.
    Each query record is then normalized using that config.

    Args:
        ref_record:     Reference GenBank record (used for virus detection)
        query_records:  List of query GenBank records to normalize
        registry_path:  Path to virus_alias_registry.json

    Returns:
        Tuple of:
            - results: list of per-record dicts with keys:
                record_id, strategy, feature_type, status,
                resolved, total, normalized_features
            - payload: aggregated LLM payload dict, or None if all resolved
    """
    from app.src.alias.alias_registry import (
        detect_alias_config_for_record,
        get_detected_virus_name,
    )
    from app.src.features.annotation_strategy import select_feature_type
    from app.src.alias.gene_alias import (
        apply_alias_to_features,
        apply_ref_naming,
        build_canonical_to_ref_map,
        load_alias_lookup,
    )
    from app.src.io.genbank_parser import parse_cds_features, parse_mat_peptides

    # Detect virus and load alias config from reference
    alias_config_path = detect_alias_config_for_record(ref_record, registry_path)
    virus_name = get_detected_virus_name(ref_record, registry_path)

    alias_lookup = load_alias_lookup(alias_config_path) if alias_config_path else {}
    canonical_names = _unique_ordered(list(alias_lookup.values())) if alias_lookup else []

    # Determine expected feature type from reference
    ref_feature_type = select_feature_type(ref_record, alias_lookup or None)

    # Build canonical -> ref name map for output naming
    if use_ref_naming and alias_lookup:
        if ref_feature_type == "mat_peptide":
            ref_features = parse_mat_peptides(ref_record)
        elif ref_feature_type == "CDS":
            ref_features = parse_cds_features(ref_record)
        else:
            ref_features = []
        canonical_to_ref = build_canonical_to_ref_map(ref_features, alias_lookup)
    else:
        canonical_to_ref = {}

    results = []
    payload_records = []

    for record in query_records:
        feature_type = select_feature_type(record, alias_lookup or None)

        needs_lifting = feature_type is None

        if feature_type == "mat_peptide":
            raw_features = parse_mat_peptides(record)
        elif feature_type == "CDS":
            raw_features = parse_cds_features(record)
        else:
            results.append({
                "record_id": record.id,
                "strategy": "tblastn",
                "feature_type": None,
                "status": "no_annotation",
                "resolved": 0,
                "total": 0,
                "normalized_features": [],
            })
            continue

        if alias_config_path is None:
            # Virus not in registry — build alias map from scratch using ref
            results.append({
                "record_id": record.id,
                "strategy": "direct" if not needs_lifting else "tblastn",
                "feature_type": feature_type,
                "status": "new_virus",
                "resolved": 0,
                "total": len(raw_features),
                "normalized_features": raw_features,
            })
            payload_records.append(build_new_virus_payload(record, raw_features, feature_type))
            continue

        normalized = apply_alias_to_features(raw_features, alias_lookup)
        if canonical_to_ref:
            normalized = apply_ref_naming(normalized, canonical_to_ref)
        unresolved = get_unresolved_features(normalized)
        resolved_count = len(normalized) - len(unresolved)

        if not unresolved:
            results.append({
                "record_id": record.id,
                "strategy": "direct" if not needs_lifting else "tblastn",
                "feature_type": feature_type,
                "status": "all_resolved",
                "resolved": resolved_count,
                "total": len(normalized),
                "normalized_features": normalized,
            })
            continue

        results.append({
            "record_id": record.id,
            "strategy": "direct" if not needs_lifting else "tblastn",
            "feature_type": feature_type,
            "status": "has_unresolved",
            "resolved": resolved_count,
            "total": len(normalized),
            "normalized_features": normalized,
        })
        payload_records.append(build_map_aliases_payload(
            record=record,
            unresolved_features=unresolved,
            existing_canonical_names=canonical_names,
            feature_type=feature_type,
        ))

    if not payload_records:
        return results, None

    aggregated_payload = {
        "virus": virus_name,
        "alias_config": str(alias_config_path) if alias_config_path else None,
        "records_with_unresolved": len(payload_records),
        "records": payload_records,
    }

    return results, aggregated_payload


# ---------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------

def save_payload(payload: Dict, output_path: Path) -> None:
    """Save a payload dictionary to a JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
