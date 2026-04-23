import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.SeqRecord import SeqRecord


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
    - Payloads are saved to disk for inspection (LLM API not yet integrated).
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


# ---------------------------------------------------------------------
# Feature helpers
# ---------------------------------------------------------------------

def get_unresolved_features(normalized_features: List[Dict]) -> List[Dict]:
    """Return features that were not resolved by alias lookup (name_source == 'raw')."""
    return [f for f in normalized_features if f.get("name_source") == "raw"]


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
    from app.src.annotation.alias_registry import (
        detect_alias_config_for_record,
        get_detected_virus_name,
    )
    from app.src.annotation.annotation_strategy import choose_strategy
    from app.src.annotation.gene_alias import (
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
    _, ref_feature_type = choose_strategy(ref_record)

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
        _, feature_type = choose_strategy(record)

        # mat_peptide always has meaningful annotation regardless of ref type.
        # CDS on a mat_peptide virus: only block if it's a single CDS (polyprotein shell).
        # CDS count > 1 means real gene-level annotation even if ref uses mat_peptide.
        if feature_type == "mat_peptide":
            strategy = "direct"
        elif feature_type == "CDS":
            if ref_feature_type == "mat_peptide" and len(parse_cds_features(record)) == 1:
                feature_type = None
                strategy = "minimap"
            else:
                strategy = "direct"
        else:
            feature_type = None
            strategy = "minimap"

        if feature_type == "mat_peptide":
            raw_features = parse_mat_peptides(record)
        elif feature_type == "CDS":
            raw_features = parse_cds_features(record)
        else:
            results.append({
                "record_id": record.id,
                "strategy": "minimap",
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
                "strategy": strategy,
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
                "strategy": strategy,
                "feature_type": feature_type,
                "status": "all_resolved",
                "resolved": resolved_count,
                "total": len(normalized),
                "normalized_features": normalized,
            })
            continue

        results.append({
            "record_id": record.id,
            "strategy": strategy,
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
