import json
from pathlib import Path
from typing import Dict, List, Optional

from Bio.SeqRecord import SeqRecord


"""
Module: alias_payload.py

Purpose:
    Extract structured feature information from GenBank records
    and build JSON payloads intended for LLM-assisted alias building.

Two payload types:
    - "map_aliases":    virus is known, some feature names are unresolved
    - "build_alias_map": virus is new, need to build alias map from scratch

Notes:
    - Payloads are saved to disk for inspection (LLM API not yet integrated).
    - Nucleotide sequences are never included.
    - Only naming-relevant fields are extracted per feature.
"""


def _extract_feature_info(feature: Dict) -> Dict:
    """
    Extract minimal naming-relevant info from one feature dict.

    Args:
        feature: Parsed feature dictionary

    Returns:
        Minimal feature info dict for LLM payload
    """
    return {
        "name": feature.get("name"),
        "gene": feature.get("gene"),
        "product": feature.get("product"),
        "length": feature.get("length"),
        "order": feature.get("order"),
    }


def _get_organism(record: SeqRecord) -> str:
    """Extract organism name from a GenBank record."""
    return record.annotations.get("organism", "") or record.description or record.id


def build_map_aliases_payload(
    record: SeqRecord,
    unresolved_features: List[Dict],
    existing_canonical_names: List[str],
    feature_type: str,
) -> Dict:
    """
    Build payload for the case where virus is known but some feature
    names are not yet in the alias map.

    Args:
        record: Query GenBank record
        unresolved_features: Features with name_source == "raw"
        existing_canonical_names: Canonical names already in alias config
        feature_type: "CDS" or "mat_peptide"

    Returns:
        Payload dictionary ready for JSON serialization
    """
    return {
        "task": "map_aliases",
        "record_id": record.id,
        "organism": _get_organism(record),
        "feature_type": feature_type,
        "existing_canonical_names": existing_canonical_names,
        "unresolved_features": [
            _extract_feature_info(f) for f in unresolved_features
        ],
    }


def build_new_virus_payload(
    record: SeqRecord,
    all_features: List[Dict],
    feature_type: str,
) -> Dict:
    """
    Build payload for the case where the virus is not yet in the registry.
    LLM needs to suggest canonical names and aliases from scratch.

    Args:
        record: Reference GenBank record (well-annotated, used as ground truth)
        all_features: All parsed features from the reference
        feature_type: "CDS" or "mat_peptide"

    Returns:
        Payload dictionary ready for JSON serialization
    """
    return {
        "task": "build_alias_map",
        "record_id": record.id,
        "organism": _get_organism(record),
        "feature_type": feature_type,
        "all_features": [
            _extract_feature_info(f) for f in all_features
        ],
    }


def get_unresolved_features(normalized_features: List[Dict]) -> List[Dict]:
    """
    Filter features that were not resolved by alias lookup.

    Args:
        normalized_features: Features after apply_alias_to_features()

    Returns:
        List of features with name_source == "raw"
    """
    return [f for f in normalized_features if f.get("name_source") == "raw"]


def save_payload(payload: Dict, output_path: Path) -> None:
    """
    Save a payload dictionary to a JSON file.

    Args:
        payload: Payload dictionary to save
        output_path: Destination file path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
