from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.SeqRecord import SeqRecord

from app.src.alias.alias_registry import (
    detect_alias_config_for_record,
    get_detected_virus_name,
)
from app.src.features.annotation_strategy import get_best_feature_type, get_feature_type
from app.src.alias.gene_alias import apply_alias_to_features, load_alias_lookup
from app.src.io.genbank_parser import parse_cds_features, parse_mat_peptides


"""
Module: ref_loader.py

Purpose:
    Load and prepare reference features for the pipeline.

    Handles feature type detection, alias config resolution (auto-detect or
    user-provided), and alias normalization of reference feature names.

Main function:
    prepare_reference_features() — called once per run before processing queries.
"""


def prepare_reference_features(
    ref_record: SeqRecord,
    alias_config_arg: Optional[str],
    alias_registry_arg: str,
) -> Tuple[List[Dict], Optional[Path], Optional[str], Dict[str, str]]:
    """
    Parse reference features and optionally normalize names using alias config.

    Priority for alias config:
        1. User-provided path (--alias-config)
        2. Auto-detected from virus alias registry
        3. No config → keep raw names

    Args:
        ref_record:         Reference SeqRecord.
        alias_config_arg:   Optional CLI value for --alias-config.
        alias_registry_arg: Path to virus_alias_registry.json.

    Returns:
        Tuple of:
            - ref_features:         Parsed (and optionally alias-normalized) features.
            - alias_config_path:    Path to alias config used, or None.
            - detected_virus_name:  Virus name from registry, or None.
            - alias_lookup:         {raw_name: canonical} dict, or {}.
    """
    preliminary_feature_type = get_feature_type(ref_record)

    if preliminary_feature_type == "mat_peptide":
        ref_features = parse_mat_peptides(ref_record)
    elif preliminary_feature_type == "CDS":
        ref_features = parse_cds_features(ref_record)
    else:
        raise ValueError("Reference record has no CDS or mat_peptide features.")

    if not ref_features:
        raise ValueError(f"Reference record has no {preliminary_feature_type} features.")

    alias_config_path: Optional[Path] = None
    detected_virus_name: Optional[str] = None
    alias_lookup: Dict[str, str] = {}

    if alias_config_arg:
        alias_config_path = Path(alias_config_arg)
    else:
        registry_path = Path(alias_registry_arg)
        try:
            alias_config_path = detect_alias_config_for_record(ref_record, registry_path)
            if alias_config_path is not None:
                detected_virus_name = get_detected_virus_name(ref_record, registry_path)
        except (FileNotFoundError, ValueError):
            alias_config_path = None

    if alias_config_path is not None:
        alias_lookup = load_alias_lookup(alias_config_path)
        feature_type = get_best_feature_type(ref_record, alias_lookup) or preliminary_feature_type
        if feature_type == "mat_peptide":
            ref_features = parse_mat_peptides(ref_record)
        elif feature_type == "CDS":
            ref_features = parse_cds_features(ref_record)
        ref_features = apply_alias_to_features(ref_features, alias_lookup)

    return ref_features, alias_config_path, detected_virus_name, alias_lookup
