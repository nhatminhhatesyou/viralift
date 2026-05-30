from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.SeqRecord import SeqRecord

from app.src.alias.alias_registry import (
    detect_alias_config_for_record,
    get_detected_virus_name,
)
from app.src.features.annotation_strategy import select_feature_type
from app.src.alias.gene_alias import apply_alias_to_features, load_alias_lookup
from app.src.io.genbank_parser import parse_cds_features, parse_mat_peptides


"""
Module: ref_loader.py

Purpose:
    Load and prepare reference features for the pipeline.

    Handles alias config resolution (auto-detect or user-provided), feature type
    selection, and alias normalization of reference feature names — all in one
    pass with no redundant parsing.

Main function:
    prepare_reference_features() — called once per run before processing queries.
    Returns a 5-tuple: (ref_features, ref_feature_type, alias_config_path,
                        detected_virus_name, alias_lookup)
"""


def prepare_reference_features(
    ref_record: SeqRecord,
    alias_config_arg: Optional[str],
    alias_registry_arg: str,
) -> Tuple[List[Dict], str, Optional[Path], Optional[str], Dict[str, str]]:
    """
    Parse reference features and optionally normalize names using alias config.

    Processing order (no redundant parsing):
        1. Resolve alias config path (CLI arg → auto-detect → none)
        2. Load alias lookup from config
        3. Select feature type once with full alias information
        4. Parse ref features using the selected type
        5. Apply alias normalization

    Priority for alias config:
        1. User-provided path (--alias-config)
        2. Auto-detected from virus alias registry
        3. No config → raw names kept, legacy feature type detection

    Args:
        ref_record:         Reference SeqRecord.
        alias_config_arg:   Optional CLI value for --alias-config.
        alias_registry_arg: Path to virus_alias_registry.json.

    Returns:
        Tuple of:
            - ref_features:         Parsed (and optionally alias-normalized) features.
            - ref_feature_type:     "CDS" or "mat_peptide" — authoritative for this run.
            - alias_config_path:    Path to alias config used, or None.
            - detected_virus_name:  Virus name from registry, or None.
            - alias_lookup:         {normalised_name: canonical} dict, or {}.

    Raises:
        ValueError: If the reference has no usable CDS or mat_peptide features.
    """
    # Step 1: resolve alias config path
    alias_config_path: Optional[Path] = None
    detected_virus_name: Optional[str] = None

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

    # Step 2: load alias lookup
    alias_lookup: Dict[str, str] = {}
    if alias_config_path is not None:
        alias_lookup = load_alias_lookup(alias_config_path)

    # Step 3: select feature type once with full information
    feature_type = select_feature_type(ref_record, alias_lookup or None)
    if feature_type is None:
        raise ValueError("Reference record has no usable CDS or mat_peptide features.")

    # Step 4: parse features using the chosen type
    if feature_type == "mat_peptide":
        ref_features = parse_mat_peptides(ref_record)
    else:
        ref_features = parse_cds_features(ref_record)

    if not ref_features:
        raise ValueError(f"Reference record has no {feature_type} features.")

    # Step 5: apply alias normalization
    if alias_lookup:
        ref_features = apply_alias_to_features(ref_features, alias_lookup)

    return ref_features, feature_type, alias_config_path, detected_virus_name, alias_lookup
