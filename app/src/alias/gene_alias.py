import json
from pathlib import Path
from typing import Dict, List, Optional


"""
Module: gene_alias.py

Purpose:
    Load alias configuration from JSON files and normalize raw gene/feature names
    into canonical names.

Design:
    - Alias data is stored outside the codebase in config/*.json
    - This module contains only logic for:
        1. text normalization
        2. config loading
        3. alias lookup construction
        4. alias resolution
        5. applying canonical naming to parsed features

Expected JSON format:
    {
      "virus": "PRRSV",
      "canonical_names": {
        "GP5": ["GP5", "ORF5", "glycoprotein 5"],
        "M": ["M", "ORF6", "membrane protein"]
      }
    }

Notes:
    - Alias matching is exact after normalization.
    - This module does NOT do fuzzy matching.
    - When no alias match is found, the original raw name is preserved.
"""


def normalize_text(text: Optional[str]) -> str:
    """
    Normalize text for stable alias lookup.

    Normalization rules:
        - strip surrounding whitespace
        - lowercase
        - remove spaces
        - remove hyphens
        - remove underscores

    Args:
        text: Raw text to normalize

    Returns:
        Normalized string, or empty string if input is None/empty
    """
    if not text:
        return ""

    return (
        text.strip()
        .lower()
        .replace("-", "")
        .replace("_", "")
        .replace(" ", "")
    )


def load_alias_config(config_path: Path) -> Dict:
    """
    Load alias configuration from a JSON file.

    Args:
        config_path: Path to alias JSON config

    Returns:
        Parsed JSON dictionary

    Raises:
        FileNotFoundError: If config file does not exist
        ValueError: If JSON structure is invalid
    """
    if not config_path.exists():
        raise FileNotFoundError(f"Alias config file not found: {config_path}")

    with open(config_path, "r", encoding="utf-8") as handle:
        config_data = json.load(handle)

    if not isinstance(config_data, dict):
        raise ValueError("Alias config must be a JSON object.")

    if "canonical_names" not in config_data:
        raise ValueError("Alias config must contain 'canonical_names'.")

    if not isinstance(config_data["canonical_names"], dict):
        raise ValueError("'canonical_names' must be a JSON object.")

    return config_data


def build_alias_lookup(config_data: Dict) -> Dict[str, str]:
    """
    Build a normalized alias lookup table from config data.

    Output format:
        normalized_alias -> canonical_name

    Example:
        "orf5" -> "GP5"
        "glycoprotein5" -> "GP5"

    Args:
        config_data: Parsed alias config dictionary

    Returns:
        Alias lookup dictionary

    Raises:
        ValueError: If one normalized alias maps to multiple canonical names
    """
    lookup: Dict[str, str] = {}
    canonical_names = config_data.get("canonical_names", {})

    for canonical_name, aliases in canonical_names.items():
        if not isinstance(aliases, list):
            raise ValueError(
                f"Aliases for canonical name '{canonical_name}' must be a list."
            )

        # Include canonical name itself in the lookup
        all_names = [canonical_name] + aliases

        for alias in all_names:
            if not isinstance(alias, str):
                raise ValueError(
                    f"Alias under '{canonical_name}' must be a string. Got: {type(alias)}"
                )

            normalized_alias = normalize_text(alias)

            if not normalized_alias:
                continue

            if normalized_alias in lookup and lookup[normalized_alias] != canonical_name:
                raise ValueError(
                    f"Alias conflict detected: '{alias}' normalizes to '{normalized_alias}' "
                    f"and maps to both '{lookup[normalized_alias]}' and '{canonical_name}'."
                )

            lookup[normalized_alias] = canonical_name

    return lookup


def resolve_alias(raw_name: Optional[str], alias_lookup: Dict[str, str]) -> Optional[str]:
    """
    Resolve one raw name to its canonical name using alias lookup.

    Args:
        raw_name: Original feature name
        alias_lookup: normalized_alias -> canonical_name

    Returns:
        Canonical name if matched, otherwise original raw_name
    """
    if not raw_name:
        return raw_name

    normalized_name = normalize_text(raw_name)
    return alias_lookup.get(normalized_name, raw_name)


def apply_alias_to_feature(feature: Dict, alias_lookup: Dict[str, str]) -> Dict:
    """
    Apply alias normalization to a single feature dictionary.

    Input feature is expected to have at least:
        - name

    Output feature will include:
        - raw_name
        - name (canonical if matched, else original)
        - name_source ("alias" or "raw")

    Args:
        feature: Feature dictionary
        alias_lookup: normalized_alias -> canonical_name

    Returns:
        Updated feature dictionary
    """
    new_feature = feature.copy()

    raw_name = feature.get("name")
    in_lookup = normalize_text(raw_name) in alias_lookup

    if in_lookup:
        canonical_name = resolve_alias(raw_name, alias_lookup)
        name_source = "alias"
    else:
        # fallback: try product field
        product = feature.get("product")
        in_product = normalize_text(product) in alias_lookup if product else False
        if in_product:
            canonical_name = resolve_alias(product, alias_lookup)
            name_source = "product_alias"
        else:
            canonical_name = raw_name
            name_source = "raw"

    new_feature["raw_name"] = raw_name
    new_feature["name"] = canonical_name
    new_feature["name_source"] = name_source

    return new_feature


def apply_alias_to_features(features: List[Dict], alias_lookup: Dict[str, str]) -> List[Dict]:
    """
    Apply alias normalization to a list of feature dictionaries.

    Args:
        features: List of parsed features
        alias_lookup: normalized_alias -> canonical_name

    Returns:
        List of updated feature dictionaries
    """
    return [apply_alias_to_feature(feature, alias_lookup) for feature in features]


def load_alias_lookup(config_path: Path) -> Dict[str, str]:
    """
    Convenience wrapper to load a config file and build alias lookup.

    Args:
        config_path: Path to alias JSON config

    Returns:
        Alias lookup dictionary
    """
    config_data = load_alias_config(config_path)
    return build_alias_lookup(config_data)


def get_config_virus_name(config_path: Path) -> Optional[str]:
    """
    Return the virus name declared in the alias config.

    Args:
        config_path: Path to alias JSON config

    Returns:
        Virus name if present, else None
    """
    config_data = load_alias_config(config_path)
    return config_data.get("virus")

def build_canonical_to_ref_map(ref_features: List[Dict], alias_lookup: Dict[str, str]) -> Dict[str, str]:
    """
    Build a reverse map from canonical name -> ref feature name.

    Normalizes each ref feature name via alias lookup to get its canonical name,
    then stores canonical -> raw ref name. Used to translate query output names
    back to ref naming convention.

    Args:
        ref_features: Parsed features from the reference record
        alias_lookup: normalized_alias -> canonical_name

    Returns:
        Dict mapping canonical_name -> ref raw name
    """
    canonical_to_ref: Dict[str, str] = {}

    for feature in ref_features:
        raw_name = feature.get("name")
        if not raw_name:
            continue
        canonical = alias_lookup.get(normalize_text(raw_name), raw_name)
        canonical_to_ref[canonical] = raw_name

    return canonical_to_ref


def apply_ref_naming(normalized_features: List[Dict], canonical_to_ref: Dict[str, str]) -> List[Dict]:
    """
    Translate canonical names in normalized features back to ref naming convention.

    For each feature, if its canonical name exists in canonical_to_ref,
    replace the name with the ref name. Records the canonical name in
    'canonical_name' field for reference.

    Args:
        normalized_features: Features after apply_alias_to_features()
        canonical_to_ref: canonical_name -> ref raw name

    Returns:
        Updated feature list with ref-convention names
    """
    result = []

    for feature in normalized_features:
        new_feature = feature.copy()
        canonical = feature.get("name")

        if canonical in canonical_to_ref:
            new_feature["canonical_name"] = canonical
            new_feature["name"] = canonical_to_ref[canonical]

        result.append(new_feature)

    return result
