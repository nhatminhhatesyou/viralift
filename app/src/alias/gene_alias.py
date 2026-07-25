import json
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from app.src.io.genbank_parser import _LOOKUP_QUALIFIER_KEYS


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


EXCLUDED_SENTINEL  = "__excluded__"
IGNORED_SENTINEL   = "__ignored__"
AMBIGUOUS_SENTINEL = "__ambiguous__"
EXCLUSION_SENTINELS = {EXCLUDED_SENTINEL, IGNORED_SENTINEL, AMBIGUOUS_SENTINEL}


def get_excluded_names(config_data: Dict) -> List[str]:
    """
    Return the unified do-not-map name list.

    `excluded_names` is the current config key. `ignored_names` and
    `ambiguous_names` are accepted as legacy input so older configs continue to
    behave the same until they are saved again by the manager.
    """
    result: List[str] = []
    seen = set()
    for field in ("excluded_names", "ignored_names", "ambiguous_names"):
        for name in config_data.get(field, []) or []:
            if not isinstance(name, str):
                continue
            key = normalize_text(name)
            if key and key not in seen:
                seen.add(key)
                result.append(name)
    return sorted(result, key=normalize_text)


def canonical_entry_aliases(canonical_name: str, entry) -> List[str]:
    """
    Read the alias list of one canonical entry, in either supported shape.

    Two shapes are accepted so configs can be migrated gradually:

        "ORF1a": ["polyprotein 1a", ...]                        # legacy
        "NSP2":  {"aliases": [...], "parent": "ORF1a"}          # current

    The second shape adds `parent`, which records that this canonical lies
    inside another one (a mature peptide within its polyprotein). `parent` is
    metadata only: it never changes what an alias resolves to, so the flat
    lookup is identical under both shapes.
    """
    if isinstance(entry, list):
        aliases = entry
    elif isinstance(entry, dict):
        aliases = entry.get("aliases", [])
        if not isinstance(aliases, list):
            raise ValueError(
                f"'aliases' for canonical name '{canonical_name}' must be a list."
            )
    else:
        raise ValueError(
            f"Entry for canonical name '{canonical_name}' must be a list of "
            f"aliases or an object with an 'aliases' key. Got: {type(entry)}"
        )
    return aliases


def canonical_entry_parent(entry) -> Optional[str]:
    """Return the declared parent of a canonical entry, if any."""
    if isinstance(entry, dict):
        parent = entry.get("parent")
        return parent if isinstance(parent, str) and parent else None
    return None


def build_parent_map(config_data: Dict) -> Dict[str, str]:
    """
    Map each canonical name to its parent, validating the hierarchy.

    Raises when a parent does not exist or when the parent chain forms a cycle,
    so a typo in a hand-edited config fails at load instead of silently
    producing a broken containment relation.
    """
    canonical_names = config_data.get("canonical_names", {})
    parents: Dict[str, str] = {}
    for name, entry in canonical_names.items():
        parent = canonical_entry_parent(entry)
        if parent is None:
            continue
        if parent not in canonical_names:
            raise ValueError(
                f"Canonical '{name}' declares parent '{parent}', which is not a "
                f"canonical name in this config."
            )
        parents[name] = parent

    for name in parents:
        seen = {name}
        cursor = parents.get(name)
        while cursor is not None:
            if cursor in seen:
                raise ValueError(
                    f"Cycle in canonical parent chain involving '{name}'."
                )
            seen.add(cursor)
            cursor = parents.get(cursor)
    return parents


def build_children_map(config_data: Dict) -> Dict[str, List[str]]:
    """
    Invert the parent map into {parent: [children]}.

    Stored as `parent` on each child so a write only ever touches one entry;
    this index exists for the read direction (UI trees, "what is inside X").
    """
    children: Dict[str, List[str]] = {}
    for child, parent in build_parent_map(config_data).items():
        children.setdefault(parent, []).append(child)
    for names in children.values():
        names.sort()
    return children


def build_alias_lookup(config_data: Dict) -> Dict[str, str]:
    """
    Build a normalized alias lookup table from config data.

    Output format:
        normalized_alias    -> canonical_name
        normalized_excluded -> EXCLUDED_SENTINEL

    Conflict handling:
        If the same normalized alias appears under two different canonical names
        (i.e. a name was accidentally duplicated), a warning is emitted and the
        alias is demoted to AMBIGUOUS_SENTINEL instead of raising an error. Add
        unsafe reusable names to "excluded_names" to suppress automatic mapping.

    Args:
        config_data: Parsed alias config dictionary

    Returns:
        Alias lookup dictionary
    """
    lookup: Dict[str, str] = {}
    canonical_names = config_data.get("canonical_names", {})

    for canonical_name, entry in canonical_names.items():
        aliases = canonical_entry_aliases(canonical_name, entry)

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
                existing = lookup[normalized_alias]
                if existing != AMBIGUOUS_SENTINEL:
                    # First conflict: demote to ambiguous and warn
                    warnings.warn(
                        f"Alias conflict: '{alias}' (normalized: '{normalized_alias}') "
                        f"maps to both '{existing}' and '{canonical_name}' — "
                        f"marked as ambiguous. Move it to 'excluded_names' in the "
                        f"alias config to suppress this warning.",
                        stacklevel=2,
                    )
                    lookup[normalized_alias] = AMBIGUOUS_SENTINEL
                # else: already ambiguous, nothing to do
                continue

            lookup[normalized_alias] = canonical_name

    # Excluded names: added last, but only if not already mapped to a real
    # canonical. These are intentionally skipped during automatic mapping.
    for excluded_name in get_excluded_names(config_data):
        normalized = normalize_text(excluded_name)
        if normalized and normalized not in lookup:
            lookup[normalized] = EXCLUDED_SENTINEL

    return lookup


def lookup_field_value(value: Optional[str], alias_lookup: Dict[str, str]) -> Optional[str]:
    """
    Look up one qualifier field value in the alias lookup.

    Tries the whole string first, then splits on ";" to handle compound
    note fields like "M; ORF6" or "GP5; ORF5; glycoprotein 5".
    Prefers a real canonical over exclusion sentinels.

    Args:
        value:        Raw qualifier string (e.g. "M; ORF6").
        alias_lookup: normalized_alias -> canonical_name mapping.

    Returns:
        Canonical string (including sentinels) if any part matched, else None.
    """
    if not value:
        return None

    # Try whole string first
    canonical = alias_lookup.get(normalize_text(value))
    if canonical is not None:
        return canonical

    # Try each semicolon-separated part
    parts = [p.strip() for p in value.split(";") if p.strip()]
    if len(parts) <= 1:
        return None

    _SPECIAL = EXCLUSION_SENTINELS
    best: Optional[str] = None
    for part in parts:
        canonical = alias_lookup.get(normalize_text(part))
        if canonical is None:
            continue
        if canonical not in _SPECIAL:
            return canonical          # real canonical — return immediately
        if best is None:
            best = canonical          # keep as fallback
    return best


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

    result = lookup_field_value(raw_name, alias_lookup)
    return result if result is not None else raw_name


def apply_alias_to_feature(feature: Dict, alias_lookup: Dict[str, str]) -> Dict:
    """
    Apply alias normalization to a single feature dictionary.

    Strategy:
        Iterate over all qualifier fields (gene, product, label, standard_name,
        locus_tag, note) in priority order and collect every field whose value
        has a hit in alias_lookup. Then decide:

        - 0 hits              → keep raw display name, name_source="raw"
        - 1 hit               → use that canonical
        - multiple hits, same canonical → unanimous, use it
        - multiple hits, different canonicals → conflict; use the hit from the
          highest-priority field (first in _LOOKUP_QUALIFIER_KEYS order)

        Exclusion sentinel hits are treated like misses for canonical resolution
        but still propagate name_source="excluded" if they are the only result.

    Output fields added to feature dict:
        - raw_name   : original display name before any resolution
        - name       : canonical name (or raw_name if unresolved)
        - name_source: one of "alias", "alias_conflict_resolved",
                       "excluded", "raw"

    Args:
        feature:      Feature dictionary produced by genbank_parser.
        alias_lookup: normalized_alias -> canonical_name mapping.

    Returns:
        Updated feature dictionary.
    """
    new_feature = feature.copy()
    raw_name = feature.get("name")

    # Collect (field, raw_value, canonical) for every qualifier that hits alias.
    # Uses lookup_field_value which handles compound note fields like "M; ORF6".
    hits: List[Tuple[str, str, str]] = []
    for field in _LOOKUP_QUALIFIER_KEYS:
        value = feature.get(field)
        if not value:
            continue
        canonical = lookup_field_value(value, alias_lookup)
        if canonical is not None:
            hits.append((field, value, canonical))

    _SPECIAL = EXCLUSION_SENTINELS

    alias_source_value = None

    if not hits:
        # Nothing matched alias at all → keep raw display name
        canonical_name = raw_name
        name_source = "raw"

    else:
        unique_canonicals = {canonical for _, _, canonical in hits}

        if len(unique_canonicals) == 1:
            # All hits agree (or only one hit) → use it
            _, alias_source_value, canonical_name = hits[0]
            if canonical_name in EXCLUSION_SENTINELS:
                canonical_name = raw_name
                alias_source_value = None
                name_source = "excluded"
            else:
                name_source = "alias"

        else:
            # Conflict: multiple different canonicals.
            # Priority: real canonical > excluded/ambiguous placeholders.
            real_hits = [
                (field, value, canonical)
                for field, value, canonical in hits
                if canonical not in _SPECIAL
            ]
            if real_hits:
                # At least one field resolves to a real canonical — use
                # the highest-priority one and note the conflict was resolved.
                _, alias_source_value, canonical_name = real_hits[0]
                name_source = "alias_conflict_resolved"
            else:
                canonical_name = raw_name
                name_source = "excluded"

    new_feature["raw_name"] = raw_name
    new_feature["alias_source_value"] = alias_source_value
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
