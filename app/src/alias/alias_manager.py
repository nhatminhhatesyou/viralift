import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from app.src.alias.gene_alias import (
    canonical_entry_aliases,
    canonical_entry_parent,
    get_excluded_names,
    normalize_text,
)


"""
Module: alias_manager.py

Purpose:
    Small CRUD helpers for editing alias config JSON files from the UI.

Design:
    - Keep file edits centralized and conservative.
    - Create a timestamped backup before every save.
    - Validate obvious alias conflicts before writing.
"""


def load_registry(registry_path: Path) -> Dict:
    with registry_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def save_registry(registry_path: Path, registry: Dict, create_backup: bool = True) -> Path:
    backup_path = None
    if create_backup and registry_path.exists():
        backup_path = backup_registry(registry_path)
    with registry_path.open("w", encoding="utf-8") as handle:
        json.dump(registry, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return backup_path


def backup_registry(registry_path: Path) -> Path:
    backup_dir = registry_path.parent / "backups"
    backup_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = backup_dir / f"{registry_path.stem}.{stamp}{registry_path.suffix}"
    shutil.copy2(registry_path, backup_path)
    return backup_path


def list_registry_entries(registry_path: Path) -> List[Dict]:
    registry = load_registry(registry_path)
    return registry.get("viruses", [])


def add_registry_keyword(registry_path: Path, alias_config: str, keyword: str) -> Path:
    registry = load_registry(registry_path)
    keyword = str(keyword or "").strip()
    if not keyword:
        return None

    for entry in registry.get("viruses", []):
        if entry.get("alias_config") != alias_config:
            continue
        keywords = entry.setdefault("keywords", [])
        if normalize_text(keyword) not in {normalize_text(item) for item in keywords}:
            keywords.append(keyword)
            return save_registry(registry_path, registry)
        return None

    raise ValueError(f"Alias config not found in registry: {alias_config}")


def update_registry_entry(
    registry_path: Path,
    alias_config: str,
    virus_name: str,
    keywords: List[str],
) -> Path:
    registry = load_registry(registry_path)
    for entry in registry.get("viruses", []):
        if entry.get("alias_config") != alias_config:
            continue
        entry["virus_name"] = virus_name
        entry["keywords"] = _dedupe_strings(keywords)
        return save_registry(registry_path, registry)
    raise ValueError(f"Alias config not found in registry: {alias_config}")


def remove_registry_entry(
    registry_path: Path,
    alias_config: str,
    root: Optional[Path] = None,
    archive_alias_config: bool = False,
) -> Tuple[Path, Optional[Path]]:
    """
    Remove one virus entry from the alias registry.

    If archive_alias_config is True, copy the active alias JSON into the normal
    backups directory, then remove the active file. The config file is only
    archived when no remaining registry entry points to the same alias_config.
    """
    registry = load_registry(registry_path)
    entries = registry.get("viruses", [])
    kept_entries = [
        entry
        for entry in entries
        if entry.get("alias_config") != alias_config
    ]
    if len(kept_entries) == len(entries):
        raise ValueError(f"Alias config not found in registry: {alias_config}")

    registry["viruses"] = kept_entries
    registry_backup = save_registry(registry_path, registry)

    config_backup = None
    if archive_alias_config:
        still_referenced = any(
            entry.get("alias_config") == alias_config
            for entry in kept_entries
        )
        if still_referenced:
            raise ValueError(
                f"Cannot archive alias config still used by another registry entry: {alias_config}"
            )
        config_root = root or registry_path.parent
        config_path = resolve_config_path(Path(alias_config), config_root)
        if config_path.exists():
            config_backup = backup_alias_config(config_path)
            config_path.unlink()

    return registry_backup, config_backup


def resolve_config_path(config_path: Path, root: Path) -> Path:
    if config_path.is_absolute():
        return config_path
    return root / config_path


def load_alias_config(config_path: Path) -> Dict:
    with config_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def backup_alias_config(config_path: Path) -> Path:
    backup_dir = config_path.parent / "backups"
    backup_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = backup_dir / f"{config_path.stem}.{stamp}{config_path.suffix}"
    shutil.copy2(config_path, backup_path)
    return backup_path


def save_alias_config(config_path: Path, config: Dict, create_backup: bool = True) -> Path:
    backup_path = None
    if create_backup and config_path.exists():
        backup_path = backup_alias_config(config_path)

    with config_path.open("w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    return backup_path


def alias_config_to_tables(config: Dict) -> Tuple[List[Dict], List[Dict], List[Dict]]:
    alias_rows = []
    for canonical, entry in sorted(config.get("canonical_names", {}).items()):
        alias_rows.append({
            "canonical_name": canonical,
            "aliases": "\n".join(canonical_entry_aliases(canonical, entry) or []),
            "parent": canonical_entry_parent(entry) or "",
        })

    excluded_rows = [
        {"excluded_name": name}
        for name in get_excluded_names(config)
    ]
    return alias_rows, excluded_rows, []


def tables_to_alias_config(
    original_config: Dict,
    alias_rows: List[Dict],
    excluded_rows: Optional[List[Dict]] = None,
    ignored_rows: Optional[List[Dict]] = None,
    ambiguous_rows: Optional[List[Dict]] = None,
) -> Dict:
    config = dict(original_config)
    canonical_names: Dict = {}

    # Parent links are not editable through the alias table, so carry them over
    # from the config being edited. Without this, every UI save would silently
    # flatten the containment hierarchy.
    inherited_parents = {
        name: canonical_entry_parent(entry)
        for name, entry in (original_config.get("canonical_names", {}) or {}).items()
    }

    for row in alias_rows:
        canonical = str(row.get("canonical_name") or "").strip()
        if not canonical:
            continue
        aliases = _split_multiline_values(row.get("aliases"))
        aliases = [alias for alias in aliases if normalize_text(alias) != normalize_text(canonical)]
        existing = canonical_entry_aliases(canonical, canonical_names.get(canonical, []))
        merged = _dedupe_strings([*existing, *aliases])

        parent = str(row.get("parent") or "").strip() or inherited_parents.get(canonical)
        canonical_names[canonical] = (
            {"aliases": merged, "parent": parent} if parent else merged
        )

    config["canonical_names"] = canonical_names
    merged_excluded_rows = [
        *(excluded_rows or []),
        *(ignored_rows or []),
        *(ambiguous_rows or []),
    ]
    config["excluded_names"] = _dedupe_strings(
        row.get("excluded_name") or row.get("ignored_name") or row.get("ambiguous_name")
        for row in merged_excluded_rows
        if row.get("excluded_name") or row.get("ignored_name") or row.get("ambiguous_name")
    )
    config.pop("ignored_names", None)
    config.pop("ambiguous_names", None)
    return config


def validate_alias_config(config: Dict) -> List[str]:
    warnings = []
    alias_hits: Dict[str, Dict[str, set]] = {}

    for canonical, aliases in config.get("canonical_names", {}).items():
        canonical_norm = normalize_text(canonical)
        if not canonical_norm:
            warnings.append("Canonical name cannot be blank.")
            continue
        names = [canonical] + list(aliases or [])
        for name in names:
            norm = normalize_text(name)
            if not norm:
                continue
            hit = alias_hits.setdefault(norm, {"names": set(), "canonicals": set()})
            hit["names"].add(name)
            hit["canonicals"].add(canonical)

    for hit in alias_hits.values():
        if len(hit["canonicals"]) <= 1:
            continue
        names = ", ".join(f"`{name}`" for name in sorted(hit["names"], key=normalize_text))
        canonicals = ", ".join(
            f"`{canonical}`"
            for canonical in sorted(hit["canonicals"], key=normalize_text)
        )
        warnings.append(f"{names} maps to multiple canonicals: {canonicals}.")

    excluded = {normalize_text(name): name for name in get_excluded_names(config)}

    for norm, name in excluded.items():
        if norm in alias_hits:
            canonicals = _format_canonical_set(alias_hits[norm]["canonicals"])
            warnings.append(f"`{name}` is both excluded and an alias for {canonicals}.")

    return warnings


def save_validated_alias_config(
    config_path: Path,
    config: Dict,
    create_backup: bool = True,
) -> Path:
    """
    Validate and save an alias config through the normal backup-aware path.

    Raises:
        ValueError: If obvious alias conflicts are found.
    """
    warnings = validate_alias_config(config)
    if warnings:
        raise ValueError("Alias config validation failed:\n" + "\n".join(warnings))
    return save_alias_config(config_path, config, create_backup=create_backup)


def move_excluded_to_alias(config: Dict, excluded_name: str, canonical_name: str) -> Dict:
    excluded_norm = normalize_text(excluded_name)
    canonical_names = config.setdefault("canonical_names", {})
    aliases = canonical_names.setdefault(canonical_name, [])

    if excluded_name and excluded_name not in aliases and normalize_text(excluded_name) != normalize_text(canonical_name):
        aliases.append(excluded_name)

    config["excluded_names"] = [
        name
        for name in get_excluded_names(config)
        if normalize_text(name) != excluded_norm
    ]
    config.pop("ignored_names", None)
    config.pop("ambiguous_names", None)
    return config


def move_ignored_to_alias(config: Dict, ignored_name: str, canonical_name: str) -> Dict:
    """Backward-compatible wrapper for the old function name."""
    return move_excluded_to_alias(config, ignored_name, canonical_name)


def _split_multiline_values(value) -> List[str]:
    if value is None:
        return []
    text = str(value)
    parts = []
    for line in text.replace(";", "\n").replace(",", "\n").splitlines():
        item = line.strip()
        if item:
            parts.append(item)
    return parts


def _dedupe_strings(values) -> List[str]:
    result = []
    seen = set()
    for value in values:
        item = str(value or "").strip()
        key = normalize_text(item)
        if not key or key in seen:
            continue
        seen.add(key)
        result.append(item)
    return sorted(result, key=normalize_text)


def _format_canonical_set(canonicals: set) -> str:
    return ", ".join(
        f"`{canonical}`"
        for canonical in sorted(canonicals, key=normalize_text)
    )
