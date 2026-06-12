import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

from app.src.alias.gene_alias import normalize_text


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
    for canonical, aliases in sorted(config.get("canonical_names", {}).items()):
        alias_rows.append({
            "canonical_name": canonical,
            "aliases": "\n".join(aliases or []),
        })

    ignored_rows = [
        {"ignored_name": name}
        for name in sorted(config.get("ignored_names", []), key=normalize_text)
    ]
    ambiguous_rows = [
        {"ambiguous_name": name}
        for name in sorted(config.get("ambiguous_names", []), key=normalize_text)
    ]
    return alias_rows, ignored_rows, ambiguous_rows


def tables_to_alias_config(
    original_config: Dict,
    alias_rows: List[Dict],
    ignored_rows: List[Dict],
    ambiguous_rows: List[Dict],
) -> Dict:
    config = dict(original_config)
    canonical_names: Dict[str, List[str]] = {}

    for row in alias_rows:
        canonical = str(row.get("canonical_name") or "").strip()
        if not canonical:
            continue
        aliases = _split_multiline_values(row.get("aliases"))
        aliases = [alias for alias in aliases if normalize_text(alias) != normalize_text(canonical)]
        existing = canonical_names.setdefault(canonical, [])
        canonical_names[canonical] = _dedupe_strings([*existing, *aliases])

    config["canonical_names"] = canonical_names
    config["ignored_names"] = _dedupe_strings(
        row.get("ignored_name")
        for row in ignored_rows
        if row.get("ignored_name")
    )
    config["ambiguous_names"] = _dedupe_strings(
        row.get("ambiguous_name")
        for row in ambiguous_rows
        if row.get("ambiguous_name")
    )
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

    ignored = {normalize_text(name): name for name in config.get("ignored_names", [])}
    ambiguous = {normalize_text(name): name for name in config.get("ambiguous_names", [])}

    for norm, name in ignored.items():
        if norm in alias_hits:
            canonicals = _format_canonical_set(alias_hits[norm]["canonicals"])
            warnings.append(f"`{name}` is both ignored and an alias for {canonicals}.")
    for norm, name in ambiguous.items():
        if norm in alias_hits:
            canonicals = _format_canonical_set(alias_hits[norm]["canonicals"])
            warnings.append(f"`{name}` is both ambiguous and an alias for {canonicals}.")
        if norm in ignored:
            warnings.append(f"`{name}` is both ignored and ambiguous.")

    return warnings


def move_ignored_to_alias(config: Dict, ignored_name: str, canonical_name: str) -> Dict:
    ignored_norm = normalize_text(ignored_name)
    canonical_names = config.setdefault("canonical_names", {})
    aliases = canonical_names.setdefault(canonical_name, [])

    if ignored_name and ignored_name not in aliases and normalize_text(ignored_name) != normalize_text(canonical_name):
        aliases.append(ignored_name)

    config["ignored_names"] = [
        name
        for name in config.get("ignored_names", [])
        if normalize_text(name) != ignored_norm
    ]
    return config


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
