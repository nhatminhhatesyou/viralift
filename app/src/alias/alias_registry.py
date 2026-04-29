import json
from pathlib import Path
from typing import Dict, List, Optional

from Bio.SeqRecord import SeqRecord


"""
Module: alias_registry.py

Purpose:
    Resolve which alias configuration file should be used for a given virus.

Design:
    - Registry data is stored in config/virus_alias_registry.json
    - Each registry entry defines:
        * virus_name
        * keywords
        * alias_config

Example registry format:
    {
      "viruses": [
        {
          "virus_name": "PRRSV",
          "keywords": [
            "porcine reproductive and respiratory syndrome virus",
            "prrsv",
            "prrs virus"
          ],
          "alias_config": "config/prrsv_alias.json"
        }
      ]
    }

Matching logic:
    - Build a lowercase searchable text from the GenBank record
    - Search registry keywords in that text
    - Return the first matched alias config path

Notes:
    - Matching is keyword-based, not fuzzy matching
    - If no match is found, return None
"""


def load_alias_registry(registry_path: Path) -> Dict:
    """
    Load the virus alias registry JSON file.

    Args:
        registry_path: Path to virus_alias_registry.json

    Returns:
        Parsed registry dictionary

    Raises:
        FileNotFoundError: If registry file does not exist
        ValueError: If registry structure is invalid
    """
    if not registry_path.exists():
        raise FileNotFoundError(f"Alias registry file not found: {registry_path}")

    with open(registry_path, "r", encoding="utf-8") as handle:
        registry_data = json.load(handle)

    if not isinstance(registry_data, dict):
        raise ValueError("Alias registry must be a JSON object.")

    if "viruses" not in registry_data:
        raise ValueError("Alias registry must contain a 'viruses' field.")

    if not isinstance(registry_data["viruses"], list):
        raise ValueError("'viruses' must be a list.")

    return registry_data


def get_record_search_text(record: SeqRecord) -> str:
    """
    Build searchable lowercase text from a GenBank record.

    Sources used:
        - record.annotations['organism']
        - record.description
        - record.name
        - record.id

    Args:
        record: Biopython SeqRecord

    Returns:
        Lowercase combined text for keyword matching
    """
    organism = record.annotations.get("organism", "")
    description = getattr(record, "description", "") or ""
    name = getattr(record, "name", "") or ""
    record_id = getattr(record, "id", "") or ""

    text_parts = [organism, description, name, record_id]
    return " ".join(str(part) for part in text_parts if part).lower()


def find_registry_entry_for_text(search_text: str, registry_data: Dict) -> Optional[Dict]:
    """
    Find the first registry entry whose keyword appears in the search text.

    Args:
        search_text: Lowercase searchable text
        registry_data: Parsed registry dictionary

    Returns:
        Matched registry entry, or None if no match is found
    """
    viruses = registry_data.get("viruses", [])

    for entry in viruses:
        keywords = entry.get("keywords", [])

        if not isinstance(keywords, list):
            continue

        for keyword in keywords:
            if not isinstance(keyword, str):
                continue

            if keyword.lower() in search_text:
                return entry

    return None


def find_registry_entry_for_record(record: SeqRecord, registry_data: Dict) -> Optional[Dict]:
    """
    Find the best registry entry for a GenBank record.

    Args:
        record: Biopython SeqRecord
        registry_data: Parsed registry dictionary

    Returns:
        Matched registry entry, or None
    """
    search_text = get_record_search_text(record)
    return find_registry_entry_for_text(search_text, registry_data)


def get_alias_config_path_from_entry(entry: Dict) -> Optional[Path]:
    """
    Extract alias config path from one registry entry.

    Args:
        entry: One registry entry

    Returns:
        Path to alias config file, or None if missing
    """
    alias_config = entry.get("alias_config")
    if not alias_config:
        return None

    return Path(alias_config)


def detect_alias_config_for_record(
    record: SeqRecord,
    registry_path: Path,
) -> Optional[Path]:
    """
    Auto-detect the alias config file for a GenBank record.

    Args:
        record: Biopython SeqRecord
        registry_path: Path to virus_alias_registry.json

    Returns:
        Path to matched alias config file, or None if no match is found
    """
    registry_data = load_alias_registry(registry_path)
    entry = find_registry_entry_for_record(record, registry_data)

    if entry is None:
        return None

    return get_alias_config_path_from_entry(entry)


def get_detected_virus_name(record: SeqRecord, registry_path: Path) -> Optional[str]:
    """
    Return the detected virus_name from the registry for a given record.

    Args:
        record: Biopython SeqRecord
        registry_path: Path to virus_alias_registry.json

    Returns:
        virus_name string if matched, else None
    """
    registry_data = load_alias_registry(registry_path)
    entry = find_registry_entry_for_record(record, registry_data)

    if entry is None:
        return None

    return entry.get("virus_name")


def list_registered_viruses(registry_path: Path) -> List[str]:
    """
    List all registered virus names in the alias registry.

    Args:
        registry_path: Path to virus_alias_registry.json

    Returns:
        List of virus_name strings
    """
    registry_data = load_alias_registry(registry_path)
    viruses = registry_data.get("viruses", [])

    names: List[str] = []
    for entry in viruses:
        virus_name = entry.get("virus_name")
        if isinstance(virus_name, str) and virus_name.strip():
            names.append(virus_name)

    return names