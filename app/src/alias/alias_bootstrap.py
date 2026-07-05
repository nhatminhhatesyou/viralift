import re
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple

from Bio.SeqRecord import SeqRecord

from app.src.alias.alias_classifier import (
    GENERIC_NAME_BLACKLIST,
    classify_alias_candidate,
)
from app.src.alias.alias_manager import (
    load_alias_config,
    load_registry,
    save_registry,
    save_validated_alias_config,
)
from app.src.alias.gene_alias import normalize_text
from app.src.features.annotation_strategy import select_feature_type
from app.src.io.genbank_parser import (
    _LOOKUP_QUALIFIER_KEYS,
    parse_cds_features,
    parse_mat_peptides,
)
from app.src.lifting.tblastn_lifter import process_one_query_record


"""
Module: alias_bootstrap.py

Purpose:
    Build alias suggestions for viruses that do not yet have an alias config.

Design:
    - Reference feature names are treated as temporary canonical names.
    - tblastn lifting provides coordinate evidence from reference -> query.
    - Query annotations are matched to lifted reference features by IoU.
    - Each raw qualifier string is classified independently as:
        save_alias / manual_review / ignore
    - User approval is still required before writing aliases permanently.
"""


def parse_features_for_type(record: SeqRecord, feature_type: Optional[str]) -> List[Dict]:
    """Parse a record using the selected feature type."""
    if feature_type == "mat_peptide":
        return parse_mat_peptides(record)
    if feature_type == "CDS":
        return parse_cds_features(record)
    return []


def select_bootstrap_feature_type(record: SeqRecord) -> Optional[str]:
    """
    Select feature type for a new-virus bootstrap run.

    This intentionally uses the no-alias selection path because a new virus
    does not have an alias lookup yet.
    """
    return select_feature_type(record, alias_lookup=None)


def build_seed_alias_config_from_ref(
    ref_record: SeqRecord,
    ref_features: List[Dict],
    virus_name: Optional[str] = None,
) -> Dict:
    """
    Build a minimal alias config using reference feature names as canonicals.

    The current alias lookup builder already maps canonical keys to themselves,
    so each new canonical starts with an empty alias list.
    """
    canonical_names: Dict[str, List[str]] = {}
    for feature in ref_features:
        name = feature.get("name")
        if not name:
            continue
        canonical_names.setdefault(str(name), [])

    return {
        "virus": virus_name or _get_organism(ref_record) or ref_record.id,
        "notes": "Bootstrapped from reference feature names. Review aliases before production use.",
        "ignored_names": sorted(GENERIC_NAME_BLACKLIST),
        "ambiguous_names": [],
        "canonical_names": canonical_names,
    }


def collect_query_name_candidates(
    query_record: SeqRecord,
    feature_type: str,
) -> List[Dict]:
    """
    Collect query annotation features and naming qualifier values.

    Nucleotide/protein sequences are not included.
    """
    features = parse_features_for_type(query_record, feature_type)
    candidates = []
    for feature in features:
        raw_values = []
        seen = set()
        for field in _LOOKUP_QUALIFIER_KEYS:
            value = feature.get(field)
            if not value:
                continue
            for part in split_compound_name(value):
                key = (field, normalize_text(part))
                if not key[1] or key in seen:
                    continue
                seen.add(key)
                raw_values.append({"field": field, "value": part})

        if not raw_values:
            continue

        candidates.append({
            "record_id": query_record.id,
            "feature_type": feature_type,
            "query_name": feature.get("name"),
            "start": feature.get("start"),
            "end": feature.get("end"),
            "strand": feature.get("strand"),
            "length": feature.get("length"),
            "order": feature.get("order"),
            "raw_values": raw_values,
        })

    return candidates


def split_compound_name(value: str) -> List[str]:
    """
    Split semicolon-delimited qualifier text while keeping the full string.

    Example:
        "M; ORF6" -> ["M; ORF6", "M", "ORF6"]
    """
    if not value:
        return []
    values = [value.strip()]
    parts = [part.strip() for part in value.split(";") if part.strip()]
    if len(parts) > 1:
        values.extend(parts)
    result = []
    seen = set()
    for item in values:
        norm = normalize_text(item)
        if norm and norm not in seen:
            seen.add(norm)
            result.append(item)
    return result


def compute_iou(
    start_a: Optional[int],
    end_a: Optional[int],
    start_b: Optional[int],
    end_b: Optional[int],
) -> float:
    """Compute inclusive interval IoU for two 1-based coordinate spans."""
    if None in (start_a, end_a, start_b, end_b):
        return 0.0
    if start_a > end_a or start_b > end_b:
        return 0.0
    intersection = max(0, min(end_a, end_b) - max(start_a, start_b) + 1)
    union = max(end_a, end_b) - min(start_a, start_b) + 1
    return intersection / union if union else 0.0


def build_coordinate_supported_alias_suggestions(
    ref_record: SeqRecord,
    query_records: List[SeqRecord],
    ref_features: List[Dict],
    ref_feature_type: str,
    min_iou: float = 0.90,
    min_coverage: float = 0.5,
    min_identity: float = 0.3,
    evalue: float = 1e-5,
    rescue_window: int = 200,
    diagnostics: Optional[Dict] = None,
    progress_callback: Optional[Callable[[int, int, str], None]] = None,
) -> List[Dict]:
    """
    Produce alias suggestions from query annotations + tblastn coordinate evidence.
    """
    suggestions: List[Dict] = []
    if diagnostics is not None:
        diagnostics.clear()
        diagnostics.update({
            "total_records": len(query_records),
            "records_without_usable_annotation": 0,
            "records_without_name_candidates": 0,
            "records_tblastn_run": 0,
            "records_with_lifted_features": 0,
            "candidate_features_total": 0,
            "lifted_features_total": 0,
            "matched_query_features_total": 0,
            "raw_suggestions_total": 0,
            "query_feature_type_counts": {},
        })

    total_records = len(query_records)
    for index, query_record in enumerate(query_records, start=1):
        if progress_callback is not None:
            progress_callback(index - 1, total_records, f"Checking annotation: {query_record.id}")

        query_feature_type = select_bootstrap_feature_type(query_record)
        if query_feature_type is None:
            if diagnostics is not None:
                diagnostics["records_without_usable_annotation"] += 1
            if progress_callback is not None:
                progress_callback(index, total_records, f"Skipped no usable annotation: {query_record.id}")
            continue
        if diagnostics is not None:
            counts = diagnostics["query_feature_type_counts"]
            counts[query_feature_type] = counts.get(query_feature_type, 0) + 1

        query_features = collect_query_name_candidates(query_record, query_feature_type)
        if not query_features:
            if diagnostics is not None:
                diagnostics["records_without_name_candidates"] += 1
            if progress_callback is not None:
                progress_callback(index, total_records, f"Skipped no name candidates: {query_record.id}")
            continue
        if diagnostics is not None:
            diagnostics["candidate_features_total"] += len(query_features)

        if progress_callback is not None:
            progress_callback(index - 1, total_records, f"Running tblastn: {query_record.id}")

        lifted_features = process_one_query_record(
            ref_record=ref_record,
            query_record=query_record,
            ref_cds=ref_features,
            ref_feature_type=ref_feature_type,
            min_coverage=min_coverage,
            min_identity=min_identity,
            evalue=evalue,
            rescue_window=rescue_window,
            quiet=True,
        )
        if diagnostics is not None:
            diagnostics["records_tblastn_run"] += 1
            diagnostics["lifted_features_total"] += len(lifted_features)
            if lifted_features:
                diagnostics["records_with_lifted_features"] += 1

        for query_feature in query_features:
            match = find_best_lifted_match(query_feature, lifted_features)
            if match is None or match["iou"] < min_iou:
                continue
            if diagnostics is not None:
                diagnostics["matched_query_features_total"] += 1

            for raw in query_feature["raw_values"]:
                classified = classify_alias_candidate(
                    raw_value=raw["value"],
                    field=raw["field"],
                    canonical_name=match["canonical_name"],
                    evidence=match,
                )
                suggestions.append({
                    "record_id": query_record.id,
                    "query_feature_type": query_feature_type,
                    "query_name": query_feature.get("query_name"),
                    "query_start": query_feature.get("start"),
                    "query_end": query_feature.get("end"),
                    "query_strand": query_feature.get("strand"),
                    "canonical_name": match["canonical_name"],
                    "tblastn_start": match["tblastn_start"],
                    "tblastn_end": match["tblastn_end"],
                    "tblastn_strand": match["tblastn_strand"],
                    "iou": round(match["iou"], 4),
                    "coverage": match.get("coverage"),
                    "identity": match.get("identity"),
                    "raw_value": raw["value"],
                    "field": raw["field"],
                    **classified,
                })
                if diagnostics is not None:
                    diagnostics["raw_suggestions_total"] += 1

    deduped = deduplicate_suggestions(suggestions)
    if diagnostics is not None:
        diagnostics["deduplicated_suggestions_total"] = len(deduped)
    if progress_callback is not None:
        progress_callback(total_records, total_records, "Suggestion generation complete")
    return deduped


def find_best_lifted_match(query_feature: Dict, lifted_features: Iterable) -> Optional[Dict]:
    """Find the lifted reference feature with highest IoU to a query annotation."""
    best: Optional[Dict] = None
    for lifted in lifted_features:
        iou = compute_iou(
            query_feature.get("start"),
            query_feature.get("end"),
            lifted.query_start,
            lifted.query_end,
        )
        if iou <= 0:
            continue
        strand_match = query_feature.get("strand") == lifted.strand
        candidate = {
            "canonical_name": getattr(lifted, "canonical_name", None) or lifted.name,
            "tblastn_start": lifted.query_start,
            "tblastn_end": lifted.query_end,
            "tblastn_strand": lifted.strand,
            "iou": iou,
            "coverage": lifted.coverage,
            "identity": lifted.identity,
            "status": lifted.status,
            "strand_match": strand_match,
        }
        if best is None or candidate["iou"] > best["iou"]:
            best = candidate
    return best


def deduplicate_suggestions(suggestions: List[Dict]) -> List[Dict]:
    """
    Merge duplicate raw_value/canonical suggestions across records.

    Keeps one row per raw_value/field/canonical_name and records how many
    unique query records support that alias decision.
    """
    grouped: Dict[Tuple[str, str, str], Dict] = {}
    for row in suggestions:
        key = (
            normalize_text(row.get("raw_value")),
            row.get("field") or "",
            row.get("canonical_name") or "",
        )
        record_id = row.get("record_id")
        existing = grouped.get(key)
        if existing is None:
            grouped[key] = {
                **row,
                "_support_records": {record_id} if record_id else set(),
            }
            continue
        if record_id:
            existing["_support_records"].add(record_id)
        if row.get("score", 0) > existing.get("score", 0):
            support_records = existing["_support_records"]
            grouped[key] = {**row, "_support_records": support_records}

    rows = []
    for row in grouped.values():
        support_records = sorted(row.pop("_support_records", set()))
        row["support_count"] = len(support_records)
        row["support_records"] = ", ".join(support_records[:10])
        if len(support_records) > 10:
            row["support_records"] += f", ... (+{len(support_records) - 10})"
        rows.append(row)

    return sorted(
        rows,
        key=lambda r: (
            r.get("suggested_action") != "save_alias",
            r.get("suggested_action") != "manual_review",
            -(r.get("score") or 0),
            r.get("canonical_name") or "",
            r.get("raw_value") or "",
        ),
    )


def write_new_alias_config(config: Dict, output_path: Path) -> Path:
    """Write a new alias config JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    save_validated_alias_config(output_path, config, create_backup=False)
    return output_path


def apply_approved_alias_suggestions(
    alias_config_path: Path,
    approved_rows: List[Dict],
    ignored_rows: Optional[List[Dict]] = None,
) -> Dict:
    """Apply user-approved alias suggestions to an alias config file."""
    config = load_alias_config(alias_config_path)

    canonical_names = config.setdefault("canonical_names", {})
    for row in approved_rows:
        raw_value = row.get("raw_value")
        canonical_name = row.get("canonical_name")
        if not raw_value or not canonical_name:
            continue
        aliases = canonical_names.setdefault(canonical_name, [])
        if raw_value != canonical_name and raw_value not in aliases:
            aliases.append(raw_value)

    ignored = config.setdefault("ignored_names", [])
    for row in ignored_rows or []:
        raw_value = row.get("raw_value")
        if raw_value and raw_value not in ignored:
            ignored.append(raw_value)

    config["ignored_names"] = sorted(set(ignored), key=normalize_text)
    save_validated_alias_config(alias_config_path, config)
    return config


def append_alias_registry_entry(
    registry_path: Path,
    virus_name: str,
    keywords: List[str],
    alias_config_path: Path,
) -> Dict:
    """Append or update one registry entry for a bootstrapped virus."""
    registry = load_registry(registry_path)

    viruses = registry.setdefault("viruses", [])
    config_str = str(alias_config_path)
    entry = {
        "virus_name": virus_name,
        "keywords": [kw for kw in keywords if kw],
        "alias_config": config_str,
    }

    for i, existing in enumerate(viruses):
        if existing.get("alias_config") == config_str or existing.get("virus_name") == virus_name:
            viruses[i] = entry
            break
    else:
        viruses.append(entry)

    save_registry(registry_path, registry)

    return entry


def safe_alias_filename(virus_name: str) -> str:
    """Create a conservative alias config filename stem from a virus name."""
    stem = re.sub(r"[^a-zA-Z0-9]+", "_", virus_name.strip().lower()).strip("_")
    return f"{stem or 'new_virus'}_alias.json"


def _get_organism(record: SeqRecord) -> str:
    return record.annotations.get("organism", "") or record.description or record.id
