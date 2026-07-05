# Auto-split from the original monolithic streamlit_app.py.
import pandas as pd
import re
import streamlit as st
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from app.src.alias.alias_manager import (
    add_registry_keyword,
    load_alias_config as manager_load_alias_config,
    resolve_config_path,
    save_validated_alias_config,
)
from app.src.alias.gene_alias import (
    AMBIGUOUS_SENTINEL,
    IGNORED_SENTINEL,
    apply_alias_to_features,
    load_alias_lookup,
    lookup_field_value,
    normalize_text,
)
from app.src.features.annotation_strategy import select_feature_type
from app.src.io.genbank_parser import _LOOKUP_QUALIFIER_KEYS
from app.src.io.run_logger import log_alias_added, log_canonical_added, log_run_complete, log_run_start
from app.src.pipeline import PipelineConfig, run_pipeline
from ui.i18n import _t
from ui.state import CONFIG_DIR, REGISTRY_PATH, ROOT


def _save_upload(uploaded_file) -> Path:
    """Save a Streamlit UploadedFile to the session temp dir. Returns path."""
    tmp_dir = Path(st.session_state.tmp.name)
    dest = tmp_dir / uploaded_file.name
    dest.write_bytes(uploaded_file.read())
    return dest


def _suggest_virus_name(record: SeqRecord) -> str:
    """Best-effort display name for a virus without an existing alias config."""
    organism = record.annotations.get("organism")
    if organism:
        return organism
    return record.description or record.id or "new virus"


def _record_metadata_candidates(record: SeqRecord) -> List[str]:
    candidates = [
        record.annotations.get("organism", ""),
        getattr(record, "description", "") or "",
        getattr(record, "name", "") or "",
        getattr(record, "id", "") or "",
    ]
    result = []
    seen = set()
    for item in candidates:
        value = str(item or "").strip()
        key = normalize_text(value)
        if key and key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _split_keywords(raw_text: str, virus_name: str) -> List[str]:
    """Parse comma/newline separated registry keywords and include virus name."""
    values = [virus_name]
    for part in re.split(r"[,\n]+", raw_text or ""):
        item = part.strip()
        if item:
            values.append(item)

    deduped = []
    seen = set()
    for item in values:
        key = normalize_text(item)
        if key and key not in seen:
            seen.add(key)
            deduped.append(item)
    return deduped


def _unique_alias_config_paths(filename: str) -> Tuple[Path, Path]:
    """Return absolute/relative alias config paths without overwriting existing files."""
    candidate = CONFIG_DIR / filename
    if not candidate.exists():
        return candidate, Path("app/config") / filename

    stem = Path(filename).stem
    suffix = Path(filename).suffix
    index = 2
    while True:
        next_name = f"{stem}_{index}{suffix}"
        candidate = CONFIG_DIR / next_name
        if not candidate.exists():
            return candidate, Path("app/config") / next_name
        index += 1


def _continue_with_existing_alias_config(entry: Dict, save_keyword: Optional[str] = None) -> None:
    alias_config = entry.get("alias_config")
    alias_config_path = resolve_config_path(Path(alias_config), ROOT)
    if save_keyword:
        add_registry_keyword(REGISTRY_PATH, alias_config, save_keyword)

    alias_lookup = load_alias_lookup(alias_config_path)
    ref_features = apply_alias_to_features(st.session_state.ref_features, alias_lookup)
    ignored = _load_ignored_names(alias_config_path)
    unknown = _scan_unknown_names(st.session_state.query_records, alias_lookup, ignored)
    unknown_ref = _scan_unknown_ref_names(ref_features, ignored)

    st.session_state.ref_features = ref_features
    st.session_state.alias_lookup = alias_lookup
    st.session_state.alias_config_path = alias_config_path
    st.session_state.virus_name = entry.get("virus_name")
    st.session_state.canonical_list = sorted(set(alias_lookup.values()))
    st.session_state.unknown_names = unknown
    st.session_state.unknown_ref_names = unknown_ref
    st.session_state.stage = "resolve" if (unknown or unknown_ref) else "running"


def _selected_dataframe_rows(state) -> List[int]:
    """Extract selected row indices from Streamlit dataframe selection state."""
    if not state:
        return []
    selection = getattr(state, "selection", None)
    if selection is None and isinstance(state, dict):
        selection = state.get("selection")
    if not selection:
        return []
    rows = getattr(selection, "rows", None)
    if rows is None and isinstance(selection, dict):
        rows = selection.get("rows")
    return list(rows or [])


def _scan_unknown_names(
    query_records: List[SeqRecord],
    alias_lookup: Dict[str, str],
    ignored_names: set,
) -> Dict[str, Dict]:
    """
    Return a dict keyed by representative name for every feature group in query
    records that needs user resolution: either the name is completely unknown
    (misses alias lookup entirely) or it is explicitly ambiguous (maps to
    AMBIGUOUS_SENTINEL, shared by multiple genes).

    Each value is:
        {
            "records":    [record_id, ...],   # records containing this feature
            "candidates": [val, ...],         # all qualifier values, priority order
            "ambiguous":  bool,               # True means known-ambiguous; False means unknown
        }

    Resolution logic mirrors apply_alias_to_feature:
        - ANY candidate hits a real canonical: feature is resolved, skip
        - ALL candidates miss OR hit AMBIGUOUS: include (unknown or ambiguous)
        - ALL candidates hit IGNORED: skip (intentionally excluded)
    """
    result: Dict[str, Dict] = {}

    for rec in query_records:
        selected_feature_type = select_feature_type(rec, alias_lookup)
        if selected_feature_type is None:
            continue

        for feat in rec.features:
            if feat.type != selected_feature_type:
                continue

            # Collect all unique qualifier values in priority order
            seen_vals: set = set()
            candidates = []
            for field in _LOOKUP_QUALIFIER_KEYS:
                val = feat.qualifiers.get(field, [None])[0]
                if val and val not in seen_vals:
                    seen_vals.add(val)
                    candidates.append(val)

            if not candidates:
                continue

            # Classify each candidate's hit using the same resolver as the core
            # pipeline, including semicolon-separated note fields like "ORF6; M".
            hits = {v: lookup_field_value(v, alias_lookup) for v in candidates}

            # If ANY candidate resolves to a real canonical, already handled, skip
            if any(
                h is not None and h not in (AMBIGUOUS_SENTINEL, IGNORED_SENTINEL)
                for h in hits.values()
            ):
                continue

            # If ALL candidates resolve to IGNORED, skip entirely.
            # Keep this in sync with alias matching: both sides use normalize_text,
            # so punctuation/case/spacing variants do not leak into Resolve.
            if all(normalize_text(v) in ignored_names for v in candidates):
                continue

            # Determine if this is ambiguous or fully unknown
            is_ambiguous = any(h == AMBIGUOUS_SENTINEL for h in hits.values())

            non_ignored_candidates = [
                v for v in candidates if normalize_text(v) not in ignored_names
            ]
            if not non_ignored_candidates:
                continue

            representative = non_ignored_candidates[0]

            if representative not in result:
                result[representative] = {
                    "records":   [],
                    "candidates": non_ignored_candidates,
                    "ambiguous": is_ambiguous,
                }
            if rec.id not in result[representative]["records"]:
                result[representative]["records"].append(rec.id)

    return result


def _scan_unknown_ref_names(ref_features: list, ignored_names: set) -> List[str]:
    """
    Return a sorted list of ref feature raw names that were not resolved by the
    alias DB (i.e. name_source == 'raw') and are not explicitly ignored.
    These will be lifted correctly (tblastn uses protein sequence), but their
    output name will just be the raw annotation name, no canonical key.
    """
    seen = []
    for f in ref_features:
        if f.get("name_source") != "raw":
            continue
        raw = f.get("raw_name") or f.get("name") or ""
        if normalize_text(raw) in ignored_names:
            continue
        if raw and raw not in seen:
            seen.append(raw)
    return sorted(seen)


def _add_new_canonicals_to_config(
    alias_config_path: Path,
    new_canonicals: List[str],
) -> int:
    """
    Add brand-new canonical entries (with empty alias list) to the alias JSON config.
    Skips any canonical key that already exists.
    Returns the number of entries actually added.
    """
    if not alias_config_path or not alias_config_path.exists():
        return 0
    cfg = manager_load_alias_config(alias_config_path)
    added = 0
    config_name = alias_config_path.name
    for name in new_canonicals:
        if name not in cfg["canonical_names"]:
            cfg["canonical_names"][name] = []
            log_canonical_added(config_name, name)
            added += 1
    if added:
        save_validated_alias_config(alias_config_path, cfg)
    return added


def _load_ignored_names(alias_config_path: Optional[Path]) -> set:
    if alias_config_path is None or not alias_config_path.exists():
        return set()
    cfg = manager_load_alias_config(alias_config_path)
    return {normalize_text(n) for n in cfg.get("ignored_names", [])}


def _build_effective_lookup(
    base_lookup: Dict[str, str],
    resolver: Dict[str, str],
) -> Dict[str, str]:
    """Merge base alias lookup with user resolver decisions."""
    effective = dict(base_lookup)
    for raw_name, canonical in resolver.items():
        if canonical and canonical != "-- ignore --":
            effective[normalize_text(raw_name)] = canonical
    return effective


def _canonical_to_ref_map(ref_features: list) -> Dict[str, str]:
    """
    Build {canonical_name: ref_raw_name} from ref features after alias normalization.
    e.g. {"Lpro": "Lab", "3Cpro": "3C", ...}
    """
    return {
        f["name"]: f["raw_name"]
        for f in ref_features
        if f.get("raw_name") and f.get("name")
    }


def _results_to_df(all_results, ref_name_map: Dict[str, str] = None) -> pd.DataFrame:
    rows = []
    for query_id, features in all_results:
        for lf in features:
            display_name = (
                ref_name_map.get(lf.name, lf.name)
                if ref_name_map else lf.name
            )
            rows.append({
                "record_id":  query_id,
                "name":       display_name,
                "source_name": lf.source_name or "",
                "start":      lf.query_start,
                "end":        lf.query_end,
                "strand":     lf.strand,
                "status":     lf.status,
                "coverage":   round(lf.coverage, 3) if lf.coverage else None,
                "identity":   lf.identity,
                "method":     lf.method,
                "sequence":   lf.sequence or "",
            })
    return pd.DataFrame(rows)


def _run_pipeline(
    ref_record, query_records, ref_features, ref_feature_type,
    effective_lookup, min_coverage, min_identity, evalue, rescue_window,
    progress_bar,
    virus_name: Optional[str] = None,
    alias_config_path=None,
    run_errors: Optional[List[Dict[str, str]]] = None,
) -> List[Tuple[str, List]]:

    log_run_start(
        ref_id=ref_record.id,
        n_queries=len(query_records),
        min_coverage=min_coverage,
        min_identity=min_identity,
        evalue=evalue,
        rescue_window=rescue_window,
        virus_name=virus_name,
        alias_config=Path(alias_config_path).name if alias_config_path else None,
    )

    all_results = []
    n = len(query_records)

    def update_progress(index: int, total: int, record_id: str) -> None:
        if index >= total:
            progress_bar.progress(1.0, text=_t("done"))
            return
        progress_bar.progress(
            index / total,
            text=f"{_t('processing')} {record_id}  ({index + 1}/{total})",
        )

    pipeline_result = run_pipeline(
        ref_record=ref_record,
        query_records=query_records,
        ref_features=ref_features,
        ref_feature_type=ref_feature_type,
        alias_lookup=effective_lookup,
        config=PipelineConfig(
            min_coverage=min_coverage,
            min_identity=min_identity,
            evalue=evalue,
            rescue_window=rescue_window,
            catch_record_errors=True,
        ),
        progress_callback=update_progress,
    )

    all_results = pipeline_result.all_results
    if run_errors is not None:
        run_errors.extend(pipeline_result.errors)

    summary = pipeline_result.summary
    log_run_complete(ref_id=ref_record.id, n_queries=n, summary=summary)

    return all_results


def _save_to_alias_config(alias_config_path: Path, mappings: Dict[str, str]) -> int:
    """
    Persist user-confirmed mappings into the alias JSON config file.

    For each (raw_name -> canonical) pair, appends raw_name to the canonical's
    alias list if not already present.

    Returns the number of new aliases written.
    """
    if not alias_config_path or not alias_config_path.exists():
        return 0

    cfg = manager_load_alias_config(alias_config_path)

    written = 0
    config_name = alias_config_path.name
    for raw_name, canonical in mappings.items():
        aliases = cfg["canonical_names"].get(canonical)
        if aliases is None:
            continue
        if raw_name not in aliases:
            aliases.append(raw_name)
            log_alias_added(config_name, raw_name, canonical)
            written += 1

    if written:
        save_validated_alias_config(alias_config_path, cfg)

    return written
