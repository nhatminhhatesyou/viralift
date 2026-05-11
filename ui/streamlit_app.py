"""
ViraLift — Streamlit Web UI

Stages:
    upload  → user uploads ref + query, pipeline is configured
    resolve → unmapped gene names found, user maps or ignores each
    results → pipeline ran, show results + export options
"""

import sys
import json
import tempfile
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import streamlit as st
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ── project root on path ─────────────────────────────────────────────
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from app.src.alias.alias_registry import (
    detect_alias_config_for_record,
    get_detected_virus_name,
)
from app.src.features.annotation_strategy import get_feature_type, get_strategy
from app.src.alias.gene_alias import (
    apply_alias_to_features,
    load_alias_lookup,
    normalize_text,
)
from app.src.io.genbank_parser import (
    load_single_genbank,
    load_genbank_records,
    parse_cds_features,
    parse_mat_peptides,
    _LOOKUP_QUALIFIER_KEYS,
)
from app.src.features.direct_extractor import direct_extract_with_alias
from app.src.features.ref_loader import prepare_reference_features
from app.src.io.result_writer import summarize_counts
from app.src.lifting.tblastn_lifter import process_one_query_record
from app.src.lifting.base import LiftedFeature
from app.src.io.run_logger import (
    log_alias_added,
    log_canonical_added,
    log_session_decisions,
    log_run_start,
    log_run_complete,
    log_error,
    log_warning,
)

# ── constants ────────────────────────────────────────────────────────
REGISTRY_PATH  = ROOT / "app/config/virus_alias_registry.json"
STATUS_ICON = {
    "ok":                  "🟢",
    "ok_rescued":          "🟡",
    "invalid_boundaries":  "🟠",
    "low_coverage":        "🟠",
    "no_hit":              "🔴",
    "translation_fail":    "🔴",
    "direct":              "🔵",
}
GOOD_STATUSES = {"ok", "ok_rescued", "direct"}


# ═══════════════════════════════════════════════════════════════════
# Session-state bootstrap
# ═══════════════════════════════════════════════════════════════════

def _init_state():
    defaults = dict(
        stage="upload",
        tmp=None,               # tempfile.TemporaryDirectory object
        ref_record=None,
        query_records=None,
        ref_features=None,
        use_ref_names=False,    # if True, output uses ref's raw names instead of alias keys
        ref_feature_type=None,
        alias_lookup={},
        alias_config_path=None,
        virus_name=None,
        canonical_list=[],
        unknown_ref_names=[],   # [raw_name] ref features not in alias DB
        unknown_names={},       # {raw_name: [record_id, ...]}
        resolver={},            # {raw_name: canonical or "-- ignore --"}
        all_results=None,       # [(query_id, [LiftedFeature])]
        min_coverage=0.5,
        min_identity=0.3,
        evalue=1e-5,
        rescue_window=50,
    )
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


def _reset():
    for k in list(st.session_state.keys()):
        del st.session_state[k]
    _init_state()


# ═══════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════

def _save_upload(uploaded_file) -> Path:
    """Save a Streamlit UploadedFile to the session temp dir. Returns path."""
    tmp_dir = Path(st.session_state.tmp.name)
    dest = tmp_dir / uploaded_file.name
    dest.write_bytes(uploaded_file.read())
    return dest


def _scan_unknown_names(
    query_records: List[SeqRecord],
    alias_lookup: Dict[str, str],
    ignored_names: set,
) -> Dict[str, Dict]:
    """
    Return a dict keyed by representative name for every feature group in query
    records where NONE of the qualifier fields (_LOOKUP_QUALIFIER_KEYS) hit the
    alias lookup and none are explicitly ignored.

    Each value is:
        {
            "records":    [record_id, ...],       # records containing this feature
            "candidates": [val, ...],             # all qualifier values, priority order
        }

    The representative key is the highest-priority non-ignored candidate — used
    as the widget key in the resolver UI. All candidates are preserved so the
    user can see them and all get saved to alias config when a mapping is confirmed.

    Mirrors the logic of apply_alias_to_feature: a feature is only considered
    unknown when every candidate field misses the alias — not just gene/product.
    """
    unknown: Dict[str, Dict] = {}
    for rec in query_records:
        for feat in rec.features:
            if feat.type not in ("CDS", "mat_peptide"):
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

            # Skip if ANY candidate hits the alias (same logic as apply_alias_to_feature)
            if any(normalize_text(v) in alias_lookup for v in candidates):
                continue

            # Skip if ALL candidates are in ignored_names
            if all(v.lower() in ignored_names for v in candidates):
                continue

            # Representative = highest-priority non-ignored name (widget key)
            representative = next(
                (v for v in candidates if v.lower() not in ignored_names),
                candidates[0],
            )

            non_ignored_candidates = [
                v for v in candidates if v.lower() not in ignored_names
            ]

            if representative not in unknown:
                unknown[representative] = {"records": [], "candidates": non_ignored_candidates}
            if rec.id not in unknown[representative]["records"]:
                unknown[representative]["records"].append(rec.id)

    return unknown


def _scan_unknown_ref_names(ref_features: list, ignored_names: set) -> List[str]:
    """
    Return a sorted list of ref feature raw names that were not resolved by the
    alias DB (i.e. name_source == 'raw') and are not explicitly ignored.
    These will be lifted correctly (tblastn uses protein sequence), but their
    output name will just be the raw annotation name — no canonical key.
    """
    seen = []
    for f in ref_features:
        if f.get("name_source") != "raw":
            continue
        raw = f.get("raw_name") or f.get("name") or ""
        if raw.lower() in ignored_names:
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
    with open(alias_config_path) as f:
        cfg = json.load(f)
    added = 0
    config_name = alias_config_path.name
    for name in new_canonicals:
        if name not in cfg["canonical_names"]:
            cfg["canonical_names"][name] = []
            log_canonical_added(config_name, name)
            added += 1
    if added:
        with open(alias_config_path, "w") as f:
            json.dump(cfg, f, indent=2)
    return added


def _load_ignored_names(alias_config_path: Optional[Path]) -> set:
    if alias_config_path is None or not alias_config_path.exists():
        return set()
    with open(alias_config_path) as f:
        cfg = json.load(f)
    return {n.lower() for n in cfg.get("ignored_names", [])}


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
                "source_name": lf.source_name or "—",
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
    for i, qrec in enumerate(query_records):
        progress_bar.progress(i / n, text=f"Processing {qrec.id}  ({i+1}/{n})")
        strategy = get_strategy(qrec, ref_feature_type)
        try:
            if strategy == "direct":
                query_feature_type = get_feature_type(qrec)
                results = direct_extract_with_alias(
                    qrec, query_feature_type, ref_features, effective_lookup
                )
            else:
                results = process_one_query_record(
                    ref_record=ref_record,
                    query_record=qrec,
                    ref_cds=ref_features,
                    ref_feature_type=ref_feature_type,
                    min_coverage=min_coverage,
                    min_identity=min_identity,
                    evalue=evalue,
                    rescue_window=rescue_window,
                    quiet=True,
                )
        except Exception as e:
            log_error(f"processing record {qrec.id}", e)
            results = []

        all_results.append((qrec.id, results))

    progress_bar.progress(1.0, text="Done!")

    summary = summarize_counts(all_results)
    log_run_complete(ref_id=ref_record.id, n_queries=n, summary=summary)

    return all_results


# ═══════════════════════════════════════════════════════════════════
# Stage: UPLOAD
# ═══════════════════════════════════════════════════════════════════

def stage_upload():
    st.header("📂 Upload files")

    col_ref, col_query = st.columns(2)
    ref_file   = col_ref.file_uploader("Reference GenBank (.gb)", type=["gb", "gbk"])
    query_file = col_query.file_uploader("Query GenBank (.gb)", type=["gb", "gbk"])

    st.divider()
    st.subheader("⚙️ Advanced options")
    adv = st.expander("Lifting thresholds", expanded=False)
    with adv:
        c1, c2, c3, c4 = st.columns(4)
        st.session_state.min_coverage   = c1.number_input("Min coverage",  0.0, 1.0, 0.5, 0.05)
        st.session_state.min_identity   = c2.number_input("Min identity",  0.0, 1.0, 0.3, 0.05)
        st.session_state.evalue         = c3.number_input("E-value",       value=1e-5, format="%.0e")
        st.session_state.rescue_window  = c4.number_input("Rescue window", 10,  200,   50,    10)

    st.divider()
    st.session_state.use_ref_names = st.toggle(
        "Use ref gene names as output",
        value=False,
        help=(
            "OFF (default): output canonical names from the alias config key (e.g. 'Lpro').\n\n"
            "ON: output the ref's original gene name instead (e.g. 'Lab' if that's what the ref says)."
        ),
    )

    ready = ref_file and query_file
    if st.button("▶ Run ViraLift", disabled=not ready, type="primary", use_container_width=True):
        # persist files to temp dir
        if st.session_state.tmp is None:
            st.session_state.tmp = tempfile.TemporaryDirectory()

        with st.spinner("Loading files…"):
            ref_path   = _save_upload(ref_file)
            query_path = _save_upload(query_file)

            ref_record    = load_single_genbank(ref_path)
            query_records = load_genbank_records(query_path)

            ref_features, alias_config_path, virus_name, alias_lookup = (
                prepare_reference_features(
                    ref_record=ref_record,
                    alias_config_arg=None,
                    alias_registry_arg=str(REGISTRY_PATH),
                )
            )
            ref_feature_type = get_feature_type(ref_record)

            # canonical list for resolver dropdowns
            canonical_list = sorted(set(alias_lookup.values())) if alias_lookup else []

            # scan for unknown names in query records AND in ref
            ignored = _load_ignored_names(alias_config_path)
            unknown     = _scan_unknown_names(query_records, alias_lookup, ignored)
            unknown_ref = _scan_unknown_ref_names(ref_features, ignored)

        # store to session state
        st.session_state.ref_record        = ref_record
        st.session_state.query_records     = query_records
        st.session_state.ref_features      = ref_features
        st.session_state.ref_feature_type  = ref_feature_type
        st.session_state.alias_lookup      = alias_lookup
        st.session_state.alias_config_path = alias_config_path
        st.session_state.virus_name        = virus_name
        st.session_state.canonical_list    = canonical_list
        st.session_state.unknown_names     = unknown
        st.session_state.unknown_ref_names = unknown_ref

        if unknown or unknown_ref:
            st.session_state.stage = "resolve"
        else:
            st.session_state.stage = "running"

        st.rerun()

    if ref_file and query_file:
        st.caption(f"Ref: `{ref_file.name}`  |  Query: `{query_file.name}`")


# ═══════════════════════════════════════════════════════════════════
# Stage: RESOLVE
# ═══════════════════════════════════════════════════════════════════

def _save_to_alias_config(alias_config_path: Path, mappings: Dict[str, str]) -> int:
    """
    Persist user-confirmed mappings into the alias JSON config file.

    For each (raw_name -> canonical) pair, appends raw_name to the canonical's
    alias list if not already present.

    Returns the number of new aliases written.
    """
    if not alias_config_path or not alias_config_path.exists():
        return 0

    with open(alias_config_path) as f:
        cfg = json.load(f)

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

    with open(alias_config_path, "w") as f:
        json.dump(cfg, f, indent=2)

    return written


def stage_resolve():
    unknown     = st.session_state.unknown_names
    unknown_ref = st.session_state.unknown_ref_names
    canonicals  = st.session_state.canonical_list
    virus       = st.session_state.virus_name or "unknown virus"

    st.header("⚠️ Unrecognised gene names")

    # ── Ref-side warning ────────────────────────────────────────────────
    if unknown_ref:
        with st.expander(
            f"📋 {len(unknown_ref)} ref gene(s) not in alias DB — click to review",
            expanded=True,
        ):
            st.markdown(
                "These names from the **reference** were not found in the "
                f"**{virus}** alias config. Lifting still works (tblastn uses protein "
                "sequence, not the name), but they will appear in output with their "
                "**raw annotation name** instead of a canonical key.\n\n"
                "You can add them as new canonical entries now so they're recognised "
                "in future runs — or just continue as-is."
            )
            st.markdown("**Unrecognised ref names:** " +
                        ", ".join(f"`{n}`" for n in unknown_ref))

            ref_save_flags: Dict[str, bool] = {}
            for name in unknown_ref:
                ref_save_flags[name] = st.checkbox(
                    f"➕ Add **`{name}`** as a new canonical to alias config",
                    key=f"ref_add_{name}",
                    value=False,
                    help=(
                        f"Creates a new entry '{name}: []' in the alias config. "
                        "You can add aliases to it later by editing the JSON file."
                    ),
                )

            if st.button("💾 Save selected ref names to alias config",
                         key="save_ref_canonicals",
                         disabled=not any(ref_save_flags.values())):
                to_add = [n for n, checked in ref_save_flags.items() if checked]
                if to_add and st.session_state.alias_config_path:
                    added = _add_new_canonicals_to_config(
                        Path(st.session_state.alias_config_path), to_add
                    )
                    if added:
                        # also add them to the canonical_list for the query resolver below
                        st.session_state.canonical_list = sorted(
                            set(st.session_state.canonical_list) | set(to_add)
                        )
                        canonicals = st.session_state.canonical_list
                        st.toast(f"➕ {added} new canonical(s) added to alias config", icon="✅")
                    else:
                        st.info("All selected names already exist in the config.")

        st.divider()

    # ── Query-side resolver ──────────────────────────────────────────────
    if unknown:
        st.markdown(
            f"The query file contains **{len(unknown)} name(s)** not found in the "
            f"**{virus}** alias config. Decide what to do with each one before running."
        )
        st.divider()
    else:
        st.markdown(
            "✅ All query gene names are already in the alias config. "
            "No query-side decisions needed."
        )
        st.divider()

    decisions  = {}
    save_flags = {}
    options    = ["-- ignore (keep raw name) --"] + canonicals

    for rep, info in unknown.items():
        record_ids = info["records"]
        candidates = info["candidates"]

        col_name, col_action, col_save = st.columns([3, 3, 1])

        # Show all candidates so user has full context
        chips = " ".join(f"`{v}`" for v in candidates)
        col_name.markdown(chips)
        col_name.caption(f"Appears in: {', '.join(record_ids[:5])}"
                         + ("…" if len(record_ids) > 5 else ""))

        choice = col_action.selectbox(
            "Map to canonical →",
            options,
            key=f"resolve_{rep}",
            label_visibility="collapsed",
        )
        mapped = None if choice.startswith("--") else choice
        decisions[rep] = mapped

        # only offer save if user actually picked a canonical
        if mapped:
            save_flags[rep] = col_save.checkbox(
                "💾 Save", key=f"save_{rep}", value=True,
                help=(
                    "Add ALL names shown above to the alias config "
                    "so they're recognised next time"
                ),
            )
        else:
            col_save.write("")   # keep layout aligned

        st.divider()

    col_back, col_run = st.columns([1, 3])
    if col_back.button("← Back"):
        _reset()
        st.rerun()

    if col_run.button("▶ Continue with these decisions", type="primary", use_container_width=True):
        # Expand all candidates for each group into flat {candidate: canonical} dicts.
        # This ensures every variant name (product, note, etc.) is covered — both
        # for the session-only effective lookup and for permanent alias config saves.

        # Session resolver: all candidates of every decided group
        resolver_expanded: Dict[str, str] = {}
        for rep, canonical in decisions.items():
            if canonical:
                for candidate in unknown[rep]["candidates"]:
                    resolver_expanded[candidate] = canonical

        # Persist to alias config: only groups where 💾 Save was checked
        to_save: Dict[str, str] = {}
        for rep, canonical in decisions.items():
            if canonical and save_flags.get(rep, False):
                for candidate in unknown[rep]["candidates"]:
                    to_save[candidate] = canonical

        if to_save and st.session_state.alias_config_path:
            written = _save_to_alias_config(
                Path(st.session_state.alias_config_path), to_save
            )
            if written:
                st.toast(f"💾 {written} alias(es) saved to config", icon="✅")

        # log ALL decisions so there is always a trace
        if decisions:
            log_session_decisions(
                decisions=resolver_expanded,
                saved_names=list(to_save.keys()),
            )

        st.session_state.resolver = resolver_expanded
        st.session_state.stage    = "running"
        st.rerun()


# ═══════════════════════════════════════════════════════════════════
# Stage: RUNNING  (transient — immediately transitions to results)
# ═══════════════════════════════════════════════════════════════════

def stage_running():
    st.header("🔬 Running ViraLift…")
    progress = st.progress(0, text="Starting…")

    effective_lookup = _build_effective_lookup(
        st.session_state.alias_lookup,
        st.session_state.resolver,
    )

    all_results = _run_pipeline(
        ref_record        = st.session_state.ref_record,
        query_records     = st.session_state.query_records,
        ref_features      = st.session_state.ref_features,
        ref_feature_type  = st.session_state.ref_feature_type,
        effective_lookup  = effective_lookup,
        min_coverage      = st.session_state.min_coverage,
        min_identity      = st.session_state.min_identity,
        evalue            = st.session_state.evalue,
        rescue_window     = st.session_state.rescue_window,
        progress_bar      = progress,
        virus_name        = st.session_state.virus_name,
        alias_config_path = st.session_state.alias_config_path,
    )

    st.session_state.all_results = all_results
    st.session_state.stage       = "results"
    st.rerun()


# ═══════════════════════════════════════════════════════════════════
# Stage: RESULTS
# ═══════════════════════════════════════════════════════════════════

def stage_results():
    all_results = st.session_state.all_results
    summary     = summarize_counts(all_results)

    # build ref name map if user wants ref names as output
    # {canonical_key: ref_raw_name}  e.g. {"Lpro": "Lab", "3Cpro": "3C", ...}
    ref_name_map: Optional[Dict[str, str]] = (
        _canonical_to_ref_map(st.session_state.ref_features)
        if st.session_state.use_ref_names else None
    )

    df = _results_to_df(all_results, ref_name_map)

    # ── header + new run button ──────────────────────────────────
    h_col, btn_col = st.columns([5, 1])
    h_col.header("📊 Results")
    if st.session_state.use_ref_names:
        h_col.caption("📌 Names shown as ref gene names (toggle off in Upload to use canonical keys)")
    if btn_col.button("🔄 New run"):
        _reset()
        st.rerun()

    # ── summary badges ───────────────────────────────────────────
    cols = st.columns(6)
    badges = [
        ("🟢 OK",               summary["ok"]),
        ("🟡 Rescued",          summary["ok_rescued"]),
        ("🟠 Invalid boundary", summary["invalid_boundaries"]),
        ("🟠 Low coverage",     summary["low_coverage"]),
        ("🔴 No hit",           summary["no_hit"]),
        ("🔴 Translation fail", summary["translation_fail"]),
    ]
    for col, (label, count) in zip(cols, badges):
        col.metric(label, count)

    st.divider()

    # ── per-record expanders ──────────────────────────────────────
    st.subheader("🗂️ Per-record breakdown")
    for query_id, features in all_results:
        ok_count     = sum(1 for f in features if f.status in GOOD_STATUSES)
        total_count  = len(features)
        icon         = "✅" if ok_count == total_count else "⚠️"

        methods      = {f.method for f in features if f.method}
        if methods == {"direct"}:
            method_tag = "📋 direct"
        elif "tblastn" in methods and "direct" not in methods:
            method_tag = "🔬 tblastn"
        else:
            method_tag = "🔬 tblastn + 📋 direct"

        label = f"{icon} {query_id}   —   {ok_count}/{total_count} mapped   ·   {method_tag}"

        with st.expander(label, expanded=False):
            rec_rows = []
            for lf in features:
                display_name = (
                    ref_name_map.get(lf.name, lf.name)
                    if ref_name_map else lf.name
                )
                rec_rows.append({
                    "status":    STATUS_ICON.get(lf.status, "❓") + " " + lf.status,
                    "canonical": display_name,
                    "raw name":  lf.source_name if lf.method == "direct" else "—",
                    "start":     lf.query_start,
                    "end":       lf.query_end,
                    "coverage":  f"{lf.coverage:.0%}" if lf.coverage else "—",
                    "identity":  f"{lf.identity:.1f}%" if lf.identity else "—",
                    "method":    lf.method,
                })
            st.dataframe(pd.DataFrame(rec_rows), use_container_width=True, hide_index=True)

    st.divider()

    # ── export section ────────────────────────────────────────────
    st.subheader("⬇️ Export")
    tab_tsv, tab_fasta = st.tabs(["📄 TSV", "🧬 FASTA extraction"])

    # TSV tab
    with tab_tsv:
        st.markdown("Download the full results table.")
        col_canon, col_raw = st.columns(2)

        tsv_canonical = df.drop(columns=["sequence"]).to_csv(sep="\t", index=False)
        col_canon.download_button(
            "⬇ Download TSV (canonical names)",
            data=tsv_canonical,
            file_name="viralift_canonical.tsv",
            mime="text/tab-separated-values",
            use_container_width=True,
        )

        df_raw = df.copy()
        df_raw["canonical"] = df_raw["name"]   # overwrite canonical col with raw
        tsv_raw = df_raw.drop(columns=["sequence"]).to_csv(sep="\t", index=False)
        col_raw.download_button(
            "⬇ Download TSV (raw names)",
            data=tsv_raw,
            file_name="viralift_raw.tsv",
            mime="text/tab-separated-values",
            use_container_width=True,
        )

    # FASTA extraction tab
    with tab_fasta:
        st.markdown(
            "Select which genes to extract. Only features with status "
            "**ok** or **ok_rescued** are included."
        )

        # gene list — display names respect the use_ref_names toggle
        if ref_name_map:
            ref_gene_names = sorted({
                ref_name_map.get(f.get("name"), f.get("name"))
                for f in st.session_state.ref_features
                if f.get("name")
            })
        else:
            ref_gene_names = sorted({
                f.get("name")
                for f in st.session_state.ref_features
                if f.get("name")
            })

        col_select, col_format = st.columns([3, 2])
        selected_genes = col_select.multiselect(
            "Genes to extract",
            options=ref_gene_names,
            default=ref_gene_names,
        )
        fasta_mode = col_format.radio(
            "Output format",
            ["One FASTA per gene", "All genes in one FASTA"],
        )

        # quality filter
        st.markdown("**Quality filter**")
        qf_col1, qf_col2, qf_col3 = st.columns(3)
        min_cov_export = qf_col1.slider("Min coverage", 0.0, 1.0, 0.5, 0.05)
        min_id_export  = qf_col2.slider("Min identity",  0.0, 100.0, 0.0, 5.0)
        include_rescued = qf_col3.checkbox("Include ok_rescued", value=True)

        accepted_statuses = {"ok", "direct"}
        if include_rescued:
            accepted_statuses.add("ok_rescued")

        if st.button("⬇ Generate & Download FASTA", type="primary"):
            # build gene → sequences dict (keyed by display name)
            gene_seqs: Dict[str, List[str]] = {g: [] for g in selected_genes}
            skipped = 0

            for query_id, features in all_results:
                for lf in features:
                    # resolve display name the same way as the gene list above
                    gene = (
                        ref_name_map.get(lf.name, lf.name)
                        if ref_name_map else lf.name
                    )
                    if gene not in gene_seqs:
                        continue
                    if lf.status not in accepted_statuses:
                        skipped += 1
                        continue
                    if lf.coverage and lf.coverage < min_cov_export:
                        skipped += 1
                        continue
                    if lf.identity and lf.identity < min_id_export:
                        skipped += 1
                        continue
                    if not lf.sequence:
                        continue
                    header = f">{query_id}|{gene}|{lf.query_start}|{lf.query_end}|{lf.strand}"
                    gene_seqs[gene].append(f"{header}\n{lf.sequence}")

            if fasta_mode == "All genes in one FASTA":
                all_seqs = "\n".join(
                    seq for gene in selected_genes for seq in gene_seqs[gene]
                )
                st.download_button(
                    "⬇ Download all_genes.fasta",
                    data=all_seqs,
                    file_name="all_genes.fasta",
                    mime="text/plain",
                    use_container_width=True,
                )
            else:
                # one download button per gene
                for gene in selected_genes:
                    seqs = gene_seqs[gene]
                    if not seqs:
                        st.warning(f"No sequences for `{gene}` passed the filter.")
                        continue
                    st.download_button(
                        f"⬇ {gene}.fasta  ({len(seqs)} sequences)",
                        data="\n".join(seqs),
                        file_name=f"{gene}.fasta",
                        mime="text/plain",
                        key=f"dl_{gene}",
                    )

            if skipped:
                st.caption(f"ℹ️ {skipped} features skipped due to quality filter or missing sequence.")


# ═══════════════════════════════════════════════════════════════════
# App entry point
# ═══════════════════════════════════════════════════════════════════

st.set_page_config(
    page_title="ViraLift",
    page_icon="🧬",
    layout="wide",
)

_init_state()

# sidebar — always visible
with st.sidebar:
    st.title("🧬 ViraLift")
    st.caption("Reference-guided viral gene name standardisation")
    st.divider()

    if st.session_state.ref_record:
        st.markdown(f"**Reference**  \n`{st.session_state.ref_record.id}`")
    if st.session_state.query_records:
        st.markdown(f"**Query records**  \n{len(st.session_state.query_records)} loaded")
    if st.session_state.virus_name:
        st.markdown(f"**Detected virus**  \n{st.session_state.virus_name}")
    if st.session_state.alias_config_path:
        st.markdown(f"**Alias config**  \n`{Path(st.session_state.alias_config_path).name}`")

    st.divider()
    stage_labels = {"upload": "1 — Upload", "resolve": "2 — Resolve names",
                    "running": "3 — Running", "results": "4 — Results"}
    st.markdown(f"**Stage:** {stage_labels.get(st.session_state.stage, '?')}")

# route to current stage
stage = st.session_state.stage
if stage == "upload":
    stage_upload()
elif stage == "resolve":
    stage_resolve()
elif stage == "running":
    stage_running()
elif stage == "results":
    stage_results()
