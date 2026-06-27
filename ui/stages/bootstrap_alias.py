# Auto-split from the original monolithic streamlit_app.py.
import pandas as pd
import streamlit as st
from app.src.alias.alias_bootstrap import append_alias_registry_entry, apply_approved_alias_suggestions, build_coordinate_supported_alias_suggestions, build_seed_alias_config_from_ref, safe_alias_filename, write_new_alias_config
from app.src.alias.gene_alias import apply_alias_to_features, load_alias_lookup
from app.src.io.run_logger import log_error
from ui.components import _render_context_panel, _render_page_intro
from ui.i18n import _t
from ui.services import _load_ignored_names, _scan_unknown_names, _scan_unknown_ref_names, _split_keywords, _unique_alias_config_paths
from ui.state import REGISTRY_PATH, _reset


def stage_bootstrap_alias():
    seed_config = st.session_state.bootstrap_alias_config
    if not seed_config:
        st.warning("No bootstrap alias state is available. Start a new run.")
        if st.button(_t("new_run")):
            _reset()
            st.rerun()
        return

    _render_page_intro(
        "Alias seed",
        "Create alias config for a new virus",
        (
            "No existing alias config matched this reference. Seed canonical "
            "names from the reference, then use tblastn coordinate evidence to "
            "suggest which query names are safe aliases."
        ),
    )

    canonical_names = sorted(seed_config.get("canonical_names", {}))
    _render_context_panel([
        (_t("reference"), st.session_state.ref_record.id),
        (_t("query_records"), len(st.session_state.query_records)),
        ("Ref feature type", st.session_state.ref_feature_type),
        ("Seed canonicals", len(canonical_names)),
    ])

    st.subheader("1. Reference canonical names")
    st.caption(
        "For a new virus, ViraLift treats the reference feature names as the first canonical names."
    )
    st.dataframe(
        pd.DataFrame({"canonical_name": canonical_names}),
        width="stretch",
        hide_index=True,
    )

    st.subheader("2. Virus registry entry")
    default_name = st.session_state.bootstrap_virus_name or seed_config.get("virus") or "new virus"
    virus_name = st.text_input(
        "Virus name",
        value=default_name,
        help="Stored in the alias config and registry.",
    ).strip() or default_name
    keyword_default = st.session_state.bootstrap_keywords or virus_name
    keywords_raw = st.text_area(
        "Registry keywords",
        value=keyword_default,
        help="Comma or newline separated strings used to auto-detect this virus next time.",
    )

    st.subheader("3. Query-supported alias suggestions")
    st.caption(
        "Suggestions are generated only when a query annotation overlaps a tblastn-lifted reference feature with IoU >= 0.90. "
        "Duplicate aliases across records are merged into one review row with a support count."
    )

    gen_col, info_col = st.columns([1, 2])
    if gen_col.button("Generate suggestions", type="primary", width="stretch"):
        progress = st.progress(0.0, text="Starting alias suggestion scan...")

        def _update_bootstrap_progress(done: int, total: int, message: str) -> None:
            fraction = done / total if total else 1.0
            progress.progress(min(1.0, max(0.0, fraction)), text=message)

        with st.spinner("Running tblastn and matching query annotations..."):
            try:
                diagnostics = {}
                suggestions = build_coordinate_supported_alias_suggestions(
                    ref_record=st.session_state.ref_record,
                    query_records=st.session_state.query_records,
                    ref_features=st.session_state.ref_features,
                    ref_feature_type=st.session_state.ref_feature_type,
                    min_iou=0.90,
                    min_coverage=st.session_state.min_coverage,
                    min_identity=st.session_state.min_identity,
                    evalue=st.session_state.evalue,
                    rescue_window=st.session_state.rescue_window,
                    diagnostics=diagnostics,
                    progress_callback=_update_bootstrap_progress,
                )
            except Exception as e:
                log_error("bootstrap alias suggestions", e)
                st.error(f"Could not generate suggestions: {e}")
                suggestions = []
                diagnostics = {"error": str(e)}
                progress.empty()
            st.session_state.bootstrap_suggestions = suggestions
            st.session_state.bootstrap_diagnostics = diagnostics
            st.rerun()

    info_col.info(
        "Review each raw query name independently. Save strong gene symbols like `GP5`; "
        "ignore broad descriptions like `major envelope glycoprotein`. You review alias patterns once, not record by record."
    )

    diagnostics = st.session_state.bootstrap_diagnostics or {}
    if diagnostics:
        if diagnostics.get("error"):
            st.error(f"Suggestion diagnostics: {diagnostics['error']}")
        else:
            d_cols = st.columns(5)
            d_cols[0].metric("Records", diagnostics.get("total_records", 0))
            d_cols[1].metric("tblastn runs", diagnostics.get("records_tblastn_run", 0))
            d_cols[2].metric("Lifted features", diagnostics.get("lifted_features_total", 0))
            d_cols[3].metric("Matched features", diagnostics.get("matched_query_features_total", 0))
            d_cols[4].metric("Suggestions", diagnostics.get("deduplicated_suggestions_total", 0))

            skipped = diagnostics.get("records_without_usable_annotation", 0)
            no_names = diagnostics.get("records_without_name_candidates", 0)
            feature_counts = diagnostics.get("query_feature_type_counts", {})
            st.caption(
                "Diagnostics: "
                f"{skipped} record(s) skipped because no usable annotation was detected; "
                f"{no_names} record(s) had usable features but no name qualifiers; "
                f"query feature types: {feature_counts or '{}'}."
            )

    suggestions = st.session_state.bootstrap_suggestions or []
    edited_rows = pd.DataFrame()
    if suggestions:
        suggestion_df = pd.DataFrame(suggestions)
        suggestion_df["save"] = suggestion_df["default_save"].fillna(False).astype(bool)
        suggestion_df["ignore"] = suggestion_df["suggested_action"].eq("ignore")
        show_cols = [
            "save",
            "ignore",
            "raw_value",
            "field",
            "canonical_name",
            "suggested_action",
            "confidence",
            "score",
            "support_count",
            "support_records",
            "iou",
            "coverage",
            "identity",
            "reason",
            "record_id",
            "query_name",
            "query_start",
            "query_end",
            "tblastn_start",
            "tblastn_end",
        ]
        show_cols = [c for c in show_cols if c in suggestion_df.columns]
        edited_rows = st.data_editor(
            suggestion_df[show_cols],
            width="stretch",
            hide_index=True,
            disabled=[c for c in show_cols if c not in {"save", "ignore"}],
            column_config={
                "save": st.column_config.CheckboxColumn("Save alias"),
                "ignore": st.column_config.CheckboxColumn("Ignore"),
                "raw_value": st.column_config.TextColumn("Raw query name"),
                "canonical_name": st.column_config.TextColumn("Canonical"),
                "suggested_action": st.column_config.TextColumn("Suggestion"),
                "support_count": st.column_config.NumberColumn("Support", step=1),
                "support_records": st.column_config.TextColumn("Supporting records"),
                "iou": st.column_config.NumberColumn("IoU", format="%.3f"),
                "coverage": st.column_config.NumberColumn("Coverage", format="%.3f"),
                "identity": st.column_config.NumberColumn("Identity", format="%.3f"),
            },
            key="bootstrap_suggestion_editor",
        )
    else:
        st.info("No suggestions yet. You can still save a seed config with only reference canonical names.")

    st.divider()
    back_col, save_col = st.columns([1, 3])
    if back_col.button(_t("back")):
        _reset()
        st.rerun()

    if save_col.button("Save alias config and continue", type="primary", width="stretch"):
        config = build_seed_alias_config_from_ref(
            ref_record=st.session_state.ref_record,
            ref_features=st.session_state.ref_features,
            virus_name=virus_name,
        )
        config["notes"] = (
            "Bootstrapped from reference feature names. Query aliases were added "
            "only after coordinate-supported user approval."
        )

        filename = safe_alias_filename(virus_name)
        absolute_config_path, relative_config_path = _unique_alias_config_paths(filename)

        try:
            write_new_alias_config(config, absolute_config_path)

            if not edited_rows.empty:
                approved_rows = edited_rows[edited_rows["save"].fillna(False)].to_dict("records")
                ignored_rows = edited_rows[
                    edited_rows["ignore"].fillna(False) & ~edited_rows["save"].fillna(False)
                ].to_dict("records")
                config = apply_approved_alias_suggestions(
                    absolute_config_path,
                    approved_rows=approved_rows,
                    ignored_rows=ignored_rows,
                )

            append_alias_registry_entry(
                REGISTRY_PATH,
                virus_name=virus_name,
                keywords=_split_keywords(keywords_raw, virus_name),
                alias_config_path=relative_config_path,
            )

            alias_lookup = load_alias_lookup(absolute_config_path)
            ref_features = apply_alias_to_features(
                st.session_state.ref_features,
                alias_lookup,
            )
            ignored = _load_ignored_names(absolute_config_path)
            unknown = _scan_unknown_names(st.session_state.query_records, alias_lookup, ignored)
            unknown_ref = _scan_unknown_ref_names(ref_features, ignored)

            st.session_state.ref_features = ref_features
            st.session_state.alias_lookup = alias_lookup
            st.session_state.alias_config_path = absolute_config_path
            st.session_state.virus_name = virus_name
            st.session_state.canonical_list = sorted(config.get("canonical_names", {}))
            st.session_state.unknown_names = unknown
            st.session_state.unknown_ref_names = unknown_ref
            st.session_state.resolver = {}
            st.session_state.bootstrap_alias_config_path = absolute_config_path

            st.toast(f"Alias config saved: {relative_config_path}")
            st.session_state.stage = "resolve" if (unknown or unknown_ref) else "running"
            st.rerun()
        except Exception as e:
            log_error("saving bootstrap alias config", e)
            st.error(f"Could not save alias config: {e}")
