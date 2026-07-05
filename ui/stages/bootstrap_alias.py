# Auto-split from the original monolithic streamlit_app.py.
import pandas as pd
import streamlit as st
from app.src.alias.alias_bootstrap import append_alias_registry_entry, apply_approved_alias_suggestions, build_coordinate_supported_alias_suggestions, build_seed_alias_config_from_ref, safe_alias_filename, write_new_alias_config
from app.src.alias.gene_alias import apply_alias_to_features, load_alias_lookup
from app.src.io.run_logger import log_error
from app.src.llm.alias_review import review_uncertain_alias_suggestions
from ui.components import _render_context_panel, _render_page_intro
from ui.i18n import _t
from ui.services import _load_ignored_names, _scan_unknown_names, _scan_unknown_ref_names, _split_keywords, _unique_alias_config_paths
from ui.state import REGISTRY_PATH, _reset


def _suggestion_flag_conflicts(rows: pd.DataFrame) -> list[str]:
    conflicts = []
    if rows.empty:
        return conflicts
    for _, row in rows.iterrows():
        selected = [
            label
            for label, column in (("Save", "save"), ("Ignore", "ignore"), ("Skip", "skip"))
            if bool(row.get(column))
        ]
        if len(selected) > 1:
            conflicts.append(f"`{row.get('raw_value', '')}`: {', '.join(selected)}")
    return conflicts


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
    add_canonical_col, add_canonical_button_col = st.columns([3, 1])
    new_seed_canonical = add_canonical_col.text_input(
        "Add canonical target",
        value="",
        placeholder="e.g. ORF1ab",
        key="bootstrap_new_canonical",
    ).strip()
    add_canonical_button_col.markdown("<div style='height: 1.78rem'></div>", unsafe_allow_html=True)
    if add_canonical_button_col.button(
        "Add canonical",
        disabled=not new_seed_canonical,
        width="stretch",
        key="bootstrap_add_canonical",
    ):
        canonical_map = seed_config.setdefault("canonical_names", {})
        already_exists = new_seed_canonical in canonical_map
        canonical_map.setdefault(new_seed_canonical, [])
        st.session_state.bootstrap_alias_config = seed_config
        if already_exists:
            st.session_state.bootstrap_canonical_notice = (
                "info",
                f"`{new_seed_canonical}` already exists.",
            )
        else:
            st.session_state.bootstrap_canonical_notice = (
                "success",
                f"Added `{new_seed_canonical}` successfully.",
            )
        st.rerun()
    notice = st.session_state.pop("bootstrap_canonical_notice", None)
    if notice:
        kind, message = notice
        getattr(st, kind)(message)

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
                suggestions, llm_diagnostics = review_uncertain_alias_suggestions(
                    suggestions,
                    virus_name=virus_name,
                    canonical_names=canonical_names,
                    ignored_names=seed_config.get("ignored_names", []),
                    ambiguous_names=seed_config.get("ambiguous_names", []),
                    cache=st.session_state.llm_alias_review_cache,
                )
                diagnostics["llm_review"] = llm_diagnostics
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
            llm_diag = diagnostics.get("llm_review") or {}
            if llm_diag:
                if llm_diag.get("status") == "reviewed":
                    cache_note = " (cached)" if llm_diag.get("cache_hit") else ""
                    st.caption(
                        "LLM review: "
                        f"{llm_diag.get('reviewed_rows', 0)} uncertain alias row(s) reviewed"
                        f"{cache_note}."
                    )
                elif llm_diag.get("status") == "error":
                    st.warning(f"LLM review failed; deterministic suggestions were kept. {llm_diag.get('error')}")
                elif llm_diag.get("status") == "missing_api_key":
                    st.caption("LLM review skipped: set OPENAI_API_KEY to review uncertain aliases.")

    suggestions = st.session_state.bootstrap_suggestions or []
    edited_rows = pd.DataFrame()
    if suggestions:
        suggestion_df = pd.DataFrame(suggestions)
        suggestion_df["save"] = False
        suggestion_df["ignore"] = suggestion_df["suggested_action"].eq("ignore")
        suggestion_df["skip"] = False
        if "llm_reviewed" in suggestion_df.columns:
            llm_save = (
                suggestion_df["llm_reviewed"].fillna(False).astype(bool)
                & suggestion_df["llm_action"].eq("save_alias")
                & suggestion_df["llm_confidence"].isin(["medium", "high"])
                & suggestion_df["llm_canonical_name"].isin(canonical_names)
            )
            suggestion_df.loc[llm_save, "canonical_name"] = suggestion_df.loc[
                llm_save,
                "llm_canonical_name",
            ]
            suggestion_df.loc[llm_save, ["save", "ignore", "skip"]] = [True, False, False]

            llm_ignore = (
                suggestion_df["llm_reviewed"].fillna(False).astype(bool)
                & suggestion_df["llm_action"].eq("ignore")
                & suggestion_df["llm_confidence"].isin(["medium", "high"])
            )
            suggestion_df.loc[llm_ignore, ["save", "ignore", "skip"]] = [False, True, False]

            llm_skip = (
                suggestion_df["llm_reviewed"].fillna(False).astype(bool)
                & suggestion_df["llm_action"].isin(["skip", "move_to_ambiguous"])
            )
            suggestion_df.loc[llm_skip, ["save", "ignore", "skip"]] = [False, False, True]

        suggestion_df.loc[
            suggestion_df["default_save"].fillna(False).astype(bool)
            & ~(suggestion_df["save"] | suggestion_df["ignore"] | suggestion_df["skip"]),
            ["save", "ignore", "skip"],
        ] = [True, False, False]
        suggestion_df.loc[
            ~(suggestion_df["save"] | suggestion_df["ignore"]),
            "skip",
        ] = True
        show_cols = [
            "save",
            "ignore",
            "skip",
            "raw_value",
            "field",
            "canonical_name",
            "confidence",
            "support_count",
            "support_records",
            "reason",
        ]
        show_cols = [c for c in show_cols if c in suggestion_df.columns]

        editor_config = {
            "save": st.column_config.CheckboxColumn(
                "Save",
                help="Add this raw name as an alias.",
                width="small",
            ),
            "ignore": st.column_config.CheckboxColumn(
                "Ignore",
                help="Store this raw name as intentionally ignored.",
                width="small",
            ),
            "skip": st.column_config.CheckboxColumn(
                "Skip",
                help="Do nothing with this raw name.",
                width="small",
            ),
            "raw_value": st.column_config.TextColumn("Raw query name", width="medium"),
            "field": st.column_config.TextColumn("Field", width="small"),
            "canonical_name": st.column_config.SelectboxColumn(
                "Canonical",
                options=canonical_names,
                required=False,
                width="small",
            ),
            "confidence": st.column_config.TextColumn("Confidence", width="small"),
            "support_count": st.column_config.NumberColumn("Support", step=1, width="small"),
            "support_records": st.column_config.TextColumn("Supporting records", width="large"),
            "reason": st.column_config.TextColumn("Reason", width="large"),
            "llm_action": st.column_config.TextColumn("LLM action", width="small"),
            "llm_confidence": st.column_config.TextColumn("LLM confidence", width="small"),
            "llm_reason": st.column_config.TextColumn("LLM reason", width="large"),
        }
        editable_cols = {"save", "ignore", "skip", "canonical_name"}
        edited_frames = []
        llm_mask = (
            suggestion_df["llm_reviewed"].fillna(False).astype(bool)
            if "llm_reviewed" in suggestion_df.columns
            else pd.Series(False, index=suggestion_df.index)
        )
        deterministic_rows = suggestion_df[~llm_mask].copy()
        llm_rows = suggestion_df[llm_mask].copy()

        if not deterministic_rows.empty:
            st.markdown("**Deterministic suggestions**")
            st.caption(
                "Tick exactly one action per raw name. Edit Canonical when a good alias should point to a different target."
            )
            edited_frames.append(st.data_editor(
                deterministic_rows[show_cols],
                width="stretch",
                hide_index=True,
                height=min(520, 88 + 36 * len(deterministic_rows)),
                disabled=[c for c in show_cols if c not in editable_cols],
                column_config=editor_config,
                key="bootstrap_suggestion_editor_deterministic",
            ))

        if not llm_rows.empty:
            llm_show_cols = [
                "save",
                "ignore",
                "skip",
                "raw_value",
                "field",
                "canonical_name",
                "llm_action",
                "llm_confidence",
                "llm_reason",
                "confidence",
                "reason",
                "support_count",
                "support_records",
            ]
            llm_show_cols = [c for c in llm_show_cols if c in llm_rows.columns]
            with st.expander("LLM reviewed uncertain aliases", expanded=True):
                st.caption(
                    "These rows were separated because deterministic scoring was uncertain or risky. "
                    "Sequences and full GenBank records are not sent to the LLM."
                )
                edited_frames.append(st.data_editor(
                    llm_rows[llm_show_cols],
                    width="stretch",
                    hide_index=True,
                    height=min(420, 88 + 36 * len(llm_rows)),
                    disabled=[c for c in llm_show_cols if c not in editable_cols],
                    column_config=editor_config,
                    key="bootstrap_suggestion_editor_llm",
                ))

        if edited_frames:
            edited_rows = pd.concat(edited_frames, ignore_index=True)
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
            conflicts = _suggestion_flag_conflicts(edited_rows)
            if conflicts:
                st.error(
                    "Each suggestion can have only one action. Fix these rows: "
                    + "; ".join(conflicts[:8])
                    + ("..." if len(conflicts) > 8 else "")
                )
                return

            write_new_alias_config(config, absolute_config_path)

            if not edited_rows.empty:
                approved_rows = edited_rows[
                    edited_rows["save"].fillna(False)
                ].to_dict("records")
                ignored_rows = edited_rows[
                    edited_rows["ignore"].fillna(False)
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
