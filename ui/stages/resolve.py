# Auto-split from the original monolithic streamlit_app.py.
import html
import streamlit as st
from pathlib import Path
from typing import Dict
from app.src.alias.gene_alias import normalize_text
from app.src.llm.alias_review import review_unresolved_names
from app.src.io.run_logger import log_session_decisions
from ui.components import _render_context_panel, _render_page_intro
from ui.i18n import _t
from ui.services import _add_new_canonicals_to_config, _load_excluded_names, _save_to_alias_config
from ui.state import _reset


def _llm_review_for_row(llm_reviews: Dict, rep: str, info: Dict) -> Dict:
    """Find an LLM review even when the representative differs from candidates."""
    if not llm_reviews:
        return {}
    lookup_values = [rep, normalize_text(rep)]
    for candidate in info.get("candidates", []) or []:
        text = str(candidate or "").strip()
        if text:
            lookup_values.extend([text, normalize_text(text)])
    for value in lookup_values:
        if value and value in llm_reviews:
            return llm_reviews[value]
    return {}


def stage_resolve():
    unknown     = st.session_state.unknown_names
    unknown_ref = st.session_state.unknown_ref_names
    canonicals  = st.session_state.canonical_list
    virus       = st.session_state.virus_name or "unknown virus"

    _render_page_intro(
        _t("alias_review"),
        _t("resolve_title"),
        _t("resolve_body"),
    )

    _render_context_panel([
        (_t("reference"), st.session_state.ref_record.id),
        (_t("query_records"), len(st.session_state.query_records)),
        (_t("detected_virus"), virus),
        (_t("alias_keys"), len(canonicals)),
    ])

    # ── Ref-side warning ────────────────────────────────────────────────
    if unknown_ref:
        with st.expander(
            _t("ref_missing_title", count=len(unknown_ref)),
            expanded=True,
        ):
            st.markdown(_t("ref_missing_body", virus=virus))
            st.markdown(f"**{_t('unrecognised_ref_names')}** " +
                        ", ".join(f"`{n}`" for n in unknown_ref))

            ref_save_flags: Dict[str, bool] = {}
            for name in unknown_ref:
                ref_save_flags[name] = st.checkbox(
                    _t("add_canonical", name=name),
                    key=f"ref_add_{name}",
                    value=False,
                    help=_t("add_canonical_help", name=name),
                )

            if st.button(_t("save_ref_names"),
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
                        st.toast(_t("canonicals_added", count=added))
                    else:
                        st.info(_t("already_exists"))

        st.divider()

    # ── Query-side resolver ──────────────────────────────────────────────
    unknown_items = dict(unknown)

    if unknown_items:
        st.markdown(_t("unknown_query_intro", count=len(unknown_items), virus=virus))
        st.divider()
    else:
        st.markdown(_t("all_known"))
        st.divider()

    if unknown_items:
        st.markdown("**Add new canonical target**")
        add_col, button_col = st.columns([3, 1])
        new_canonical = add_col.text_input(
            "New canonical name",
            value="",
            placeholder="e.g. ORF1ab",
            key="resolve_new_canonical_name",
            label_visibility="collapsed",
        ).strip()
        button_col.markdown("<div style='height: 1.78rem'></div>", unsafe_allow_html=True)
        if button_col.button(
            "Add canonical",
            disabled=not new_canonical,
            width="stretch",
            key="resolve_add_canonical",
        ):
            already_exists = new_canonical in st.session_state.canonical_list
            if st.session_state.alias_config_path:
                added = _add_new_canonicals_to_config(
                    Path(st.session_state.alias_config_path),
                    [new_canonical],
                )
            else:
                added = 0
            st.session_state.canonical_list = sorted(
                set(st.session_state.canonical_list) | {new_canonical}
            )
            if added and not already_exists:
                st.session_state.resolve_canonical_notice = (
                    "success",
                    f"Added `{new_canonical}` successfully.",
                )
            else:
                st.session_state.resolve_canonical_notice = (
                    "info",
                    f"`{new_canonical}` already exists.",
                )
            st.session_state.resolve_llm_reviews = {}
            st.session_state.resolve_llm_diagnostics = {}
            st.rerun()
        notice = st.session_state.pop("resolve_canonical_notice", None)
        if notice:
            kind, message = notice
            getattr(st, kind)(message)
        st.caption("Add targets like `ORF1ab` here before mapping unknown names below.")
        st.divider()
        canonicals = st.session_state.canonical_list

    llm_reviews = st.session_state.get("resolve_llm_reviews", {}) or {}
    llm_diagnostics = st.session_state.get("resolve_llm_diagnostics", {}) or {}

    if unknown_items:
        st.markdown("**LLM assist for unresolved names**")
        st.caption(
            "Use this for names that did not become coordinate-backed tblastn suggestions. "
            "Only raw names, candidate qualifier values, record IDs, and available canonicals are sent. "
            "The model pre-fills suggestions; you still approve each dropdown before continuing."
        )
        assist_cols = st.columns([1.2, 2])
        assist_cols[0].metric("Rows to review", len(unknown_items))
        if assist_cols[1].button(
            "Ask LLM to suggest mappings",
            type="primary",
            width="stretch",
            key="resolve_run_llm_review",
        ):
            with st.spinner("Reviewing unresolved names with LLM..."):
                excluded = (
                    _load_excluded_names(Path(st.session_state.alias_config_path))
                    if st.session_state.alias_config_path
                    else set()
                )
                reviews, diagnostics = review_unresolved_names(
                    unknown_items=unknown_items,
                    ambiguous_items={},
                    virus_name=virus,
                    canonical_names=canonicals,
                    excluded_names=excluded,
                    cache=st.session_state.llm_alias_review_cache,
                )
            st.session_state.resolve_llm_reviews = reviews
            st.session_state.resolve_llm_diagnostics = diagnostics
            ignore_option = _t("ignore_option")
            for representative, info in unknown_items.items():
                review = _llm_review_for_row(reviews, representative, info)
                if not review:
                    continue
                action = review.get("action")
                canonical = review.get("canonical_name")
                key = f"resolve_{representative}"
                if action == "save_alias" and canonical in canonicals:
                    st.session_state[key] = canonical
                elif action in {"ignore", "skip", "move_to_ambiguous"}:
                    st.session_state[key] = ignore_option
            st.rerun()

        if llm_diagnostics:
            status = llm_diagnostics.get("status")
            if status == "reviewed":
                cache_note = " cached" if llm_diagnostics.get("cache_hit") else ""
                st.success(
                    "LLM reviewed "
                    f"{llm_diagnostics.get('reviewed_rows', 0)} / "
                    f"{llm_diagnostics.get('submitted_rows', 0)} submitted row(s){cache_note}."
                )
            elif status == "missing_api_key":
                st.warning("LLM review is enabled but no API key is configured.")
            elif status == "disabled":
                st.info("LLM review is disabled. Set `VIRALIFT_LLM_ENABLED=true` to use this assist.")
            elif status == "error":
                st.error(f"LLM review failed: {llm_diagnostics.get('error')}")
            elif status:
                st.caption(f"LLM review status: {status}")
        st.divider()

    decisions = {}
    save_candidates = {}
    options    = [_t("ignore_option")] + canonicals

    def _default_save_candidate(candidate: str, representative: str) -> bool:
        """Persist only the clean representative by default."""
        text = str(candidate or "").strip()
        if not text or ";" in text:
            return False
        return text == representative

    def _render_resolver_row(rep: str, info: Dict) -> None:
        record_ids = info["records"]
        candidates = info["candidates"]

        col_name, col_action, col_save = st.columns([3, 2, 2])

        chips = " ".join(
            f"<span class='vl-pill'>{html.escape(str(v))}</span>"
            for v in candidates
        )
        col_name.markdown(chips, unsafe_allow_html=True)
        col_name.caption(
            _t("unknown_prefix") + " · "
            + f"{_t('appears_in')}: {', '.join(record_ids[:5])}"
            + ("..." if len(record_ids) > 5 else "")
        )

        choice = col_action.selectbox(
            _t("map_to_canonical"),
            options,
            key=f"resolve_{rep}",
            label_visibility="collapsed",
        )
        mapped = None if choice.startswith("--") else choice
        decisions[rep] = mapped

        review = _llm_review_for_row(llm_reviews, rep, info)
        if review:
            action = review.get("action")
            canonical = review.get("canonical_name") or "none"
            confidence = review.get("confidence") or "unknown"
            reason = review.get("reason") or "No reason returned."
            col_action.caption(
                f"LLM: `{action}` -> `{canonical}` ({confidence}). {reason}"
            )

        if mapped:
            col_save.caption(_t("save_aliases"))
            selected = []
            for idx, candidate in enumerate(candidates):
                checked = col_save.checkbox(
                    str(candidate),
                    key=f"save_{rep}_{idx}",
                    value=_default_save_candidate(str(candidate), rep),
                    help=_t("save_candidate_help"),
                )
                if checked:
                    selected.append(candidate)
            save_candidates[rep] = selected
        else:
            col_save.write("")

        st.divider()

    # Unknown names (completely unrecognised)
    for rep, info in unknown_items.items():
        _render_resolver_row(rep, info)

    col_back, col_run = st.columns([1, 3])
    if col_back.button(_t("back")):
        _reset()
        st.rerun()

    if col_run.button(_t("continue"), type="primary", width="stretch"):
        # Expand all candidates for each group into flat {candidate: canonical} dicts.
        # This ensures every variant name (product, note, etc.) is covered, both
        # for the session-only effective lookup and for permanent alias config saves.

        # Session resolver: all candidates of every decided group
        resolver_expanded: Dict[str, str] = {}
        for rep, canonical in decisions.items():
            if canonical:
                for candidate in unknown[rep]["candidates"]:
                    resolver_expanded[candidate] = canonical

        # Persist to alias config: only candidate aliases explicitly selected
        # by the user. Session resolution remains broader than permanent saves.
        to_save: Dict[str, str] = {}
        for rep, canonical in decisions.items():
            if canonical:
                for candidate in save_candidates.get(rep, []):
                    to_save[candidate] = canonical

        if to_save and st.session_state.alias_config_path:
            written = _save_to_alias_config(
                Path(st.session_state.alias_config_path), to_save
            )
            if written:
                st.toast(_t("aliases_saved", count=written))

        # log ALL decisions so there is always a trace
        if decisions:
            log_session_decisions(
                decisions=resolver_expanded,
                saved_names=list(to_save.keys()),
            )

        st.session_state.resolver = resolver_expanded
        st.session_state.stage    = "running"
        st.rerun()
