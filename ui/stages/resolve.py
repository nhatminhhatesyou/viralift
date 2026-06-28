# Auto-split from the original monolithic streamlit_app.py.
import html
import streamlit as st
from pathlib import Path
from typing import Dict
from app.src.io.run_logger import log_session_decisions
from ui.components import _render_context_panel, _render_page_intro
from ui.i18n import _t
from ui.services import _add_new_canonicals_to_config, _save_to_alias_config
from ui.state import _reset


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
    unknown_items   = {k: v for k, v in unknown.items() if not v.get("ambiguous")}
    ambiguous_items = {k: v for k, v in unknown.items() if v.get("ambiguous")}

    if unknown_items:
        st.markdown(_t("unknown_query_intro", count=len(unknown_items), virus=virus))
        st.divider()
    elif ambiguous_items:
        pass  # header shown below in ambiguous section
    else:
        st.markdown(_t("all_known"))
        st.divider()

    if unknown_items or ambiguous_items:
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
            st.rerun()
        notice = st.session_state.pop("resolve_canonical_notice", None)
        if notice:
            kind, message = notice
            getattr(st, kind)(message)
        st.caption("Add targets like `ORF1ab` here before mapping unknown names below.")
        st.divider()
        canonicals = st.session_state.canonical_list

    decisions  = {}
    save_flags = {}
    options    = [_t("ignore_option")] + canonicals

    def _render_resolver_row(rep: str, info: Dict, is_ambiguous: bool) -> None:
        record_ids = info["records"]
        candidates = info["candidates"]

        col_name, col_action, col_save = st.columns([3, 3, 1])

        chips = " ".join(
            f"<span class='vl-pill'>{html.escape(str(v))}</span>"
            for v in candidates
        )
        col_name.markdown(chips, unsafe_allow_html=True)
        col_name.caption(
            (_t("ambiguous_prefix") + " · " if is_ambiguous else _t("unknown_prefix") + " · ")
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

        if mapped:
            save_flags[rep] = col_save.checkbox(
                _t("save"), key=f"save_{rep}", value=True,
                help=_t("save_help"),
            )
        else:
            col_save.write("")

        st.divider()

    # Unknown names (completely unrecognised)
    for rep, info in unknown_items.items():
        _render_resolver_row(rep, info, is_ambiguous=False)

    # Ambiguous names (known to map to multiple genes, user must pick which one)
    if ambiguous_items:
        st.markdown(_t("ambiguous_intro", count=len(ambiguous_items), virus=virus))
        st.divider()
        for rep, info in ambiguous_items.items():
            _render_resolver_row(rep, info, is_ambiguous=True)

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

        # Persist to alias config: only groups where Save was checked
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
