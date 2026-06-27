# Auto-split from the original monolithic streamlit_app.py.
import streamlit as st
from ui.components import _render_page_intro
from ui.i18n import _t
from ui.services import _build_effective_lookup, _run_pipeline


def stage_running():
    _render_page_intro(
        _t("processing"),
        _t("running_title"),
        _t("running_body"),
    )
    progress = st.progress(0, text=_t("starting"))
    st.session_state.run_errors = []

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
        run_errors        = st.session_state.run_errors,
    )

    st.session_state.all_results = all_results
    st.session_state.stage       = "results"
    st.rerun()
