"""
ViraLift Streamlit Web UI — application entry point.

This file is intentionally thin: it wires together the ui/ package
(theme, i18n, state, components, services, stages, alias manager) and
routes to the active stage. All real logic lives in the submodules.

Stages:
    upload: user uploads ref + query, pipeline is configured
    bootstrap_alias: unknown virus, user creates a first alias config
    resolve: unmapped gene names found, user maps or ignores each
    results: pipeline ran, show results + export options
"""

import sys
from pathlib import Path

import streamlit as st

# project root on path (before importing ui.* which pull in app.src)
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ui.state import _init_state
from ui.theme import _inject_css
from ui.i18n import _t
from ui.components import _sidebar_item
from ui.stages.upload import stage_upload
from ui.stages.virus_review import stage_virus_review
from ui.stages.bootstrap_alias import stage_bootstrap_alias
from ui.stages.resolve import stage_resolve
from ui.stages.running import stage_running
from ui.stages.results import stage_results
from ui.alias_manager_page import page_alias_manager
from ui.data_crawler_page import page_data_crawler


APP_MODES = ("Run pipeline", "Data crawler", "Alias manager")


def _mode_key(label: str) -> str:
    return label.lower().replace(" ", "_")


def _set_app_mode(mode: str) -> None:
    _sync_crawler_credentials()
    st.session_state.app_mode = mode


def _sync_crawler_credentials() -> None:
    st.session_state.crawler_email_value = st.session_state.get("crawler_email_input", "")
    st.session_state.crawler_api_key_value = st.session_state.get("crawler_api_key_input", "")


# ═══════════════════════════════════════════════════════════════════
# App entry point
# ═══════════════════════════════════════════════════════════════════

st.set_page_config(
    page_title="ViraLift",
    page_icon="V",
    layout="wide",
)

_init_state()

# sidebar: always visible
with st.sidebar:
    st.title("ViraLift")
    st.caption(_t("sidebar_subtitle"))
    st.divider()

    st.markdown("**Workspace**")
    for mode in APP_MODES:
        st.button(
            mode,
            key=f"mode_{_mode_key(mode)}",
            type="primary" if st.session_state.app_mode == mode else "secondary",
            width="stretch",
            on_click=_set_app_mode,
            args=(mode,),
        )

    if st.session_state.app_mode == "Data crawler":
        st.divider()
        st.text_input(
            "NCBI email",
            value=st.session_state.get("crawler_email_value", ""),
            key="crawler_email_input",
            placeholder="you@lab.org",
            help="Required by NCBI Entrez requests.",
            on_change=_sync_crawler_credentials,
        )
        st.text_input(
            "NCBI API key",
            value=st.session_state.get("crawler_api_key_value", ""),
            key="crawler_api_key_input",
            type="password",
            help="Optional. Increases Entrez rate limit.",
            on_change=_sync_crawler_credentials,
        )
    st.divider()

_inject_css()

with st.sidebar:
    if st.session_state.ref_record:
        _sidebar_item(_t("reference"), st.session_state.ref_record.id)
    if st.session_state.query_records:
        _sidebar_item(_t("query_records"), _t("loaded", count=len(st.session_state.query_records)))
    if st.session_state.virus_name:
        _sidebar_item(_t("detected_virus"), st.session_state.virus_name)
    if st.session_state.alias_config_path:
        _sidebar_item(_t("alias_config"), Path(st.session_state.alias_config_path).name)

    if st.session_state.app_mode == "Run pipeline":
        st.divider()
        stage_labels = {
            "upload": _t("stage_upload"),
            "virus_review": _t("stage_virus_review"),
            "bootstrap_alias": _t("stage_bootstrap"),
            "resolve": _t("stage_resolve"),
            "running": _t("stage_run"),
            "results": _t("stage_review"),
        }
        _sidebar_item(_t("stage"), stage_labels.get(st.session_state.stage, "?"))

# route to current stage
if st.session_state.app_mode == "Alias manager":
    page_alias_manager()
elif st.session_state.app_mode == "Data crawler":
    page_data_crawler()
else:
    stage = st.session_state.stage
    if stage == "upload":
        stage_upload()
    elif stage == "virus_review":
        stage_virus_review()
    elif stage == "bootstrap_alias":
        stage_bootstrap_alias()
    elif stage == "resolve":
        stage_resolve()
    elif stage == "running":
        stage_running()
    elif stage == "results":
        stage_results()
