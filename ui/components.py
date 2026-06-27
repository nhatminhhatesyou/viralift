# Auto-split from the original monolithic streamlit_app.py.
import html
import streamlit as st
from typing import Dict, List, Optional, Tuple
from ui.i18n import _t
from ui.state import STATUS_TONE


def _render_page_intro(kicker: str, title: str, body: str, show_stages: bool = True):
    stages = [
        ("upload", _t("stage_upload")),
        ("virus_review", _t("stage_virus_review")),
        ("bootstrap_alias", _t("stage_bootstrap")),
        ("resolve", _t("stage_resolve")),
        ("running", _t("stage_run")),
        ("results", _t("stage_review")),
    ]
    rail = ""
    if show_stages:
        rail = "".join(
            "<div class='vl-stage {active}'>{label}</div>".format(
                active="vl-stage-active" if st.session_state.stage == key else "",
                label=label,
            )
            for key, label in stages
        )
        rail = f"<div class='vl-stage-rail'>{rail}</div>"
    st.markdown(
        f"""
        <section class="vl-hero">
            <div class="vl-kicker">{html.escape(kicker)}</div>
            <h1>{html.escape(title)}</h1>
            <div class="vl-help">{html.escape(body)}</div>
            {rail}
        </section>
        """,
        unsafe_allow_html=True,
    )


def _render_context_panel(items: List[Tuple[str, object]]):
    cells = "".join(
        (
            "<div class='vl-context'>"
            "<div class='vl-context-label'>{label}</div>"
            "<div class='vl-context-value'>{value}</div>"
            "</div>"
        ).format(
            label=html.escape(str(label)),
            value=html.escape(str(value)),
        )
        for label, value in items
    )
    st.markdown(f"<div class='vl-context-grid'>{cells}</div>", unsafe_allow_html=True)


def _status_text(status: str) -> str:
    tone = _t(f"tone_{STATUS_TONE.get(status, 'review')}")
    return f"{tone} · {_t(f'status_{status}')}"


def _ordered_methods(methods: List[str]) -> List[str]:
    priority = {"direct": 0, "tblastn": 1}
    return sorted(
        {method or "unknown" for method in methods},
        key=lambda method: (priority.get(method, 99), method),
    )


def _method_section_label(method: str) -> str:
    labels = {
        "direct": "Direct extraction",
        "tblastn": "tblastn lifting",
    }
    return labels.get(method or "unknown", method or "Unknown method")


def _sidebar_item(label: str, value: object):
    st.markdown(
        """
        <div class="vl-sidebar-card">
            <div class="vl-sidebar-label">{label}</div>
            <div class="vl-sidebar-value">{value}</div>
        </div>
        """.format(
            label=html.escape(str(label)),
            value=html.escape(str(value)),
        ),
        unsafe_allow_html=True,
    )


def _format_percent(value: Optional[float], digits: int = 1) -> str:
    if value is None:
        return "NA"
    return f"{value:.{digits}f}%"


def _format_fraction_percent(value: Optional[float], digits: int = 0) -> str:
    if value is None:
        return "NA"
    return f"{value * 100:.{digits}f}%"


def _error_by_record(errors: List[Dict[str, str]]) -> Dict[str, str]:
    return {err["record_id"]: err["error"] for err in errors}
