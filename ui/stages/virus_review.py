# Auto-split from the original monolithic streamlit_app.py.
import pandas as pd
import streamlit as st
from app.src.alias.alias_manager import list_registry_entries
from ui.components import _render_context_panel, _render_page_intro
from ui.i18n import _t
from ui.services import _continue_with_existing_alias_config
from ui.state import REGISTRY_PATH


def stage_virus_review():
    _render_page_intro(
        "Virus review",
        "Virus not recognised",
        (
            "The reference metadata did not match any registered virus keyword. "
            "Use an existing alias config if this is a known virus with new metadata, "
            "or create a new virus config."
        ),
    )

    metadata = st.session_state.virus_review_metadata or []
    _render_context_panel([
        (_t("reference"), st.session_state.ref_record.id),
        ("Suggested name", st.session_state.virus_review_guess),
        ("Metadata candidates", len(metadata)),
    ])
    if metadata:
        st.markdown("**Reference metadata found**")
        st.dataframe(pd.DataFrame({"candidate_keyword": metadata}), hide_index=True, width="stretch")

    entries = list_registry_entries(REGISTRY_PATH)
    choice = st.radio(
        "What do you want to do?",
        ["Use existing virus config", "Create new virus config"],
        horizontal=True,
    )

    if choice == "Use existing virus config":
        options = [
            f"{entry.get('virus_name', 'unknown')} — {entry.get('alias_config', '')}"
            for entry in entries
        ]
        selected = st.selectbox("Existing virus", options)
        entry = entries[options.index(selected)]

        keyword_options = [""] + metadata
        keyword_to_save = st.selectbox(
            "Save metadata as new keyword for this virus",
            keyword_options,
            help="Optional. Saving a keyword helps auto-detect this metadata next time.",
        )
        custom_keyword = st.text_input("Or custom keyword", value="")
        final_keyword = custom_keyword.strip() or keyword_to_save.strip()

        if st.button("Use this config", type="primary", width="stretch"):
            try:
                _continue_with_existing_alias_config(entry, save_keyword=final_keyword or None)
                st.rerun()
            except Exception as e:
                st.error(f"Could not continue with selected config: {e}")
    else:
        st.info("Continue to seed a new alias config from the reference feature names.")
        if st.button("Create new virus config", type="primary", width="stretch"):
            st.session_state.stage = "bootstrap_alias"
            st.rerun()
