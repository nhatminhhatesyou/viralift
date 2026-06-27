# Alias-manager page for the Streamlit UI.
#
# This page intentionally avoids st.dataframe/st.data_editor. Those widgets are
# strongly tied to Streamlit's browser-level theme setting, so they can stay dark
# even when ViraLift's own theme is light. Checkbox lists and text areas are less
# fancy, but they obey the app CSS and are easier to review/edit.

import json
from pathlib import Path
from typing import Dict, Iterable, List

import streamlit as st

from app.src.alias.alias_manager import (
    alias_config_to_tables,
    list_registry_entries,
    load_alias_config as manager_load_alias_config,
    move_ignored_to_alias,
    resolve_config_path,
    save_alias_config as manager_save_alias_config,
    tables_to_alias_config,
    update_registry_entry,
    validate_alias_config,
)
from app.src.alias.gene_alias import normalize_text
from ui.components import _render_context_panel, _render_page_intro
from ui.state import REGISTRY_PATH, ROOT


def _split_lines(raw_text: str) -> List[str]:
    values: List[str] = []
    seen = set()
    for line in (raw_text or "").splitlines():
        value = line.strip()
        key = normalize_text(value)
        if key and key not in seen:
            seen.add(key)
            values.append(value)
    return values


def _format_values(values: Iterable[str]) -> str:
    return "\n".join(str(value) for value in values if str(value).strip())


def _deepcopy_config(config: Dict) -> Dict:
    return json.loads(json.dumps(config))


def _checkbox_list(
    values: List[str],
    key_prefix: str,
    empty_message: str,
) -> List[str]:
    selected: List[str] = []
    if not values:
        st.caption(empty_message)
        return selected

    with st.container(border=True):
        for index, value in enumerate(values):
            checked = st.checkbox(
                value,
                value=False,
                key=f"{key_prefix}_{index}_{normalize_text(value)}",
            )
            if checked:
                selected.append(value)
    return selected


def _delete_aliases(config_path: Path, config: Dict, canonical: str, aliases: List[str]) -> None:
    st.session_state.alias_manager_open_canonical = canonical
    updated = _deepcopy_config(config)
    selected_norms = {normalize_text(alias) for alias in aliases}
    updated["canonical_names"][canonical] = [
        alias
        for alias in updated.get("canonical_names", {}).get(canonical, [])
        if normalize_text(alias) not in selected_norms
    ]
    backup = manager_save_alias_config(config_path, updated)
    st.success(
        f"Deleted {len(aliases)} alias(es) from `{canonical}`. "
        f"Backup: `{backup.name if backup else 'none'}`"
    )
    st.rerun()


def _add_aliases(config_path: Path, config: Dict, canonical: str, aliases: List[str]) -> None:
    st.session_state.alias_manager_open_canonical = canonical
    updated = _deepcopy_config(config)
    existing = updated.setdefault("canonical_names", {}).setdefault(canonical, [])
    existing_norms = {normalize_text(alias) for alias in existing}
    added = []
    for alias in aliases:
        key = normalize_text(alias)
        if key and key not in existing_norms:
            existing.append(alias)
            existing_norms.add(key)
            added.append(alias)

    if not added:
        st.info("No new alias to add.")
        return

    backup = manager_save_alias_config(config_path, updated)
    st.success(
        f"Added {len(added)} alias(es) to `{canonical}`. "
        f"Backup: `{backup.name if backup else 'none'}`"
    )
    st.rerun()


def _delete_names(config_path: Path, config: Dict, field: str, names: List[str], label: str) -> None:
    updated = _deepcopy_config(config)
    selected_norms = {normalize_text(name) for name in names}
    updated[field] = [
        name
        for name in updated.get(field, [])
        if normalize_text(name) not in selected_norms
    ]
    backup = manager_save_alias_config(config_path, updated)
    st.success(
        f"Deleted {len(names)} {label}(s). "
        f"Backup: `{backup.name if backup else 'none'}`"
    )
    st.rerun()


def page_alias_manager():
    _render_page_intro(
        "Alias manager",
        "Review and edit alias maps",
        (
            "Manage canonical names, aliases, ignored names, and ambiguous names. "
            "Every save creates a timestamped backup beside the config file."
        ),
        show_stages=False,
    )

    try:
        entries = list_registry_entries(REGISTRY_PATH)
    except Exception as e:
        st.error(f"Could not load alias registry: {e}")
        return

    if not entries:
        st.warning("No registered virus alias configs found.")
        return

    options = [
        f"{entry.get('virus_name', 'unknown')} - {entry.get('alias_config', '')}"
        for entry in entries
    ]
    selected = st.selectbox("Alias config", options, index=0)
    entry = entries[options.index(selected)]
    config_path = resolve_config_path(Path(entry.get("alias_config", "")), ROOT)
    st.session_state.alias_manager_config_path = config_path

    try:
        config = manager_load_alias_config(config_path)
    except Exception as e:
        st.error(f"Could not load alias config `{config_path}`: {e}")
        return

    _render_context_panel([
        ("Virus", entry.get("virus_name", config.get("virus", "unknown"))),
        ("Config", config_path),
        ("Canonical names", len(config.get("canonical_names", {}))),
    ])

    _, ignored_rows, ambiguous_rows = alias_config_to_tables(config)
    ignored_existing = [row["ignored_name"] for row in ignored_rows]
    ambiguous_existing = [row["ambiguous_name"] for row in ambiguous_rows]

    tab_registry, tab_alias, tab_ignored, tab_ambiguous, tab_raw = st.tabs([
        "Registry",
        "Canonical aliases",
        "Ignored names",
        "Ambiguous names",
        "Raw JSON",
    ])

    alias_records_for_save: List[Dict[str, str]] = []
    edited_ignored_records = [{"ignored_name": name} for name in ignored_existing]
    edited_ambiguous_records = [{"ambiguous_name": name} for name in ambiguous_existing]

    with tab_registry:
        st.caption("Auto-detection uses keywords. Virus name is the display label.")
        registry_virus_name = st.text_input(
            "Virus name",
            value=entry.get("virus_name", config.get("virus", "")),
            key=f"registry_virus_name_{config_path}",
        )
        keywords_text = st.text_area(
            "Detection keywords",
            value=_format_values(entry.get("keywords", [])),
            height=150,
            help="One keyword per line. These strings are used to auto-detect this virus.",
            key=f"registry_keywords_{config_path}",
        )
        if st.button("Save registry entry", type="primary", width="stretch"):
            try:
                backup = update_registry_entry(
                    REGISTRY_PATH,
                    alias_config=entry.get("alias_config", ""),
                    virus_name=registry_virus_name.strip(),
                    keywords=_split_lines(keywords_text),
                )
                st.success(f"Registry saved. Backup: `{backup.name if backup else 'none'}`")
                st.rerun()
            except Exception as e:
                st.error(f"Could not save registry entry: {e}")

    with tab_alias:
        st.caption(
            "Each canonical has its own editor. Tick aliases to delete them immediately, "
            "or add new aliases and click Save alias config at the bottom."
        )
        search_text = st.text_input(
            "Search canonical or alias",
            value="",
            key=f"alias_manager_search_{config_path}",
        )

        st.markdown("**Add new canonical**")
        add_col1, add_col2 = st.columns([2, 3])
        new_canonical = add_col1.text_input(
            "Canonical name",
            value="",
            key=f"alias_manager_new_canonical_{config_path}",
        ).strip()
        new_canonical_aliases = add_col2.text_area(
            "Initial aliases",
            value="",
            height=96,
            help="Optional. One alias per line.",
            key=f"alias_manager_new_canonical_aliases_{config_path}",
        )
        if new_canonical:
            alias_records_for_save.append({
                "canonical_name": new_canonical,
                "aliases": _format_values(_split_lines(new_canonical_aliases)),
            })

        canonical_names = config.get("canonical_names", {})
        search_norm = normalize_text(search_text)
        for canonical in sorted(canonical_names, key=normalize_text):
            aliases = canonical_names.get(canonical, []) or []
            searchable = normalize_text(" ".join([canonical] + aliases))
            if search_norm and search_norm not in searchable:
                continue

            should_expand = (
                bool(search_norm)
                or st.session_state.get("alias_manager_open_canonical") == canonical
            )
            with st.expander(f"{canonical} - {len(aliases)} alias(es)", expanded=should_expand):
                delete_canonical = st.checkbox(
                    f"Delete canonical `{canonical}`",
                    value=False,
                    key=f"delete_canonical_{config_path}_{normalize_text(canonical)}",
                    help="Removes this canonical and all aliases when you save.",
                )

                selected_aliases = _checkbox_list(
                    aliases,
                    key_prefix=f"alias_select_{config_path}_{normalize_text(canonical)}",
                    empty_message="No aliases yet.",
                )
                if selected_aliases:
                    st.caption(
                        "Selected: " + ", ".join(f"`{alias}`" for alias in selected_aliases)
                    )

                if st.button(
                    f"Delete selected ({len(selected_aliases)})",
                    disabled=not selected_aliases,
                    key=f"delete_selected_aliases_{config_path}_{normalize_text(canonical)}",
                ):
                    _delete_aliases(config_path, config, canonical, selected_aliases)

                if delete_canonical:
                    st.warning(f"`{canonical}` will be deleted on save.")
                    continue

                with st.form(
                    key=f"alias_add_form_{config_path}_{normalize_text(canonical)}",
                    clear_on_submit=True,
                ):
                    new_alias_text = st.text_area(
                        f"Add aliases for {canonical}",
                        value="",
                        height=96,
                        placeholder="One alias per line",
                        key=f"alias_add_{config_path}_{normalize_text(canonical)}",
                    )
                    add_submitted = st.form_submit_button(
                        "Add aliases",
                    )

                if add_submitted:
                    new_aliases = _split_lines(new_alias_text)
                    if not new_aliases:
                        st.warning("Add at least one alias first.")
                        st.stop()
                    _add_aliases(config_path, config, canonical, new_aliases)

                alias_records_for_save.append({
                    "canonical_name": canonical,
                    "aliases": _format_values(aliases),
                })

    with tab_ignored:
        st.caption("Names here are intentionally ignored by automatic alias matching.")
        selected_ignored = _checkbox_list(
            ignored_existing,
            key_prefix=f"ignored_select_{config_path}",
            empty_message="No ignored names yet.",
        )
        if selected_ignored:
            st.caption("Selected: " + ", ".join(f"`{name}`" for name in selected_ignored))
        if st.button(
            f"Delete selected ({len(selected_ignored)})",
            disabled=not selected_ignored,
            key=f"delete_selected_ignored_{config_path}",
        ):
            _delete_names(config_path, config, "ignored_names", selected_ignored, "ignored name")

        new_ignored_text = st.text_area(
            "Add ignored names",
            value="",
            height=120,
            placeholder="One ignored name per line",
            key=f"alias_manager_ignored_add_{config_path}",
        )
        edited_ignored_records = (
            [{"ignored_name": name} for name in ignored_existing]
            + [{"ignored_name": name} for name in _split_lines(new_ignored_text)]
        )

        st.divider()
        st.markdown("**Move ignored name to alias**")
        canonical_options = sorted(config.get("canonical_names", {}))
        c1, c2, c3 = st.columns([2, 2, 1])
        ignored_choice = c1.selectbox(
            "Ignored name",
            ignored_existing,
            key=f"move_ignored_name_{config_path}",
        ) if ignored_existing else None
        canonical_choice = c2.selectbox(
            "Canonical",
            canonical_options,
            key=f"move_ignored_canonical_{config_path}",
        ) if canonical_options else None
        if c3.button(
            "Move",
            disabled=not (ignored_choice and canonical_choice),
            key=f"move_ignored_btn_{config_path}",
        ):
            updated = move_ignored_to_alias(config, ignored_choice, canonical_choice)
            warnings = validate_alias_config(updated)
            if warnings:
                st.warning("\n".join(warnings))
            else:
                backup = manager_save_alias_config(config_path, updated)
                st.success(
                    f"Moved `{ignored_choice}` to `{canonical_choice}`. "
                    f"Backup: {backup.name if backup else 'none'}"
                )
                st.rerun()

    with tab_ambiguous:
        st.caption("Ambiguous names are known shared terms that require user decision per dataset.")
        selected_ambiguous = _checkbox_list(
            ambiguous_existing,
            key_prefix=f"ambiguous_select_{config_path}",
            empty_message="No ambiguous names yet.",
        )
        if selected_ambiguous:
            st.caption("Selected: " + ", ".join(f"`{name}`" for name in selected_ambiguous))
        if st.button(
            f"Delete selected ({len(selected_ambiguous)})",
            disabled=not selected_ambiguous,
            key=f"delete_selected_ambiguous_{config_path}",
        ):
            _delete_names(
                config_path,
                config,
                "ambiguous_names",
                selected_ambiguous,
                "ambiguous name",
            )

        new_ambiguous_text = st.text_area(
            "Add ambiguous names",
            value="",
            height=120,
            placeholder="One ambiguous name per line",
            key=f"alias_manager_ambiguous_add_{config_path}",
        )
        edited_ambiguous_records = (
            [{"ambiguous_name": name} for name in ambiguous_existing]
            + [{"ambiguous_name": name} for name in _split_lines(new_ambiguous_text)]
        )

    with tab_raw:
        st.json(config)

    st.divider()
    updated_config = tables_to_alias_config(
        config,
        alias_records_for_save,
        edited_ignored_records,
        edited_ambiguous_records,
    )
    warnings = validate_alias_config(updated_config)
    if warnings:
        st.warning("Review before saving:\n\n" + "\n".join(f"- {w}" for w in warnings))

    save_col, reload_col = st.columns([3, 1])
    if save_col.button(
        "Save alias config",
        type="primary",
        width="stretch",
        disabled=bool(warnings),
    ):
        try:
            backup = manager_save_alias_config(config_path, updated_config)
            st.success(f"Saved `{config_path.name}`. Backup created: `{backup.name if backup else 'none'}`")
            st.rerun()
        except Exception as e:
            st.error(f"Could not save alias config: {e}")

    if reload_col.button("Reload", width="stretch"):
        st.rerun()
