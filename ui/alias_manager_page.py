# Alias-manager page for the Streamlit UI.
#
# This page intentionally avoids st.dataframe/st.data_editor. Those widgets are
# strongly tied to Streamlit's browser-level theme setting, so they can stay dark
# even when ViraLift's own theme is light. Checkbox lists and text areas are less
# fancy, but they obey the app CSS and are easier to review/edit.
#
# Save model (single, predictable rule): every edit persists immediately and writes
# a timestamped backup. Destructive actions (delete canonical / delete virus) ask
# for confirmation first. Success is surfaced through a flash toast that survives the
# rerun, so an action never looks like "nothing happened".

import json
from pathlib import Path
from typing import Dict, Iterable, List

import streamlit as st
import streamlit.components.v1 as components

from app.src.alias.alias_manager import (
    alias_config_to_tables,
    list_registry_entries,
    load_alias_config as manager_load_alias_config,
    move_excluded_to_alias,
    remove_registry_entry,
    resolve_config_path,
    save_alias_config as manager_save_alias_config,
    update_registry_entry,
    validate_alias_config,
)
from app.src.alias.gene_alias import get_excluded_names, normalize_text
from ui.components import _render_context_panel, _render_page_intro
from ui.state import REGISTRY_PATH, ROOT


# ── Flash feedback ──────────────────────────────────────────────────────
# st.success() written right before st.rerun() is discarded by the rerun, so the
# user never sees it. Stash the message and show it as a toast on the next run.

def _flash(message: str, icon: str = "✅") -> None:
    st.session_state["alias_manager_flash"] = (message, icon)


def _render_flash() -> None:
    payload = st.session_state.pop("alias_manager_flash", None)
    if payload:
        message, icon = payload
        st.toast(message, icon=icon)


def _scroll_to_top_once() -> None:
    if not st.session_state.pop("alias_manager_scroll_top", False):
        return
    components.html(
        """
        <script>
        const root = window.parent;
        root.scrollTo({ top: 0, behavior: "smooth" });
        </script>
        """,
        height=0,
    )


@st.dialog("Virus deleted")
def _delete_success_dialog(message: str) -> None:
    st.success(message)
    if st.button("OK", type="primary", width="stretch"):
        st.session_state.pop("alias_manager_delete_success", None)
        st.rerun()


@st.dialog("Delete virus?")
def _confirm_delete_virus_dialog(entry: Dict, archive_alias_config: bool) -> None:
    virus_name = entry.get("virus_name") or entry.get("alias_config") or "this virus"
    st.warning(f"Are you sure you want to delete `{virus_name}`?")
    st.caption(
        "This removes the virus from the Alias Manager registry. "
        "Auto-detection will no longer select it."
    )
    if archive_alias_config:
        st.caption("The active alias config file will also be archived and removed.")

    cancel_col, delete_col = st.columns(2)
    if cancel_col.button("Cancel", width="stretch"):
        st.rerun()
    if delete_col.button("Yes, delete virus", type="primary", width="stretch"):
        _delete_virus_entry(entry, archive_alias_config)


@st.dialog("Delete canonical?")
def _confirm_delete_canonical_dialog(config_path: Path, config: Dict, canonical: str) -> None:
    aliases = config.get("canonical_names", {}).get(canonical, []) or []
    st.warning(f"Delete canonical `{canonical}` and its {len(aliases)} alias(es)?")
    st.caption("This cannot be undone from the UI, but a backup is written first.")
    cancel_col, delete_col = st.columns(2)
    if cancel_col.button("Cancel", width="stretch"):
        st.rerun()
    if delete_col.button("Yes, delete canonical", type="primary", width="stretch"):
        _delete_canonical_now(config_path, config, canonical)


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


def _dedupe_values(values: Iterable[str]) -> List[str]:
    result: List[str] = []
    seen = set()
    for value in values:
        item = str(value or "").strip()
        key = normalize_text(item)
        if not key or key in seen:
            continue
        seen.add(key)
        result.append(item)
    return sorted(result, key=normalize_text)


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
    manager_save_alias_config(config_path, updated)
    _flash(f"Removed {len(aliases)} alias(es) from {canonical}.", icon="🗑️")
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

    manager_save_alias_config(config_path, updated)
    _flash(f"Added {len(added)} alias(es) to {canonical}.")
    st.rerun()


def _add_canonical(
    config_path: Path,
    config: Dict,
    canonical: str,
    aliases: List[str],
) -> None:
    canonical = canonical.strip()
    if not canonical:
        st.warning("Add a canonical name first.")
        return

    updated = _deepcopy_config(config)
    canonical_names = updated.setdefault("canonical_names", {})
    existing_key = next(
        (
            name
            for name in canonical_names
            if normalize_text(name) == normalize_text(canonical)
        ),
        None,
    )
    target = existing_key or canonical
    existing = canonical_names.setdefault(target, [])
    existing_norms = {normalize_text(alias) for alias in existing}

    for alias in aliases:
        key = normalize_text(alias)
        if key and key not in existing_norms:
            existing.append(alias)
            existing_norms.add(key)

    st.session_state.alias_manager_open_canonical = target
    manager_save_alias_config(config_path, updated)
    _flash(f"{'Updated' if existing_key else 'Added'} canonical {target}.")
    st.rerun()


def _delete_canonical_now(config_path: Path, config: Dict, canonical: str) -> None:
    updated = _deepcopy_config(config)
    updated.setdefault("canonical_names", {}).pop(canonical, None)
    manager_save_alias_config(config_path, updated)
    st.session_state.pop("alias_manager_open_canonical", None)
    _flash(f"Deleted canonical {canonical}.", icon="🗑️")
    st.rerun()


def _delete_names(config_path: Path, config: Dict, field: str, names: List[str], label: str) -> None:
    updated = _deepcopy_config(config)
    selected_norms = {normalize_text(name) for name in names}
    source_names = get_excluded_names(updated) if field == "excluded_names" else updated.get(field, [])
    updated[field] = [
        name
        for name in source_names
        if normalize_text(name) not in selected_norms
    ]
    if field == "excluded_names":
        updated.pop("ignored_names", None)
        updated.pop("ambiguous_names", None)
    manager_save_alias_config(config_path, updated)
    _flash(f"Removed {len(names)} {label}(s).", icon="🗑️")
    st.rerun()


def _add_excluded(config_path: Path, config: Dict, names: List[str]) -> None:
    updated = _deepcopy_config(config)
    existing = get_excluded_names(updated)
    existing_norms = {normalize_text(name) for name in existing}
    added = []
    for name in names:
        key = normalize_text(name)
        if key and key not in existing_norms:
            existing.append(name)
            existing_norms.add(key)
            added.append(name)

    if not added:
        st.info("No new excluded name to add.")
        return

    updated["excluded_names"] = sorted(existing, key=normalize_text)
    updated.pop("ignored_names", None)
    updated.pop("ambiguous_names", None)
    manager_save_alias_config(config_path, updated)
    _flash(f"Added {len(added)} excluded name(s).")
    st.rerun()


def _save_registry_entry(
    alias_config: str,
    virus_name: str,
    keywords: List[str],
    message: str,
) -> None:
    try:
        update_registry_entry(
            REGISTRY_PATH,
            alias_config=alias_config,
            virus_name=virus_name.strip(),
            keywords=keywords,
        )
        _flash(message)
        st.rerun()
    except Exception as e:
        st.error(f"Could not save registry entry: {e}")


def _delete_registry_keywords(entry: Dict, virus_name: str, keywords: List[str]) -> None:
    current_keywords = entry.get("keywords", [])
    selected_norms = {normalize_text(keyword) for keyword in keywords}
    kept_keywords = [
        keyword
        for keyword in current_keywords
        if normalize_text(keyword) not in selected_norms
    ]
    _save_registry_entry(
        alias_config=entry.get("alias_config", ""),
        virus_name=virus_name,
        keywords=kept_keywords,
        message=f"Removed {len(keywords)} keyword(s).",
    )


def _add_registry_keywords(entry: Dict, virus_name: str, keywords: List[str]) -> None:
    if not keywords:
        st.warning("Add at least one keyword first.")
        return

    existing = entry.get("keywords", [])
    merged = _dedupe_values([*existing, *keywords])
    added_count = len(merged) - len(_dedupe_values(existing))
    if added_count <= 0:
        st.info("No new keyword to add.")
        return

    _save_registry_entry(
        alias_config=entry.get("alias_config", ""),
        virus_name=virus_name,
        keywords=merged,
        message=f"Added {added_count} keyword(s).",
    )


def _delete_virus_entry(entry: Dict, archive_alias_config: bool) -> None:
    try:
        registry_backup, config_backup = remove_registry_entry(
            REGISTRY_PATH,
            alias_config=entry.get("alias_config", ""),
            root=ROOT,
            archive_alias_config=archive_alias_config,
        )
        if st.session_state.get("alias_manager_config_path"):
            del st.session_state["alias_manager_config_path"]
        message = f"Deleted virus from registry. Registry backup: `{registry_backup.name}`"
        if config_backup:
            message += f". Alias config archived as `{config_backup.name}`"
        elif archive_alias_config:
            message += ". Alias config file was already missing"
        st.session_state.alias_manager_delete_success = message
        st.session_state.alias_manager_scroll_top = True
        st.rerun()
    except Exception as e:
        st.error(f"Could not delete virus: {e}")


def _virus_picker(entries: List[Dict]) -> Dict:
    """Pick a virus by name (not by config file path)."""
    def _label(index: int) -> str:
        entry = entries[index]
        name = entry.get("virus_name") or Path(entry.get("alias_config", "")).stem or "unknown"
        keywords = entry.get("keywords", []) or []
        return f"{name}  ·  {len(keywords)} keyword(s)"

    selected_index = st.selectbox(
        "Virus",
        list(range(len(entries))),
        index=0,
        format_func=_label,
        help="Select the virus alias map you want to review or edit.",
    )
    entry = entries[selected_index]
    st.caption(f"Config file: `{Path(entry.get('alias_config', '')).name}`")
    return entry


def page_alias_manager():
    _scroll_to_top_once()
    _render_flash()
    _render_page_intro(
        "Alias manager",
        "Review and edit alias maps",
        (
            "Manage canonical names, aliases, and excluded names. "
            "Every edit saves immediately and writes a timestamped backup."
        ),
        show_stages=False,
    )

    if st.session_state.get("alias_manager_delete_success"):
        _delete_success_dialog(st.session_state.alias_manager_delete_success)

    try:
        entries = list_registry_entries(REGISTRY_PATH)
    except Exception as e:
        st.error(f"Could not load alias registry: {e}")
        return

    if not entries:
        st.warning("No registered virus alias configs found.")
        return

    entry = _virus_picker(entries)
    config_path = resolve_config_path(Path(entry.get("alias_config", "")), ROOT)
    st.session_state.alias_manager_config_path = config_path

    try:
        config = manager_load_alias_config(config_path)
    except Exception as e:
        st.error(f"Could not load alias config `{config_path}`: {e}")
        return

    _, excluded_rows, _ = alias_config_to_tables(config)
    excluded_existing = [row["excluded_name"] for row in excluded_rows]

    canonical_names_map = config.get("canonical_names", {})
    total_aliases = sum(len(aliases or []) for aliases in canonical_names_map.values())
    _render_context_panel([
        ("Virus", entry.get("virus_name", config.get("virus", "unknown"))),
        ("Canonical", len(canonical_names_map)),
        ("Aliases", total_aliases),
        ("Excluded", len(excluded_existing)),
    ])

    # Always-on health check so config issues are visible without hunting for a Save button.
    config_warnings = validate_alias_config(config)
    if config_warnings:
        with st.expander(f"⚠ {len(config_warnings)} config warning(s)", expanded=False):
            for warning in config_warnings:
                st.markdown(f"- {warning}")

    tab_registry, tab_alias, tab_excluded = st.tabs([
        "Registry",
        "Canonical aliases",
        "Excluded names",
    ])

    with tab_registry:
        st.caption("Auto-detection uses keywords. Virus name is the display label.")

        # Virus name with an inline, right-sized Save (not a heavy full-width bar).
        vn_input_col, vn_btn_col = st.columns([4, 1])
        registry_virus_name = vn_input_col.text_input(
            "Virus name",
            value=entry.get("virus_name", config.get("virus", "")),
            key=f"registry_virus_name_{config_path}",
        )
        vn_btn_col.markdown("<div style='height: 1.75rem'></div>", unsafe_allow_html=True)
        if vn_btn_col.button("Save", key=f"save_virus_name_{config_path}", width="stretch"):
            _save_registry_entry(
                alias_config=entry.get("alias_config", ""),
                virus_name=registry_virus_name,
                keywords=entry.get("keywords", []),
                message="Virus name saved.",
            )

        st.markdown("**Detection keywords**")
        st.caption("Tick keywords to remove them, or add one keyword per line below.")
        selected_keywords = _checkbox_list(
            entry.get("keywords", []),
            key_prefix=f"registry_keyword_select_{config_path}",
            empty_message="No detection keywords yet.",
        )
        if selected_keywords:
            st.caption(
                "Selected: " + ", ".join(f"`{keyword}`" for keyword in selected_keywords)
            )

        if st.button(
            f"Remove selected ({len(selected_keywords)})",
            disabled=not selected_keywords,
            key=f"delete_selected_registry_keywords_{config_path}",
        ):
            _delete_registry_keywords(entry, registry_virus_name, selected_keywords)

        with st.form(f"registry_add_keywords_form_{config_path}", clear_on_submit=True):
            new_keywords_text = st.text_area(
                "Add detection keywords",
                value="",
                height=80,
                placeholder="One keyword per line",
                key=f"registry_keywords_add_{config_path}",
            )
            add_keywords_submitted = st.form_submit_button(
                "Add keywords",
                type="primary",
            )

        if add_keywords_submitted:
            _add_registry_keywords(
                entry,
                registry_virus_name,
                _split_lines(new_keywords_text),
            )

        # Destructive action tucked away in a collapsed expander so it can't be
        # hit by accident and doesn't weigh down the tab.
        with st.expander("⚠ Danger zone — delete this virus", expanded=False):
            st.caption(
                "Delete this virus from the alias registry so ViraLift will no longer "
                "auto-detect or list it."
            )
            archive_alias_config = st.checkbox(
                "Also archive and remove the active alias config file",
                value=False,
                help=(
                    "Copies the alias JSON into app/config/backups/ before removing "
                    "the active config file."
                ),
                key=f"delete_virus_archive_config_{config_path}",
            )
            if st.button(
                "Delete virus from Alias Manager",
                key="delete_virus_entry",
            ):
                _confirm_delete_virus_dialog(entry, archive_alias_config)

    with tab_alias:
        st.caption(
            "Each canonical has its own editor. Ticking aliases and adding aliases "
            "save immediately; deleting a canonical asks for confirmation first."
        )
        search_text = st.text_input(
            "Search canonical or alias",
            value="",
            key=f"alias_manager_search_{config_path}",
        )

        with st.expander("➕ Add new canonical", expanded=False):
            with st.form(f"alias_manager_add_canonical_form_{config_path}", clear_on_submit=True):
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
                add_canonical_submitted = st.form_submit_button(
                    "Add canonical",
                    type="primary",
                )

            if add_canonical_submitted:
                _add_canonical(
                    config_path,
                    config,
                    new_canonical,
                    _split_lines(new_canonical_aliases),
                )

        canonical_names = config.get("canonical_names", {})
        search_norm = normalize_text(search_text)
        matches = 0
        for canonical in sorted(canonical_names, key=normalize_text):
            aliases = canonical_names.get(canonical, []) or []
            searchable = normalize_text(" ".join([canonical] + aliases))
            if search_norm and search_norm not in searchable:
                continue
            matches += 1

            should_expand = (
                bool(search_norm)
                or st.session_state.get("alias_manager_open_canonical") == canonical
            )
            with st.expander(f"{canonical} — {len(aliases)} alias(es)", expanded=should_expand):
                selected_aliases = _checkbox_list(
                    aliases,
                    key_prefix=f"alias_select_{config_path}_{normalize_text(canonical)}",
                    empty_message="No aliases yet.",
                )
                if selected_aliases:
                    st.caption(
                        "Selected: " + ", ".join(f"`{alias}`" for alias in selected_aliases)
                    )

                action_col, delete_col = st.columns([1, 1])
                if action_col.button(
                    f"Remove selected ({len(selected_aliases)})",
                    disabled=not selected_aliases,
                    key=f"delete_selected_aliases_{config_path}_{normalize_text(canonical)}",
                    width="stretch",
                ):
                    _delete_aliases(config_path, config, canonical, selected_aliases)

                if delete_col.button(
                    "Delete canonical",
                    key=f"delete_canonical_{config_path}_{normalize_text(canonical)}",
                    help="Removes this canonical and all its aliases (asks for confirmation).",
                    width="stretch",
                ):
                    _confirm_delete_canonical_dialog(config_path, config, canonical)

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
                    add_submitted = st.form_submit_button("Add aliases")

                if add_submitted:
                    new_aliases = _split_lines(new_alias_text)
                    if not new_aliases:
                        st.warning("Add at least one alias first.")
                        st.stop()
                    _add_aliases(config_path, config, canonical, new_aliases)

        if search_norm and matches == 0:
            st.caption("No canonical or alias matches your search.")

    with tab_excluded:
        st.caption(
            "Names here are skipped during automatic alias matching. Use this for "
            "generic, ambiguous, or unsafe names that should not become reusable aliases."
        )
        selected_excluded = _checkbox_list(
            excluded_existing,
            key_prefix=f"excluded_select_{config_path}",
            empty_message="No excluded names yet.",
        )
        if selected_excluded:
            st.caption("Selected: " + ", ".join(f"`{name}`" for name in selected_excluded))
        if st.button(
            f"Remove selected ({len(selected_excluded)})",
            disabled=not selected_excluded,
            key=f"delete_selected_excluded_{config_path}",
        ):
            _delete_names(config_path, config, "excluded_names", selected_excluded, "excluded name")

        with st.form(f"alias_manager_excluded_add_form_{config_path}", clear_on_submit=True):
            new_excluded_text = st.text_area(
                "Add excluded names",
                value="",
                height=120,
                placeholder="One excluded name per line",
                key=f"alias_manager_excluded_add_{config_path}",
            )
            add_excluded_submitted = st.form_submit_button("Add excluded names", type="primary")

        if add_excluded_submitted:
            names = _split_lines(new_excluded_text)
            if not names:
                st.warning("Add at least one excluded name first.")
                st.stop()
            _add_excluded(config_path, config, names)

        st.divider()
        st.markdown("**Move excluded name to alias**")
        canonical_options = sorted(config.get("canonical_names", {}))
        c1, c2, c3 = st.columns([2, 2, 1])
        excluded_choice = c1.selectbox(
            "Excluded name",
            excluded_existing,
            key=f"move_excluded_name_{config_path}",
        ) if excluded_existing else None
        canonical_choice = c2.selectbox(
            "Canonical",
            canonical_options,
            key=f"move_excluded_canonical_{config_path}",
        ) if canonical_options else None
        c3.markdown("<div style='height: 1.75rem'></div>", unsafe_allow_html=True)
        if c3.button(
            "Move",
            disabled=not (excluded_choice and canonical_choice),
            key=f"move_excluded_btn_{config_path}",
            width="stretch",
        ):
            updated = move_excluded_to_alias(config, excluded_choice, canonical_choice)
            warnings = validate_alias_config(updated)
            if warnings:
                st.warning("\n".join(warnings))
            else:
                manager_save_alias_config(config_path, updated)
                _flash(f"Moved {excluded_choice} to {canonical_choice}.")
                st.rerun()

    # Advanced: raw config view + export. Not a daily-use surface, so it's collapsed.
    with st.expander("Advanced — raw JSON & export", expanded=False):
        st.download_button(
            "Download alias config (JSON)",
            data=json.dumps(config, ensure_ascii=False, indent=2),
            file_name=config_path.name,
            mime="application/json",
        )
        if st.button("Reload from disk", key=f"alias_manager_reload_{config_path}"):
            st.rerun()
        st.json(config)
