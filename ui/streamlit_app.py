"""
ViraLift Streamlit Web UI

Stages:
    upload: user uploads ref + query, pipeline is configured
    bootstrap_alias: unknown virus, user creates a first alias config
    resolve: unmapped gene names found, user maps or ignores each
    results: pipeline ran, show results + export options
"""

import sys
import json
import html
import re
import tempfile
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import streamlit as st
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ── project root on path ─────────────────────────────────────────────
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from app.src.alias.alias_registry import (
    detect_alias_config_for_record,
    get_detected_virus_name,
)
from app.src.alias.alias_bootstrap import (
    append_alias_registry_entry,
    apply_approved_alias_suggestions,
    build_coordinate_supported_alias_suggestions,
    build_seed_alias_config_from_ref,
    safe_alias_filename,
    write_new_alias_config,
)
from app.src.alias.alias_manager import (
    add_registry_keyword,
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
from app.src.features.annotation_strategy import get_strategy, select_feature_type
from app.src.alias.gene_alias import (
    apply_alias_to_features,
    load_alias_lookup,
    normalize_text,
    AMBIGUOUS_SENTINEL,
    IGNORED_SENTINEL,
)
from app.src.io.genbank_parser import (
    load_single_genbank,
    load_genbank_records,
    parse_cds_features,
    parse_mat_peptides,
    _LOOKUP_QUALIFIER_KEYS,
)
from app.src.features.direct_extractor import direct_extract_with_alias
from app.src.features.ref_loader import prepare_reference_features
from app.src.io.result_writer import summarize_counts
from app.src.lifting.tblastn_lifter import process_one_query_record
from app.src.lifting.base import LiftedFeature
from app.src.io.run_logger import (
    log_alias_added,
    log_canonical_added,
    log_session_decisions,
    log_run_start,
    log_run_complete,
    log_error,
    log_warning,
)

# ── constants ────────────────────────────────────────────────────────
REGISTRY_PATH  = ROOT / "app/config/virus_alias_registry.json"
CONFIG_DIR = ROOT / "app/config"
GOOD_STATUSES = {"ok", "ok_rescued", "direct"}
REVIEW_STATUSES = {
    "invalid_boundaries",
    "low_coverage",
    "no_hit",
    "translation_fail",
    "unresolved_name",
    "ambiguous_name",
    "not_in_reference",
}
STATUS_LABEL = {
    "ok": "OK",
    "ok_rescued": "Rescued",
    "direct": "Direct",
    "invalid_boundaries": "Invalid boundary",
    "low_coverage": "Low coverage",
    "no_hit": "No hit",
    "translation_fail": "Translation fail",
    "unresolved_name": "Unresolved name",
    "ambiguous_name": "Ambiguous name",
    "not_in_reference": "Not in reference",
}
STATUS_TONE = {
    "ok": "pass",
    "ok_rescued": "pass",
    "direct": "pass",
    "invalid_boundaries": "review",
    "low_coverage": "review",
    "no_hit": "fail",
    "translation_fail": "fail",
    "unresolved_name": "review",
    "ambiguous_name": "review",
    "not_in_reference": "review",
}

UI_TEXT = {
    "en": {
        "theme": "Theme",
        "light": "Light",
        "dark": "Dark",
        "run_setup": "Run setup",
        "upload_title": "Upload GenBank files",
        "upload_body": (
            "Use one well-annotated reference record and one query file containing "
            "one or more genomes. ViraLift will detect the virus alias config, "
            "choose direct extraction or tblastn per record, then ask only for "
            "gene-name decisions it cannot resolve."
        ),
        "reference_record": "Reference record",
        "reference_caption": "Single annotated GenBank file",
        "query_records": "Query records",
        "query_caption": "One or many genomes in GenBank format",
        "file_requirements": "File requirements",
        "file_requirements_body": (
            "- Reference file: exactly one GenBank record with usable CDS or mat_peptide annotations.\n"
            "- Query file: one or more GenBank records; annotations are optional.\n"
            "- Supported extensions: `.gb` and `.gbk`.\n"
            "- Alias config is auto-detected from `app/config/virus_alias_registry.json`."
        ),
        "advanced_options": "Advanced options",
        "lifting_thresholds": "Lifting thresholds",
        "threshold_caption": "Defaults are permissive for divergent viral genomes. Tighten them when you want higher-confidence extraction.",
        "min_coverage": "Min coverage",
        "min_identity": "Min identity",
        "evalue": "E-value",
        "rescue_window": "Rescue window",
        "use_ref_names": "Use ref gene names as output",
        "use_ref_names_help": (
            "OFF (default): output canonical names from the alias config key (e.g. 'Lpro').\n\n"
            "ON: output the ref's original gene name instead (e.g. 'Lab' if that's what the ref says)."
        ),
        "run_button": "Run ViraLift",
        "loading_files": "Loading files...",
        "load_error": "We could not load these files",
        "reference_file": "Reference file",
        "query_file": "Query file",
        "name_mode": "Name mode",
        "reference_names": "reference names",
        "canonical_names": "canonical names",
        "alias_review": "Alias review",
        "resolve_title": "Resolve gene names",
        "resolve_body": (
            "ViraLift found names that are not cleanly covered by the alias database. "
            "Map only the names you trust; ignored names stay raw for this run."
        ),
        "reference": "Reference",
        "detected_virus": "Detected virus",
        "alias_keys": "Alias keys",
        "ref_missing_title": "{count} ref gene(s) not in alias DB. Review",
        "ref_missing_body": (
            "These names from the **reference** were not found in the **{virus}** alias config. "
            "Lifting still works (tblastn uses protein sequence, not the name), but they will appear "
            "in output with their **raw annotation name** instead of a canonical key.\n\n"
            "You can add them as new canonical entries now so they're recognised in future runs, or just continue as-is."
        ),
        "unrecognised_ref_names": "Unrecognised ref names:",
        "add_canonical": "Add **`{name}`** as a new canonical to alias config",
        "add_canonical_help": "Creates a new entry '{name}: []' in the alias config. You can add aliases to it later by editing the JSON file.",
        "save_ref_names": "Save selected ref names to alias config",
        "canonicals_added": "{count} new canonical(s) added to alias config",
        "already_exists": "All selected names already exist in the config.",
        "unknown_query_intro": "The query file contains **{count} unrecognised name(s)** not found in the **{virus}** alias config. Decide what to do with each one before running.",
        "all_known": "All query gene names are already in the alias config. No query-side decisions needed.",
        "ignore_option": "-- ignore (keep raw name) --",
        "ambiguous_prefix": "Ambiguous",
        "unknown_prefix": "Unknown",
        "appears_in": "Appears in",
        "map_to_canonical": "Map to canonical",
        "save": "Save",
        "save_help": "Add ALL names shown above to the alias config so they're recognised next time",
        "ambiguous_intro": "**{count} ambiguous name(s)**. These names appear in the alias config but are shared across multiple genes in **{virus}**. Select which gene each one refers to in this dataset.",
        "back": "Back",
        "continue": "Continue with these decisions",
        "aliases_saved": "{count} alias(es) saved to config",
        "processing": "Processing",
        "running_title": "Running ViraLift",
        "running_body": "Each query record is routed independently. Records that fail are logged and shown in the results review.",
        "starting": "Starting...",
        "done": "Done.",
        "run_review": "Run review",
        "results_title": "Results",
        "results_body": "Review run health first, then inspect records that need attention before exporting TSV or FASTA.",
        "ref_names_caption": "Names shown as ref gene names. Start a new run to use canonical keys.",
        "new_run": "New run",
        "records_processed": "Records processed",
        "features_found": "Features found",
        "pass_rate": "Pass rate",
        "needs_review": "Needs review",
        "failed_delta": "{count} failed",
        "processing_error": "{count} query record(s) failed during processing. They are excluded from pass-rate calculations and listed below.",
        "processing_errors": "Processing errors",
        "status_breakdown": "Status breakdown",
        "record_overview": "Record overview",
        "health_filter": "Health filter",
        "all": "all",
        "yes": "yes",
        "no": "no",
        "search_record": "Search record ID",
        "open_all": "Open all details",
        "record_details": "Record details",
        "failed": "Failed",
        "empty": "Empty",
        "passed": "Passed",
        "no_features": "No features were returned for this record.",
        "export": "Export",
        "tsv": "TSV",
        "fasta": "FASTA extraction",
        "download_table": "Download the full results table.",
        "download_display": "Download TSV (display names)",
        "download_raw": "Download TSV (raw/source names)",
        "fasta_intro": "Select genes and quality filters before generating FASTA downloads.",
        "genes_to_extract": "Genes to extract",
        "output_format": "Output format",
        "one_fasta": "One FASTA per gene",
        "all_fasta": "All genes in one FASTA",
        "quality_filter": "Quality filter",
        "include_rescued": "Include ok_rescued",
        "candidate_count": "{count} sequence(s) currently pass these FASTA filters.",
        "generate_fasta": "Generate FASTA downloads",
        "download_all_fasta": "Download all_genes.fasta",
        "no_sequences": "No sequences for `{gene}` passed the filter.",
        "gene_download": "{gene}.fasta  ({count} sequences)",
        "skipped": "{count} features skipped due to quality filter or missing sequence.",
        "stage_upload": "1 Upload",
        "stage_virus_review": "2 Virus",
        "stage_bootstrap": "3 Alias seed",
        "stage_resolve": "4 Resolve",
        "stage_run": "5 Run",
        "stage_review": "6 Review",
        "sidebar_subtitle": "Reference-guided viral gene name standardisation",
        "alias_config": "Alias config",
        "stage": "Stage",
        "loaded": "{count} loaded",
        "status_ok": "OK",
        "status_ok_rescued": "Rescued",
        "status_direct": "Direct",
        "status_invalid_boundaries": "Invalid boundary",
        "status_low_coverage": "Low coverage",
        "status_no_hit": "No hit",
        "status_translation_fail": "Translation fail",
        "status_unresolved_name": "Unresolved name",
        "status_ambiguous_name": "Ambiguous name",
        "status_not_in_reference": "Not in reference",
        "tone_pass": "PASS",
        "tone_review": "REVIEW",
        "tone_fail": "FAIL",
    },
}


# ═══════════════════════════════════════════════════════════════════
# Session-state bootstrap
# ═══════════════════════════════════════════════════════════════════

def _init_state():
    defaults = dict(
        stage="upload",
        app_mode="Run pipeline",
        tmp=None,               # tempfile.TemporaryDirectory object
        ref_record=None,
        query_records=None,
        ref_features=None,
        use_ref_names=False,    # if True, output uses ref's raw names instead of alias keys
        ref_feature_type=None,
        alias_lookup={},
        alias_config_path=None,
        virus_name=None,
        canonical_list=[],
        unknown_ref_names=[],   # [raw_name] ref features not in alias DB
        unknown_names={},       # {raw_name: [record_id, ...]}
        resolver={},            # {raw_name: canonical or "-- ignore --"}
        all_results=None,       # [(query_id, [LiftedFeature])]
        run_errors=[],          # [{record_id, error}]
        min_coverage=0.5,
        min_identity=0.3,
        evalue=1e-5,
        rescue_window=50,
        ui_theme="dark",
        bootstrap_alias_config=None,
        bootstrap_alias_config_path=None,
        bootstrap_suggestions=[],
        bootstrap_diagnostics={},
        bootstrap_virus_name="",
        bootstrap_keywords="",
        alias_manager_config_path=None,
        virus_review_metadata=[],
        virus_review_guess="",
    )
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


def _reset():
    ui_theme = st.session_state.get("ui_theme", "dark")
    for k in list(st.session_state.keys()):
        del st.session_state[k]
    _init_state()
    st.session_state.ui_theme = ui_theme


# ═══════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════

def _t(key: str, **kwargs) -> str:
    template = UI_TEXT["en"].get(key, key)
    return template.format(**kwargs) if kwargs else template


def _theme_overrides() -> str:
    if st.session_state.get("ui_theme", "dark") == "light":
        return """
        <style>
            :root {
                --vl-bg: #f4f7f3;
                --vl-bg-2: #e8efe8;
                --vl-surface: #fbfcf8;
                --vl-surface-strong: #ffffff;
                --vl-surface-muted: #eef4ee;
                --vl-border: #d6ded3;
                --vl-border-strong: #bac9b7;
                --vl-text: #17221d;
                --vl-muted: #637269;
                --vl-faint: #859188;
                --vl-accent: #1d6c63;
                --vl-accent-2: #174f49;
                --vl-accent-soft: #dceee9;
                --vl-danger: #a63a3a;
                --vl-shadow: 0 22px 70px rgba(42, 70, 58, 0.12);
                --vl-app-bg:
                    radial-gradient(circle at 6% 4%, rgba(29, 108, 99, 0.14), transparent 30rem),
                    radial-gradient(circle at 86% 0%, rgba(71, 105, 86, 0.12), transparent 28rem),
                    linear-gradient(180deg, var(--vl-bg), #fbfcf8 46%, #f6f8f4);
                --vl-sidebar-bg:
                    linear-gradient(180deg, rgba(251, 252, 248, 0.96), rgba(236, 243, 235, 0.98)),
                    var(--vl-surface);
                --vl-card-bg: linear-gradient(180deg, rgba(255, 255, 255, 0.96), rgba(247, 250, 246, 0.96));
                --vl-hero-bg:
                    radial-gradient(circle at 95% 10%, rgba(29, 108, 99, 0.15), transparent 16rem),
                    linear-gradient(135deg, rgba(255, 255, 255, 0.92), rgba(239, 246, 238, 0.88));
                --vl-soft-panel: rgba(255, 255, 255, 0.72);
                --vl-soft-panel-2: rgba(255, 255, 255, 0.74);
                --vl-upload-zone: rgba(238, 244, 238, 0.58);
                --vl-upload-zone-hover: rgba(220, 238, 233, 0.68);
                --vl-orb: rgba(29, 108, 99, 0.08);
            }
        </style>
        """
    return """
    <style>
        :root {
            --vl-bg: #07110f;
            --vl-bg-2: #0b1916;
            --vl-surface: #0d1b18;
            --vl-surface-strong: #10231f;
            --vl-surface-muted: #132d28;
            --vl-border: #24413b;
            --vl-border-strong: #38625a;
            --vl-text: #eef8f1;
            --vl-muted: #a1b4ab;
            --vl-faint: #6f877d;
            --vl-accent: #4fe0c6;
            --vl-accent-2: #21a894;
            --vl-accent-soft: rgba(79, 224, 198, 0.16);
            --vl-danger: #ff8f7f;
            --vl-shadow: 0 24px 80px rgba(0, 0, 0, 0.34);
            --vl-app-bg:
                radial-gradient(circle at 8% 0%, rgba(79, 224, 198, 0.16), transparent 31rem),
                radial-gradient(circle at 88% 3%, rgba(61, 126, 111, 0.18), transparent 30rem),
                linear-gradient(180deg, #06100e, #091614 42%, #0b1110);
            --vl-sidebar-bg:
                linear-gradient(180deg, rgba(12, 28, 24, 0.98), rgba(7, 17, 15, 0.98)),
                var(--vl-surface);
            --vl-card-bg: linear-gradient(180deg, rgba(19, 42, 37, 0.94), rgba(10, 26, 23, 0.96));
            --vl-hero-bg:
                radial-gradient(circle at 90% 7%, rgba(79, 224, 198, 0.18), transparent 18rem),
                linear-gradient(135deg, rgba(18, 43, 38, 0.96), rgba(7, 18, 16, 0.93));
            --vl-soft-panel: rgba(16, 35, 31, 0.72);
            --vl-soft-panel-2: rgba(14, 31, 28, 0.82);
            --vl-upload-zone: rgba(17, 42, 37, 0.68);
            --vl-upload-zone-hover: rgba(27, 65, 58, 0.78);
            --vl-orb: rgba(79, 224, 198, 0.1);
        }

        div[data-testid="stDataFrame"] * {
            color-scheme: dark;
        }
    </style>
    """


def _inject_css():
    st.markdown(
        """
        <style>
            :root {
                --vl-bg: #f4f7f3;
                --vl-bg-2: #e8efe8;
                --vl-surface: #fbfcf8;
                --vl-surface-strong: #ffffff;
                --vl-surface-muted: #eef4ee;
                --vl-border: #d6ded3;
                --vl-border-strong: #bac9b7;
                --vl-text: #17221d;
                --vl-muted: #637269;
                --vl-faint: #859188;
                --vl-accent: #1d6c63;
                --vl-accent-2: #174f49;
                --vl-accent-soft: #dceee9;
                --vl-warn: #9a6a1f;
                --vl-danger: #a63a3a;
                --vl-danger-soft: #f5dfdc;
                --vl-ok-soft: #dceee3;
                --vl-shadow: 0 22px 70px rgba(42, 70, 58, 0.12);
            }

            html, body, [class*="css"], .stApp {
                font-family: "Geist", ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
                font-variant-numeric: tabular-nums;
            }

            .stApp {
                color: var(--vl-text);
                background: var(--vl-app-bg);
            }

            .main .block-container {
                max-width: 1320px;
                padding: 2.1rem 2.5rem 4rem;
            }

            section[data-testid="stSidebar"] {
                background: var(--vl-sidebar-bg);
                border-right: 1px solid var(--vl-border);
            }

            section[data-testid="stSidebar"] h1 {
                font-size: 1.28rem;
                letter-spacing: 0;
                font-weight: 800;
            }

            section[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p {
                color: var(--vl-muted);
            }

            h1, h2, h3 {
                color: var(--vl-text);
                letter-spacing: 0;
            }

            h1 {
                font-size: clamp(2.35rem, 4vw, 4.4rem);
                line-height: 0.96;
                margin-bottom: 0.65rem;
                font-weight: 800;
            }

            h2, h3 {
                font-weight: 750;
            }

            hr {
                margin: 1.45rem 0;
                border-color: var(--vl-border);
            }

            div[data-testid="stFileUploader"] {
                background: var(--vl-card-bg);
                border: 1px solid var(--vl-border);
                border-radius: 14px;
                padding: 1rem;
                box-shadow: 0 16px 45px rgba(45, 74, 63, 0.07);
            }

            div[data-testid="stFileUploader"] section {
                border: 1.4px dashed var(--vl-border-strong);
                border-radius: 12px;
                background: var(--vl-upload-zone);
                transition: border-color 160ms ease, background 160ms ease, transform 160ms ease;
            }

            div[data-testid="stFileUploader"] section:hover {
                border-color: var(--vl-accent);
                background: var(--vl-upload-zone-hover);
                transform: translateY(-1px);
            }

            div[data-testid="stFileUploader"] button,
            div[data-testid="stFileUploader"] button[kind],
            div[data-testid="stFileUploader"] [data-testid="stBaseButton-secondary"] {
                background: linear-gradient(180deg, var(--vl-surface-muted), var(--vl-surface-strong));
                border: 1px solid var(--vl-border-strong);
                color: var(--vl-text);
                border-radius: 10px;
                box-shadow: none;
                font-weight: 750;
            }

            div[data-testid="stFileUploader"] button:hover,
            div[data-testid="stFileUploader"] [data-testid="stBaseButton-secondary"]:hover {
                background: var(--vl-accent-soft);
                border-color: var(--vl-accent);
                color: var(--vl-text);
                box-shadow: 0 10px 22px rgba(79, 224, 198, 0.13);
            }

            div[data-testid="stFileUploader"] button *,
            div[data-testid="stFileUploader"] [data-testid="stBaseButton-secondary"] * {
                color: var(--vl-text);
                fill: var(--vl-text);
                stroke: var(--vl-text);
            }

            div[data-testid="stFileUploader"] small,
            div[data-testid="stFileUploader"] span,
            div[data-testid="stFileUploader"] p {
                color: var(--vl-muted);
            }

            label, .stMarkdown p, .stCaptionContainer {
                color: var(--vl-muted);
            }

            div[data-testid="stMetric"] {
                background: var(--vl-card-bg), var(--vl-surface);
                border: 1px solid var(--vl-border);
                border-radius: 14px;
                padding: 1rem 1.05rem;
                box-shadow: var(--vl-shadow);
                min-height: 7rem;
            }

            div[data-testid="stMetricLabel"] p {
                color: var(--vl-muted);
                font-size: 0.76rem;
                font-weight: 700;
            }

            div[data-testid="stMetricValue"] {
                color: var(--vl-text);
                font-weight: 800;
            }

            div[data-testid="stMetricDelta"] {
                color: var(--vl-danger);
            }

            .vl-kicker {
                color: var(--vl-accent);
                font-size: 0.72rem;
                font-weight: 800;
                margin-bottom: 0.55rem;
                text-transform: uppercase;
                letter-spacing: 0.12em;
            }

            .vl-help {
                color: var(--vl-muted);
                font-size: 1.02rem;
                max-width: 820px;
                line-height: 1.62;
                margin-bottom: 1.1rem;
            }

            .vl-hero {
                position: relative;
                overflow: hidden;
                border: 1px solid var(--vl-border);
                border-radius: 18px;
                padding: 1.25rem 1.35rem 1.45rem;
                margin-bottom: 1.35rem;
                background: var(--vl-hero-bg);
                box-shadow: var(--vl-shadow);
            }

            .vl-hero::after {
                content: "";
                position: absolute;
                inset: auto -6rem -7rem auto;
                width: 18rem;
                height: 18rem;
                border-radius: 999px;
                background: var(--vl-orb);
                pointer-events: none;
            }

            .vl-stage-rail {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(8.5rem, 1fr));
                gap: 0.55rem;
                margin: 0.95rem 0 0.2rem;
            }

            .vl-stage {
                border: 1px solid var(--vl-border);
                border-radius: 999px;
                padding: 0.45rem 0.65rem;
                color: var(--vl-muted);
                background: var(--vl-soft-panel);
                font-size: 0.78rem;
                font-weight: 700;
                white-space: nowrap;
                text-align: center;
            }

            .vl-stage-active {
                border-color: rgba(29, 108, 99, 0.42);
                background: var(--vl-accent-soft);
                color: var(--vl-accent-2);
            }

            .vl-panel {
                background: var(--vl-card-bg);
                border: 1px solid var(--vl-border);
                border-radius: 14px;
                padding: 1.05rem 1.15rem;
                margin: 0.6rem 0 1.1rem;
                box-shadow: 0 14px 40px rgba(45, 74, 63, 0.06);
            }

            .vl-panel strong {
                color: var(--vl-text);
            }

            .vl-pill {
                display: inline-block;
                border: 1px solid var(--vl-border);
                border-radius: 999px;
                padding: 0.22rem 0.55rem;
                margin: 0.1rem;
                background: var(--vl-surface-muted);
                color: var(--vl-text);
                font-size: 0.82rem;
                font-family: "Geist Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
            }

            .vl-step {
                color: var(--vl-muted);
                font-size: 0.86rem;
                line-height: 1.55;
            }

            .vl-context-grid {
                display: grid;
                grid-template-columns: repeat(3, minmax(0, 1fr));
                gap: 0.75rem;
                margin: 1rem 0 1.25rem;
            }

            .vl-context {
                background: var(--vl-soft-panel);
                border: 1px solid var(--vl-border);
                border-radius: 12px;
                padding: 0.85rem 0.95rem;
            }

            .vl-context-label {
                color: var(--vl-faint);
                font-size: 0.7rem;
                font-weight: 800;
                letter-spacing: 0.1em;
                text-transform: uppercase;
            }

            .vl-context-value {
                color: var(--vl-text);
                font-family: "Geist Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
                font-size: 0.88rem;
                font-weight: 700;
                margin-top: 0.28rem;
                overflow-wrap: anywhere;
            }

            .vl-sidebar-card {
                border: 1px solid var(--vl-border);
                border-radius: 14px;
                padding: 0.9rem;
                margin: 0.75rem 0;
                background: var(--vl-soft-panel);
            }

            .vl-sidebar-label {
                color: var(--vl-faint);
                font-size: 0.68rem;
                font-weight: 800;
                letter-spacing: 0.1em;
                text-transform: uppercase;
            }

            .vl-sidebar-value {
                color: var(--vl-text);
                font-family: "Geist Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
                font-size: 0.82rem;
                font-weight: 700;
                margin-top: 0.22rem;
                overflow-wrap: anywhere;
            }

            .stButton > button,
            .stDownloadButton > button {
                border-radius: 10px;
                border: 1px solid var(--vl-border-strong);
                font-weight: 750;
                transition: transform 160ms ease, border-color 160ms ease, background 160ms ease, box-shadow 160ms ease;
            }

            .stButton > button:hover,
            .stDownloadButton > button:hover {
                transform: translateY(-1px);
                border-color: var(--vl-accent);
                box-shadow: 0 12px 24px rgba(29, 108, 99, 0.12);
            }

            .stButton > button:active,
            .stDownloadButton > button:active {
                transform: translateY(0);
            }

            button[kind="primary"],
            .stButton > button[kind="primary"] {
                border-color: var(--vl-accent-2);
                background: linear-gradient(180deg, var(--vl-accent), var(--vl-accent-2));
                color: #ffffff;
            }

            div[data-testid="stExpander"] {
                border: 1px solid var(--vl-border);
                border-radius: 13px;
                background: var(--vl-soft-panel-2);
                box-shadow: 0 12px 34px rgba(45, 74, 63, 0.05);
            }

            div[data-testid="stExpander"] summary {
                font-weight: 750;
                color: var(--vl-text);
            }

            div[data-testid="stDataFrame"] {
                border: 1px solid var(--vl-border);
                border-radius: 13px;
                overflow: hidden;
                box-shadow: 0 16px 40px rgba(45, 74, 63, 0.06);
            }

            div[data-baseweb="tab-list"] {
                gap: 0.45rem;
                border-bottom: 1px solid var(--vl-border);
            }

            button[data-baseweb="tab"] {
                border-radius: 10px 10px 0 0;
                font-weight: 750;
            }

            div[data-testid="stAlert"] {
                border-radius: 12px;
                border: 1px solid var(--vl-border);
            }

            .stProgress > div > div > div {
                background-color: var(--vl-accent);
            }

            input, textarea, div[data-baseweb="select"] > div {
                border-radius: 10px;
            }

            div[role="radiogroup"] {
                gap: 0.35rem;
            }

            div[role="radiogroup"] label {
                border: 1px solid var(--vl-border);
                border-radius: 999px;
                padding: 0.25rem 0.65rem;
                background: var(--vl-soft-panel);
                margin-right: 0.25rem;
            }

            div[role="radiogroup"] label:has(input:checked) {
                border-color: var(--vl-accent);
                background: var(--vl-accent-soft);
                color: var(--vl-text);
            }

            div[data-baseweb="input"] > div,
            div[data-baseweb="select"] > div {
                background: var(--vl-soft-panel);
                border-color: var(--vl-border);
                color: var(--vl-text);
            }

            div[data-baseweb="popover"] {
                color: var(--vl-text);
            }

            @media (max-width: 760px) {
                .main .block-container {
                    padding: 1.15rem 1rem 2rem;
                }

                .vl-stage-rail,
                .vl-context-grid {
                    grid-template-columns: 1fr;
                }

                h1 {
                    font-size: 2.35rem;
                }
            }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.markdown(_theme_overrides(), unsafe_allow_html=True)


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

def _save_upload(uploaded_file) -> Path:
    """Save a Streamlit UploadedFile to the session temp dir. Returns path."""
    tmp_dir = Path(st.session_state.tmp.name)
    dest = tmp_dir / uploaded_file.name
    dest.write_bytes(uploaded_file.read())
    return dest


def _suggest_virus_name(record: SeqRecord) -> str:
    """Best-effort display name for a virus without an existing alias config."""
    organism = record.annotations.get("organism")
    if organism:
        return organism
    return record.description or record.id or "new virus"


def _record_metadata_candidates(record: SeqRecord) -> List[str]:
    candidates = [
        record.annotations.get("organism", ""),
        getattr(record, "description", "") or "",
        getattr(record, "name", "") or "",
        getattr(record, "id", "") or "",
    ]
    result = []
    seen = set()
    for item in candidates:
        value = str(item or "").strip()
        key = normalize_text(value)
        if key and key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _split_keywords(raw_text: str, virus_name: str) -> List[str]:
    """Parse comma/newline separated registry keywords and include virus name."""
    values = [virus_name]
    for part in re.split(r"[,\n]+", raw_text or ""):
        item = part.strip()
        if item:
            values.append(item)

    deduped = []
    seen = set()
    for item in values:
        key = normalize_text(item)
        if key and key not in seen:
            seen.add(key)
            deduped.append(item)
    return deduped


def _unique_alias_config_paths(filename: str) -> Tuple[Path, Path]:
    """Return absolute/relative alias config paths without overwriting existing files."""
    candidate = CONFIG_DIR / filename
    if not candidate.exists():
        return candidate, Path("app/config") / filename

    stem = Path(filename).stem
    suffix = Path(filename).suffix
    index = 2
    while True:
        next_name = f"{stem}_{index}{suffix}"
        candidate = CONFIG_DIR / next_name
        if not candidate.exists():
            return candidate, Path("app/config") / next_name
        index += 1


def _continue_with_existing_alias_config(entry: Dict, save_keyword: Optional[str] = None) -> None:
    alias_config = entry.get("alias_config")
    alias_config_path = resolve_config_path(Path(alias_config), ROOT)
    if save_keyword:
        add_registry_keyword(REGISTRY_PATH, alias_config, save_keyword)

    alias_lookup = load_alias_lookup(alias_config_path)
    ref_features = apply_alias_to_features(st.session_state.ref_features, alias_lookup)
    ignored = _load_ignored_names(alias_config_path)
    unknown = _scan_unknown_names(st.session_state.query_records, alias_lookup, ignored)
    unknown_ref = _scan_unknown_ref_names(ref_features, ignored)

    st.session_state.ref_features = ref_features
    st.session_state.alias_lookup = alias_lookup
    st.session_state.alias_config_path = alias_config_path
    st.session_state.virus_name = entry.get("virus_name")
    st.session_state.canonical_list = sorted(set(alias_lookup.values()))
    st.session_state.unknown_names = unknown
    st.session_state.unknown_ref_names = unknown_ref
    st.session_state.stage = "resolve" if (unknown or unknown_ref) else "running"


def _selected_dataframe_rows(state) -> List[int]:
    """Extract selected row indices from Streamlit dataframe selection state."""
    if not state:
        return []
    selection = getattr(state, "selection", None)
    if selection is None and isinstance(state, dict):
        selection = state.get("selection")
    if not selection:
        return []
    rows = getattr(selection, "rows", None)
    if rows is None and isinstance(selection, dict):
        rows = selection.get("rows")
    return list(rows or [])


def _scan_unknown_names(
    query_records: List[SeqRecord],
    alias_lookup: Dict[str, str],
    ignored_names: set,
) -> Dict[str, Dict]:
    """
    Return a dict keyed by representative name for every feature group in query
    records that needs user resolution: either the name is completely unknown
    (misses alias lookup entirely) or it is explicitly ambiguous (maps to
    AMBIGUOUS_SENTINEL, shared by multiple genes).

    Each value is:
        {
            "records":    [record_id, ...],   # records containing this feature
            "candidates": [val, ...],         # all qualifier values, priority order
            "ambiguous":  bool,               # True means known-ambiguous; False means unknown
        }

    Resolution logic mirrors apply_alias_to_feature:
        - ANY candidate hits a real canonical: feature is resolved, skip
        - ALL candidates miss OR hit AMBIGUOUS: include (unknown or ambiguous)
        - ALL candidates hit IGNORED: skip (intentionally excluded)
    """
    result: Dict[str, Dict] = {}

    for rec in query_records:
        selected_feature_type = select_feature_type(rec, alias_lookup)
        if selected_feature_type is None:
            continue

        for feat in rec.features:
            if feat.type != selected_feature_type:
                continue

            # Collect all unique qualifier values in priority order
            seen_vals: set = set()
            candidates = []
            for field in _LOOKUP_QUALIFIER_KEYS:
                val = feat.qualifiers.get(field, [None])[0]
                if val and val not in seen_vals:
                    seen_vals.add(val)
                    candidates.append(val)

            if not candidates:
                continue

            # Classify each candidate's hit in the alias lookup
            hits = {v: alias_lookup.get(normalize_text(v)) for v in candidates}

            # If ANY candidate resolves to a real canonical, already handled, skip
            if any(
                h is not None and h not in (AMBIGUOUS_SENTINEL, IGNORED_SENTINEL)
                for h in hits.values()
            ):
                continue

            # If ALL candidates resolve to IGNORED, skip entirely
            if all(v.lower() in ignored_names for v in candidates):
                continue

            # Determine if this is ambiguous or fully unknown
            is_ambiguous = any(h == AMBIGUOUS_SENTINEL for h in hits.values())

            non_ignored_candidates = [
                v for v in candidates if v.lower() not in ignored_names
            ]
            if not non_ignored_candidates:
                continue

            representative = non_ignored_candidates[0]

            if representative not in result:
                result[representative] = {
                    "records":   [],
                    "candidates": non_ignored_candidates,
                    "ambiguous": is_ambiguous,
                }
            if rec.id not in result[representative]["records"]:
                result[representative]["records"].append(rec.id)

    return result


def _scan_unknown_ref_names(ref_features: list, ignored_names: set) -> List[str]:
    """
    Return a sorted list of ref feature raw names that were not resolved by the
    alias DB (i.e. name_source == 'raw') and are not explicitly ignored.
    These will be lifted correctly (tblastn uses protein sequence), but their
    output name will just be the raw annotation name, no canonical key.
    """
    seen = []
    for f in ref_features:
        if f.get("name_source") != "raw":
            continue
        raw = f.get("raw_name") or f.get("name") or ""
        if raw.lower() in ignored_names:
            continue
        if raw and raw not in seen:
            seen.append(raw)
    return sorted(seen)


def _add_new_canonicals_to_config(
    alias_config_path: Path,
    new_canonicals: List[str],
) -> int:
    """
    Add brand-new canonical entries (with empty alias list) to the alias JSON config.
    Skips any canonical key that already exists.
    Returns the number of entries actually added.
    """
    if not alias_config_path or not alias_config_path.exists():
        return 0
    with open(alias_config_path) as f:
        cfg = json.load(f)
    added = 0
    config_name = alias_config_path.name
    for name in new_canonicals:
        if name not in cfg["canonical_names"]:
            cfg["canonical_names"][name] = []
            log_canonical_added(config_name, name)
            added += 1
    if added:
        with open(alias_config_path, "w") as f:
            json.dump(cfg, f, indent=2)
    return added


def _load_ignored_names(alias_config_path: Optional[Path]) -> set:
    if alias_config_path is None or not alias_config_path.exists():
        return set()
    with open(alias_config_path) as f:
        cfg = json.load(f)
    return {n.lower() for n in cfg.get("ignored_names", [])}


def _build_effective_lookup(
    base_lookup: Dict[str, str],
    resolver: Dict[str, str],
) -> Dict[str, str]:
    """Merge base alias lookup with user resolver decisions."""
    effective = dict(base_lookup)
    for raw_name, canonical in resolver.items():
        if canonical and canonical != "-- ignore --":
            effective[normalize_text(raw_name)] = canonical
    return effective


def _canonical_to_ref_map(ref_features: list) -> Dict[str, str]:
    """
    Build {canonical_name: ref_raw_name} from ref features after alias normalization.
    e.g. {"Lpro": "Lab", "3Cpro": "3C", ...}
    """
    return {
        f["name"]: f["raw_name"]
        for f in ref_features
        if f.get("raw_name") and f.get("name")
    }


def _results_to_df(all_results, ref_name_map: Dict[str, str] = None) -> pd.DataFrame:
    rows = []
    for query_id, features in all_results:
        for lf in features:
            display_name = (
                ref_name_map.get(lf.name, lf.name)
                if ref_name_map else lf.name
            )
            rows.append({
                "record_id":  query_id,
                "name":       display_name,
                "source_name": lf.source_name or "",
                "start":      lf.query_start,
                "end":        lf.query_end,
                "strand":     lf.strand,
                "status":     lf.status,
                "coverage":   round(lf.coverage, 3) if lf.coverage else None,
                "identity":   lf.identity,
                "method":     lf.method,
                "sequence":   lf.sequence or "",
            })
    return pd.DataFrame(rows)


def _run_pipeline(
    ref_record, query_records, ref_features, ref_feature_type,
    effective_lookup, min_coverage, min_identity, evalue, rescue_window,
    progress_bar,
    virus_name: Optional[str] = None,
    alias_config_path=None,
    run_errors: Optional[List[Dict[str, str]]] = None,
) -> List[Tuple[str, List]]:

    log_run_start(
        ref_id=ref_record.id,
        n_queries=len(query_records),
        min_coverage=min_coverage,
        min_identity=min_identity,
        evalue=evalue,
        rescue_window=rescue_window,
        virus_name=virus_name,
        alias_config=Path(alias_config_path).name if alias_config_path else None,
    )

    all_results = []
    n = len(query_records)
    for i, qrec in enumerate(query_records):
        progress_bar.progress(i / n, text=f"{_t('processing')} {qrec.id}  ({i+1}/{n})")
        strategy, query_feature_type = get_strategy(qrec, effective_lookup)
        try:
            if strategy == "direct":
                results = direct_extract_with_alias(
                    qrec, query_feature_type, ref_features, effective_lookup
                )
            else:
                results = process_one_query_record(
                    ref_record=ref_record,
                    query_record=qrec,
                    ref_cds=ref_features,
                    ref_feature_type=ref_feature_type,
                    min_coverage=min_coverage,
                    min_identity=min_identity,
                    evalue=evalue,
                    rescue_window=rescue_window,
                    quiet=True,
                )
        except Exception as e:
            log_error(f"processing record {qrec.id}", e)
            if run_errors is not None:
                run_errors.append({"record_id": qrec.id, "error": str(e)})
            results = []

        all_results.append((qrec.id, results))

    progress_bar.progress(1.0, text=_t("done"))

    summary = summarize_counts(all_results)
    log_run_complete(ref_id=ref_record.id, n_queries=n, summary=summary)

    return all_results


# ═══════════════════════════════════════════════════════════════════
# Stage: UPLOAD
# ═══════════════════════════════════════════════════════════════════

def stage_upload():
    _render_page_intro(
        _t("run_setup"),
        _t("upload_title"),
        _t("upload_body"),
    )

    col_ref, col_query = st.columns(2)
    with col_ref:
        st.markdown(f"**{_t('reference_record')}**")
        st.caption(_t("reference_caption"))
        ref_file = st.file_uploader("Reference GenBank (.gb)", type=["gb", "gbk"], label_visibility="collapsed")
    with col_query:
        st.markdown(f"**{_t('query_records')}**")
        st.caption(_t("query_caption"))
        query_file = st.file_uploader("Query GenBank (.gb)", type=["gb", "gbk"], label_visibility="collapsed")

    with st.expander(_t("file_requirements"), expanded=False):
        st.markdown(_t("file_requirements_body"))

    st.divider()
    st.subheader(_t("advanced_options"))
    adv = st.expander(_t("lifting_thresholds"), expanded=False)
    with adv:
        st.caption(_t("threshold_caption"))
        c1, c2, c3, c4 = st.columns(4)
        st.session_state.min_coverage   = c1.number_input(_t("min_coverage"),  0.0, 1.0, 0.5, 0.05)
        st.session_state.min_identity   = c2.number_input(_t("min_identity"),  0.0, 1.0, 0.3, 0.05)
        st.session_state.evalue         = c3.number_input(_t("evalue"),       value=1e-5, format="%.0e")
        st.session_state.rescue_window  = c4.number_input(_t("rescue_window"), 10,  200,   50,    10)

    st.divider()
    st.session_state.use_ref_names = st.toggle(
        _t("use_ref_names"),
        value=False,
        help=_t("use_ref_names_help"),
    )

    ready = ref_file and query_file
    if st.button(_t("run_button"), disabled=not ready, type="primary", width="stretch"):
        # persist files to temp dir
        if st.session_state.tmp is None:
            st.session_state.tmp = tempfile.TemporaryDirectory()

        with st.spinner(_t("loading_files")):
            try:
                ref_path   = _save_upload(ref_file)
                query_path = _save_upload(query_file)

                ref_record    = load_single_genbank(ref_path)
                query_records = load_genbank_records(query_path)

                ref_features, ref_feature_type, alias_config_path, virus_name, alias_lookup = (
                    prepare_reference_features(
                        ref_record=ref_record,
                        alias_config_arg=None,
                        alias_registry_arg=str(REGISTRY_PATH),
                    )
                )

                if alias_config_path is None:
                    virus_guess = _suggest_virus_name(ref_record)
                    seed_config = build_seed_alias_config_from_ref(
                        ref_record=ref_record,
                        ref_features=ref_features,
                        virus_name=virus_guess,
                    )

                    st.session_state.ref_record        = ref_record
                    st.session_state.query_records     = query_records
                    st.session_state.ref_features      = ref_features
                    st.session_state.ref_feature_type  = ref_feature_type
                    st.session_state.alias_lookup      = {}
                    st.session_state.alias_config_path = None
                    st.session_state.virus_name        = virus_guess
                    st.session_state.canonical_list    = sorted(seed_config["canonical_names"])
                    st.session_state.unknown_names     = {}
                    st.session_state.unknown_ref_names = []
                    st.session_state.bootstrap_alias_config = seed_config
                    st.session_state.bootstrap_virus_name = virus_guess
                    st.session_state.bootstrap_keywords = virus_guess
                    st.session_state.bootstrap_suggestions = []
                    st.session_state.bootstrap_diagnostics = {}
                    st.session_state.virus_review_metadata = _record_metadata_candidates(ref_record)
                    st.session_state.virus_review_guess = virus_guess
                    st.session_state.stage = "virus_review"
                    st.rerun()

                # canonical list for resolver dropdowns
                canonical_list = sorted(set(alias_lookup.values())) if alias_lookup else []

                # scan for unknown names in query records AND in ref
                ignored = _load_ignored_names(alias_config_path)
                unknown     = _scan_unknown_names(query_records, alias_lookup, ignored)
                unknown_ref = _scan_unknown_ref_names(ref_features, ignored)
            except Exception as e:
                log_error("loading UI inputs", e)
                st.error(f"{_t('load_error')}: {e}")
                return

        # store to session state
        st.session_state.ref_record        = ref_record
        st.session_state.query_records     = query_records
        st.session_state.ref_features      = ref_features
        st.session_state.ref_feature_type  = ref_feature_type
        st.session_state.alias_lookup      = alias_lookup
        st.session_state.alias_config_path = alias_config_path
        st.session_state.virus_name        = virus_name
        st.session_state.canonical_list    = canonical_list
        st.session_state.unknown_names     = unknown
        st.session_state.unknown_ref_names = unknown_ref

        if unknown or unknown_ref:
            st.session_state.stage = "resolve"
        else:
            st.session_state.stage = "running"

        st.rerun()

    if ref_file and query_file:
        _render_context_panel([
            (_t("reference_file"), ref_file.name),
            (_t("query_file"), query_file.name),
            (_t("name_mode"), _t("reference_names") if st.session_state.use_ref_names else _t("canonical_names")),
        ])


# ═══════════════════════════════════════════════════════════════════
# Stage: VIRUS REVIEW
# ═══════════════════════════════════════════════════════════════════

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


# ═══════════════════════════════════════════════════════════════════
# Stage: BOOTSTRAP ALIAS CONFIG
# ═══════════════════════════════════════════════════════════════════

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


# ═══════════════════════════════════════════════════════════════════
# Stage: RESOLVE
# ═══════════════════════════════════════════════════════════════════

def _save_to_alias_config(alias_config_path: Path, mappings: Dict[str, str]) -> int:
    """
    Persist user-confirmed mappings into the alias JSON config file.

    For each (raw_name -> canonical) pair, appends raw_name to the canonical's
    alias list if not already present.

    Returns the number of new aliases written.
    """
    if not alias_config_path or not alias_config_path.exists():
        return 0

    with open(alias_config_path) as f:
        cfg = json.load(f)

    written = 0
    config_name = alias_config_path.name
    for raw_name, canonical in mappings.items():
        aliases = cfg["canonical_names"].get(canonical)
        if aliases is None:
            continue
        if raw_name not in aliases:
            aliases.append(raw_name)
            log_alias_added(config_name, raw_name, canonical)
            written += 1

    with open(alias_config_path, "w") as f:
        json.dump(cfg, f, indent=2)

    return written


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


# ═══════════════════════════════════════════════════════════════════
# Stage: RUNNING  (transient, immediately transitions to results)
# ═══════════════════════════════════════════════════════════════════

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


# ═══════════════════════════════════════════════════════════════════
# Stage: RESULTS
# ═══════════════════════════════════════════════════════════════════

def stage_results():
    all_results = st.session_state.all_results
    summary     = summarize_counts(all_results)
    run_errors  = st.session_state.run_errors or []
    error_map   = _error_by_record(run_errors)

    # build ref name map if user wants ref names as output
    # {canonical_key: ref_raw_name}  e.g. {"Lpro": "Lab", "3Cpro": "3C", ...}
    ref_name_map: Optional[Dict[str, str]] = (
        _canonical_to_ref_map(st.session_state.ref_features)
        if st.session_state.use_ref_names else None
    )

    df = _results_to_df(all_results, ref_name_map)

    # ── header + new run button ──────────────────────────────────
    h_col, btn_col = st.columns([5, 1])
    with h_col:
        _render_page_intro(
            _t("run_review"),
            _t("results_title"),
            _t("results_body"),
        )
    if st.session_state.use_ref_names:
        h_col.caption(_t("ref_names_caption"))
    if btn_col.button(_t("new_run")):
        _reset()
        st.rerun()

    total_records = len(all_results)
    failed_records = len(run_errors)
    total_features = sum(len(features) for _, features in all_results)
    passing_features = sum(
        1
        for _, features in all_results
        for feature in features
        if feature.status in GOOD_STATUSES
    )
    review_features = sum(
        1
        for _, features in all_results
        for feature in features
        if feature.status in REVIEW_STATUSES
    )
    pass_rate = (passing_features / total_features * 100) if total_features else 0.0

    # ── run health ────────────────────────────────────────────────
    health_cols = st.columns(4)
    health_cols[0].metric(_t("records_processed"), total_records - failed_records, delta=_t("failed_delta", count=failed_records) if failed_records else None)
    health_cols[1].metric(_t("features_found"), total_features)
    health_cols[2].metric(_t("pass_rate"), f"{pass_rate:.1f}%")
    health_cols[3].metric(_t("needs_review"), review_features)

    if run_errors:
        st.error(_t("processing_error", count=failed_records))
        with st.expander(_t("processing_errors"), expanded=True):
            st.dataframe(pd.DataFrame(run_errors), width="stretch", hide_index=True)

    status_rows = [
        {
            "status": _t(f"status_{status}"),
            "count": count,
            "review": _t("yes") if status in REVIEW_STATUSES else _t("no"),
        }
        for status, count in summary.items()
        if count
    ]
    if status_rows:
        with st.expander(_t("status_breakdown"), expanded=review_features > 0):
            st.dataframe(pd.DataFrame(status_rows), width="stretch", hide_index=True)

    st.divider()

    # ── per-record overview ───────────────────────────────────────
    st.subheader(_t("record_overview"))
    overview_rows = []
    for query_id, features in all_results:
        ok_count     = sum(1 for f in features if f.status in GOOD_STATUSES)
        total_count  = len(features)
        needs_review = sum(1 for f in features if f.status in REVIEW_STATUSES)
        methods      = {f.method for f in features if f.method}
        if methods == {"direct"}:
            method_tag = "direct"
        elif "tblastn" in methods and "direct" not in methods:
            method_tag = "tblastn"
        elif methods:
            method_tag = "tblastn + direct"
        else:
            method_tag = "none"

        if query_id in error_map:
            health = _t("failed").lower()
        elif total_count == 0:
            health = _t("empty").lower()
        elif needs_review:
            health = _t("needs_review").lower()
        else:
            health = _t("passed").lower()

        overview_rows.append({
            "record_id": query_id,
            "health": health,
            "mapped": ok_count,
            "total": total_count,
            "needs_review": needs_review,
            "method": method_tag,
            "error": error_map.get(query_id, ""),
        })

    overview_df = pd.DataFrame(overview_rows)
    f_col1, f_col2, f_col3 = st.columns([2, 2, 1])
    all_filter = _t("all")
    health_options = [all_filter] + sorted(overview_df["health"].unique().tolist())
    selected_health = f_col1.selectbox(_t("health_filter"), health_options)
    search_text = f_col2.text_input(_t("search_record"))
    show_details = f_col3.checkbox(_t("open_all"), value=False)

    visible_overview = overview_df.copy()
    if selected_health != all_filter:
        visible_overview = visible_overview[visible_overview["health"] == selected_health]
    if search_text:
        visible_overview = visible_overview[
            visible_overview["record_id"].str.contains(search_text, case=False, na=False)
        ]

    st.dataframe(visible_overview, width="stretch", hide_index=True)

    st.subheader(_t("record_details"))
    visible_ids = set(visible_overview["record_id"].tolist())
    for query_id, features in all_results:
        if query_id not in visible_ids:
            continue

        ok_count     = sum(1 for f in features if f.status in GOOD_STATUSES)
        total_count  = len(features)
        needs_review = sum(1 for f in features if f.status in REVIEW_STATUSES)

        if query_id in error_map:
            label = f"{_t('failed')} · {query_id}"
        elif total_count == 0:
            label = f"{_t('empty')} · {query_id}"
        elif needs_review:
            label = f"{_t('needs_review')} · {query_id} · {ok_count}/{total_count} mapped"
        else:
            label = f"{_t('passed')} · {query_id} · {ok_count}/{total_count} mapped"

        with st.expander(label, expanded=show_details or query_id in error_map):
            if query_id in error_map:
                st.error(error_map[query_id])
                continue
            if not features:
                st.warning(_t("no_features"))
                continue
            rec_rows = []
            for lf in features:
                display_name = (
                    ref_name_map.get(lf.name, lf.name)
                    if ref_name_map else lf.name
                )
                rec_rows.append({
                    "status":    _status_text(lf.status),
                    "name":      display_name,
                    "raw_name":  lf.source_name or "",
                    "start":     lf.query_start,
                    "end":       lf.query_end,
                    "coverage":  _format_fraction_percent(lf.coverage),
                    "identity":  _format_percent(lf.identity),
                    "method":    lf.method,
                })
            st.dataframe(pd.DataFrame(rec_rows), width="stretch", hide_index=True)

    st.divider()

    # ── export section ────────────────────────────────────────────
    st.subheader(_t("export"))
    tab_tsv, tab_fasta = st.tabs([_t("tsv"), _t("fasta")])

    # TSV tab
    with tab_tsv:
        st.markdown(_t("download_table"))
        col_canon, col_raw = st.columns(2)

        tsv_canonical = df.drop(columns=["sequence"]).to_csv(sep="\t", index=False)
        col_canon.download_button(
            _t("download_display"),
            data=tsv_canonical,
            file_name="viralift_canonical.tsv",
            mime="text/tab-separated-values",
            width="stretch",
        )

        df_raw = df.copy()
        has_source = df_raw["source_name"].notna() & (df_raw["source_name"] != "")
        df_raw.loc[has_source, "name"] = df_raw.loc[has_source, "source_name"]
        tsv_raw = df_raw.drop(columns=["sequence"]).to_csv(sep="\t", index=False)
        col_raw.download_button(
            _t("download_raw"),
            data=tsv_raw,
            file_name="viralift_raw.tsv",
            mime="text/tab-separated-values",
            width="stretch",
        )

    # FASTA extraction tab
    with tab_fasta:
        st.markdown(_t("fasta_intro"))

        # gene list: display names respect the use_ref_names toggle
        if ref_name_map:
            ref_gene_names = sorted({
                ref_name_map.get(f.get("name"), f.get("name"))
                for f in st.session_state.ref_features
                if f.get("name")
            })
        else:
            ref_gene_names = sorted({
                f.get("name")
                for f in st.session_state.ref_features
                if f.get("name")
            })

        col_select, col_format = st.columns([3, 2])
        selected_genes = col_select.multiselect(
            _t("genes_to_extract"),
            options=ref_gene_names,
            default=ref_gene_names,
        )
        fasta_mode = col_format.radio(
            _t("output_format"),
            [_t("one_fasta"), _t("all_fasta")],
        )

        # quality filter
        st.markdown(f"**{_t('quality_filter')}**")
        qf_col1, qf_col2, qf_col3 = st.columns(3)
        min_cov_export = qf_col1.slider(_t("min_coverage"), 0.0, 1.0, 0.5, 0.05)
        min_id_export  = qf_col2.slider(_t("min_identity"),  0.0, 100.0, 0.0, 5.0)
        include_rescued = qf_col3.checkbox(_t("include_rescued"), value=True)

        accepted_statuses = {"ok", "direct"}
        if include_rescued:
            accepted_statuses.add("ok_rescued")

        candidate_count = 0
        for _, features in all_results:
            for lf in features:
                gene = ref_name_map.get(lf.name, lf.name) if ref_name_map else lf.name
                if gene not in selected_genes:
                    continue
                if lf.status not in accepted_statuses:
                    continue
                if lf.coverage is not None and lf.coverage < min_cov_export:
                    continue
                if lf.identity is not None and lf.identity < min_id_export:
                    continue
                if lf.sequence:
                    candidate_count += 1
        st.caption(_t("candidate_count", count=candidate_count))

        if st.button(_t("generate_fasta"), type="primary"):
            # build gene to sequences dict (keyed by display name)
            gene_seqs: Dict[str, List[str]] = {g: [] for g in selected_genes}
            skipped = 0

            for query_id, features in all_results:
                for lf in features:
                    # resolve display name the same way as the gene list above
                    gene = (
                        ref_name_map.get(lf.name, lf.name)
                        if ref_name_map else lf.name
                    )
                    if gene not in gene_seqs:
                        continue
                    if lf.status not in accepted_statuses:
                        skipped += 1
                        continue
                    if lf.coverage is not None and lf.coverage < min_cov_export:
                        skipped += 1
                        continue
                    if lf.identity is not None and lf.identity < min_id_export:
                        skipped += 1
                        continue
                    if not lf.sequence:
                        continue
                    header = f">{query_id}|{gene}|{lf.query_start}|{lf.query_end}|{lf.strand}"
                    gene_seqs[gene].append(f"{header}\n{lf.sequence}")

            if fasta_mode == _t("all_fasta"):
                all_seqs = "\n".join(
                    seq for gene in selected_genes for seq in gene_seqs[gene]
                )
                st.download_button(
                    _t("download_all_fasta"),
                    data=all_seqs,
                    file_name="all_genes.fasta",
                    mime="text/plain",
                    width="stretch",
                )
            else:
                # one download button per gene
                for gene in selected_genes:
                    seqs = gene_seqs[gene]
                    if not seqs:
                        st.warning(_t("no_sequences", gene=gene))
                        continue
                    st.download_button(
                        _t("gene_download", gene=gene, count=len(seqs)),
                        data="\n".join(seqs),
                        file_name=f"{gene}.fasta",
                        mime="text/plain",
                        key=f"dl_{gene}",
                    )

            if skipped:
                st.caption(_t("skipped", count=skipped))


# ═══════════════════════════════════════════════════════════════════
# Page: ALIAS MANAGER
# ═══════════════════════════════════════════════════════════════════

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
        f"{entry.get('virus_name', 'unknown')} — {entry.get('alias_config', '')}"
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

    tab_registry, tab_alias, tab_ignored, tab_ambiguous, tab_raw = st.tabs([
        "Registry",
        "Canonical aliases",
        "Ignored names",
        "Ambiguous names",
        "Raw JSON",
    ])

    with tab_registry:
        st.caption("Auto-detection uses keywords. Virus name is the display label.")
        registry_virus_name = st.text_input(
            "Virus name",
            value=entry.get("virus_name", config.get("virus", "")),
            key=f"registry_virus_name_{config_path}",
        )
        keyword_rows = [{"keyword": keyword} for keyword in entry.get("keywords", [])]
        edited_keywords = st.data_editor(
            pd.DataFrame(keyword_rows, columns=["keyword"]),
            num_rows="dynamic",
            width="stretch",
            hide_index=True,
            column_config={"keyword": st.column_config.TextColumn("Detection keyword")},
            key=f"registry_keywords_{config_path}",
        )
        if st.button("Save registry entry", type="primary", width="stretch"):
            try:
                keywords = [
                    str(row.get("keyword") or "").strip()
                    for row in edited_keywords.to_dict("records")
                    if str(row.get("keyword") or "").strip()
                ]
                backup = update_registry_entry(
                    REGISTRY_PATH,
                    alias_config=entry.get("alias_config", ""),
                    virus_name=registry_virus_name.strip(),
                    keywords=keywords,
                )
                st.success(f"Registry saved. Backup: `{backup.name if backup else 'none'}`")
                st.rerun()
            except Exception as e:
                st.error(f"Could not save registry entry: {e}")

    with tab_alias:
        st.caption("Each canonical has its own table. Select alias row(s), then use the delete button that appears below the table.")
        search_text = st.text_input(
            "Search canonical or alias",
            value="",
            key=f"alias_manager_search_{config_path}",
        )

        st.markdown("**Add canonical / alias**")
        new_alias_rows = st.data_editor(
            pd.DataFrame(columns=["canonical_name", "alias"]),
            num_rows="dynamic",
            width="stretch",
            hide_index=True,
            column_config={
                "canonical_name": st.column_config.TextColumn("Canonical name", required=True),
                "alias": st.column_config.TextColumn("Alias"),
            },
            key=f"alias_manager_new_aliases_{config_path}",
        )

        alias_records_for_save = []
        canonical_names = config.get("canonical_names", {})
        search_norm = normalize_text(search_text)
        for canonical in sorted(canonical_names, key=normalize_text):
            aliases = canonical_names.get(canonical, []) or []
            searchable = normalize_text(" ".join([canonical] + aliases))
            if search_norm and search_norm not in searchable:
                continue

            label = f"{canonical} · {len(aliases)} alias(es)"
            with st.expander(label, expanded=bool(search_norm)):
                delete_canonical = st.checkbox(
                    f"Delete canonical `{canonical}`",
                    value=False,
                    key=f"delete_canonical_{config_path}_{normalize_text(canonical)}",
                    help="Removes this canonical and all aliases when you save.",
                )
                alias_df = pd.DataFrame(
                    [{"alias": alias} for alias in aliases],
                    columns=["alias"],
                )
                alias_table_state = st.dataframe(
                    alias_df,
                    width="stretch",
                    hide_index=True,
                    on_select="rerun",
                    selection_mode="multi-row",
                    column_config={
                        "alias": st.column_config.TextColumn("Alias"),
                    },
                    key=f"alias_editor_{config_path}_{normalize_text(canonical)}",
                )
                selected_rows = _selected_dataframe_rows(alias_table_state)
                selected_aliases = [
                    alias_df.iloc[row]["alias"]
                    for row in selected_rows
                    if 0 <= row < len(alias_df)
                ]

                if selected_aliases:
                    st.caption(
                        "Selected: " + ", ".join(f"`{alias}`" for alias in selected_aliases)
                    )
                delete_col, edit_hint_col = st.columns([1, 3])
                if delete_col.button(
                    f"Delete selected ({len(selected_aliases)})",
                    disabled=not selected_aliases,
                    key=f"delete_selected_aliases_{config_path}_{normalize_text(canonical)}",
                ):
                    updated = json.loads(json.dumps(config))
                    selected_norms = {normalize_text(alias) for alias in selected_aliases}
                    updated["canonical_names"][canonical] = [
                        alias
                        for alias in updated.get("canonical_names", {}).get(canonical, [])
                        if normalize_text(alias) not in selected_norms
                    ]
                    backup = manager_save_alias_config(config_path, updated)
                    st.success(
                        f"Deleted {len(selected_aliases)} alias(es) from `{canonical}`. "
                        f"Backup: `{backup.name if backup else 'none'}`"
                    )
                    remaining_warnings = validate_alias_config(updated)
                    if remaining_warnings:
                        st.info(
                            "Alias was deleted. Remaining config warnings still need review:\n\n"
                            + "\n".join(f"- {w}" for w in remaining_warnings)
                        )
                    st.rerun()
                edit_hint_col.caption("To edit an alias value, delete the old row and add the corrected alias in the Add canonical / alias table.")

                if delete_canonical:
                    st.warning(f"`{canonical}` will be deleted on save.")
                    continue

                alias_records_for_save.append({
                    "canonical_name": canonical,
                    "aliases": "\n".join(aliases),
                })

        for row in new_alias_rows.to_dict("records"):
            canonical = str(row.get("canonical_name") or "").strip()
            alias = str(row.get("alias") or "").strip()
            if not canonical:
                continue
            alias_records_for_save.append({
                "canonical_name": canonical,
                "aliases": alias,
            })

    with tab_ignored:
        st.caption("Names here are intentionally ignored by automatic alias matching.")
        ignored_df = pd.DataFrame(ignored_rows, columns=["ignored_name"])
        ignored_table_state = st.dataframe(
            ignored_df,
            width="stretch",
            hide_index=True,
            on_select="rerun",
            selection_mode="multi-row",
            column_config={
                "ignored_name": st.column_config.TextColumn("Ignored name"),
            },
            key=f"alias_manager_ignored_table_{config_path}",
        )
        selected_ignored_rows = _selected_dataframe_rows(ignored_table_state)
        selected_ignored = [
            ignored_df.iloc[row]["ignored_name"]
            for row in selected_ignored_rows
            if 0 <= row < len(ignored_df)
        ]
        if selected_ignored:
            st.caption("Selected: " + ", ".join(f"`{name}`" for name in selected_ignored))
        del_ignored_col, add_ignored_col = st.columns([1, 3])
        if del_ignored_col.button(
            f"Delete selected ({len(selected_ignored)})",
            disabled=not selected_ignored,
            key=f"delete_selected_ignored_{config_path}",
        ):
            updated = json.loads(json.dumps(config))
            selected_norms = {normalize_text(name) for name in selected_ignored}
            updated["ignored_names"] = [
                name
                for name in updated.get("ignored_names", [])
                if normalize_text(name) not in selected_norms
            ]
            backup = manager_save_alias_config(config_path, updated)
            st.success(
                f"Deleted {len(selected_ignored)} ignored name(s). "
                f"Backup: `{backup.name if backup else 'none'}`"
            )
            st.rerun()
        add_ignored_col.caption("Use the table below to add new ignored names, then save.")

        edited_ignored_new = st.data_editor(
            pd.DataFrame(columns=["ignored_name"]),
            num_rows="dynamic",
            width="stretch",
            hide_index=True,
            column_config={
                "ignored_name": st.column_config.TextColumn("Add ignored name"),
            },
            key=f"alias_manager_ignored_add_{config_path}",
        )
        edited_ignored = pd.concat(
            [
                ignored_df,
                edited_ignored_new,
            ],
            ignore_index=True,
        )

        st.divider()
        st.markdown("**Move ignored name to alias**")
        ignored_options = [row["ignored_name"] for row in ignored_rows]
        canonical_options = sorted(config.get("canonical_names", {}))
        c1, c2, c3 = st.columns([2, 2, 1])
        ignored_choice = c1.selectbox("Ignored name", ignored_options, key=f"move_ignored_name_{config_path}") if ignored_options else None
        canonical_choice = c2.selectbox("Canonical", canonical_options, key=f"move_ignored_canonical_{config_path}") if canonical_options else None
        if c3.button("Move", disabled=not (ignored_choice and canonical_choice), key=f"move_ignored_btn_{config_path}"):
            updated = move_ignored_to_alias(config, ignored_choice, canonical_choice)
            warnings = validate_alias_config(updated)
            if warnings:
                st.warning("\n".join(warnings))
            else:
                backup = manager_save_alias_config(config_path, updated)
                st.success(f"Moved `{ignored_choice}` to `{canonical_choice}`. Backup: {backup.name if backup else 'none'}")
                st.rerun()

    with tab_ambiguous:
        st.caption("Ambiguous names are known shared terms that require user decision per dataset.")
        ambiguous_df = pd.DataFrame(ambiguous_rows, columns=["ambiguous_name"])
        ambiguous_table_state = st.dataframe(
            ambiguous_df,
            width="stretch",
            hide_index=True,
            on_select="rerun",
            selection_mode="multi-row",
            column_config={
                "ambiguous_name": st.column_config.TextColumn("Ambiguous name"),
            },
            key=f"alias_manager_ambiguous_table_{config_path}",
        )
        selected_ambiguous_rows = _selected_dataframe_rows(ambiguous_table_state)
        selected_ambiguous = [
            ambiguous_df.iloc[row]["ambiguous_name"]
            for row in selected_ambiguous_rows
            if 0 <= row < len(ambiguous_df)
        ]
        if selected_ambiguous:
            st.caption("Selected: " + ", ".join(f"`{name}`" for name in selected_ambiguous))
        del_ambiguous_col, add_ambiguous_col = st.columns([1, 3])
        if del_ambiguous_col.button(
            f"Delete selected ({len(selected_ambiguous)})",
            disabled=not selected_ambiguous,
            key=f"delete_selected_ambiguous_{config_path}",
        ):
            updated = json.loads(json.dumps(config))
            selected_norms = {normalize_text(name) for name in selected_ambiguous}
            updated["ambiguous_names"] = [
                name
                for name in updated.get("ambiguous_names", [])
                if normalize_text(name) not in selected_norms
            ]
            backup = manager_save_alias_config(config_path, updated)
            st.success(
                f"Deleted {len(selected_ambiguous)} ambiguous name(s). "
                f"Backup: `{backup.name if backup else 'none'}`"
            )
            st.rerun()
        add_ambiguous_col.caption("Use the table below to add new ambiguous names, then save.")

        edited_ambiguous_new = st.data_editor(
            pd.DataFrame(columns=["ambiguous_name"]),
            num_rows="dynamic",
            width="stretch",
            hide_index=True,
            column_config={
                "ambiguous_name": st.column_config.TextColumn("Add ambiguous name"),
            },
            key=f"alias_manager_ambiguous_add_{config_path}",
        )
        edited_ambiguous = pd.concat(
            [
                ambiguous_df,
                edited_ambiguous_new,
            ],
            ignore_index=True,
        )

    with tab_raw:
        st.json(config)

    st.divider()
    updated_config = tables_to_alias_config(
        config,
        alias_records_for_save,
        edited_ignored.to_dict("records"),
        edited_ambiguous.to_dict("records"),
    )
    warnings = validate_alias_config(updated_config)
    if warnings:
        st.warning("Review before saving:\n\n" + "\n".join(f"- {w}" for w in warnings))

    save_col, reload_col = st.columns([3, 1])
    if save_col.button("Save alias config", type="primary", width="stretch", disabled=bool(warnings)):
        try:
            backup = manager_save_alias_config(config_path, updated_config)
            st.success(f"Saved `{config_path.name}`. Backup created: `{backup.name if backup else 'none'}`")
            st.rerun()
        except Exception as e:
            st.error(f"Could not save alias config: {e}")

    if reload_col.button("Reload", width="stretch"):
        st.rerun()


# ═══════════════════════════════════════════════════════════════════
# App entry point
# ═══════════════════════════════════════════════════════════════════

st.set_page_config(
    page_title="ViraLift",
    page_icon="V",
    layout="wide",
)

_init_state()
_inject_css()

# sidebar: always visible
with st.sidebar:
    st.title("ViraLift")
    st.caption(_t("sidebar_subtitle"))
    st.divider()

    st.radio(
        "Mode",
        options=["Run pipeline", "Alias manager"],
        key="app_mode",
    )
    st.divider()

    st.radio(
        _t("theme"),
        options=["dark", "light"],
        format_func=lambda value: _t(value),
        key="ui_theme",
        horizontal=True,
    )
    st.divider()

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
