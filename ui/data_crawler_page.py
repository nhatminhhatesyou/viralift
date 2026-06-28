from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import streamlit as st

from gbcrawler import fetch as fetch_mod
from gbcrawler import split as split_mod
from ui.components import _render_context_panel, _render_page_intro
from ui.state import REGISTRY_PATH, ROOT


DEFAULT_OUT_DIR = ROOT / "output" / "gbcrawler"
DEFAULT_LEDGER = ROOT / "output" / "gbcrawler_ledger.txt"


FALLBACK_VIRUSES = {
    "PRRSV": "Porcine reproductive and respiratory syndrome virus",
    "FMDV": "Foot-and-mouth disease virus",
    "PEDV": "Porcine epidemic diarrhea virus",
}


def _virus_options(registry_path: Path) -> Dict[str, str]:
    try:
        registry = split_mod.load_registry(registry_path)
    except Exception:
        return dict(FALLBACK_VIRUSES)

    options = {}
    for name, keywords in registry:
        organism = max(keywords, key=len) if keywords else name
        options[name] = organism.strip()
    return options or dict(FALLBACK_VIRUSES)


def _parse_accessions(raw: str) -> List[str]:
    accessions = []
    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        accessions.extend(part.strip() for part in line.replace(",", " ").split())
    return [acc for acc in accessions if acc]


def _dataset_type_clause(dataset_type: str) -> Optional[str]:
    if dataset_type == "Complete genomes":
        return "complete genome[Title]"
    if dataset_type == "Gene records / fragments":
        return "NOT complete genome[Title]"
    return None


def build_ncbi_query(
    organism: str,
    dataset_type: str,
    length: Optional[Tuple[int, int]],
    dates,
    region: str,
    extra: str,
) -> str:
    parts = []
    if organism:
        parts.append(f'"{organism}"[Organism]')

    dataset_clause = _dataset_type_clause(dataset_type)
    if dataset_clause:
        parts.append(dataset_clause)

    if length:
        lo, hi = length
        parts.append(f'("{lo}"[SLEN] : "{hi}"[SLEN])')

    if dates:
        d0, d1 = dates
        parts.append(f'("{d0:%Y/%m/%d}"[PDAT] : "{d1:%Y/%m/%d}"[PDAT])')

    if region:
        parts.append(f'"{region}"[All Fields]')

    if extra and extra.strip():
        parts.append(extra.strip())

    return " AND ".join(parts)


def _need_email(email: str) -> bool:
    if not email.strip():
        st.error("Enter an email address first. NCBI requires it for Entrez requests.")
        return True
    return False


def _crawl_signature(mode: str, query: str, accessions: List[str], db: str, retmax: int) -> str:
    payload = {
        "mode": mode,
        "query": query.strip(),
        "accessions": sorted(accessions),
        "db": db,
        "retmax": int(retmax),
    }
    raw = json.dumps(payload, sort_keys=True, ensure_ascii=True)
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()


def _use_species_file_as_query(file_path: Path, label: str) -> None:
    st.session_state.crawler_query_path = str(file_path)
    st.session_state.crawler_query_label = label
    st.session_state.app_mode = "Run pipeline"
    st.session_state.stage = "upload"
    st.toast(f"Selected `{file_path.name}` as query input.")


def page_data_crawler() -> None:
    _render_page_intro(
        "Data preparation",
        "Crawl GenBank records from NCBI",
        (
            "Build a GenBank dataset from NCBI, split it by virus, then send "
            "one species file into the ViraLift pipeline as query input."
        ),
        show_stages=False,
    )

    registry_path = REGISTRY_PATH
    viruses = _virus_options(registry_path)

    email = st.session_state.get("crawler_email_value") or st.session_state.get("crawler_email_input", "")
    api_key = st.session_state.get("crawler_api_key_value") or st.session_state.get("crawler_api_key_input", "")

    with st.expander("Advanced NCBI settings", expanded=False):
        db = st.text_input("Database", value="nuccore", key="crawler_db")
        out_dir = st.text_input(
            "Output folder",
            value=str(DEFAULT_OUT_DIR),
            key="crawler_out_dir_input",
        )
        use_ledger = st.checkbox("Use dedup ledger", value=True, key="crawler_use_ledger")
        ledger_path = st.text_input(
            "Ledger file",
            value=str(DEFAULT_LEDGER),
            disabled=not use_ledger,
            key="crawler_ledger_path",
        )

    st.subheader("1. Build input dataset")
    mode = st.radio(
        "Input source",
        ["By NCBI query", "By accession list"],
        horizontal=True,
        label_visibility="collapsed",
        key="crawler_mode",
    )

    query = ""
    accessions: List[str] = []

    with st.container(border=True):
        if mode == "By NCBI query":
            top_left, top_right = st.columns([2, 1])
            with top_left:
                labels = list(viruses.keys()) + ["Other / custom"]
                selected_virus = st.selectbox("Virus", labels, key="crawler_virus")
                if selected_virus == "Other / custom":
                    organism = st.text_input(
                        "NCBI organism name",
                        placeholder="e.g. Classical swine fever virus",
                        key="crawler_custom_organism",
                    )
                else:
                    organism = viruses[selected_virus]
                    st.caption(f'NCBI organism query: `"{organism}"[Organism]`')
            with top_right:
                dataset_type = st.selectbox(
                    "Dataset type",
                    ["Complete genomes", "Gene records / fragments", "Both"],
                    key="crawler_dataset_type",
                )

            st.markdown("**Optional filters**")
            f1, f2, f3 = st.columns(3)
            with f1:
                use_len = st.checkbox("Length range", key="crawler_use_length")
                length = None
                if use_len:
                    length = st.slider(
                        "Sequence length",
                        0,
                        35000,
                        (14000, 16000),
                        step=500,
                        key="crawler_length",
                    )
            with f2:
                use_date = st.checkbox("Publication date", key="crawler_use_date")
                dates = None
                if use_date:
                    picked_dates = st.date_input(
                        "Date range",
                        value=(date(2020, 1, 1), date.today()),
                        max_value=date.today(),
                        key="crawler_dates",
                    )
                    if isinstance(picked_dates, (tuple, list)) and len(picked_dates) == 2:
                        dates = picked_dates
            with f3:
                region_choice = st.selectbox(
                    "Country",
                    ["Any (global)", "Vietnam", "Custom country…"],
                    key="crawler_region_choice",
                )
                region = ""
                if region_choice == "Vietnam":
                    region = "Vietnam"
                elif region_choice == "Custom country…":
                    region = st.text_input("Country name", placeholder="e.g. China")
                if region:
                    st.caption(
                        "⚠️ Country filter is text-based and approximate on nuccore. "
                        "For exact geographic subsets, filter by the /country "
                        "qualifier after download (or use NCBI Virus)."
                    )

            extra = st.text_input(
                "Extra NCBI condition",
                placeholder='e.g. NOT partial[Title]',
                key="crawler_extra",
            )
            query = build_ncbi_query(organism, dataset_type, length, dates, region, extra)

            with st.expander("Advanced: edit query manually"):
                use_manual_query = st.checkbox(
                    "Use manual query override",
                    value=False,
                    key="crawler_use_manual_query",
                    help="Leave this off to use the query built from the virus and filters above.",
                )
                manual_query = st.text_area(
                    "Override query",
                    value=query,
                    height=90,
                    key="crawler_manual_query",
                    disabled=not use_manual_query,
                )
                if use_manual_query and manual_query.strip():
                    query = manual_query.strip()

            st.markdown("**NCBI query preview**")
            st.code(query or "No query yet.", language="text")
        else:
            uploaded = st.file_uploader("Accession file (.txt/.csv)", type=["txt", "csv"])
            pasted = st.text_area(
                "Or paste accessions",
                height=140,
                placeholder="PP209408.1\nAF176348.2",
                key="crawler_accessions_text",
            )
            raw = pasted
            if uploaded is not None:
                raw = uploaded.getvalue().decode("utf-8", "ignore") + "\n" + raw
            accessions = _parse_accessions(raw)
            st.caption(f"{len(accessions)} accession(s) detected.")

        retmax = st.number_input(
            "Max records to fetch",
            min_value=1,
            max_value=100000,
            value=1000,
            step=100,
            key="crawler_retmax",
        )

    current_signature = _crawl_signature(mode, query, accessions, db, int(retmax))

    st.subheader("2. Preview and crawl")
    count_col, crawl_col = st.columns(2)
    count_disabled = mode != "By NCBI query"
    if count_col.button("Count records", disabled=count_disabled, width="stretch"):
        if not _need_email(email):
            if not query.strip():
                st.error("Build or enter an NCBI query first.")
            else:
                fetch_mod.configure(email, api_key or None)
                with st.spinner("Counting matching records on NCBI..."):
                    try:
                        st.session_state.crawler_count_result = fetch_mod.count_query(
                            query,
                            db=db,
                            api_key=api_key or None,
                        )
                        st.session_state.crawler_count_signature = current_signature
                    except Exception as exc:
                        st.session_state.crawler_count_result = None
                        st.session_state.crawler_count_signature = None
                        st.error(f"NCBI count failed: {exc}")

    if (
        st.session_state.crawler_count_result is not None
        and st.session_state.get("crawler_count_signature") == current_signature
    ):
        count = st.session_state.crawler_count_result
        tone = st.warning if count > 50000 else st.success
        tone(f"{count:,} record(s) match this query.")

    if crawl_col.button("Crawl GenBank and split", type="primary", width="stretch"):
        if not _need_email(email):
            if mode == "By NCBI query" and not query.strip():
                st.error("Build or enter an NCBI query first.")
            elif mode == "By accession list" and not accessions:
                st.error("Provide at least one accession.")
            else:
                out = Path(out_dir)
                out.mkdir(parents=True, exist_ok=True)
                raw_path = out / "raw_combined.gb"
                try:
                    fetch_mod.configure(email, api_key or None)
                    progress = st.progress(0.0, text="Starting NCBI download...")
                    with open(raw_path, "w", encoding="utf-8") as handle:
                        if mode == "By NCBI query":
                            fetched = fetch_mod.fetch_from_query(
                                query,
                                handle,
                                db=db,
                                retmax=int(retmax),
                                api_key=api_key or None,
                                quiet=True,
                            )
                        else:
                            fetched = fetch_mod.fetch_from_accessions(
                                accessions,
                                handle,
                                db=db,
                                api_key=api_key or None,
                                quiet=True,
                            )
                    progress.progress(0.75, text=f"Downloaded {fetched} record(s). Splitting...")

                    registry = split_mod.load_registry(registry_path)
                    seen = (
                        split_mod.load_ledger(Path(ledger_path))
                        if use_ledger and ledger_path
                        else set()
                    )
                    summary = split_mod.split_genbank(raw_path, registry, out, seen=seen)
                    if use_ledger and ledger_path:
                        split_mod.append_ledger(Path(ledger_path), summary["new_accessions"])
                    progress.progress(1.0, text="Done.")
                    st.session_state.crawler_summary = summary
                    st.session_state.crawler_out_dir = str(out)
                    st.session_state.crawler_summary_signature = current_signature
                except Exception as exc:
                    st.error(f"Crawl failed: {exc}")

    summary = st.session_state.crawler_summary
    if not summary:
        return
    if st.session_state.get("crawler_summary_signature") != current_signature:
        st.info("The crawler output below belongs to a previous query. Run `Crawl GenBank and split` again for the current query.")
        return

    st.subheader("3. Use crawler output")
    _render_context_panel([
        ("Total fetched", summary["total"]),
        ("New records", summary["new"]),
        ("Duplicate in batch", summary["dup_in_batch"]),
        ("Duplicate in ledger", summary["dup_in_ledger"]),
    ])

    if summary["counts"]:
        st.markdown("**Species files**")
        for key, count in sorted(summary["counts"].items()):
            file_path = Path(summary["files"][key])
            label = key if key != split_mod.UNMATCHED else "Unmatched records"
            c1, c2, c3 = st.columns([3, 1, 1])
            c1.write(f"**{label}** · {count} record(s) · `{file_path.name}`")
            if file_path.exists():
                c2.download_button(
                    "Download .gb",
                    data=file_path.read_bytes(),
                    file_name=file_path.name,
                    key=f"crawler_download_{file_path.name}",
                    width="stretch",
                )
                c3.button(
                    "Use as query",
                    key=f"crawler_use_{file_path.name}",
                    width="stretch",
                    on_click=_use_species_file_as_query,
                    args=(file_path, label),
                )

    manifest_path = Path(summary["manifest"])
    if manifest_path.exists():
        with st.expander("manifest.csv", expanded=False):
            manifest = pd.read_csv(manifest_path)
            st.dataframe(manifest, width="stretch", height=320)
            st.download_button(
                "Download manifest.csv",
                data=manifest_path.read_bytes(),
                file_name="manifest.csv",
            )
