# Auto-split from the original monolithic streamlit_app.py.
import pandas as pd
import streamlit as st
from typing import Dict, List, Optional
from app.src.io.result_writer import summarize_counts
from app.src.lifting.base import LiftedFeature
from ui.components import _error_by_record, _format_fraction_percent, _format_percent, _method_section_label, _ordered_methods, _render_page_intro, _status_text
from ui.i18n import _t
from ui.services import _canonical_to_ref_map, _results_to_df
from ui.state import GOOD_STATUSES, REVIEW_STATUSES, _reset


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

    # ── header ───────────────────────────────────────────────────
    _render_page_intro(
        _t("run_review"),
        _t("results_title"),
        _t("results_body"),
    )
    if st.session_state.use_ref_names:
        st.caption(_t("ref_names_caption"))
    st.button(_t("new_run"), key="results_new_run", on_click=_reset)

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
            "passed": ok_count,
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
    visible_results = [
        (query_id, features)
        for query_id, features in all_results
        if query_id in visible_ids
    ]

    def _record_method_bucket(query_id: str, features: List[LiftedFeature]) -> str:
        if query_id in error_map or not features:
            return "other"
        methods = {feature.method for feature in features if feature.method}
        if methods == {"direct"}:
            return "direct"
        if methods == {"tblastn"}:
            return "tblastn"
        if "direct" in methods and "tblastn" in methods:
            return "mixed"
        return "other"

    def _boundary_check(feature: LiftedFeature) -> str:
        checks = []
        for label, value in (
            ("start", feature.has_start_codon),
            ("stop", feature.has_stop_codon),
            ("frame", feature.in_frame),
        ):
            if value is not None:
                checks.append(f"{label}:{'yes' if value else 'no'}")
        return ", ".join(checks)

    def _rescue_action(feature: LiftedFeature) -> str:
        if feature.rescue_action:
            return feature.rescue_action
        if feature.rescue_target and feature.rescue_offset is not None:
            return f"{feature.rescue_target} {feature.rescue_offset:+d} bp"
        if feature.rescue_offset is not None:
            return f"start {feature.rescue_offset:+d} bp"
        return ""

    def _display_gene(feature: LiftedFeature) -> str:
        return ref_name_map.get(feature.name, feature.name) if ref_name_map else feature.name

    def _export_key(query_id: str, feature: LiftedFeature, gene: str) -> str:
        return "|".join(
            str(value)
            for value in (
                query_id,
                gene,
                feature.status,
                feature.query_start,
                feature.query_end,
                feature.strand,
            )
        )

    manually_approvable_statuses = {"invalid_boundaries", "low_coverage"}

    def _approval_label(query_id: str, feature: LiftedFeature, gene: str) -> str:
        return (
            f"{query_id} | {gene} | {_t(f'status_{feature.status}')} | "
            f"{feature.query_start}-{feature.query_end} | "
            f"cov {_format_fraction_percent(feature.coverage)} | "
            f"id {_format_percent(feature.identity)}"
        )

    def _approved_export_keys() -> set:
        return set(st.session_state.get("approved_export_features", []))

    def _save_approved_export_keys(keys: set) -> None:
        st.session_state.approved_export_features = sorted(keys)

    def _render_record_expander(query_id: str, features: List[LiftedFeature]) -> None:
        ok_count     = sum(1 for f in features if f.status in GOOD_STATUSES)
        total_count  = len(features)
        needs_review = sum(1 for f in features if f.status in REVIEW_STATUSES)

        if query_id in error_map:
            label = f"{_t('failed')} · {query_id}"
        elif total_count == 0:
            label = f"{_t('empty')} · {query_id}"
        elif needs_review:
            label = f"{_t('needs_review')} · {query_id} · {ok_count}/{total_count} passed"
        else:
            label = f"{_t('passed')} · {query_id} · {ok_count}/{total_count} passed"

        with st.expander(label, expanded=show_details or query_id in error_map):
            if query_id in error_map:
                st.error(error_map[query_id])
                return
            if not features:
                st.warning(_t("no_features"))
                return
            rec_rows = []
            for lf in features:
                display_name = _display_gene(lf)
                rec_rows.append({
                    "status":    _status_text(lf.status),
                    "name":      display_name,
                    "raw_name":  lf.source_name or "",
                    "start":     lf.query_start,
                    "end":       lf.query_end,
                    "coverage":  _format_fraction_percent(lf.coverage),
                    "identity":  _format_percent(lf.identity),
                    "method":    lf.method,
                    "boundary_check": _boundary_check(lf),
                    "rescue_action": _rescue_action(lf),
                    "approved_for_fasta": (
                        "yes"
                        if _export_key(query_id, lf, display_name) in _approved_export_keys()
                        else ""
                    ),
                })
            rec_df = pd.DataFrame(rec_rows)
            for method in _ordered_methods(rec_df["method"].dropna().unique().tolist()):
                method_df = rec_df[rec_df["method"] == method]
                if method_df.empty:
                    continue

                st.markdown(f"##### {_method_section_label(method)} ({len(method_df)})")
                if method == "tblastn" and method_df["status"].str.contains("Invalid boundary", na=False).any():
                    st.caption(
                        "Invalid boundary means tblastn found a strong protein match, "
                        "but the extracted CDS did not pass start/stop/frame checks. "
                        "In boundary_check, start means ATG at the beginning, "
                        "stop means TAA/TAG/TGA at the end, and frame means the "
                        "CDS length is divisible by 3. rescue_action shows which "
                        "boundary was moved, for example `start -12 bp` or `stop +6 bp`."
                    )
                review_df = method_df[
                    method_df["status"].str.startswith(f"{_t('tone_review')} ·", na=False)
                ]
                pass_df = method_df[
                    method_df["status"].str.startswith(f"{_t('tone_pass')} ·", na=False)
                ]
                other_df = method_df.drop(review_df.index.union(pass_df.index))

                if not review_df.empty:
                    st.caption(_t("needs_review"))
                    review_candidates = [
                        lf for lf in features
                        if lf.method == method
                        and lf.status in manually_approvable_statuses
                        and lf.sequence
                    ]
                    if review_candidates:
                        approved_keys = _approved_export_keys()
                        selected_keys = set()
                        with st.container(border=True):
                            for idx, lf in enumerate(review_candidates):
                                gene = _display_gene(lf)
                                key = _export_key(query_id, lf, gene)
                                checked = st.checkbox(
                                    _approval_label(query_id, lf, gene),
                                    value=key in approved_keys,
                                    key=(
                                        f"approve_export_{query_id}_{method}_{idx}_"
                                        f"{lf.status}_{lf.query_start}_{lf.query_end}"
                                    ),
                                )
                                if checked:
                                    selected_keys.add(key)

                        approve_col, clear_col = st.columns([2, 1])
                        if approve_col.button(
                            _t("approve_for_fasta"),
                            key=f"approve_selected_export_{query_id}_{method}",
                            width="stretch",
                        ):
                            record_candidate_keys = {
                                _export_key(query_id, lf, _display_gene(lf))
                                for lf in review_candidates
                            }
                            approved_keys.difference_update(record_candidate_keys)
                            approved_keys.update(selected_keys)
                            _save_approved_export_keys(approved_keys)
                            st.rerun()
                        record_candidate_keys = {
                            _export_key(query_id, lf, _display_gene(lf))
                            for lf in review_candidates
                        }
                        if clear_col.button(
                            _t("clear_approved_for_record"),
                            key=f"clear_approved_export_{query_id}_{method}",
                            width="stretch",
                            disabled=not (approved_keys & record_candidate_keys),
                        ):
                            approved_keys.difference_update(record_candidate_keys)
                            _save_approved_export_keys(approved_keys)
                            st.rerun()
                    st.dataframe(review_df, width="stretch", hide_index=True)
                if not pass_df.empty:
                    st.caption(_t("passed"))
                    st.dataframe(pass_df, width="stretch", hide_index=True)
                if not other_df.empty:
                    st.caption("Other")
                    st.dataframe(other_df, width="stretch", hide_index=True)

    detail_sections = [
        ("direct", "Direct extraction records"),
        ("tblastn", "tblastn lifting records"),
        ("mixed", "Mixed direct + tblastn records"),
        ("other", "Other / errors / empty records"),
    ]
    for bucket, title in detail_sections:
        bucket_results = [
            (query_id, features)
            for query_id, features in visible_results
            if _record_method_bucket(query_id, features) == bucket
        ]
        st.markdown(f"### {title} ({len(bucket_results)})")
        if not bucket_results:
            st.caption("No records in this section for the current filter.")
            continue
        for query_id, features in bucket_results:
            _render_record_expander(query_id, features)

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
        qf_col1, qf_col2, qf_col3, qf_col4 = st.columns(4)
        min_cov_export = qf_col1.slider(_t("min_coverage"), 0.0, 1.0, 0.5, 0.05)
        min_id_export  = qf_col2.slider(_t("min_identity"),  0.0, 100.0, 0.0, 5.0)
        include_rescued = qf_col3.checkbox(_t("include_rescued"), value=True)
        include_extrapolated = qf_col4.checkbox(_t("include_extrapolated"), value=True)

        accepted_statuses = {"ok", "direct"}
        if include_rescued:
            accepted_statuses.add("ok_rescued")
        if include_extrapolated:
            accepted_statuses.add("ok_extrapolated")

        approved_review_keys = _approved_export_keys()
        if approved_review_keys:
            st.caption(_t("approved_export_count", count=len(approved_review_keys)))

        def _passes_fasta_filters(query_id: str, feature: LiftedFeature, gene: str) -> bool:
            if gene not in selected_genes or not feature.sequence:
                return False
            if _export_key(query_id, feature, gene) in approved_review_keys:
                return True
            if feature.status not in accepted_statuses:
                return False
            if feature.coverage is not None and feature.coverage < min_cov_export:
                return False
            if feature.identity is not None and feature.identity < min_id_export:
                return False
            return True

        candidate_count = 0
        for query_id, features in all_results:
            for lf in features:
                gene = _display_gene(lf)
                if _passes_fasta_filters(query_id, lf, gene):
                    candidate_count += 1
        st.caption(_t("candidate_count", count=candidate_count))

        if st.button(_t("generate_fasta"), type="primary"):
            # build gene to sequences dict (keyed by display name)
            gene_seqs: Dict[str, List[str]] = {g: [] for g in selected_genes}
            skipped = 0

            for query_id, features in all_results:
                for lf in features:
                    # resolve display name the same way as the gene list above
                    gene = _display_gene(lf)
                    if gene not in gene_seqs:
                        continue
                    if not _passes_fasta_filters(query_id, lf, gene):
                        skipped += 1
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
