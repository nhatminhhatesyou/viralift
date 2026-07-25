# ViraLift — Validation Notebook Rebuild Plan (Q1 paper)

> **Index chuẩn = `app/validation/README.md`** (đã tổ chức lại theo output-list: descriptive
> 1–6 + analysis). Folder đánh số lại, lifting accuracy tách mỗi virus 1 notebook, denominator
> dùng `summarize_evaluable` (lọc `name_match`). File này giữ để tra lịch sử gap/quyết định;
> tên NB-xx bên dưới là số cũ trước khi reorg.

Bối cảnh: toàn bộ notebook trong `app/validation/**` chưa từng được commit (`.gitignore` từng ignore cả cây, notebook local đã mất) — chỉ còn output archive. Để đạt reproducibility cho Q1, mỗi figure/table phải tái tạo được từ notebook. Các notebook dưới đây **scaffold khung chuẩn** (Purpose → Inputs → Run module thật → Metrics → Figures → Interpretation) và **gọi trực tiếp** helper cấp cao trong `app/validation/_shared/validation_utils.py` + module trong `app/src`.

## Đã scaffold (mới)

| Notebook | Unit | Tái tạo | Paper |
|---|---|---|---|
| `00_dataset_curation.ipynb` | 00 | dataset_table, gene_presence | Table 1 |
| `engine_lifting_minimap_vs_tblastn.ipynb` ♻️ *restored from git* | 01 | engine_summary, engine_accuracy | 5.2 baseline |
| `02_alias_coverage.ipynb` | 02 | alias_summary_overall, alias_coverage | 5.4 |
| `03_end_to_end.ipynb` | 03 | e2e_summary, e2e_accuracy | 5.1 |
| `04a_accuracy.ipynb` | 04 | overall/per_gene accuracy | Fig 3 / Table 2 |
| `04b_ablation.ipynb` | 04 | start-rescue + terminal-extrapolation | Fig 5 |
| `04c_failure_attribution.ipynb` | 04 | failure_attribution | Fig 6 |
| `05_gatu_compare.ipynb` | 05 | ViraLift vs GATU | 5.2 baseline |
| `08_baselines_and_runtime.ipynb` | 08 (mới) | without-alias, direct/tblastn/hybrid, runtime | 5.2/5.3 gaps |
| `09_unannotated_recovery.ipynb` | 09 (mới) | strip-recovery accuracy + noAnno demo | 5.x use case chính |
| `10_routing_audit_orf5_demo.ipynb` | 10 (mới) | routing crosstab + ORF5 coverage | novelty routing + demo |

| `07_llm_alias_validation.ipynb` | 07 (mới) | wrap script → track A/B accuracy | 5.x LLM review |

Đã có sẵn, chỉ cần **rerun**: `06_ped_validation/*.ipynb` (3 notebook — PEDV, 5.5 + ablation pre/post). Logic LLM vẫn nằm ở `07_llm_alias_validation/*.py` (notebook chỉ wrap để render figure; real run cần `OPENAI_API_KEY`).

## Thứ tự chạy đề xuất

1. **06** (rerun, gần như free) — quick win, PEDV.
2. **04a** — con số accuracy chính của bài (nặng nhất).
3. **01** → **03** — baseline engine + end-to-end.
4. **02** — alias coverage (nhẹ).
5. **04b** + **08** — ablation & gap.
6. **05** — GATU compare.
7. **04c** — failure attribution.
8. **07** — LLM review (nếu cần).

Mọi notebook có toggle `RUN_FULL = False` (smoke test 10 records) → đổi `True` cho full 100.

## Gap còn phải bịt (đánh dấu `TODO` trong notebook)

- **04b**: `run_tblastn_against_truth` chạy rescue ON mặc định — cần thêm đường **rescue OFF** (kwarg qua `run_tblastn_batch`/`process_one_query_record`, hoặc bỏ `extrapolate_terminal_boundaries` + start-rescue) để so sánh. Output `terminal_extrapolation` cũ **đang mất**, phải regenerate.
- ~~**01**: cạnh minimap2 chưa nối~~ ✅ **Đã khôi phục từ git** (commit `68ecbd1`): notebook gốc `engine_lifting_minimap_vs_tblastn.ipynb` (1462 dòng, có sẵn logic minimap-vs-tblastn) + `app/src/alignment/minimap_runner.py` (`run_minimap2_with_fallback`), `app/src/alignment/sam_lifter.py`, `app/src/annotation/extractor.py`. Cả 3 chỉ phụ thuộc `app/src/lifting/validator.py` (đã có). Yêu cầu: `pip install -e ".[validation]"` (thêm `pysam`) + binary `minimap2` trên PATH.
- **05**: cần chạy GATU bên ngoài, lưu `gatu_results.tsv` rồi join.
- **08**: without-alias (empty lookup) + direct-only cần force `strategy`; runtime đo bằng `time.perf_counter`.

## Faithfulness fix (validation phải chạy đúng code tool)

- ✅ `run_production_pipeline_against_truth` giờ gọi thẳng `app/src/pipeline.run_pipeline` với `PipelineConfig()` default — trước đây tự dựng lại routing và truyền `rescue_window=50` (lệch default **200** của tool). Giờ NB-03 đo đúng pipeline ship ra.
- ✅ `run_tblastn_against_truth` dùng `PipelineConfig()` tường minh (coverage/identity/evalue/rescue_window) thay vì default ngầm.
- ✅ Bỏ hardcode `if virus == PRRS: add ORF1ab` trong `truth_target_names` + `should_use_target_truth_filter`. Mẫu số accuracy = **R ∩ query-truth** (gene có trong ref ∩ có trong truth); ca granularity (ORF1ab phủ ORF1a/ORF1b) để NB-04c phát hiện **generic bằng overlap**, không nhắc tên virus/gene.
- ⏳ **Còn lại (scoring layer)**: NB-04a đang tính `summarize_comparison(mọi prediction)` → vẫn phạt oan ca name-gap. Cần đổi mẫu số sang truth-based R∩T (giống hàm `summarize_truth_available` của NB-06). Chưa chốt: prediction-based lọc `name_match` hay truth-based join đầy đủ.

## Cần verify trước khi chạy

- Ref records đang trỏ `FMD_ref_test.gb` / `PRRS_ref_test.gb` / `PED_ref_1|2.gb` — xác nhận đúng ref dùng cho bài (có bản `*_Anno.gb` thay thế).
- Env chạy notebook cần: Biopython, BLAST+ (`tblastn`, `makeblastdb`) trong PATH — dùng `.venv` của repo.
