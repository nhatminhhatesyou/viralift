# ViraLift Validation — Notebook Index

Tổ chức theo **output list của paper**: nhóm *descriptive* (mục 1–6) và nhóm *analysis*
(ablation + failure + runtime + phụ trợ). Mỗi folder single-purpose; output ghi vào
`<folder>/outputs/` (git-ignored, tái tạo được).

Nguyên tắc: notebook **gọi lại code tool thật** (`app/src/pipeline.run_pipeline`,
`lift_all_tblastn`, …) qua helper trong `_shared/validation_utils.py` — không chép lại
logic, không magic-number; phần truth + scoring tách riêng và **generic** (không hard-code
virus). Toggle `RUN_FULL=False` để smoke-test, `True` để chạy full 100 record.

## Descriptive

| # | Folder | Notebook | Trả lời |
|---|---|---|---|
| 1 | `01_alias_coverage/` | `alias_coverage.ipynb` | Alias resolve được bao nhiêu tên GenBank (PRRS/FMD/PED) |
| 2 | `02_lifting_accuracy/` | `fmdv_accuracy.ipynb`, `prrsv_accuracy.ipynb`, `pedv_accuracy.ipynb` | Độ chính xác tblastn lift — **mỗi virus 1 notebook** |
| 3 | `03_alias_suggestion/` | `alias_suggestion_tool_llm.ipynb` (+ scripts) | Alias suggestion accuracy: classifier (tool) + LLM |
| 4 | `04_engine_comparison/` | `engine_lifting_minimap_vs_tblastn.ipynb` | tblastn vs minimap2 |
| 5 | `05_tool_comparison/` | `gatu_compare.ipynb` | So với GATU / AEGIS |
| 6 | `06_run_tool_extract/` | `06a_mixed_records_extract.ipynb`, `06b_crawled_records_extract.ipynb` | Chạy tool: (a) mixed annotated+unannotated, (b) record tự crawl → trích gene |

## Analysis (ablation + failure + runtime)

| Folder | Notebook | Trả lời |
|---|---|---|
| `07_ablation_runtime/` | `rescue_ablation.ipynb`, `ablation_baselines_runtime.ipynb` | rescue on/off, without-alias, direct/tblastn/hybrid, runtime |
| `08_failure_attribution/` | `failure_attribution.ipynb` | Fail do tool vs ref-gap vs annotation-disagreement (generic) |
| `09_recovery_unannotated/` | `recovery_unannotated.ipynb` | Recovery accuracy khi strip annotation + record noAnno thật |
| `10_pipeline_end_to_end/` | `pipeline_end_to_end.ipynb` | Full pipeline (routing) cross-check |

Phụ trợ: `00_dataset/` (Table 1 dataset). Routing audit nằm trong `06a` (Part A).

## Định nghĩa accuracy (denominator)

Mẫu số = **gene ∈ ref ∩ có trong truth của query** (`summarize_evaluable`, lọc `name_match`).
Prediction không có same-name truth (name-gap/granularity, vd ORF1a lift lên record chỉ có
ORF1ab) **bị loại khỏi mẫu số**, không tính sai — phân tích riêng ở `08_failure_attribution/`.

## Archives

`_archive_tblastn_breakdown/`, `_archive_ped/`: output + notebook từ các lần chạy trước
(pre-refactor). Giữ để tham chiếu, **không phải nguồn số chính thức** — số chính thức lấy từ
notebook chạy lại ở trên.

## Chạy

```bash
pip install -e ".[validation]"     # pysam, pandas, matplotlib, jupyter
# + BLAST+ (tblastn, makeblastdb) và minimap2 (cho unit 04) trên PATH
```
