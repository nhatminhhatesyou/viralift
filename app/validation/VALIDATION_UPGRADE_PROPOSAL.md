# Đề xuất update bộ notebook (tham khảo notebook cũ trong git)

Nguồn tham khảo (git blob, không nằm trên working tree):
- `fmd_prrsv_accuracy_summary.ipynb` @`b40188ec` — accuracy summary + bucket.
- `prrsv_tblastn_breakdown_strict_clean.ipynb` @`b10900fe` — breakdown PRRSV (granularity/ORF1ab, ORF7, ORF2b, frameshift).
- `tblastn_truth_breakdown_strict_clean.ipynb` @`7b1648e4` — adjudication FMD.
- `prrsv_start_rescue_full_comparison.ipynb` @`b40188ec` — ablation before/after.
- `alias_coverage.ipynb` @`68ecbd16` — coverage by field/feature-type + heatmap.

Ưu tiên từ cao xuống thấp.

## 1. Nâng cách chấm accuracy → `02_lifting_accuracy/*` + `_shared/validation_utils.py`

**Notebook cũ làm gì:** dùng `build_truth_presence` (mỗi record × mỗi ref-gene → `has_truth_gene`) rồi `summarize_per_gene` align theo **`ref_name`** (không phải pred_name), denominator = **số record có gene đó trong truth**. Buckets rời nhau: `exact / coord_only(=coord−exact) / failed(=total−coord)`, `accuracy_pct=(exact+coord_only)/total`. Có đếm riêng `no_hit` và `extra_predictions_without_truth`.

**Update:** thay `summarize_evaluable` (đang lọc `name_match` — prediction-based) bằng cặp helper `build_truth_presence` + `summarize_per_gene`/`summarize_overall` (truth-based, align `ref_name`). Chặt hơn vì:
- bắt được ca **tool không lift ra prediction nào** cho gene có trong truth (name_match bỏ sót),
- tách `no_hit` và `extra_predictions_without_truth` (ca ORF1a lift lên record ORF1ab) thành cột riêng thay vì âm thầm loại.

Port 3 hàm này vào `validation_utils.py` để mọi notebook dùng chung.

## 2. Thêm cột chẩn đoán vào bảng so sánh → `compare_predictions_to_truth`

**Cũ có:** `ref_len, pred_len, truth_len, pred_minus_ref_len, truth_minus_ref_len, pred_minus_truth_len, delta_start, delta_end`. Đây là xương sống để chẩn đoán lệch biên (convention vs tool).

**Update:** thêm các cột này vào output của `compare_predictions_to_truth` (hiện mới có iou/best_iou/failure_mode). Rẻ, và mọi notebook breakdown/ablation đều cần.

## 3. Failure attribution generic → `08_failure_attribution/`

**Cũ làm:** `adjudicate_failure` → bucket `likely_validation_artifact / likely_tool_error / manual_review`; `_fmd_final_reason` → `final_blame = ref_truth` (truth_feature_absent, ref/query boundary convention mismatch) vs tool_error. Dùng length-delta bucket để quyết.

**Update + phải sửa:** bản cũ **hardcode tên gene** (`pred in {'ORF2b','VP3','VP4'}`, nhánh riêng ORF1ab). Port sang `08` nhưng **làm generic**:
- `truth absent` (pred có, truth cùng tên không có) → suy từ `has_truth_gene`, không kể tên gene.
- `prediction≈ref nhưng truth khác` (boundary convention) → suy từ `pred_minus_ref_len≈0` & `delta_*` với truth, không hardcode.
- granularity (pred phủ union nhiều truth khác tên, vd ORF1a↔ORF1ab) → suy từ `best_overlap_name != pred_name` + overlap, không gọi tên ORF1ab.

## 4. Ablation theo cấu trúc baseline-vs-experiment → `07_ablation_runtime/rescue_ablation.ipynb`

**Cũ làm:** load baseline predictions → rerun với config thí nghiệm → merge per-gene `baseline_* / experiment_*` → cột `*_delta` → bảng `changed cases` + summary theo gene → chart baseline (bar) vs experiment (scatter).

**Update:** đổi `rescue_ablation` sang đúng khung này (load baseline → rerun rescue-OFF/ON → delta → changed cases). Vẫn cần wire đường rescue-OFF (đã ghi ở REBUILD_PLAN).

## 5. Surface "extra predictions without truth" → `02` + `08`

**Cũ:** bảng riêng các prediction mà gene không có same-name truth (chính là name-gap/granularity). Không tính vào accuracy nhưng **hiện ra** để phân tích.

**Update:** thêm bảng này vào `02` (đếm) và `08` (phân tích) — thay cho việc `summarize_evaluable` lặng lẽ loại chúng.

## 6. Alias coverage đầy đủ hơn → `01_alias_coverage/`

**Cũ có:** overall coverage, **coverage theo feature-type × field**, non-canonical names cần review, top raw names, chart outcome-mix + **coverage heatmap**. Bản hiện tại mỏng hơn (chỉ overall + stacked).

**Update:** thêm breakdown theo field/feature-type + heatmap + bảng top raw names + non-canonical để review (nối sang `03_alias_suggestion`).

## 7. Case-study evidence (cho Section 6 của bài) → `08` hoặc notebook case-study riêng

**Cũ có sẵn evidence sections:** ORF2b (truth gap), ORF1b (frameshift/start boundary), ORF1a↔ORF1ab (granularity), ORF7 (tool-side boundary), FMD boundary convention. Đây là phần interpretation Q1 thích. Gene-specific ở đây **chấp nhận được** vì là case study (diễn giải), không phải scoring — nhưng scoring (mục 1–3) thì phải generic.

---

### Thứ tự làm đề xuất
1 (accuracy scoring) → 2 (diagnostic cols) → 3 (failure adjudication generic) → 5 (extra-without-truth) → 4 (ablation) → 6 (alias) → 7 (case study).

Mục 1–2 chạm `validation_utils.py` (dùng chung) nên làm trước để các notebook sau kế thừa.
