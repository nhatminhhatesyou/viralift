# FMD Breakdown Notes

Ghi chú ngày hôm nay cho phần `strict_clean` FMD validation.

## Kết Quả Tổng Quan

- Tổng FMD predictions: `1152`
- Exact match: `1088 / 1152 = 94.44%`
- Coordinate-correct: `1131 / 1152 = 98.18%`
- Non-exact predictions: `64 / 1152 = 5.56%`
- Failed coord/name theo raw scoring: `21 / 1152 = 1.82%`

## Final Adjudication Cho 64 Non-Exact Cases

- Do ref/truth annotation model: `28 / 64 = 43.75%`
- Do tool-side prediction/boundary: `36 / 64 = 56.25%`

Tính trên toàn bộ FMD predictions:

- Ref/truth artifact: `28 / 1152 = 2.43%`
- Tool-side issue: `36 / 1152 = 3.12%`

## Nhóm Do Ref/Truth

### `ref_query_boundary_convention_mismatch`: 25 cases

- Tool prediction match reference model, nhưng query truth dùng boundary khác.
- Không nên xem đây là lỗi tool trực tiếp.

Các pattern chính:

- `2A`: 15 cases
  - `ref_len = 54 bp`
  - `pred_len = 54 bp`
  - `truth_len = 48 bp`
  - Tool đang lift đúng theo ref, nhưng query truth annotate `2A` ngắn hơn 6 bp.
  - Vì `2A` rất ngắn, lệch 6 bp làm IoU tụt xuống `0.8889`.

- `VP1`: 9 cases
  - Nhiều case `pred_len = ref_len = 633 bp`
  - Query truth thường dài hơn vài bp, ví dụ `639 bp`.
  - Đây là boundary convention mismatch giữa ref và query truth.

- `VP2`: 1 case
  - `ref_len = pred_len = 654 bp`
  - `truth_len = 954 bp`
  - Query truth rõ ràng dùng feature boundary/convention khác.

### `truth_feature_absent`: 3 cases

- Query truth thiếu same-name feature nên không thể score công bằng.
- Gồm các case như `VP3`, `VP4`, và một case `VP1`.

## Nhóm Do Tool

### `n_terminal_truncation_12bp`: 15 cases

- Chủ yếu ở `Lpro` và `3Cpro`.
- Pattern điển hình:
  - `ref_len == truth_len`
  - `pred_len = ref_len - 12 bp`
  - `delta_start = +12`
  - `delta_end = 0`
- Ví dụ `Lpro`:
  - `ref_len = 603 bp`
  - `truth_len = 603 bp`
  - `pred_len = 591 bp`
- Kết luận:
  - Ref và truth đồng ý với nhau.
  - Prediction bị thiếu 12 bp ở đầu, tức thiếu 4 aa N-terminal.
  - Đây là tool-side boundary inference issue.

### `short_peptide_boundary_offset`: 9 cases

- Chủ yếu ở `2A`.
- Pattern:
  - Lệch boundary 3 bp ở peptide rất ngắn.
  - IoU vẫn khoảng `0.94`, nhưng do `2A` ngắn nên offset nhỏ bị phạt mạnh.
- Đây là boundary precision issue, không phải lift sai vùng lớn.

### `minor_vp1_boundary_offset`: 8 cases

- `VP1` có một số case không match cả ref và truth boundary.
- Offset nhỏ, IoU cao khoảng `0.99`.
- Tính là tool-side under strict scoring, nhưng severity thấp.

### Các Tool Issues Lẻ

- `no_hit`: 1 case
  - `2A` không lift được.
  - Đây là tool miss rõ.

- `major_c_terminal_truncation`: 1 case
  - `3A`
  - Start đúng nhưng end thiếu `111 bp`.
  - Prediction không match ref cũng không match truth.
  - Đây là lỗi tool rõ nhất trong FMD breakdown.

- `minor_c_terminal_truncation_3bp`: 1 case
  - `3B`
  - Prediction thiếu 3 bp ở cuối.
  - Minor boundary issue.

- `boundary_offset_matches_neither_ref_nor_truth`: 1 case
  - Mixed boundary case, currently counted as tool-side under strict scoring.

## Kết Luận FMD

- FMD tblastn localization nhìn chung rất tốt.
- Raw exact accuracy là `94.44%`.
- Coordinate-correct accuracy là `98.18%`.
- Nhiều raw failures không phải tool lift sai vùng mà do ref/query truth annotation boundary khác nhau.
- Tool-side issues chủ yếu là boundary precision, không phải localization failure lớn.
- Severe tool errors rõ nhất:
  - `2A no_hit`: 1 case
  - `3A` bị thiếu C-terminal `111 bp`: 1 case

Một câu report có thể dùng:

> For FMDV, tblastn lifting is highly reliable at the localization level. Most non-exact predictions are caused by reference-vs-query annotation boundary differences or small terminal boundary offsets. Clear severe tool-side failures are rare.

## Insight Quan Trọng Về Tool

FMD dùng `mat_peptide`, nên code không chạy codon validation/rescue:

- Không dùng `rescue_start_codon`
- Không dùng `rescue_stop_codon`

Vì vậy các case thiếu `12 bp` ở `Lpro` / `3Cpro` không phải do rescue codon.

Nguyên nhân khả dĩ:

- `tblastn` là local alignment.
- Nếu HSP không cover 4 aa đầu của protein, prediction sẽ thiếu:
  - `4 aa * 3 bp = 12 bp`
- Code hiện tại lấy boundary trực tiếp từ HSP span:
  - `pred_start = min HSP subject coordinate`
  - `pred_end = max HSP subject coordinate`
- Do đó terminal amino acids không align sẽ bị mất khỏi predicted feature.

## Hướng Cải Thiện Tool

### Ý tưởng: terminal extrapolation cho tblastn HSPs

Thay vì lấy HSP span làm full feature boundary, dùng query protein coordinates để suy ra terminal missing aa.

Ví dụ:

- `protein_length = 201 aa`
- HSP bắt đầu ở `query_start = 5`
- Missing N-terminal:
  - `5 - 1 = 4 aa`
  - `4 * 3 = 12 bp`
- Với strand `+`, extend prediction start upstream 12 bp.

Tương tự C-terminal:

- Nếu `hsp.query_end < protein_length`
- Missing C-terminal:
  - `protein_length - hsp.query_end`
  - extend downstream `missing_aa * 3 bp`

### Rule Nên Conservative

Chỉ extrapolate nếu:

- coverage cao, ví dụ `>= 0.90`
- terminal missing nhỏ, ví dụ `<= 10 aa`
- strand rõ ràng
- extension không vượt genome bounds
- không tạo tọa độ invalid

### Features Có Thể Được Cải Thiện

- `Lpro`: thiếu 12 bp N-terminal
- `3Cpro`: thiếu 12 bp N-terminal
- `3B`: thiếu 3 bp C-terminal
- Một phần `2A`: boundary offset 3-6 bp

## Việc Nên Làm Tiếp

- Implement thử terminal extrapolation trong một branch riêng.
- Chạy lại `strict_clean` FMD validation.
- So sánh:
  - exact_pct
  - coord_pct
  - số `n_terminal_truncation_12bp`
  - số `short_peptide_boundary_offset`
  - kiểm tra có làm xấu PRRS không.
- Nếu sợ ảnh hưởng PRRS gene overlap, có thể bật extrapolation trước cho `mat_peptide` only.
