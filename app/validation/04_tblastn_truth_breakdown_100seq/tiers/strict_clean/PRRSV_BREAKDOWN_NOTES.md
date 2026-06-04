# PRRSV Breakdown Notes

Ghi chú cho phần `strict_clean` PRRSV validation, gồm baseline và kết quả sau khi thử cải thiện start-codon rescue.

## Cách Tính Accuracy

PRRSV annotation không đồng nhất giữa các record, nên không nên tính accuracy trên toàn bộ prediction rows một cách mù.

Metric chính dùng trong breakdown:

- `truth_available_records`: số record có same-gene truth cho gene đó.
- `exact_correct`: pred name, start, end khớp hoàn toàn với truth.
- `coord_correct`: tọa độ overlap tốt với same-gene truth theo IoU threshold.
- `extra_predictions_without_truth`: tool predict gene nhưng query truth không annotate gene đó; không tính vào accuracy chính.

Ví dụ `ORF2b`:

- Reference có `ORF2b`.
- Tool predict `ORF2b` cho 95 records.
- Nhưng truth chỉ annotate `ORF2b` trong 48 records.
- Accuracy chính của `ORF2b` phải tính trên 48 records đó.

## Baseline: Per-Gene Accuracy

Baseline là PRRSV strict-clean validation trước khi sửa start rescue.

| Gene | Truth records | Exact | Coord | Nhận xét |
|---|---:|---:|---:|---|
| ORF2a | 95 | 95 | 95 | Perfect |
| ORF3 | 95 | 95 | 95 | Perfect |
| ORF2b | 48 | 45 | 48 | Tốt; 3 case chỉ lệch boundary nhỏ |
| ORF4 | 95 | 94 | 95 | 1 case boundary offset |
| ORF5 | 95 | 91 | 95 | 4 case boundary offset |
| ORF6 | 94 | 91 | 94 | 3 case boundary offset |
| ORF7 | 94 | 66 | 66 | 28 case wrong coords do start bị lệch |
| ORF1a | 62 | 61 | 62 | Tốt; 1 case boundary offset |
| ORF1b | 83 | 0 | 83 | Coord đúng hết, exact fail do start/frameshift convention |

## Finding Theo Gene Ở Baseline

### `ORF2a` và `ORF3`

- `95/95 exact`
- `95/95 coord`
- Không có vấn đề đáng kể.

### `ORF2b`

- `45/48 exact`
- `48/48 coord`
- 3 non-exact cases vẫn đúng tọa độ, chỉ lệch boundary nhỏ.
- 47 predictions còn lại không có same-gene truth, nên để riêng là `extra_predictions_without_truth`, không tính vào accuracy chính.

Kết luận:

- Tool lift `ORF2b` tốt trên các records có truth để validate.
- Các case không có truth chủ yếu phản ánh annotation gap/convention của PRRSV overlap gene.

### `ORF4`, `ORF5`, `ORF6`

- Coordinate accuracy đều `100%`.
- Exact fail là boundary offset nhỏ.
- Đây là nhóm có thể cải thiện bằng rescue chọn start tốt hơn.

### `ORF1a`

- `61/62 exact`
- `62/62 coord`
- Gần như ổn.
- Một số record PRRSV dùng convention `ORF1ab` thay vì tách `ORF1a/ORF1b`; nhóm này cần phân tích annotation convention riêng.

### `ORF1b`

- `0/83 exact`
- `83/83 coord`
- Pattern chính:
  - `delta_end = 0` nhiều/ổn định.
  - Start boundary lệch.
- Không nên xem là localization failure.

Kết luận:

- Tool tìm đúng vùng `ORF1b`.
- Exact fail chủ yếu do start/frameshift boundary convention, không phải do tblastn lift sai vùng.

### `ORF7`

Baseline là vấn đề rõ nhất:

- `66/94 exact`
- `66/94 coord`
- 28 failed cases đều có `delta_end = 0`, tức end boundary đúng.
- Pattern:
  - 25 cases: `truth_len = ref_len = 372`, nhưng `pred_len = 329`
  - 3 cases: `truth_len = 387`, `ref_len = 372`, `pred_len = 339`

Với 25 cases chính:

- Ref và truth đồng ý với nhau (`372 bp`).
- Tool prediction thiếu N-terminal.
- `pred_len = 329`, mà `329 % 3 = 2`, tức CDS length không hợp lệ theo frame.
- Sequence check cho thấy tool chọn ATG nội bộ gần HSP hơn thay vì ATG thật ở `truth_start`.

Kết luận:

- ORF7 lỗi chủ yếu do `rescue_start_codon` chọn nhầm internal ATG.
- Không phải do tblastn không tìm được gene.

## Cải Thiện Đã Thử: Ref-Length/Frame-Aware Start Rescue

Branch thử nghiệm:

`experiment/prrsv-orf7-start-rescue`

Thay đổi logic:

- `validate_cds_boundaries` check thêm `len(CDS) % 3 == 0`.
- `rescue_start_codon` không chọn ATG gần nhất một cách mù.
- Nếu biết ref protein length:
  - expected CDS length = `protein_length * 3 + 3`
  - ưu tiên ATG tạo CDS đúng frame và gần ref length.
- Nếu end boundary đã đúng, có thể suy ra start hợp lý từ:
  - `expected_start = pred_end - expected_CDS_length + 1`

Ý nghĩa:

- Với ORF7, expected length là `372 bp`.
- Nếu pred_end đúng, start hợp lý là vị trí tạo ra đoạn `372 bp`.
- Logic mới chọn ATG đó thay vì internal ATG ngắn hơn.

## Kết Quả Sau Cải Thiện

Full PRRSV validation được chạy lại từ đầu bằng notebook:

`prrsv_start_rescue_full_comparison.ipynb`

| Gene | Baseline exact | Experiment exact | Delta exact | Baseline coord | Experiment coord | Delta coord |
|---|---:|---:|---:|---:|---:|---:|
| ORF7 | 66/94 | 91/94 | +25 | 66/94 | 91/94 | +25 |
| ORF5 | 91/95 | 95/95 | +4 | 95/95 | 95/95 | 0 |
| ORF2b | 45/48 | 48/48 | +3 | 48/48 | 48/48 | 0 |
| ORF6 | 91/94 | 94/94 | +3 | 94/94 | 94/94 | 0 |
| ORF4 | 94/95 | 95/95 | +1 | 95/95 | 95/95 | 0 |
| ORF1a | 61/62 | 61/62 | 0 | 62/62 | 62/62 | 0 |
| ORF1b | 0/83 | 0/83 | 0 | 83/83 | 83/83 | 0 |
| ORF2a | 95/95 | 95/95 | 0 | 95/95 | 95/95 | 0 |
| ORF3 | 95/95 | 95/95 | 0 | 95/95 | 95/95 | 0 |

Changed-case summary:

- `ORF7`: 25 cases improved from wrong coords to correct.
- `ORF5`: 4 boundary-only cases became exact.
- `ORF6`: 3 boundary-only cases became exact.
- `ORF2b`: 3 truth-available boundary-only cases became exact.
- `ORF4`: 1 boundary-only case became exact.
- `ORF1b`: 1 status changed, but no accuracy improvement/regression.

## Case ORF7 Còn Lại

Sau cải thiện, ORF7 còn 3 non-exact cases:

| Record | Pred len | Truth len | Ref len | Nhận xét |
|---|---:|---:|---:|---|
| AY366525.1 | 339 | 387 | 372 | truth dài hơn ref 15 bp |
| DQ489311.1 | 339 | 387 | 372 | truth dài hơn ref 15 bp |
| DQ864705.1 | 339 | 387 | 372 | truth dài hơn ref 15 bp |

Không nên auto-fix 3 cases này bằng ref-length rescue vì truth dài hơn reference. Đây có thể là annotation/reference convention mismatch hoặc cần manual review.

## Kết Luận PRRSV

- PRRSV tblastn lifting có localization rất tốt.
- Baseline coordinate accuracy theo per-gene truth-available gần như perfect, trừ ORF7.
- ORF7 baseline fail do start rescue chọn internal ATG, không phải do tblastn không tìm được vùng gene.
- Ref-length/frame-aware start rescue cải thiện rõ:
  - ORF7: `66/94 -> 91/94`
  - ORF2b/ORF4/ORF5/ORF6: boundary-only cases được nâng lên exact.
- ORF1b vẫn cần xử lý riêng vì liên quan frameshift/start-boundary convention.

Một câu report có thể dùng:

> For PRRSV, tblastn reliably localizes most genes, including overlapping ORFs. The main baseline error came from start-codon rescue selecting internal ATGs, especially for ORF7. A frame- and reference-length-aware rescue strategy substantially improved exact boundary accuracy without reducing coordinate accuracy.

## Việc Nên Làm Tiếp

- Merge start-rescue experiment nếu full validation không phát hiện regression ngoài bảng changed cases.
- Giữ ORF1b như một nhóm riêng: coordinate-correct nhưng exact fail do frameshift/start convention.
- Manual review 3 ORF7 cases có `truth_len = 387` vì truth dài hơn ref.
- Sau khi merge, rerun notebook PRRSV chính để update final report/table.
