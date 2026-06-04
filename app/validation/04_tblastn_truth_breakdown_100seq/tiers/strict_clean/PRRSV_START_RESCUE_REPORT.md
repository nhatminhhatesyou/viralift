# PRRSV Start Rescue Validation Report

## Mục tiêu

Validate PRRSV trên strict-clean dataset để đánh giá độ chính xác của `tblastn` lifting sau khi sửa logic start-codon rescue. Trọng tâm là kiểm tra pattern lỗi: tool lift đúng vùng/end boundary nhưng start bị chọn nhầm ATG nội bộ.

## Metric dùng trong report

- `truth_available_records`: số record có annotation truth cho gene tương ứng với reference. Đây là denominator chính.
- `exact_correct`: pred_name, start, end khớp hoàn toàn với truth.
- `coord_correct`: tọa độ overlap tốt với truth theo IoU threshold, nhưng có thể lệch boundary nhỏ.
- `extra_predictions_without_truth`: tool predict gene nhưng record truth không annotate gene đó; không tính vào accuracy chính.
- `delta_start = pred_start - truth_start`; `delta_end = pred_end - truth_end`.

## Baseline

Baseline là PRRSV strict-clean validation trước khi sửa start rescue.

| Gene | Truth records | Baseline exact | Baseline coord | Nhận xét |
|---|---:|---:|---:|---|
| ORF2a | 95 | 95 | 95 | Perfect |
| ORF3 | 95 | 95 | 95 | Perfect |
| ORF2b | 48 | 45 | 48 | 3 case lệch boundary nhỏ |
| ORF4 | 95 | 94 | 95 | 1 case lệch boundary nhỏ |
| ORF5 | 95 | 91 | 95 | 4 case lệch boundary nhỏ |
| ORF6 | 94 | 91 | 94 | 3 case lệch boundary nhỏ |
| ORF7 | 94 | 66 | 66 | 28 case sai start rõ rệt |
| ORF1a | 62 | 61 | 62 | 1 case lệch boundary nhỏ |
| ORF1b | 83 | 0 | 83 | Coord đúng hết, exact fail do start/frameshift convention |

## Finding chính

ORF7 có pattern rất rõ:

- 28 failed cases đều có `delta_end = 0`: end boundary đúng.
- 25/28 cases có `truth_len = ref_len = 372`, nhưng `pred_len = 329`.
- `329 % 3 = 2`, tức CDS length không hợp lệ theo frame.
- Sequence check cho thấy có ATG thật ở `truth_start`, nhưng rescue cũ chọn ATG nội bộ gần HSP hơn.

Kết luận: lỗi chính không phải do `tblastn` tìm sai vùng gene, mà do post-processing `rescue_start_codon` chọn ATG gần nhất thay vì ATG tạo CDS hợp lệ/gần ref length.

## Cải thiện đã thử

Logic mới:

- `validate_cds_boundaries` check thêm `len(CDS) % 3 == 0`.
- `rescue_start_codon` không chọn ATG gần nhất một cách mù.
- Khi có ref protein length, suy ra expected CDS length = `protein_length * 3 + 3`.
- Ưu tiên ATG candidate tạo CDS:
  - đúng frame,
  - length gần ref CDS length,
  - sau đó mới xét khoảng cách tới HSP start.

## Kết quả sau cải thiện

Full PRRSV validation được chạy lại từ đầu bằng notebook:

`prrsv_start_rescue_full_comparison.ipynb`

Figure so sánh:

![PRRSV start rescue exact comparison](prrsv_start_rescue_full_outputs/prrsv_start_rescue_exact_comparison.png)

| Gene | Baseline exact | Experiment exact | Delta | Baseline coord | Experiment coord | Delta |
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

Tổng raw PRRSV prediction-level result sau experiment:

- `674/855` exact = `78.83%`
- `758/855` coordinate correct = `88.65%`

## Bằng chứng changed cases

Changed-case summary:

| Gene | Changed cases | Exact delta | Coord delta |
|---|---:|---:|---:|
| ORF7 | 25 | +25 | +25 |
| ORF5 | 4 | +4 | 0 |
| ORF6 | 3 | +3 | 0 |
| ORF2b | 4 | +3 | 0 |
| ORF4 | 1 | +1 | 0 |
| ORF1b | 1 | 0 | 0 |

ORF7 là cải thiện lớn nhất: 25 cases chuyển từ `Wrong coords` sang `Correct`.

Các gene ORF2b/ORF4/ORF5/ORF6 vốn đã coord-correct, sau logic mới chỉ được chỉnh boundary để thành exact.

## Case còn lại

ORF7 còn 3 non-exact cases:

| Record | Pred len | Truth len | Ref len | Pattern |
|---|---:|---:|---:|---|
| AY366525.1 | 339 | 387 | 372 | truth dài hơn ref 15 bp |
| DQ489311.1 | 339 | 387 | 372 | truth dài hơn ref 15 bp |
| DQ864705.1 | 339 | 387 | 372 | truth dài hơn ref 15 bp |

Các case này không nên auto-fix bằng ref-length rescue vì truth dài hơn reference. Cần manual review hoặc xem đây là annotation/reference convention mismatch.

## Kết luận

Start-rescue cải thiện PRRSV rõ rệt và không gây regression lớn trong full validation. Lỗi ORF7 chủ yếu đến từ rescue chọn nhầm internal ATG, không phải do `tblastn` không tìm được gene. Sau khi dùng frame + ref-length-aware start rescue, ORF7 tăng từ `66/94` lên `91/94` exact/coord, đồng thời các boundary-only cases ở ORF2b/ORF4/ORF5/ORF6 cũng được nâng lên exact.
