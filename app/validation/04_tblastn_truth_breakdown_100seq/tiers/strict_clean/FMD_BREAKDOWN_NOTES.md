# FMD Breakdown Notes

Ghi chú cho phần `strict_clean` FMD validation, gồm baseline và kết quả sau khi thử terminal extrapolation.

## Cách Tính Accuracy

FMD dùng `mat_peptide`, nên validation so từng peptide/gene trên các query records có same-gene truth.

Metric chính:

- `total`: số record có truth gene đó.
- `exact`: pred name, start, end khớp hoàn toàn với truth.
- `coord_only`: tọa độ đúng theo IoU threshold nhưng không exact.
- `failed`: không exact và cũng không coordinate-correct.
- `accuracy_pct = (exact + coord_only) / total`.

## Baseline Tổng Quan

Baseline trước terminal extrapolation:

- Tổng predictions: `1152`
- Exact match: `1088 / 1152 = 94.44%`
- Coordinate-correct: `1131 / 1152 = 98.18%`
- Non-exact: `64 / 1152 = 5.56%`

Trong 64 non-exact cases:

- Ref/truth annotation issue: `28 / 64`
- Tool-side boundary issue: `36 / 64`

## Finding Chính Ở Baseline

### Ref/Truth Annotation Issues

Một số case fail không phải do tool lift sai, mà do reference và query truth dùng boundary khác nhau.

Các nhóm chính:

- `2A`: nhiều query truth annotate `2A` ngắn hơn reference khoảng `6 bp`.
- `VP1`: một số query truth dài/ngắn khác reference vài bp.
- `VP2`: có case truth dùng boundary dài hơn reference rất nhiều.
- Một số gene truth absent nên không thể score công bằng.

Kết luận:

- Không nên quy toàn bộ non-exact baseline thành lỗi tool.
- FMD có vài boundary convention mismatch giữa ref và query truth.

### Tool-Side Boundary Issues

Pattern tool-side rõ nhất:

- `Lpro` và `3Cpro`: thiếu `12 bp` ở N-terminal.
- `3B`: thiếu `3 bp` ở C-terminal.
- `3A`: một case thiếu C-terminal lớn `111 bp`.
- `2A`: peptide rất ngắn nên lệch `3-6 bp` làm IoU tụt mạnh.

Pattern quan trọng:

- Nhiều case có `delta_end = 0`, nhưng `delta_start > 0`.
- Nghĩa là tool kết thúc đúng, nhưng bắt đầu muộn, thiếu N-terminal amino acids.

## Vì Sao FMD Bị Thiếu Terminal?

FMD dùng `mat_peptide`, không phải CDS riêng.

Vì vậy code không dùng:

- `rescue_start_codon`
- `rescue_stop_codon`
- start/stop codon validation

Baseline lấy boundary trực tiếp từ HSP span của `tblastn`.

Nhưng `tblastn` là local alignment:

- Nếu HSP không cover vài amino acids đầu/cuối protein,
- prediction sẽ thiếu đoạn terminal đó.

Ví dụ:

- HSP bắt đầu ở amino acid 5.
- Thiếu `4 aa` đầu.
- Tọa độ nucleotide bị thiếu `4 * 3 = 12 bp`.

## Cải Thiện Đã Thử: Terminal Extrapolation

Ý tưởng:

- Dùng HSP query protein coordinate để biết HSP thiếu bao nhiêu aa ở đầu/cuối.
- Nếu coverage cao và thiếu ít aa, extend boundary thêm `missing_aa * 3 bp`.

Rule conservative:

- Chỉ apply khi coverage đủ cao.
- Chỉ extend terminal missing nhỏ.
- Không dùng cho CDS/codon rescue path.
- Hiện chủ yếu áp dụng cho `mat_peptide` như FMD.

## Kết Quả Sau Cải Thiện

Sau terminal extrapolation:

- Exact match tăng từ `94.44%` lên `96.70%`.
- Tool-side failures giảm từ `36` xuống `9`.
- Nhiều case thiếu N-terminal `12 bp` được fix.

Các nhóm được cải thiện rõ:

- `Lpro`: fix nhiều case thiếu N-terminal.
- `3Cpro`: fix nhiều case thiếu N-terminal.
- `3B`: fix case thiếu C-terminal nhỏ.
- Một phần boundary offset ngắn cũng được cải thiện.

## Kết Quả Theo Gene Sau Cải Thiện

Theo notebook summary mới:

| Gene | Total | Exact | Coord only | Failed | Accuracy |
|---|---:|---:|---:|---:|---:|
| Lpro | 96 | 96 | 0 | 0 | 100.00% |
| VP4 | 95 | 95 | 0 | 0 | 100.00% |
| VP2 | 96 | 95 | 0 | 1 | 98.96% |
| VP3 | 95 | 95 | 0 | 0 | 100.00% |
| VP1 | 95 | 79 | 16 | 0 | 100.00% |
| 2A | 96 | 79 | 0 | 17 | 82.29% |
| 2B | 96 | 96 | 0 | 0 | 100.00% |
| 2C | 96 | 96 | 0 | 0 | 100.00% |
| 3A | 96 | 95 | 0 | 1 | 98.96% |
| 3B | 96 | 96 | 0 | 0 | 100.00% |
| 3Cpro | 96 | 96 | 0 | 0 | 100.00% |
| 3Dpol | 96 | 96 | 0 | 0 | 100.00% |

Overall FMD sau cải thiện:

- `1130 / 1149` coordinate-correct-or-exact
- Accuracy: `98.35%`
- Exact: `1114 / 1149 = 96.95%`
- Failed: `19 / 1149 = 1.65%`

## Case Còn Lại

Các failed cases còn lại chủ yếu nằm ở:

- `2A`: peptide rất ngắn, lệch `6 bp` làm IoU chỉ khoảng `0.8889`, dưới threshold `0.90`.
- `VP2`: 1 case boundary lệch lớn.
- `3A`: 1 case thiếu C-terminal lớn.

Interpretation:

- `2A` failures có thể bị strict IoU phạt mạnh vì gene quá ngắn.
- `VP2` và `3A` là case cần review riêng.

## Kết Luận FMD

- FMD tblastn lifting rất tốt ở mức localization.
- Baseline fail chủ yếu do terminal boundary thiếu vài amino acids hoặc ref/truth boundary convention mismatch.
- Terminal extrapolation cải thiện rõ exact accuracy.
- Sau cải thiện, FMD đạt `98.35%` accuracy theo truth-available per-gene scoring.
- Severe tool-side failures còn lại ít, nổi bật nhất là một số case `2A`, `VP2`, `3A`.

Một câu report có thể dùng:

> For FMDV, tblastn lifting is highly reliable at the localization level. Most baseline errors were small terminal boundary offsets caused by local HSP truncation or reference-vs-query boundary convention differences. Conservative terminal extrapolation improved exact boundary recovery while keeping coordinate accuracy high.

## Việc Nên Làm Tiếp

- Giữ terminal extrapolation cho `mat_peptide` path.
- Review riêng `2A` vì gene quá ngắn và IoU threshold dễ phạt nặng.
- Manual review `VP2` và `3A` failed cases còn lại.
- Sau khi merge PRRSV start rescue, rerun notebook summary chung để cập nhật final numbers.
