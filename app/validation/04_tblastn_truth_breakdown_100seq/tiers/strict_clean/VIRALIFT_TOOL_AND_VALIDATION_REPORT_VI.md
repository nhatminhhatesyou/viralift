# ViraLift: Tổng Quan Tool Và Báo Cáo Validation

## Mục Lục

- [1. ViraLift Là Gì?](#1-viralift-là-gì)
- [2. Tool Chạy Như Thế Nào?](#2-tool-chạy-như-thế-nào)
- [3. Ví Dụ Sử Dụng](#3-ví-dụ-sử-dụng)
- [4. Validation Dataset](#4-validation-dataset)
- [5. Cách Tính Accuracy](#5-cách-tính-accuracy)
- [6. Kết Quả Tổng Quan](#6-kết-quả-tổng-quan)
- [7. Breakdown Theo Gene: FMD](#7-breakdown-theo-gene-fmd)
- [8. Breakdown Theo Gene: PRRSV](#8-breakdown-theo-gene-prrsv)
- [9. Finding Chính](#9-finding-chính)
- [10. Kết Luận](#10-kết-luận)

## 1. ViraLift Là Gì?

ViraLift là tool hỗ trợ chuẩn hóa và chuyển annotation gene/peptide cho virus genome.

Mục tiêu chính:

- Chuẩn hóa tên gene từ nhiều cách đặt tên khác nhau trong GenBank.
- Nếu query genome đã có annotation đủ tốt, tool trích xuất trực tiếp từ annotation.
- Nếu query genome chưa có annotation hoặc annotation thiếu, tool dùng reference chuẩn do user cung cấp để lift gene bằng `tblastn`.
- Xuất kết quả thành bảng TSV và FASTA để user kiểm tra hoặc dùng tiếp.

Tool hiện được validate chính trên 2 nhóm virus:

- FMDV: dùng `mat_peptide`, gene/peptide thường liên tiếp.
- PRRSV: dùng `CDS`/ORF, có nhiều gene chồng lấp và convention annotation phức tạp hơn.

## 2. Tool Chạy Như Thế Nào?

Pipeline chính:

1. User cung cấp reference GenBank đã annotation chuẩn.
2. Tool đọc reference và xác định feature type hữu ích, ví dụ `CDS` hoặc `mat_peptide`.
3. Tool chuẩn hóa tên gene bằng alias map.
4. Với mỗi query genome:
   - Nếu query có annotation đủ hữu ích, tool direct extract.
   - Nếu query thiếu annotation, tool dùng `tblastn` để map protein từ reference sang genome query.
5. Tool validate boundary:
   - Với `CDS`: kiểm tra start/stop codon và frame.
   - Với `mat_peptide`: không dùng start/stop codon vì peptide cleavage product không nhất thiết có ATG/stop riêng.
6. Tool xuất:
   - bảng prediction,
   - FASTA sequence,
   - run summary,
   - status cho từng gene.

Với `tblastn` lifting:

- Reference feature được dịch sang protein.
- Protein được search trên query genome bằng `tblastn`.
- HSPs được merge để suy ra tọa độ nucleotide.
- Tool chỉnh boundary nếu cần, ví dụ terminal extrapolation cho FMD hoặc start-codon rescue cho PRRSV.

## 3. Ví Dụ Sử Dụng

Ví dụ chạy bằng CLI:

```bash
python -m app.src.main \
  --ref app/data/PRRS_ref_test.gb \
  --query app/data/PRRS_PP946131_noAnno.gb \
  --out output/prrs_example
```

Ví dụ logic kết quả:

```text
ORF5  -> lifted bằng tblastn, tọa độ đúng, status ok
ORF7  -> lifted bằng tblastn, start được rescue, status ok_rescued
ORF2b -> lifted đúng nhưng một số query truth có thể không annotate ORF2b riêng
```

Ý nghĩa status thường gặp:

- `ok`: lift thành công, boundary hợp lệ.
- `ok_rescued`: boundary ban đầu chưa ổn, tool đã rescue lại.
- `ok_extrapolated`: boundary được mở rộng bằng terminal extrapolation.
- `no_hit`: không tìm được hit phù hợp.
- `invalid_boundaries`: tìm được vùng nhưng boundary/codon validation chưa đạt.

## 4. Validation Dataset

Validation dùng 2 bộ virus annotated records:

- FMDV strict-clean records.
- PRRSV strict-clean records.

Nguyên tắc lọc strict-clean:

- Chỉ giữ record có annotation đủ để làm ground truth.
- Tên gene phải map được về canonical name qua alias map.
- Với PRRSV, không dùng nested-feature filter mù vì `ORF2b` có thể overlap/nest trong `ORF2a` nhưng vẫn là gene thật.
- Với các gene không xuất hiện trong truth của một query record, không dùng record đó để chấm accuracy cho gene đó.

Ví dụ:

- Reference có `ORF7`.
- Trong 95 PRRSV records, truth có `ORF7` ở 94 records.
- Vậy accuracy của `ORF7` chỉ tính trên 94 records đó.

Điều này giúp tránh đánh giá sai tool vì annotation truth thiếu hoặc dùng convention khác.

## 5. Cách Tính Accuracy

Validation được tính theo từng gene.

Các metric chính:

- `total`: số query records thật sự có gene đó trong truth.
- `exact`: prediction khớp hoàn toàn tên gene + start + end.
- `coord_only`: tọa độ đúng theo IoU threshold nhưng không exact tuyệt đối.
- `failed`: không exact và cũng không coordinate-correct.
- `accuracy_pct = (exact + coord_only) / total`.
- `IoU`: chỉ số đo mức overlap giữa tọa độ prediction và tọa độ truth.

### IoU Là Gì?

`IoU` là viết tắt của `Intersection over Union`.

Nói đơn giản:

```text
IoU = độ dài phần prediction và truth chồng lên nhau
      / độ dài vùng bao phủ bởi cả prediction và truth
```

Ví dụ 1: prediction khớp hoàn toàn truth

```text
Truth:      100 - 199
Prediction: 100 - 199
IoU = 100 / 100 = 1.00
```

Ví dụ 2: prediction lệch một chút nhưng vẫn gần đúng

```text
Truth:      100 - 199
Prediction: 103 - 199
Overlap = 97 bp
Union = 100 bp
IoU = 0.97
```

Ví dụ 3: prediction lệch nhiều hơn

```text
Truth:      100 - 199
Prediction: 130 - 199
Overlap = 70 bp
Union = 100 bp
IoU = 0.70
```

Trong validation này, một prediction được tính là `coord_correct` nếu:

```text
same-gene truth tồn tại
và IoU >= 0.90
```

Vì vậy:

- `IoU = 1.00`: tọa độ khớp hoàn toàn.
- `IoU >= 0.90`: tọa độ đủ gần, tính là coordinate-correct.
- `IoU < 0.90`: tính là failed về tọa độ.

Lưu ý quan trọng:

- Với gene dài, lệch vài bp thường IoU vẫn cao.
- Với gene rất ngắn như FMD `2A`, chỉ lệch `6 bp` cũng có thể làm IoU tụt dưới `0.90`.
- Vì vậy một số failed cases ở gene ngắn không nhất thiết là lift sai vùng lớn, mà là boundary precision bị strict scoring phạt mạnh.

Ví dụ:

```text
ORF7 total = 94
exact = 91
coord_only = 0
failed = 3
accuracy = 91 / 94 = 96.81%
```

Lý do tách `exact` và `coord_only`:

- `exact` đo boundary precision rất nghiêm ngặt.
- `coord_only` cho biết tool đã tìm đúng vùng gene, dù start/end lệch nhỏ hoặc annotation convention khác.
- Với virus annotation, nhiều case lệch vài bp có thể do convention, không nhất thiết là localization failure.

## 6. Kết Quả Tổng Quan

Kết quả sau các cải thiện hiện tại:

![Overall accuracy by virus](summary_outputs/fmd_prrsv_overall_accuracy.png)

| Virus | Total gene-record cases | Exact | Coord only | Failed | Accuracy |
|---|---:|---:|---:|---:|---:|
| FMD | 1149 | 1114 | 16 | 19 | 98.35% |
| PRRSV | 761 | 674 | 84 | 3 | 99.61% |
| All | 1910 | 1788 | 100 | 22 | 98.85% |

Diễn giải:

- FMD đạt `98.35%` khi tính exact + coordinate-correct.
- PRRSV đạt `99.61%` khi tính trên các gene thật sự có trong truth.
- Phần lớn lỗi còn lại là boundary issue nhỏ hoặc annotation/ref-truth mismatch.

Figure dưới đây breakdown theo từng gene. Màu xanh đậm là exact, xanh nhạt là coord-only, đỏ là failed.

![Per-gene accuracy](summary_outputs/fmd_prrsv_per_gene_accuracy.png)

## 7. Breakdown Theo Gene: FMD

FMD dùng `mat_peptide`, nên tool không dùng codon rescue mà dùng terminal extrapolation cho các case HSP thiếu amino acid ở đầu/cuối.

Evidence chính cho cải thiện FMD:

![FMD terminal extrapolation accuracy comparison](terminal_extrapolation_outputs/fmd_accuracy_comparison.png)

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

Finding chính của FMD:

- Hầu hết gene đạt gần hoặc đúng `100%`.
- `VP1` có nhiều `coord_only`: tool tìm đúng vùng nhưng boundary không exact.
- `2A` là gene ngắn, chỉ lệch vài bp cũng làm IoU tụt dưới threshold.
- Một số lỗi còn lại ở `VP2` và `3A` cần manual review riêng.

Cải thiện đã làm:

- Terminal extrapolation dùng HSP query protein coordinates để bù amino acids thiếu ở terminal.
- Exact accuracy FMD tăng từ `94.44%` lên khoảng `96.95%` theo per-gene truth-available scoring.
- Tool-side failures giảm rõ, đặc biệt ở các case thiếu N-terminal `12 bp`.

## 8. Breakdown Theo Gene: PRRSV

PRRSV dùng `CDS`/ORF và có overlapping genes, ví dụ `ORF2a/ORF2b`.

Evidence chính cho cải thiện PRRSV start rescue:

![PRRSV start rescue exact comparison](prrsv_start_rescue_full_outputs/prrsv_start_rescue_exact_comparison.png)

| Gene | Total | Exact | Coord only | Failed | Accuracy |
|---|---:|---:|---:|---:|---:|
| ORF1a | 62 | 61 | 1 | 0 | 100.00% |
| ORF1b | 83 | 0 | 83 | 0 | 100.00% |
| ORF2a | 95 | 95 | 0 | 0 | 100.00% |
| ORF2b | 48 | 48 | 0 | 0 | 100.00% |
| ORF3 | 95 | 95 | 0 | 0 | 100.00% |
| ORF4 | 95 | 95 | 0 | 0 | 100.00% |
| ORF5 | 95 | 95 | 0 | 0 | 100.00% |
| ORF6 | 94 | 94 | 0 | 0 | 100.00% |
| ORF7 | 94 | 91 | 0 | 3 | 96.81% |

Finding chính của PRRSV:

- `ORF2a`, `ORF2b`, `ORF3`, `ORF4`, `ORF5`, `ORF6` đạt `100%`.
- `ORF1b` có `0 exact` nhưng `83 coord_only`, nghĩa là tool tìm đúng vùng nhưng exact boundary khác truth. Đây là start/frameshift convention, không phải localization failure.
- `ORF7` là lỗi chính ở baseline, do start rescue chọn nhầm internal ATG.

Cải thiện đã làm:

- Thêm check CDS frame: `len(CDS) % 3 == 0`.
- Start rescue không chọn ATG gần nhất một cách mù.
- Ưu tiên ATG tạo CDS đúng frame và có length gần reference.

Kết quả cải thiện PRRSV:

- `ORF7` tăng từ `66/94` lên `91/94`.
- `ORF2b`, `ORF4`, `ORF5`, `ORF6` cũng tăng exact vì boundary-only cases được sửa.
- 3 case ORF7 còn lại có truth dài hơn reference, nên không nên auto-fix bằng ref-length rescue.

ORF7-only experiment cũng xác nhận đúng cơ chế lỗi:

![ORF7 rescue accuracy comparison](orf7_start_rescue_outputs/orf7_rescue_accuracy_comparison.png)

## 9. Finding Chính

Các finding quan trọng:

1. Accuracy phải tính trên gene thật sự có trong truth.
   - Nếu query truth không annotate gene đó, không thể dùng case đó để kết luận tool sai.

2. FMD và PRRSV cần xử lý khác nhau.
   - FMD dùng `mat_peptide`: phù hợp với terminal extrapolation.
   - PRRSV dùng `CDS`: phù hợp với codon/frame-aware rescue.

3. `tblastn` thường tìm đúng vùng gene.
   - Lỗi chính không phải search sai vùng lớn.
   - Lỗi thường nằm ở boundary start/end.

4. ORF1b PRRSV cần interpretation riêng.
   - Coord đúng `83/83`.
   - Exact fail do frameshift/start-boundary convention.

5. Gene ngắn như FMD `2A` dễ bị IoU phạt mạnh.
   - Lệch `6 bp` ở gene rất ngắn có thể làm IoU dưới `0.90`.

## 10. Kết Luận

ViraLift cho kết quả tốt trên cả FMD và PRRSV khi validation được tính đúng theo same-gene truth availability.

Kết quả cuối:

- FMD: `98.35%`
- PRRSV: `99.61%`
- Tổng 2 tập: `98.85%`

Kết luận kỹ thuật:

- `tblastn` là hướng phù hợp cho annotation transfer khi query thiếu annotation.
- Các lỗi còn lại chủ yếu là boundary precision hoặc annotation convention mismatch.
- Terminal extrapolation cải thiện FMD.
- Frame/ref-length-aware start rescue cải thiện PRRSV, đặc biệt ORF7.

Một câu tóm tắt:

> ViraLift reliably transfers viral gene annotations using reference-guided tblastn lifting. On strict-clean FMDV and PRRSV validation datasets, the tool achieves high localization accuracy, with remaining failures mostly caused by short-feature boundary sensitivity or reference/query annotation convention differences.
