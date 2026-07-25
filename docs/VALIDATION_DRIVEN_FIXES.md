# Validation-driven tool fixes (PRRSV case)

Validation không chỉ đo accuracy — nó **chỉ ra lỗi tool cụ thể**, dẫn tới 4 bản vá.
Điểm mấu chốt (phải nhấn trong paper): **mọi fix đều generic (kích hoạt bởi tính chất
gene / annotation của REF, không hard-code tên gene/virus) và leakage-free (chỉ dùng
reference alignment, KHÔNG bao giờ dùng truth annotation)** — nên không phải "chỉnh cho
khớp đáp án", và tự tổng quát sang virus/ref khác.

## Bug → Fix → Impact

| # | Gene lộ ra | Lỗi tool (root cause) | Fix (generic, leakage-free) | Impact đo được |
|---|---|---|---|---|
| 1 | **ORF7** | Start-rescue tóm **ATG nội bộ** khi HSP hụt đầu N. 3 lỗi con: (a) tìm ATG quanh điểm HSP thay vì nới ngược lên `missing_n_aa×3`; (b) vòng quét bỏ qua offset 0 (không xét chính vị trí neo); (c) tính length dựa trên end đã bị stop-rescue kéo dài | Nới anchor lên theo số residue N-terminal ref chưa align (từ HSP, không phải truth); quét từ offset 0; dùng raw HSP end cho length | `delta_start +72 → 0`, IoU `0.69 → >0.99`; **~29 record ORF7 failed → exact** |
| 2 | **ORF1a** | **Đọc lố qua stop đầu tiên**: lift bám ref dài hơn, chạy quá stop thật của query → CDS chứa stop nội bộ | Trim đuôi về **stop in-frame đầu tiên** (khi start đã ATG) | `delta_end +15/+66 → 0`; **~40 record ORF1a coord-only → exact** |
| 3 | **validator** | `validate_cds_boundaries` không quét **stop nội bộ** → cho qua ORF vô lý | Thêm cờ `has_internal_stop`; `valid=False` nếu có stop giữa | Chặn ORF sai; hỗ trợ fix #2 |
| 4 | **ORF1b** (frameshift) | Ép **tìm ATG cho gene không có start** (−1 PRF) → tóm ATG bừa, start lệch ~100bp, status `invalid_boundaries` | Nếu **ref CDS không bắt đầu bằng ATG** (partial/frameshift) → **bỏ start-rescue**, giữ start từ HSP, status `ok_no_start_codon` | start lệch `~100bp → ±3–12bp`, IoU `~0.98 → 0.997`, status sạch; **coord-only đúng bản chất (không ép exact)** |

## Nguyên tắc thiết kế (cho phần Methods/Discussion)

- **Generic, không hard-code:** mỗi fix kích hoạt bởi *tính chất đo được* — số aa N-terminal
  chưa align (#1), có stop nội bộ hay không (#2, #3), ref CDS có ATG hay không (#4) — chứ
  không phải `if virus == PRRSV` hay `if gene == ORF1b`. Tự áp cho bất kỳ gene/virus nào có
  cùng đặc điểm (vd frameshift-fusion ORF phổ biến ở corona/arteri/retro/astro/toti).
- **Leakage-free:** fix chỉ dùng **reference protein + tblastn alignment coverage**, không hề
  đọc vị trí start/end của truth. Validation harness chạy tool (ref+query) rồi mới so với
  truth *sau đó*. Đây là điểm bảo vệ trước reviewer: cải thiện **không** do overfit đáp án.
- **Phân biệt bug vs convention:** validation phân tách được (i) **lỗi tool thật** (ORF7, ORF1a
  — đã sửa) khỏi (ii) **annotation convention** (start ORF1b — frameshift, không có ATG, mỗi
  lab đặt khác vài bp). Cái (ii) tool **cố ý không ép** exact — báo coord-only + giải thích.

## Đưa vào paper thế nào

- **Results:** bảng before/after per-gene (ablation "validation-guided refinement"). Con số
  ORF7 (+72→0), ORF1a (+66→0), ORF1b (~100bp→±10bp, IoU 0.997).
- **Discussion:** validation framework vừa là *thước đo* vừa là *công cụ debug*; và nó phân
  biệt được lỗi tool (sửa được) với convention (không nên ép) — điểm khoa học mạnh.
- **Caveat cần nêu (trung thực):** các fix được dẫn dắt bởi validation trên chính bộ annotated
  này. Giảm nhẹ nhờ: (a) fix leakage-free (không dùng truth); (b) generic; (c) kiểm chứng thêm
  trên **ref khác / held-out** (đang làm). Nên trình bày là *phát triển phương pháp có validation
  dẫn đường*, không phải tuning trên test set.

## Trạng thái PRRSV sau fix

| Gene | Kết quả | Ghi chú |
|---|---|---|
| ORF3,4,5,6,7, ORF2a | exact | ORF7 nhờ fix #1 |
| ORF1a | exact | nhờ fix #2 (1 ca lẻ là truth ẩu) |
| ORF2b | ~exact | vài outlier truth |
| **ORF1b** | **coord-only IoU 0.997** | `ok_no_start_codon`; frameshift convention, cố ý không ép |
| AF331831.1 | loại | record truth-lỗi (ORF1b gắn vào vùng ORF1a) |
