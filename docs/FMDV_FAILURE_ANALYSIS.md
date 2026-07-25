# FMDV failure analysis (convention vs tool)

Khác PRRSV: FMDV là **polyprotein cắt thành mature peptides (`mat_peptide`)**. Các gene
(Lpro, VP4, VP2, VP3, VP1, 2A, 2B, 2C, 3A, 3B, 3Cpro, 3Dpol) **không có ATG/stop riêng** —
biên của chúng là **vị trí protease cleavage**, nên tool đi nhánh `mat_peptide` (terminal
extrapolation), không phải CDS start/stop rescue. → **các fix generic của PRRSV không tác
động FMD.**

## Key finding

**FMD gần như KHÔNG có lỗi tool hệ thống.** Phần lớn "failure" là **annotation convention ở
vị trí cleavage** (mỗi lab đặt biên mature peptide lệch vài codon) và **granularity** — tool
đã lift đúng vùng, chỉ khác truth ở quy ước biên. Đây là kết quả *tốt*: tool xử lý cấu trúc
polyprotein chính xác; sai số là do người annotate, không phải tool.

> Con số dưới đây từ lần phân tích trước (bảng `failed_cases` FMD). Chạy lại
> `02_lifting_accuracy/fmdv_accuracy.ipynb` để chốt số cuối và điền vào paper.

## Per-gene (từ failure adjudication)

| Gene | Pattern | Quy về | Tool sai? |
|---|---|---|---|
| **2A** (~16 ca) | `delta_start=-6, pred≈ref`; 2A chỉ **54bp** nên lệch 6bp kéo IoU=0.889<0.90 → rơi failed | **ref/query convention + short-peptide IoU artifact** | Không |
| **VP1** (coord-only) | `delta_end=-6` ở mối VP1/2A; nửa HSP-trim + nửa ref-vs-truth | **cleavage boundary convention** | Không |
| **3A** (1 ca, MG372730) | `delta_end=-111, coverage 0.76`; pred ngắn hơn cả ref lẫn truth | **low-coverage** → terminal extrapolation *cố ý không nới* (ngưỡng 0.90) | Borderline, 1 ca |
| **3B** (2 ca) | `pred≈ref (~213bp), truth 72bp` | **granularity**: FMDV 3B (VPg) có **3 bản tandem**; ref chú thích cả cụm, truth 1 bản | Không |
| **VP2** (1 ca, AY687334) | `delta_start=+300` | **boundary convention** mối VP4/VP2 | Không |
| **2A** (1 ca, FJ175666) | `no_hit` | tblastn không ra hit | Tool (1 ca) |

## Chốt

- **Tool-side thực sự:** ~2 ca (3A low-coverage borderline + 1 no_hit) trên ~21 failed.
- **Còn lại = annotation convention / granularity** ở biên cleavage của polyprotein.
- → FMD chứng minh tool **lift đúng vùng cho mature peptides** (IoU cao), residual là *cách
  người annotate biên*, không phải lỗi tool.

## Contrast với PRRSV (điểm hay cho paper)

| | PRRSV (CDS/frameshift) | FMDV (polyprotein/mat_peptide) |
|---|---|---|
| Lỗi tool lộ ra | ORF7 (start N-term), ORF1a (stop over-read) → **đã sửa** (`VALIDATION_DRIVEN_FIXES.md`) | Gần như không có |
| Residual | ORF1b = frameshift start convention | mature-peptide cleavage boundary convention |
| Kết luận | validation **debug** được lỗi tool | validation xác nhận tool **chính xác trên polyprotein**, tách được convention |

→ Hai virus cho hai câu chuyện bổ trợ: PRRSV cho thấy validation **tìm & sửa** lỗi tool;
FMDV cho thấy validation **phân biệt convention khỏi lỗi** trên cấu trúc polyprotein khó.

## Nhắc nhở cho Results/Discussion
- Nêu rõ **short-peptide IoU artifact** (2A 54bp): cùng lệch 6bp mà VP1 (633bp) chỉ coord-only,
  2A thành failed → cân nhắc **ngưỡng theo bp thay vì chỉ IoU** cho peptide ngắn.
- FMD 3B granularity (VPg 3 bản) là ví dụ đẹp về **annotation disagreement ≠ tool error**.
