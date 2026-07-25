# FMDV failure analysis — convention vs tool

> Viết từ **lần chạy thật** (`02_lifting_accuracy/fmdv_accuracy.ipynb`) + verify bằng sequence,
> không phải suy đoán. Bổ trợ cho `VALIDATION_DRIVEN_FIXES.md` (câu chuyện PRRSV).

FMDV là **polyprotein cắt thành mature peptides** — biên các gene là **vị trí protease cleavage**,
không có ATG/stop riêng. Nên tool đi nhánh `mat_peptide` (terminal extrapolation), khác hẳn nhánh
CDS của PRRSV. **Không có bug tool hệ thống ở FMDV**; residual gần như toàn bộ là annotation convention.

## Kết quả

| | exact | coord-only | failed | no_hit | accuracy |
|---|---|---|---|---|---|
| Trước gap-fill | 1136 | 16 | 21 | 1 | 98.21% |
| **Sau gap-fill** | **1137** | 16 | 20 | **0** | **98.29%** |
| **Sau + bp-tolerance** | 1137 | **32** | **4** | 0 | **~99.66%** |

(bp-tolerance là thay đổi **metric**, không phải tool — xem mục dưới.)

## Phát hiện chính: GenBank tự mâu thuẫn ở biên VP1/2A

Ở mối giáp VP1↔2A có **2 amino acid `MM`**. Đếm trên 50 record có annotate 2A:

- **38 record** xếp `MM` vào **2A** (2A = 18aa) — giống reference
- **12 record** xếp `MM` vào **VP1** (2A = 16aa)

→ **Không kiểu nào "đúng tuyệt đối"** — hai trường phái annotate. Tool áp **nhất quán theo ref**:
38 record cùng kiểu → exact; số còn lại lệch đúng 2 codon (`delta_start=-6` ở 2A, `delta_end=-6` ở VP1).

**Đây chính là standardization đang hoạt động** — gom các record annotate lệch về **một** convention.
Bị chấm "sai" chỉ vì truth dùng convention khác, không phải tool lift sai.

**Cùng một boundary, hai kết cục khác nhau** (minh họa đẹp cho vấn đề metric):

| gene | dài | lệch | IoU | bucket |
|---|---|---|---|---|
| VP1 | 633bp | 6bp | 0.991 | coord-only |
| 2A | 54bp | 6bp | **0.889** | **failed** |

→ **IoU phụ thuộc kích thước**: cùng 6bp mà gene dài thì qua, peptide 18aa thì rớt.

## Per-gene residual

| Gene | Ca | Nguyên nhân | Blame |
|---|---|---|---|
| **2A** | 16 | boundary convention VP1/2A (`MM`) + IoU khắt khe trên peptide ngắn | ref_truth |
| **VP1** | 16 | cùng boundary đó, nhìn từ phía VP1 | ref_truth (coord-only) |
| **3B** | 2 | VPg có **3 bản tandem**; ref chú thích cả cụm (213bp), truth 1 bản (72bp) → granularity | ref_truth |
| **VP2** | 1 | boundary convention mối VP4/VP2 (truth dài thêm 300bp 5') | ref_truth |
| **3A** | 1 | strain phân kỳ, coverage 0.76 → HSP cắt cụt C-term; extrapolation cố ý không nới (ngưỡng 0.90) | ca khó, 1/98 |
| 2A | 1 | `no_hit` — **đã cứu bằng gap-fill** | — |

## Hai cơ chế rescue cho FMDV (đều generic, leakage-free)

1. **Terminal extrapolation** *(có sẵn từ branch `experiment/fmd-terminal-extrapolation`, đã merge)* —
   tblastn là local aligner nên cắt cụt vài aa ở đầu/cuối; nới lại `missing_aa × 3` khi coverage ≥ 0.90.
2. **Polyprotein gap-fill** *(mới)* — mature peptide **liền kề nhau, không hở**, nên peptide mà tblastn
   bỏ sót (2A chỉ **18aa**, quá ngắn để ra hit có ý nghĩa) được suy từ **khe giữa 2 hàng xóm đã lift**:
   `2A = [VP1_end+1, 2B_start−1]`. Verify trên FJ175666: điền ra `3986..4039` = **đúng truth, exact**.
   Dùng tọa độ hàng xóm + thứ tự ref, **không đụng truth**.

## bp-tolerance (thay đổi METRIC, không phải tool)

IoU **phụ thuộc kích thước feature**. Cùng một lệch tuyệt đối vài bp: gene dài → IoU ~0.99 (coord-only),
peptide ngắn → IoU < 0.90 (failed). Nên thêm: prediction có **cả hai biên trong `bp_tolerance` (mặc định
6bp = 2 codon)** thì tính coord-correct, bất kể độ dài.

Kiểm chứng **không nới lỏng bừa** — các ca lệch lớn **vẫn failed**:

| ca | lệch | trong tolerance? | bucket |
|---|---|---|---|
| 2A | −6 / 0 | ✅ | coord-only |
| 3B | −69 / +72 | ❌ | failed |
| VP2 | +300 / 0 | ❌ | failed |
| 3A | 0 / −111 | ❌ | failed |

**Phải nêu rõ trong paper:** phần accuracy tăng (98.29% → ~99.66%) đến từ **metric công bằng hơn**,
KHÔNG phải lifting tốt hơn. Lý do chính đáng: không nên phạt tool vì nó **chuẩn hóa** một biên mà bản
thân GenBank chú thích theo 2 kiểu (38 vs 12).

## Chốt & contrast với PRRSV

| | PRRSV (CDS / frameshift) | FMDV (polyprotein / mat_peptide) |
|---|---|---|
| Bug tool | ORF7 (start N-term), ORF1a (stop over-read) → **đã sửa** | **Không có** |
| Cải thiện | 3 fix rescue | gap-fill (cứu peptide ngắn) |
| Residual | ORF1b = frameshift start convention | boundary/granularity convention ở cleavage |
| Câu chuyện | validation **tìm & sửa** bug | validation **phân biệt convention khỏi lỗi**, cho thấy tool **chuẩn hóa đúng** |

Hai virus bổ trợ nhau: PRRSV chứng minh validation debug được tool; FMDV chứng minh tool xử lý đúng
cấu trúc polyprotein và **vấn đề còn lại nằm ở sự không nhất quán của GenBank** — đúng thứ tool sinh ra để giải quyết.
