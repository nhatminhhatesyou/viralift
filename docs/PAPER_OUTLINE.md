# ViraLift — khung paper (IMRaD) + trạng thái từng phần

Ký hiệu trạng thái:
- ✅ **CÓ** — đã có số liệu/tài liệu, chỉ cần viết lại thành prose học thuật
- 🟡 **MỘT PHẦN** — có chất liệu nhưng chưa đủ / số cũ / cần bổ sung
- ❌ **CHƯA** — chưa có, cần làm

Nguồn viết tắt: `LV3` = `LIFTING_VALIDATION_3_VIRUSES.md`, `VDF` = `VALIDATION_DRIVEN_FIXES.md`,
`FMD` = `FMDV_FAILURE_ANALYSIS.md`, `PED` = `PED_VALIDATION_REPORT.md`, `LLM` = `LLM_ALIAS_VALIDATION_REPORT.md`,
`AEGIS` = `AEGIS_VIRALIFT_COMPARISON.md`, `ROAD` = `ROADMAP_gene_hierarchy.md`.

---

## Title / Abstract

- **Title** ❌ — chưa chốt. Gợi ý: nhấn "coordinate-driven standardisation of inconsistent viral gene
  annotations for primer design".
- **Abstract** ❌ — viết sau cùng khi đã có toàn bộ số.

---

## 1. Introduction

**1.1 Bài toán** ✅ (README, AEGIS §8)
- GenBank cho cùng một gene virus bị đặt tên rất khác nhau giữa các lab (`/gene`, `/product`, `/note`);
  và nhiều record **không có annotation**. Hai vấn đề này chặn bước hạ nguồn (VD trích CDS ORF5 cho cả
  bộ input trong pipeline thiết kế primer).

**1.2 Vì sao công cụ chung không giải quyết được** 🟡 (AEGIS — có so sánh, cần rút gọn thành 1–2 đoạn)
- Công cụ annotation phổ thông (vd AEGIS) làm gì, và vì sao không hợp với bài toán "chuẩn hóa tên +
  gán annotation theo một reference cho một loài virus".

**1.3 Đóng góp của bài** ❌ — cần viết, nhưng nội dung đã rõ:
- (i) một pipeline lift annotation dựa trên tọa độ (tblastn) + chuẩn hóa tên theo alias map per-virus;
- (ii) một khung **validation-driven** phát hiện và sửa được bug thật của chính công cụ;
- (iii) đánh giá trên **3 virus cấu trúc genome khác nhau** (arterivirus, picornavirus, coronavirus).

---

## 2. Methods

**2.1 Tổng quan pipeline** 🟡 (README, PIPELINE_RUNNER_GUIDE, CODEBASE_GUIDE — có mô tả, chưa thành hình paper)
- Routing: record có annotation → direct extraction; record không có → tblastn lifting từ reference.
- Reference quyết định mức làm việc của cả run (CDS vs mature peptide).

**2.2 Chuẩn hóa tên (alias standardisation)** 🟡 (ALIAS_MANAGER_GUIDE, code)
- Alias map per-virus: canonical + alias; `excluded_names`. Cần một hình minh họa cấu trúc.
- Cơ chế gợi ý alias 2 tầng: tọa độ (tblastn + IoU) → LLM review cho ca chưa đủ điểm. **← đang hoàn thiện**

**2.3 Annotation lifting bằng tblastn** ✅ (LV3 §Methods, VDF)
- Protein-guided lift; HSP merging; ngưỡng coverage/identity; terminal extrapolation; gap-fill polyprotein.
- −1 ribosomal frameshift, mature-peptide, các cơ chế rescue leakage-free.

**2.4 Thiết kế validation** ✅ (LV3 §"Thiết kế validation")
- Mẫu số truth-based; 4 bucket exact/coord-only/failed/no-hit; IoU ≥ 0.90 + bp-tolerance + codon check.
- Nguyên tắc **leakage-free** (chỉ dùng ref protein + tọa độ, không dùng truth) và **generic** (không
  hardcode theo virus/gene).
- **Nguyên tắc quan trọng:** validation notebook gọi đúng hàm pipeline thật, không reimplement.

**2.5 Dữ liệu** 🟡 (rải rác) — ❌ cần một bảng gọn:
- 3 virus × 100 record query; reference dùng; nguồn (GenBank accession ranges); PEDV chạy 2 reference.

---

## 3. Results

### 3.1 Annotation lifting accuracy ✅ **(phần mạnh nhất, gần xong)** — nguồn LV3

- **Bảng chính:** 3 virus, 2917 prediction, **accuracy 99.79%**, 6 failed, chỉ **1 ca blame tool**.
  (Fig 1 overall, Fig 2 per-gene — đã có trong `docs/figures/`.)
- **3.1.1 Bốn bug tool validation tìm ra & sửa** ✅ (LV3, VDF): ORF7 start rescue; ORF1a over-read;
  validator internal-stop; **HSP-merge span** (Fig 4 trước/sau). Đều generic + leakage-free.
- **3.1.2 Phần dư không phải lỗi lifting** ✅ (LV3):
  - Gene không có start codon (ORF1b): `delta_end ≡ 0` trên ~96 ca của 2 virus độc lập (Fig 3).
  - Biên cắt polyprotein FMDV (2A/VP1): standardisation, không phải bug (Fig 5).
  - Lỗi annotation trong ground truth (AF331831 tên trùng; KX550281 nhãn sai).
- **3.1.3 Chủng tái tổ hợp nhận diện đúng** ✅ (LV3) — MF577027 Belgorod, S 60% mà IoU 0.99.

> Trạng thái: **✅ gần hoàn chỉnh.** Chỉ cần chuyển prose sang tiếng Anh học thuật + đồng bộ số FMDV 99.83%.

### 3.2 Alias suggestion / config reconstruction 🟡 **(đang làm dở — mảng khuyết lớn nhất)**

- **Ý tưởng:** đưa tool virus mới (chỉ canonical từ ref) → nó dựng lại config tới đâu so với gold.
- **Đo hai chiều per-canonical:** precision (tool save đúng?) + recall (gate corpus; bắt bug kiểu sM).
- **Đã có:**
  - Harness gọi đúng hàm pipeline thật (coordinate → LLM → apply). ✅
  - FMDV: **precision 100% / recall 100%** sau khi sửa. ✅
  - Cơ chế phụ đã tìm ra: cross-canonical auto-exclude (HNZK1), gate "canonical không trong ref". ✅
- **Chưa xong:** 🟡
  - PRRSV/PEDV số chưa chốt — đang truy 3 vấn đề: (1) gate precision cho canonical-không-trong-seed;
    (2) classifier deterministic-ignore tên mô tả trước khi LLM thấy; (3) compound `;` variants.
  - Chưa có bảng kết quả 3 virus sạch cho mục này.
- **Finding phụ đã rõ (đáng đưa vào):** production cap `max_rows=20` → virus nhiều tên cần nhiều lượt;
  classifier ưu tiên hình-dạng-tên hơn tọa độ (pattern sM/mp/L/leader).

> Trạng thái: **🟡 số chưa chốt.** Đây là việc cần làm tiếp trước khi viết Results 3.2.

### 3.3 (Tùy chọn) LLM alias review — track A/B 🟡 (LLM report)
- Có report riêng nhưng N nhỏ (10 mục Track A). Có thể gộp vào 3.2 hoặc bỏ. Cân nhắc.

---

## 4. Discussion

**4.1 Ba luận điểm phương pháp** ✅ (LV3 §"Ba luận điểm")
- `exact_pct` không đo chất lượng lifting với gene không start codon (đo đồng thuận quy ước).
- IoU phụ thuộc kích thước → cần bp-tolerance; accuracy tăng do **metric công bằng hơn**, không phải
  lifting tốt hơn (phải nêu rõ, gồm cả vụ tách 3B).
- Nhiều reference khác quy ước = một thiết kế validation (đo trực tiếp mức bất nhất của GenBank).

**4.2 Tọa độ mạnh hơn tên** 🟡 (chất liệu rải, cần gom)
- Lifting 99.79% (không dùng tên) + việc tọa độ tự tìm ra sM→E, mp→ORF3 mà validation cũ bỏ lọt.
- Ranh giới: khi nào tọa độ đủ, khi nào cần LLM/người (tên mô tả mơ hồ).

**4.3 Giới hạn** 🟡/❌
- Gold standard kế thừa giới hạn của lifting (gene phân kỳ, peptide ngắn).
- Bước người không tự động được → số reconstruction là cận trên.
- Chưa lift được sub-feature ref không chứa (nsp) — future work.

**4.4 Future work** ✅ (ROAD)
- Hệ phân cấp gene (parent/child) phát hiện bằng tọa độ; điều kiện: reference có annotate mức mịn.

---

## 5. Conclusion ❌ — viết sau.

---

## Figures / Tables — trạng thái

| # | Nội dung | Trạng thái |
|---|---|---|
| Fig 1 | Accuracy overall 3 virus | ✅ `docs/figures/fig1_overall.png` |
| Fig 2 | Accuracy per-gene | ✅ `fig2_per_gene.png` |
| Fig 3 | ORF1b delta_end≡0 (2 virus) | ✅ `fig3_orf1b.png` |
| Fig 4 | Bug HSP-merge trước/sau | ✅ `fig4_hsp_bug.png` |
| Fig 5 | FMDV biên VP1/2A | ✅ `fig5_fmdv_boundary.png` |
| Fig | Sơ đồ pipeline (routing + 2 tầng alias) | ❌ chưa có |
| Fig | Cấu trúc alias map / ví dụ chuẩn hóa | ❌ chưa có |
| Fig | Kết quả alias reconstruction 3 virus | ❌ chờ số 3.2 |
| Bảng | Dataset (virus × record × reference × accession) | ❌ chưa gom |

---

## Tóm tắt: cái gì xong, cái gì thiếu

**Gần hoàn chỉnh (viết prose là được):**
- Toàn bộ **Results 3.1 (lifting accuracy)** + 5 figure.
- **Discussion 4.1** (ba luận điểm phương pháp).
- **Methods 2.3, 2.4** (lifting + thiết kế validation).
- Future work (4.4).

**Còn thiếu / phải làm trước khi viết:**
1. **Results 3.2 (alias reconstruction)** — chốt số PRRSV/PEDV (đang sửa 3 vấn đề). ← chặn lớn nhất.
2. **Sơ đồ pipeline** + hình cấu trúc alias (Methods).
3. **Bảng dataset** (accession, số record, reference).
4. **Intro 1.3** (đóng góp) + **Abstract** + **Conclusion** — viết khi số đã đủ.
5. Đồng bộ số FMDV 99.83% ở mọi tài liệu cũ.
6. Quyết định: có gộp LLM track A/B (3.3) không.

**Đường tới hạn:** phần lifting đủ để viết ngay; phần alias là nửa còn lại của câu chuyện và **đang dở** —
nên hoàn thiện validation alias trước, rồi mới ráp Results hoàn chỉnh.
