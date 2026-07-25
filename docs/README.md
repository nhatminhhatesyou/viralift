# Tài liệu nội bộ — validation & paper

Thư mục này chứa **ghi chú nội bộ** phục vụ bài báo Q1 và việc validation. Viết bằng tiếng Việt, dành cho nhóm phát triển — **không phải** hướng dẫn sử dụng.

> Hướng dẫn sử dụng sản phẩm (tiếng Anh) nằm ở thư mục gốc repo: [README](../README.md), [Data Crawler](../DATA_CRAWLER_GUIDE.md), [Pipeline Runner](../PIPELINE_RUNNER_GUIDE.md), [Alias Manager](../ALIAS_MANAGER_GUIDE.md), [Codebase Guide](../CODEBASE_GUIDE.md).

## Bài báo (paper)

| File | Nội dung |
|---|---|
| `VIRALIFT_PAPER_OUTLINE.md` (nằm ở thư mục cha `Phase_2/`, ngoài repo) | Outline bài Q1: framing, novelty, cấu trúc mục, asset dùng cho từng phần |
| [REBUILD_PLAN.md](REBUILD_PLAN.md) | Kế hoạch dựng lại notebook validation cho paper (mapping notebook → figure/table). Index notebook chuẩn: `app/validation/README.md` |

## Kết quả validation

| File | Virus / phạm vi | Nội dung |
|---|---|---|
| [PED_VALIDATION_REPORT.md](PED_VALIDATION_REPORT.md) | PEDV | Accuracy alias + tblastn trên 100 record (2 reference) |
| [LLM_ALIAS_VALIDATION_REPORT.md](LLM_ALIAS_VALIDATION_REPORT.md) | LLM alias review | Validation lớp LLM review (chồng lên kết quả PED) |
| [VALIDATION_DRIVEN_FIXES.md](VALIDATION_DRIVEN_FIXES.md) | PRRSV | 4 bug tool mà validation phát hiện → fix generic, leakage-free |
| [FMDV_FAILURE_ANALYSIS.md](FMDV_FAILURE_ANALYSIS.md) | FMDV | Phân biệt "lỗi tool" vs "convention biên cleavage" trên polyprotein |

## Nền tảng / so sánh

| File | Nội dung |
|---|---|
| [AEGIS_VIRALIFT_COMPARISON.md](AEGIS_VIRALIFT_COMPARISON.md) | Giải thích FASTA/GFF3/GenBank + so sánh use case AEGIS vs ViraLift |
