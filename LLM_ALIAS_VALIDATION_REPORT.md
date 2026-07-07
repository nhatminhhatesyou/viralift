# LLM Alias-Review Validation Report

Report này tóm tắt validation cho tính năng **LLM-assisted alias review**
(`app/src/llm/`), dựa trên script tại `app/validation/07_llm_alias_validation/`.

Dữ liệu dùng: tái sử dụng data PED đã validate ở `06_ped_validation/` (100
query records, alias map đã được validate actionable coverage 99.39%) làm
ground truth, thay vì chạy lại tblastn từ đầu — vì độ chính xác tblastn đã
được validate riêng ở mức 99.01% (xem `PED_VALIDATION_REPORT.md`). Validation
này cô lập và chỉ kiểm tra lớp **LLM naming/classification review** nằm trên
kết quả tblastn, không re-test lại phần coordinate lifting.

## 1. Mục tiêu validation

Tách thành 2 track, giống tinh thần 2-lớp của PED validation:

1. **Track B — Uncertain coordinate-supported suggestions**: các raw
   qualifier value thật từ 100 query record PED, đã qua `classify_alias_candidate`
   (bộ phân loại rule-based có sẵn) và lọc bằng đúng hàm `needs_llm_review()`
   thật trong code — chỉ những dòng "khó xử" mới được gửi LLM, giống hệt luồng
   thật trong `bootstrap_alias.py`.
2. **Track A — Known unresolved/ambiguous names**: 5 tên đã biết rõ đáp án
   đúng, lấy từ `ped_alias_names_to_review.tsv` (đã có trong PED validation
   trước), chạy qua `review_unresolved_names()` — đúng hàm mà `resolve.py`
   dùng thật.

Cả 2 track đều gọi **API OpenAI thật** (model theo cấu hình `.env`,
`gpt-5.4-mini` với fallback `gpt-5.4`), không dùng mock.

## 2. Track A — Known unresolved/ambiguous names

| Raw name | Đáp án đúng | LLM action | LLM canonical | Confidence | Đúng? |
|---|---|---|---|---|---|
| `101-bp deletion results in frameshift...` | ignore | ignore | — | high | ✅ |
| `1a polyprotein and 1b polyprotein` | save_alias → ORF1ab | save_alias | ORF1ab | high | ✅ |
| `contains ORF1a and ORF1b` | save_alias → ORF1ab | save_alias | ORF1ab | high | ✅ |
| `truncated hypothetical protein` | ignore | ignore | — | high | ✅ |
| `mp` | move_to_ambiguous | move_to_ambiguous | — | high | ✅ |

**Kết quả: 5/5 (100%) đúng, toàn bộ ở confidence "high".** Đáng chú ý nhất:
LLM nhận diện đúng cả 2 cách diễn đạt khác nhau của "ORF1a gộp ORF1b" đều nên
map về canonical riêng `ORF1ab` (đúng theo annotation granularity mismatch đã
ghi nhận ở PED report), và nhận diện đúng `mp` là tên dùng chung nhiều gene
cần giữ ambiguous — đúng như quyết định thủ công đã ghi trong
`PED_VALIDATION_REPORT.md`.

## 3. Track B — Uncertain coordinate-supported suggestions

Trên 59 suggestion row có coordinate support (IoU ≥ 0.90), 20 row bị heuristic
đánh dấu "cần LLM review". Overall:

| Metric | Giá trị |
|---|---:|
| Action accuracy (toàn bộ 20 row) | 80% (16/20) |
| Canonical accuracy khi LLM nói `save_alias` | 100% (14/14) |
| Dangerous false positive (LLM bảo save nhưng sai) | 0 |

### Điểm quan trọng nhất: accuracy trên phần thực sự được áp dụng

App hiện tại (`ui/stages/bootstrap_alias.py`) **chỉ tự động điền sẵn Action
(lúc chạy validation này còn là checkbox, giờ đã đổi thành dropdown — xem
update ở mục 4) khi LLM confidence là `medium` hoặc `high`** — recommendation
ở confidence `low` không bao giờ được auto-apply. Khi lọc theo đúng ngưỡng
này:

| Metric | Giá trị |
|---|---:|
| Recommendation ở confidence medium/high | 16/20 |
| Trong số đó, action đúng | **16/16 (100%)** |
| Trong số đó, canonical đúng (khi save_alias) | 14/14 (100%) |

**Toàn bộ 4 trường hợp LLM đoán sai đều tự gắn cờ confidence `low`** — đúng
ngưỡng mà app đã có sẵn để không auto-apply. Nói cách khác: với cấu hình
confidence-gating hiện tại, LLM assist **chưa từng đưa ra một auto-decision
sai nào trong lần test này** (0 dangerous false positive, 0 wrong save, 0
wrong ignore ở mức confidence được tin dùng).

### 4 trường hợp LLM đoán sai (đều ở confidence low)

| Raw value | Field | Đáp án đúng | LLM action | Lý do LLM đưa ra | Đánh giá |
|---|---|---|---|---|---|
| `HNZK1` | gene | ignore | skip | "not a recognized PED gene/ORF synonym... despite strong coordinate overlap" | Hedge an toàn — skip = không làm gì, không gây hại |
| `mp` | gene | move_to_ambiguous | skip | "too abbreviated... evidence is weak" | Hedge an toàn — cùng effect với ambiguous trong UI hiện tại (xem mục 4) |
| `small membrane protein` | note | save_alias → E | ignore | "generic descriptive phrase, should not be saved" | Bỏ lỡ alias hợp lệ — quá thận trọng với từ "membrane" |
| `accessory membrane protein` | product | save_alias → ORF3 | ignore | "generic/descriptive, not a clear specific alias" | Bỏ lỡ alias hợp lệ — cùng nguyên nhân |

2/4 case là hedge an toàn (không đổi gì so với hành vi deterministic gốc,
không gây hại). 2/4 case còn lại là **false negative** (bỏ lỡ alias đúng),
không phải false positive — nghĩa là chi phí của sai số này là "vẫn cần user
tự thêm alias đó," không phải "config bị lưu sai."

### Case đáng chú ý nhất: LLM sửa đúng một lỗi của bộ coordinate-matching

Raw value `ORF1ab` (field `gene`, 5 record support) — do đặc thù annotation
(một số record ghi `ORF1ab` nhưng tọa độ CDS chỉ trùng vùng `ORF1a`), bộ
matching theo tọa độ tự động gán ứng viên canonical **sai** thành `ORF1a`.
Deterministic classifier dựa vào tọa độ này cho confidence cao (`score=13`).
LLM, dựa trên `matching_available_canonical` (raw text khớp thẳng với
canonical `ORF1ab` có sẵn), đã **ghi đè đúng thành `ORF1ab`** thay vì tin theo
tọa độ, ở confidence `high`. Đây là đúng loại lỗi "annotation granularity
mismatch" mà PED report đã cảnh báo (mục 5.2) — và là ví dụ cụ thể cho thấy
LLM review bổ sung giá trị thật, không chỉ lặp lại deterministic.

## 4. Phát hiện quan trọng, độc lập với số liệu accuracy

> **Update (sau khi validation này chạy):** gap mô tả bên dưới ở
> `bootstrap_alias.py` **đã được fix**. UI giờ có 1 dropdown "Action" với 4
> lựa chọn (`Save alias` / `Save ambiguous` / `Save ignored` / `Skip this
> run`) thay vì 4 checkbox độc lập; `Save ambiguous` ghi thẳng vào
> `ambiguous_names` qua `apply_approved_alias_suggestions(...,
> ambiguous_rows=...)`. Gap ở `resolve.py` (mô tả ở đoạn dưới) **vẫn còn mở**
> — dropdown ở đó vẫn chỉ có canonical hoặc "ignore", chưa có lựa chọn đánh
> dấu ambiguous. Giữ nguyên phần dưới làm bản ghi lại phát hiện gốc tại thời
> điểm validation.

`ui/stages/bootstrap_alias.py` hiện gộp cả `skip` và `move_to_ambiguous` từ
LLM vào **cùng một checkbox "Skip"**, mà "Skip" thì không ghi gì vào
`ignored_names` hay `ambiguous_names` cả (là no-op thuần túy). Tương tự,
`ui/stages/resolve.py` không có lựa chọn "đánh dấu ambiguous" trong dropdown —
chỉ có canonical hoặc "ignore" (chọn ignore ở đây cũng không ghi gì vào
config, decisions[rep] = None).

Hệ quả: dù Track A cho thấy LLM nhận diện `mp` là ambiguous với confidence
`high` và 100% đúng, **hiện tại không có đường nào trong UI để quyết định đó
được lưu lại thành `ambiguous_names`** — nó chỉ im lặng biến mất, và tên đó
sẽ lại hiện ra "chưa xác định" ở lần chạy sau. Đây là gap kỹ thuật đáng sửa,
tách biệt hoàn toàn khỏi câu hỏi "LLM có chính xác không."

## 5. Kết luận

1. **LLM alias review hoạt động rất tốt trên PED**: 100% đúng trên Track A;
   trên Track B, 100% đúng ở đúng ngưỡng confidence mà app auto-apply
   (medium/high), 0 dangerous false positive.
2. **Sai số duy nhất là false negative ở confidence thấp** — an toàn hơn là
   có hại: chi phí là vẫn cần user tự xử lý thủ công, không phải config bị
   lưu sai.
3. **LLM có thể sửa đúng lỗi mà bộ coordinate-matching mắc phải** (case
   `ORF1ab`), cho thấy giá trị gia tăng thật, không chỉ trùng lặp deterministic.
4. **Gap cần sửa**: cả `bootstrap_alias.py` và `resolve.py` hiện không có
   đường lưu quyết định "ambiguous" xuống config, dù LLM (và cả user) có xác
   định đúng. Nên thêm: (a) ở `bootstrap_alias.py`, tách `move_to_ambiguous`
   ra khỏi bucket "skip" thành một hành động ghi vào `ambiguous_names`; (b) ở
   `resolve.py`, thêm lựa chọn "Mark as ambiguous" vào dropdown thay vì chỉ có
   canonical/ignore.
5. **Mẫu còn nhỏ** (20 + 5 = 25 case, chỉ trên 1 virus/PED) — nên xem đây là
   kết quả khả quan ban đầu, không phải kết luận tổng quát. Bước tiếp theo nên
   lặp lại trên virus thứ 2 (giống cách PED dùng cả ref_1/ref_2) trước khi coi
   tính năng này là production-ready cho mọi virus.
