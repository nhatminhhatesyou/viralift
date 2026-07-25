# Kế hoạch: hệ thống phân cấp gene cho alias

Tài liệu này ghi lại thiết kế cho việc mô hình hóa **quan hệ cha–con giữa các gene** trong ViraLift,
để phát triển về sau. Hiện tại **hoãn lại** theo phản hồi của lab (xem mục "Vì sao hoãn").

## Bối cảnh

Annotation của virus tồn tại ở nhiều **mức granularity** cùng lúc:

```
ORF1ab                     (cả replicase, một CDS)
 ├─ ORF1a                  (nửa 5')
 │   └─ nsp1 … nsp8        (mature peptide, protease cắt ra)
 └─ ORF1b                  (nửa 3')
     └─ nsp9 … nsp12
```

Cùng một vùng genome được các lab khác nhau khai ở các mức khác nhau: có record ghi `ORF1ab`, có record
tách `ORF1a`+`ORF1b`, một số ít ghi hẳn `nsp1…nsp12`. Reference thì thường chỉ ở **một mức** (PEDV/PRRSV
ref là CDS; FMDV ref là mature peptide).

Alias config phẳng không diễn đạt được quan hệ này: `NSP2` và `ORF1a` chỉ là hai tên ngang hàng, không có
gì nói `NSP2` nằm *bên trong* `ORF1a`. Hệ quả:

- Không kiểm chứng được bằng tọa độ (một feature tên `nsp2` lẽ ra *phải* nằm trong vùng `ORF1a` đã lift).
- Bật đồng thời hai mức làm `select_feature_type` phân vân, và metric lifting bị nhập nhằng vì cùng một
  vùng khớp ở hai tầng.

## Ý tưởng cốt lõi: phát hiện quan hệ cha–con bằng tọa độ

Reference thô (chỉ ORF1a/ORF1b) vẫn có thể **tự phát hiện** cấu trúc con, nếu record query có annotate:

1. Lift các feature reference lên genome query (đã có, độ chính xác ~99.7%).
2. Với mỗi feature **mịn** của query (VD `nsp2`), kiểm xem tọa độ của nó có **nằm trọn bên trong** một
   feature reference đã lift (`ORF1a`) không.
3. Nếu có → `nsp2` là con của `ORF1a`. Tự tạo canonical `NSP2` với `parent = ORF1a`, không cần người gán.

Đây là suy luận thuần hình học: `is_subfeature_of` theo quan hệ chứa đựng. Generic (không hardcode theo
virus/gene) và leakage-free (chỉ dùng tọa độ + alignment, không đụng truth). Cùng loại logic đã dùng để
phân biệt `RNA-dependent RNA polymerase` (mat_peptide nằm trong ORF1ab → sub-part) khỏi một gene thật.

## Giới hạn cứng — phải nêu rõ

Canonical con phát hiện kiểu này **chỉ tồn tại ở record nào tự annotate nó**. Vì reference không chứa nsp,
**không có protein tham chiếu để lift `NSP2` sang record khác**. Do đó:

- ✅ **Thu hoạch được**: gom nsp từ các record có annotate, gắn đúng parent, và trích xuất chúng ra.
- ❌ **Không lan được**: không suy được nsp cho record *không* annotate — không có gì để lift.

Muốn lan thật thì cần một reference **có annotate mature peptide**. Đó là quyết định về **dữ liệu**, không
phải về code.

Giá trị thực tế của cơ chế: *"tool tự phát hiện annotation trong corpus có cấu trúc phân cấp và ghi lại đúng
quan hệ, dù reference chỉ ở mức thô"* — một capability đáng đưa vào paper, nhưng không phải công cụ để lan
annotation nsp.

## Vì sao hoãn

**Lab xác nhận: với mục tiêu thiết kế primer, không cần các mảnh nsp.** Primer nhắm vào gene cấu trúc
(S, ORF5/GP5, M, N…) và các ORF trọn, không nhắm vào từng mature peptide của replicase. Nên toàn bộ tầng
nsp — và do đó phần lớn nhu cầu về hệ phân cấp — **nằm ngoài phạm vi hiện tại**.

Hệ quả cho config đang có: các canonical `NSP*` t đã thêm (PEDV, PRRSV) hiện là **thông tin thừa so với mục
tiêu**. Chúng không gây hại (lookup phẳng không đổi, không lift được nên không sinh prediction), nhưng cũng
không phục vụ gì cho primer. Hai lựa chọn khi quay lại:

- **Giữ**: coi như catalog sẵn cho tương lai; đánh dấu là "harvest-only, không dùng cho primer".
- **Rút gọn**: bỏ canonical NSP khỏi config, đưa nsp về lại `excluded_names`, giữ config tinh gọn đúng
  phạm vi primer. Các fix thật (nsp1β, accessory membrane protein → ORF3, tách 3B, gộp ORF2, khóa mức theo
  ref) **không phụ thuộc nsp** nên vẫn nguyên.

Khuyến nghị: **rút gọn** cho tới khi có use case thật cần nsp, để config phản ánh đúng phạm vi.

## Đã có sẵn (scaffold, chưa kích hoạt)

- **Schema `parent`**: parser đọc được cả hai dạng (list cũ / object `{aliases, parent}`);
  `build_parent_map` kiểm tra cha tồn tại + không chu trình; `build_children_map` cho chiều đọc ngược.
  Config cũ vẫn hợp lệ. `parent` hiện là **metadata trơ** — chưa code lift nào đọc.
- **Logic chứa đựng**: `is_subfeature_of` trong `alias_payload.build_position_context` — đã dùng cho alias
  suggestion, tái dùng được cho auto-detect.
- **Khóa mức theo reference**: `select_feature_type(..., allowed_types=(ref_feature_type,))` — reference
  quyết định mức của cả run, nên "record có cả CDS lẫn mat_peptide" không còn mơ hồ (mức của ref thắng).

## Cần làm (khi kích hoạt)

1. **Auto-detect cha–con**: sau khi lift, quét feature query mịn nằm trong feature reference thô → đề xuất
   canonical con + `parent`. Nối `is_subfeature_of` vào luồng suggestion.
2. **Lift đa mức**: cho một run xuất cả hai tầng, mỗi output gắn nhãn mức + parent; consumer chọn tầng
   (primer → mức CDS; phân tích nsp → mức mat_peptide).
3. **Metric theo tầng**: chấm accuracy riêng từng mức, tuyệt đối không gộp (cùng vùng khớp ở hai tầng sẽ
   làm accuracy vô nghĩa) — giống cách xử lý ORF1ab vs ORF1a/ORF1b hiện tại.
4. **UI set/hiển thị parent**: trường để gán cha, và sơ đồ cây quan hệ. Hiện chưa có.
5. **Reference có nsp** (quyết định dữ liệu): điều kiện tiên quyết để *lan* nsp thay vì chỉ thu hoạch.

## Tóm tắt một dòng

Hạ tầng cho hệ phân cấp gene đã dựng sẵn ở dạng trơ (`parent`, logic chứa đựng, khóa mức); việc kích hoạt
bị hoãn vì primer không cần nsp; khi cần, cơ chế phát hiện là **quan hệ chứa đựng theo tọa độ**, với giới
hạn cứng là không lan được cái mà reference không chứa.
