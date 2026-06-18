# Hướng dẫn Alias Manager trong ViraLift

Tài liệu này giải thích module Alias Manager: nó dùng để làm gì, khi nào cần dùng, và cách đọc/sửa alias map cho virus mới hoặc virus đã có config.

## Mục lục

- [Alias map là gì?](#alias-map-là-gì)
- [Vì sao cần Alias Manager?](#vì-sao-cần-alias-manager)
- [Các file chính](#các-file-chính)
- [Flow khi tool gặp virus đã biết](#flow-khi-tool-gặp-virus-đã-biết)
- [Flow khi tool gặp virus mới](#flow-khi-tool-gặp-virus-mới)
- [Cách đọc Alias Manager UI](#cách-đọc-alias-manager-ui)
- [Ý nghĩa các loại tên](#ý-nghĩa-các-loại-tên)
- [Cách tool gợi ý alias cho virus mới](#cách-tool-gợi-ý-alias-cho-virus-mới)
- [Granularity mismatch](#granularity-mismatch)
- [Ví dụ thực tế](#ví-dụ-thực-tế)
- [Các cảnh báo thường gặp](#các-cảnh-báo-thường-gặp)
- [Best practices](#best-practices)

## Alias map là gì?

Trong GenBank, cùng một gene có thể được ghi bằng nhiều tên khác nhau:

| Canonical name | Alias có thể gặp |
|---|---|
| `ORF5` | `GP5`, `ORF5 protein`, `major envelope glycoprotein` |
| `N` | `nucleoprotein`, `N protein` |
| `S` | `spike protein`, `S protein` |
| `ORF1ab` | `ORF1a/1b`, `ORF1a/b`, `Pol1`, `polyprotein 1ab` |

ViraLift cần chuẩn hóa các tên này về một tên chuẩn duy nhất gọi là `canonical name`.

Ví dụ:

```text
GP5 -> ORF5
ORF5 protein -> ORF5
major envelope glycoprotein -> có thể ignore nếu quá mô tả/generic
```

Alias map là theo từng virus, không phải dùng chung cho mọi virus. Một tên như `envelope protein` có thể quá mơ hồ ở virus này, nhưng lại map ổn định vào `E` ở virus khác nếu reference/query evidence chứng minh rõ.

## Vì sao cần Alias Manager?

Alias Manager giúp user:

- Xem virus hiện đang có alias config nào.
- Sửa canonical name và alias nếu tool map sai.
- Xóa alias bị thêm nhầm.
- Quản lý `ignored_names` và `ambiguous_names`.
- Sửa keyword dùng để tự động nhận diện virus.
- Tạo alias config cho virus mới từ reference + query annotation.

Nói ngắn gọn: đây là nơi quản lý "từ điển tên gene" của từng virus.

## Các file chính

| File | Vai trò |
|---|---|
| `app/config/virus_alias_registry.json` | Registry liệt kê virus, keyword nhận diện, và đường dẫn alias config |
| `app/config/*_alias.json` | Alias config riêng cho từng virus |
| `app/src/alias/alias_manager.py` | Hàm đọc/sửa/validate alias config |
| `app/src/alias/alias_bootstrap.py` | Tạo alias config và gợi ý alias cho virus mới |
| `ui/streamlit_app.py` | UI Alias Manager và flow tạo alias seed |

Mỗi lần lưu alias config qua UI, tool tạo backup trong:

```text
app/config/backups/
```

Folder backup này là artifact runtime và đã được ignore khỏi Git. Nó chỉ dùng để khôi phục nếu user sửa nhầm alias config.

## Flow khi tool gặp virus đã biết

Khi user upload reference, tool đọc metadata của reference, ví dụ:

```text
organism = Porcine reproductive and respiratory syndrome virus
```

Sau đó tool so với `keywords` trong `virus_alias_registry.json`.

Ví dụ registry:

```json
{
  "virus_name": "PRRSV",
  "keywords": [
    "porcine reproductive and respiratory syndrome virus",
    "prrsv",
    "prrs virus"
  ],
  "alias_config": "app/config/prrsv_alias.json"
}
```

Nếu metadata chứa keyword tương ứng, tool sẽ tự chọn alias config đó.

Sau đó pipeline chạy bình thường:

1. Parse feature từ reference.
2. Chuẩn hóa tên bằng alias config.
3. Parse query.
4. Nếu query thiếu annotation hữu ích thì dùng tblastn lifting.
5. Xuất kết quả đã chuẩn hóa tên.

## Flow khi tool gặp virus mới

Nếu reference không match virus nào trong registry, tool sẽ vào flow virus mới:

1. **Virus review**  
   Tool hiển thị metadata lấy từ reference như `organism`, `description`, `record id`.

2. **User chọn một trong hai hướng**

   - Chọn alias config đã có nếu đây thật ra là virus cũ nhưng keyword chưa có.
   - Tạo alias config mới nếu đây là virus mới.

3. **Alias seed**

   Với virus mới, ViraLift lấy tên feature trong reference làm canonical name ban đầu.

   Ví dụ reference PED có:

   ```text
   ORF1a, ORF1b, S, ORF3, E, M, N
   ```

   Thì alias config mới sẽ có canonical:

   ```json
   {
     "canonical_names": {
       "ORF1a": [],
       "ORF1b": [],
       "S": [],
       "ORF3": [],
       "E": [],
       "M": [],
       "N": []
     }
   }
   ```

4. **Generate suggestions**

   Tool chạy tblastn từ reference sang từng query record, sau đó so tọa độ với annotation thật trong query để tìm tên nào nên thêm làm alias.

5. **User approve**

   Gợi ý nào hợp lý thì tick `save`. Gợi ý nào generic hoặc sai thì để `ignore`.

6. **Save config**

   Tool lưu file alias config mới vào `app/config/` và thêm virus vào registry.

## Cách đọc Alias Manager UI

Trong sidebar chọn:

```text
Alias manager
```

Các tab chính:

### Registry

Dùng để sửa:

- `virus_name`: tên hiển thị cho user.
- `keywords`: từ khóa dùng để tự động nhận diện virus.

Lưu ý: auto-detect dựa chủ yếu vào `keywords`, không chỉ dựa vào `virus_name`.

### Canonical aliases

Mỗi canonical name có một khung riêng.

Ví dụ:

```text
E · 4 alias(es)
```

Bên trong là danh sách alias của `E`.

User có thể:

- Tick một hoặc nhiều alias.
- Bấm `Delete selected` để xóa alias đó.
- Tick `Delete canonical` để xóa cả canonical name.
- Thêm alias mới ở bảng `Add canonical / alias`.

### Ignored names

Chứa các tên không nên dùng để map gene.

Ví dụ:

```text
protein
glycoprotein
unknown protein
replicase polyprotein
```

Các tên này thường quá chung chung. Nếu đưa vào alias map, tool có thể map nhầm nhiều gene khác nhau.

Lưu ý: `envelope protein`, `membrane protein`, `nucleocapsid protein` không phải lúc nào cũng phải ignore. Với PED, các tên này có evidence rõ và được map lần lượt vào `E`, `M`, `N`.

### Ambiguous names

Chứa các tên có thể map sang nhiều gene khác nhau.

Ví dụ:

```text
envelope protein
glycosylated membrane protein
```

Nếu một tên có thể vừa giống ORF2, ORF5, ORF6 tùy virus/record, nên để ambiguous hoặc manual review.

Ví dụ PED có raw gene `mp` trong một số record. `mp` có thể hiểu là ORF3 accessory membrane protein, nhưng cũng dễ bị nhầm với membrane protein. Vì vậy để `mp` trong `ambiguous_names` sẽ an toàn hơn nếu nó xuất hiện một mình.

### Raw JSON

Hiển thị nguyên file alias config. Dùng để debug nhanh.

## Ý nghĩa các loại tên

### Canonical name

Tên chuẩn cuối cùng mà ViraLift muốn output.

Ví dụ:

```text
ORF5
```

### Alias

Tên khác nhưng đủ cụ thể để map về canonical.

Ví dụ:

```text
GP5 -> ORF5
ORF5 protein -> ORF5
```

### Ignored name

Tên nên bỏ qua vì không đủ thông tin.

Ví dụ:

```text
protein
glycoprotein
unknown protein
replicase polyprotein
```

### Ambiguous name

Tên có thông tin nhưng không đủ chắc để map vào một canonical duy nhất.

Ví dụ:

```text
glycosylated membrane protein
```

Tên này có thể chỉ các gene khác nhau tùy virus hoặc annotation convention.

## Cách tool gợi ý alias cho virus mới

Khi bấm `Generate suggestions`, tool làm như sau:

1. Chọn feature type hữu ích trong query, ví dụ `CDS` hoặc `mat_peptide`.
2. Lấy các field có thể chứa tên gene:

   ```text
   gene, product, note, label, standard_name, locus_tag
   ```

3. Chạy tblastn reference -> query.
4. So tọa độ annotation của query với tọa độ tblastn lift bằng IoU.
5. Nếu IoU đủ cao, tool coi đây là bằng chứng rằng query feature đó tương ứng với canonical từ reference.
6. Chấm điểm từng raw name độc lập.

Ví dụ cùng một feature ORF5 có thể có:

```json
{
  "gene": "GP5",
  "product": "major envelope glycoprotein",
  "note": "ORF5 protein"
}
```

Nếu tọa độ query trùng với tblastn ORF5, tool có thể gợi ý:

| Raw name | Field | Canonical | Action |
|---|---|---|---|
| `GP5` | `gene` | `ORF5` | `save_alias` |
| `ORF5 protein` | `note` | `ORF5` | `save_alias` |
| `major envelope glycoprotein` | `product` | `ORF5` | `ignore` hoặc `manual_review` |

Lý do: `GP5` và `ORF5 protein` có tên gene cụ thể. `major envelope glycoprotein` mô tả protein nhưng không nhất thiết là alias an toàn cho mọi record.

## Granularity mismatch

Một số virus có annotation convention khác nhau giữa reference và query. Alias Manager chỉ chuẩn hóa tên, không tự tách/gộp gene.

Ví dụ với PED:

```text
Reference: ORF1a + ORF1b tách riêng
Query:     ORF1ab là một feature gộp
```

Trong trường hợp này không nên map:

```text
ORF1ab -> ORF1a
ORF1a/1b -> ORF1a
Pol1 -> ORF1a
```

Vì như vậy query đang có vùng gộp nhưng bị gọi nhầm thành vùng lẻ `ORF1a`.

Cách đúng hơn là tạo canonical riêng:

```json
{
  "canonical_names": {
    "ORF1a": ["ORF1A", "ORF1a protein"],
    "ORF1b": ["ORF1B", "ORF1b polyprotein"],
    "ORF1ab": [
      "ORF1",
      "ORF1a/1b",
      "ORF1a/b",
      "ORF 1a/1b",
      "ORF1ab polyprotein",
      "polyprotein 1ab",
      "Pol1",
      "POL1"
    ]
  }
}
```

Khi reference không có `ORF1ab`, output `ORF1ab` có thể được đánh dấu `not_in_reference`. Đây là tín hiệu đúng: query và reference khác mức annotation, không phải alias sai.

## Ví dụ thực tế

Giả sử query có annotation:

```json
{
  "raw_query_names": {
    "gene": "GP5",
    "product": "major envelope glycoprotein",
    "note": "ORF5 protein"
  },
  "query_coords": {
    "start": 13788,
    "end": 14390,
    "strand": "+"
  },
  "best_tblastn_match": {
    "canonical_name": "ORF5",
    "start": 13788,
    "end": 14390,
    "strand": "+",
    "iou": 1.0,
    "coverage": 1.0,
    "identity": 0.94
  }
}
```

Đọc kết quả:

- `IoU = 1.0`: tọa độ query annotation và tblastn prediction trùng hoàn toàn.
- `coverage = 1.0`: tblastn cover đủ protein reference.
- `identity = 0.94`: trình tự protein rất giống.

Kết luận hợp lý:

```text
GP5 -> ORF5: nên save_alias
ORF5 protein -> ORF5: nên save_alias
major envelope glycoprotein: nên cân nhắc ignore/manual_review
```

Ví dụ PED sau khi review 100 records:

```text
envelope protein       -> E
membrane protein       -> M
nucleocapsid protein   -> N
accessory protein 3a   -> ORF3
ORF1a/1b, Pol1, ORF1ab -> ORF1ab
HNZK1                  -> ignore
mp                     -> ambiguous
```

Trong đó `HNZK1` là strain/isolate prefix xuất hiện ở nhiều gene khác nhau, nên không được đưa vào alias.

## Các cảnh báo thường gặp

### `X maps to multiple canonicals`

Nghĩa là cùng một alias đang được map vào nhiều canonical.

Ví dụ:

```text
HNZK1 maps to multiple canonicals: M, N, ORF3, S.
```

Trường hợp này thường là do tên `HNZK1` không phải tên gene, mà là strain/isolate prefix. Nên xóa khỏi alias.

### `X is both ignored and an alias`

Nghĩa là một tên vừa nằm trong `ignored_names`, vừa nằm trong alias của canonical.

Ví dụ:

```text
ORF3 is both ignored and an alias for ORF3.
```

Cách xử lý:

- Nếu `ORF3` là tên gene thật: xóa khỏi ignored.
- Nếu tên đó quá chung chung: xóa khỏi alias.

### Không có suggestion sau khi bấm Generate suggestions

Có thể do:

- Query không có annotation hữu ích.
- Query annotation không có field tên gene/product/note.
- tblastn không lift được feature.
- IoU giữa query annotation và tblastn prediction dưới threshold.
- Feature type trong query không phù hợp.

Xem phần diagnostics trong UI để biết record bị skip ở bước nào.

## Best practices

- Chỉ đưa vào alias những tên đủ cụ thể, ví dụ `GP5`, `ORF5 protein`, `N protein`.
- Không đưa tên quá chung như `protein`, `glycoprotein`, `replicase polyprotein` vào alias nếu nó có thể xuất hiện ở nhiều gene.
- Nếu query dùng gene gộp như `ORF1ab` nhưng reference tách `ORF1a`/`ORF1b`, hãy tạo canonical riêng `ORF1ab` thay vì ép nó vào `ORF1a`.
- Với virus mới, nên review suggestions trước khi lưu config.
- Nếu một alias làm tool map sai, vào Alias Manager xóa alias đó ngay.
- Nếu auto-detect virus chọn sai config, vào tab `Registry` sửa keyword.
- Sau khi sửa alias config, nên chạy lại vài query đại diện để kiểm tra output.
- Không cần sợ sửa nhầm quá mức: mỗi lần save qua UI đều có backup trong `app/config/backups/`.
