# Hướng dẫn Pipeline Runner trong ViraLift

Tài liệu này giải thích phần chạy pipeline chính của ViraLift: từ lúc user đưa reference + query GenBank vào tool cho tới khi nhận TSV/FASTA output.

Trong UI hiện tại, phần này được gọi là **Run pipeline**. Tên này ổn cho nút/tab thao tác. Trong tài liệu kỹ thuật, có thể gọi rõ hơn là **Pipeline Runner** hoặc **ViraLift Pipeline**.

## Mục lục

- [Pipeline Runner là gì?](#pipeline-runner-là-gì)
- [Input và output](#input-và-output)
- [Flow tổng quát](#flow-tổng-quát)
- [Các stage trong Web UI](#các-stage-trong-web-ui)
- [Tool quyết định direct hay tblastn như thế nào?](#tool-quyết-định-direct-hay-tblastn-như-thế-nào)
- [Alias config ảnh hưởng gì tới pipeline?](#alias-config-ảnh-hưởng-gì-tới-pipeline)
- [Ý nghĩa các threshold](#ý-nghĩa-các-threshold)
- [Cơ chế boundary rescue](#cơ-chế-boundary-rescue)
- [Ý nghĩa các status](#ý-nghĩa-các-status)
- [Cách đọc kết quả](#cách-đọc-kết-quả)
- [Ví dụ chạy bằng CLI](#ví-dụ-chạy-bằng-cli)
- [Các case thường gặp](#các-case-thường-gặp)
- [Best practices](#best-practices)

## Pipeline Runner là gì?

Pipeline Runner là phần xử lý chính của ViraLift. Nó nhận:

```text
1 reference GenBank + nhiều query GenBank records
```

rồi tạo output gồm:

```text
gene name chuẩn hóa + tọa độ trên query + sequence + trạng thái mapping
```

Pipeline tự xử lý hai loại query:

```text
Query đã có annotation hữu ích  -> direct extraction
Query thiếu annotation hữu ích  -> tblastn lifting từ reference
```

## Input và output

### Input

| Input | Ý nghĩa |
|---|---|
| Reference GenBank | Một record chuẩn, có annotation tốt, làm nguồn gene chuẩn |
| Query GenBank | Một hoặc nhiều records cần chuẩn hóa/lift annotation |
| Alias config | File JSON chuẩn hóa tên gene cho virus tương ứng |
| Thresholds | Điều kiện nhận tblastn hit: coverage, identity, e-value, rescue window |

### Output

| Output | Ý nghĩa |
|---|---|
| TSV canonical | Bảng kết quả dùng tên canonical hoặc tên ref tùy setting |
| TSV raw | Bảng kết quả ưu tiên tên gốc nếu query có annotation |
| FASTA | Sequence gene đã extract, có filter theo coverage/status |
| Run summary | Thống kê số feature OK, rescued, no hit, cần review |

## Flow tổng quát

```text
Upload reference + query
        |
        v
Detect virus bằng registry keyword
        |
        v
Load alias config tương ứng
        |
        v
Parse reference features
        |
        v
Chuẩn hóa tên reference bằng alias config
        |
        v
Với từng query record:
    |
    +-- Có annotation hữu ích?
    |       |
    |       +-- Có -> direct extraction
    |       |
    |       +-- Không -> tblastn lifting
    |
    v
Validate coordinate / codon / coverage
        |
        v
Show results + export TSV/FASTA
```

## Các stage trong Web UI

### 1. Upload

User upload:

- Reference GenBank file: một record.
- Query GenBank file: một hoặc nhiều records.

Ở đây cũng có advanced options:

- `min_coverage`
- `min_identity`
- `evalue`
- `rescue_window`
- `Use ref gene names as output`

Nếu tool nhận diện được virus, pipeline đi tiếp. Nếu không nhận diện được virus, UI chuyển sang stage review/tạo alias config mới.

### 2. Virus review

Chỉ hiện khi reference không match keyword nào trong registry.

User chọn:

- Dùng alias config đã có nếu đây là virus cũ nhưng metadata lạ.
- Tạo alias config mới nếu đây là virus mới.

### 3. Alias seed

Chỉ hiện khi tạo virus mới.

Tool lấy tên feature trong reference làm canonical ban đầu, rồi có thể chạy tblastn để gợi ý alias từ query annotation.

### 4. Resolve

Chỉ hiện khi query hoặc reference có tên chưa nằm trong alias config.

User có thể:

- Map tên lạ vào canonical có sẵn.
- Add tên mới làm canonical.
- Save mapping vào alias config để lần sau tự nhận.
- Ignore nếu tên quá chung hoặc không đáng tin.

### 5. Run

Tool chạy thật trên từng query record:

- record có annotation tốt -> direct
- record thiếu annotation -> tblastn

UI hiển thị progress theo record.

### 6. Review

Hiển thị:

- Tổng số records xử lý.
- Tổng số features tìm được.
- Pass rate.
- Số feature cần review.
- Bảng chi tiết từng record.
- Export TSV/FASTA.

## Tool quyết định direct hay tblastn như thế nào?

ViraLift dùng hàm strategy để chọn đường xử lý cho từng query record.

### Direct extraction

Dùng khi query đã có annotation gene-level hữu ích.

Ví dụ query có CDS:

```text
gene = ORF5
product = major envelope glycoprotein
```

Nếu alias config resolve được tên này về canonical, tool không cần alignment. Nó chỉ:

1. Parse tọa độ annotation có sẵn trong query.
2. Chuẩn hóa tên bằng alias config.
3. Extract sequence trực tiếp từ query.

Ưu điểm:

- Nhanh.
- Giữ đúng annotation gốc của query.
- Phù hợp khi query đã được annotate tốt.

### tblastn lifting

Dùng khi query không có annotation hữu ích.

Tool sẽ:

1. Lấy gene từ reference.
2. Dịch gene reference sang protein.
3. Chạy `tblastn` protein reference against query genome.
4. Merge HSP để suy ra tọa độ gene trên query.
5. Extract sequence.
6. Validate start/stop codon nếu phù hợp.

Ưu điểm:

- Dùng được cho query chưa annotate.
- Protein conserved hơn nucleotide, phù hợp virus khác dòng/serotype.
- Xử lý tốt hơn minimap2 với gene ngắn hoặc biến dị.

## Alias config ảnh hưởng gì tới pipeline?

Alias config quyết định tên raw trong GenBank sẽ được chuẩn hóa thành canonical nào.

Ví dụ PED:

```text
envelope protein       -> E
membrane protein       -> M
nucleocapsid protein   -> N
accessory protein 3a   -> ORF3
ORF1a/1b, Pol1, ORF1ab -> ORF1ab
```

Nếu alias config thiếu tên:

- Query vẫn có thể được extract.
- Nhưng tên có thể bị `unresolved_name`.
- UI sẽ hỏi user map tên đó vào canonical nào.

Nếu alias config sai:

- Direct extraction có thể gọi sai gene.
- Validation/pipeline result sẽ bị nhiễu.
- Nên sửa trong Alias Manager.

## Ý nghĩa các threshold

| Threshold | Ý nghĩa | Default |
|---|---|---|
| `min_coverage` | Tỉ lệ protein reference phải được tblastn cover | `0.5` |
| `min_identity` | Protein identity tối thiểu của hit | `0.3` |
| `evalue` | Ngưỡng ý nghĩa thống kê của tblastn hit | `1e-5` |
| `rescue_window` | Vùng bp quanh start để tìm ATG nếu boundary bị lệch | `50` |

Gợi ý:

- Nếu virus rất gần nhau: có thể giữ threshold mặc định.
- Nếu query xa reference: có thể cần giảm `min_identity`.
- Nếu nhiều hit nhiễu: tăng `min_coverage` hoặc kiểm tra reference.

## Cơ chế boundary rescue

Sau khi `tblastn` tìm được vùng gene trên query, ViraLift vẫn phải kiểm tra boundary vì HSP của tblastn là local alignment, không phải lúc nào cũng bắt đúng đầu/cuối CDS.

Với feature type `CDS`, pipeline làm thêm bước codon validation:

```text
sequence phải:
1. bắt đầu bằng ATG
2. kết thúc bằng stop codon: TAA/TAG/TGA
3. có length chia hết cho 3
```

Nếu check này không đạt, tool có thể rescue boundary.

### Start codon rescue

Start rescue chỉ trigger khi sequence hiện tại **không bắt đầu bằng ATG**.

Khi đó tool tìm các ATG quanh `q_start` trong phạm vi `rescue_window`.

Điểm quan trọng: tool không đơn giản chọn ATG gần nhất. Nó chấm điểm candidate theo thứ tự:

```text
1. ưu tiên sequence in-frame
2. ưu tiên length gần reference CDS length nhất
3. nếu vẫn ngang nhau, ưu tiên ATG gần q_start hơn
4. nếu vẫn ngang nhau, ưu tiên upstream trước downstream
```

Reference expected length được tính như sau:

```text
expected_length = protein_length * 3 + 3
```

Vì protein được translate không gồm stop codon, nên cần cộng thêm `+3` cho stop codon.

Ví dụ:

```text
Reference ORF5 length: 603 bp
tblastn lifted span: 600 bp, không bắt đầu bằng ATG
rescue_window: 50 bp
```

Tool sẽ tìm các ATG quanh start. Nếu có một ATG làm cho sequence dài gần `603 bp` và vẫn in-frame, candidate đó sẽ được ưu tiên hơn một ATG gần hơn nhưng làm length lệch nhiều.

Nếu rescue thành công và sequence sau rescue hợp lệ:

```text
status = ok_rescued
rescue_offset = vị trí start mới lệch bao nhiêu bp so với start ban đầu
```

### Stop codon rescue

Stop rescue chạy trước start validation. Nó dùng khi HSP/tblastn span bị thiếu stop codon ở cuối.

Tool scan forward theo frame từ `q_end`:

```text
q_end + 3
q_end + 6
q_end + 9
...
```

tối đa 30 codons. Nếu gặp `TAA`, `TAG`, hoặc `TGA`, tool update `q_end`.

Khác với start rescue, stop rescue hiện tại **không chọn theo length gần ref nhất**. Nó chọn stop codon in-frame gần nhất phía sau.

### Terminal extrapolation

Terminal extrapolation là cơ chế khác với codon rescue.

Nó dùng HSP query coordinates để biết protein alignment thiếu bao nhiêu amino acid ở đầu/cuối, rồi extend boundary:

```text
missing_aa * 3 bp
```

Cơ chế này dùng cho nhánh không validate codon, ví dụ `mat_peptide`, vì mat_peptide không nhất thiết có ATG/stop codon riêng như CDS.

### Vì sao ORF1b vẫn có thể lệch start?

Với các gene như PED/PRRSV `ORF1b`, start boundary có thể liên quan đến vùng `ORF1a/ORF1b` và frameshift/annotation convention. Vì vậy:

```text
end boundary thường ổn định hơn
start boundary dễ lệch hơn
```

Nếu tblastn span đã bắt đầu bằng ATG hợp lệ, start rescue sẽ **không trigger**, kể cả khi GenBank truth dùng một start coordinate khác.

Nếu start rescue có trigger, nó vẫn chọn ATG dựa trên in-frame + length gần ref. Với ORF1b, exact GenBank start đôi khi phản ánh annotation convention hơn là một start codon sinh học rõ ràng. Vì vậy có thể thấy:

```text
coord_correct = True
exact_match = False
delta_end = 0
delta_start lệch vài bp hoặc nhiều bp
```

Pattern này nên được đọc là boundary ambiguity/granularity issue trước khi kết luận tool lift sai.

## Ý nghĩa các status

| Status | Ý nghĩa | Cần review? |
|---|---|---|
| `ok` | Feature được tìm thấy và boundary hợp lệ | Không |
| `ok_rescued` | Feature được tìm thấy, start/stop được rescue | Thường không, nhưng nên chú ý nếu nhiều |
| `direct` | Feature lấy trực tiếp từ query annotation | Không |
| `invalid_boundaries` | Có hit nhưng boundary/codon không hợp lệ | Có |
| `low_coverage` | Hit có coverage thấp hơn threshold | Có |
| `no_hit` | Không tìm thấy hit bằng tblastn | Có |
| `translation_fail` | Reference feature không translate được sang protein | Có |
| `unresolved_name` | Query có tên chưa map được bằng alias config | Có |
| `ambiguous_name` | Tên nằm trong ambiguous list, cần user quyết định | Có |
| `not_in_reference` | Query resolve được tên nhưng reference không có gene đó | Có |

Lưu ý: `not_in_reference` không nhất thiết là tool sai. Ví dụ query có `ORF1ab` nhưng reference chỉ có `ORF1a` và `ORF1b`, đây là mismatch về annotation granularity.

## Cách đọc kết quả

Kết quả chính gồm các cột:

| Cột | Ý nghĩa |
|---|---|
| `query_id` | Record query |
| `name` | Tên gene sau chuẩn hóa |
| `source_name` | Tên raw ban đầu nếu có |
| `ref_start`, `ref_end` | Tọa độ feature trên reference |
| `start`, `end` | Tọa độ feature trên query |
| `strand` | Chiều gene |
| `method` | `direct` hoặc `tblastn` |
| `status` | Trạng thái mapping |
| `coverage` | Tỉ lệ protein reference được cover |
| `identity` | Protein identity của tblastn hit |
| `has_start_codon`, `has_stop_codon` | Check codon boundary |
| `rescue_offset` | Offset nếu start/stop được rescue |

### TSV canonical vs TSV raw

TSV canonical dùng tên chuẩn sau alias mapping.

TSV raw ưu tiên tên gốc từ query nếu có, hữu ích khi muốn audit annotation ban đầu.

### FASTA export

FASTA export cho phép:

- Chọn gene cần xuất.
- Chọn một file chung hoặc một file mỗi gene.
- Lọc theo coverage/identity.
- Include hoặc exclude `ok_rescued`.

Header FASTA:

```text
>{record_id}|{gene}|{start}|{end}|{strand}
```

## Ví dụ chạy bằng CLI

### Auto-detect alias config

```bash
venv/bin/python -m app.src.main \
  --reference app/data/PED/PED_ref_1.gb \
  --query app/data/PED/PED_100seqs.gb \
  --output output/ped_run
```

### Chỉ định alias config thủ công

```bash
venv/bin/python -m app.src.main \
  --reference app/data/PED/PED_ref_1.gb \
  --query app/data/PED/PED_100seqs.gb \
  --output output/ped_run \
  --alias-config app/config/porcine_epidemic_diarrhea_virus_alias.json
```

### Tăng threshold coverage

```bash
venv/bin/python -m app.src.main \
  --reference app/data/PED/PED_ref_1.gb \
  --query app/data/PED/PED_100seqs.gb \
  --output output/ped_strict \
  --min-coverage 0.9
```

## Các case thường gặp

### Query đã có annotation nhưng tên chưa chuẩn

Ví dụ:

```text
gene = GP5
product = major envelope glycoprotein
```

Nếu alias config có `GP5 -> ORF5`, pipeline sẽ đi direct và output `ORF5`.

### Query không có annotation

Pipeline dùng tblastn để lift toàn bộ gene từ reference sang query.

Output name lấy từ reference canonical.

### Query có gene gộp nhưng reference tách gene

Ví dụ:

```text
query: ORF1ab
ref:   ORF1a + ORF1b
```

Không nên ép `ORF1ab` thành `ORF1a`. Nên giữ canonical riêng `ORF1ab`. Nếu reference không có `ORF1ab`, status `not_in_reference` là hợp lý.

### Virus không có trong registry

UI sẽ chuyển sang Virus Review:

- Chọn config có sẵn nếu auto-detect thiếu keyword.
- Hoặc tạo alias config mới bằng Alias Seed.

### Nhiều `invalid_boundaries`

Có thể do:

- Reference boundary khác query truth.
- Start codon khó xác định.
- Gene bị truncated.
- tblastn local hit thiếu terminal amino acids.
- Rescue window chưa đủ hoặc rescue logic chọn boundary chưa tốt.

Nên xem theo từng gene, không chỉ nhìn tổng count.

## Best practices

- Dùng reference có annotation đầy đủ và đáng tin.
- Với virus mới, tạo alias config trước khi chạy production.
- Không map gene gộp vào gene lẻ chỉ để làm score đẹp.
- Khi pass rate thấp, đọc status breakdown trước rồi mới sửa threshold.
- `ok_rescued` thường chấp nhận được, nhưng nếu xuất hiện hàng loạt ở cùng một gene thì nên review reference/query boundary.
- Export FASTA nên lọc status `ok`, `direct`, và tùy mục tiêu có thể include `ok_rescued`.
- Nếu UI auto-detect sai virus, sửa keyword trong Alias Manager thay vì chọn thủ công mỗi lần.
