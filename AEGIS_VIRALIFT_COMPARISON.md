# AEGIS và ViraLift: giải thích và so sánh use case

## Mục lục

- [1. FASTA là gì?](#1-fasta-là-gì)
- [2. Annotation file là gì?](#2-annotation-file-là-gì)
- [3. GFF3/GTF là gì?](#3-gff3gtf-là-gì)
- [4. GenBank khác gì FASTA/GFF3?](#4-genbank-khác-gì-fastagff3)
- [5. AEGIS là gì?](#5-aegis-là-gì)
- [6. AEGIS hoạt động như thế nào?](#6-aegis-hoạt-động-như-thế-nào)
- [7. Use case phù hợp với AEGIS](#7-use-case-phù-hợp-với-aegis)
- [8. ViraLift giải quyết bài toán gì?](#8-viralift-giải-quyết-bài-toán-gì)
- [9. Nếu dùng AEGIS cho bài toán của ViraLift thì sao?](#9-nếu-dùng-aegis-cho-bài-toán-của-viralift-thì-sao)
- [10. Ví dụ: nếu dùng AEGIS để extract ORF5 thì sao?](#10-ví-dụ-nếu-dùng-aegis-để-extract-orf5-thì-sao)
- [11. So sánh nhanh AEGIS và ViraLift](#11-so-sánh-nhanh-aegis-và-viralift)
- [12. Kết luận](#12-kết-luận)
- [13. Mức độ kiểm chứng thông tin](#13-mức-độ-kiểm-chứng-thông-tin)

## 1. FASTA là gì?

FASTA là format đơn giản để lưu **trình tự nucleotide hoặc protein**.

Ví dụ FASTA nucleotide:

```text
>AF184212.1 Porcine reproductive and respiratory syndrome virus
ATGACGTAGCTAGCTAGCTAGCTAGCTAGCTAA
```

FASTA chỉ cho biết:

- ID hoặc tên sequence nằm ở dòng bắt đầu bằng `>`
- Trình tự sequence nằm ở các dòng bên dưới

FASTA **không tự nói rõ đoạn nào là gene, CDS, ORF, protein, exon...**. Muốn biết đoạn nào là gene thì cần thêm annotation.

## 2. Annotation file là gì?

Annotation file là file ghi chú các feature trên genome, ví dụ:

- gene nào nằm từ vị trí nào đến vị trí nào
- feature đó nằm trên strand `+` hay `-`
- feature đó là `gene`, `CDS`, `exon`, `mat_peptide`, `UTR`...
- tên gene/product/note là gì

Ví dụ ý nghĩa:

```text
ORF5 nằm từ 13788 đến 14390 trên strand +
product là major envelope glycoprotein
```

Trong bioinformatics, sequence và annotation thường đi thành cặp:

- FASTA: chứa trình tự
- GFF3/GTF/GenBank: chứa thông tin feature/tọa độ/tên

## 3. GFF3/GTF là gì?

`GFF3` và `GTF` là hai format phổ biến để lưu genome annotation dạng bảng.

Một dòng GFF3 thường có 9 cột:

```text
seqid  source  type  start  end  score  strand  phase  attributes
```

Ví dụ:

```text
chr1  GenBank  gene  1000  2000  .  +  .  ID=gene1;Name=ORF5
chr1  GenBank  CDS   1200  1800  .  +  0  Parent=gene1;product=major envelope glycoprotein
```

Nghĩa là:

- Trên `chr1`
- Có một `gene` từ `1000` đến `2000`
- Có một `CDS` từ `1200` đến `1800`
- Gene tên `ORF5`
- Product là `major envelope glycoprotein`

GFF3/GTF rất phổ biến trong genome lớn như thực vật, động vật, người, nấm... vì các genome này thường có cấu trúc gene phức tạp:

```text
gene
  transcript
    exon
    CDS
    UTR
```

Với virus, cấu trúc thường đơn giản hơn. Nhiều viral records chỉ có `CDS`, `gene`, `mat_peptide`, hoặc thậm chí thiếu annotation hữu ích.

## 4. GenBank khác gì FASTA/GFF3?

GenBank là format giàu thông tin hơn FASTA. Một file GenBank thường chứa cả:

- metadata: organism, accession, source
- sequence
- feature annotation
- qualifier như `/gene`, `/product`, `/note`, `/translation`

Ví dụ đơn giản:

```text
FEATURES             Location/Qualifiers
     CDS             13788..14390
                     /gene="ORF5"
                     /product="major envelope glycoprotein"
                     /translation="..."
ORIGIN
        1 atgacg...
```

Trong ViraLift, dữ liệu viral thường là GenBank. Vấn đề là GenBank viral records không đồng nhất:

- record này dùng `/gene="ORF5"`
- record khác dùng `/gene="GP5"`
- record khác chỉ có `/product="major envelope glycoprotein"`
- record khác thiếu annotation hoàn toàn

Đây là lý do cần alias map và annotation lifting.

## 5. AEGIS là gì?

AEGIS là viết tắt của **Annotation Extraction and Genomic Integration Suite**.

Paper: **AEGIS: an annotation extraction and genomic integration resource**, Bioinformatics, 2026.

Nguồn:

- Paper: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag363/8704544
- GitHub: https://github.com/Tomsbiolab/aegis

AEGIS là một Python toolkit để làm việc với genome annotation, đặc biệt là annotation dạng `GFF3/GTF`.

Theo abstract của paper và README trên GitHub, mục tiêu/chức năng chính của AEGIS gồm:

- parse annotation files
- clean/tidy annotation
- standardize format
- extract feature sequence
- merge/reformat annotation
- compare annotation giữa các genome hoặc annotation versions
- hỗ trợ comparative genomics như orthology/synteny

## 6. AEGIS hoạt động như thế nào?

Nói đơn giản, flow của AEGIS là:

1. Người dùng đưa vào genome annotation file, thường là `GFF3/GTF`, và genome FASTA nếu cần.
2. AEGIS parse annotation thành object có cấu trúc, ví dụ:

```text
Genome
  Scaffold / Chromosome
    Gene
      Transcript
        Exon
        CDS
        UTR
```

3. AEGIS kiểm tra và chuẩn hóa annotation:

- sửa format không chuẩn
- xử lý ID/Parent relationship
- tidy annotation để downstream tools đọc được
- rename/reformat nếu cần

4. AEGIS có thể extract sequence:

- gene sequence
- CDS
- protein
- promoter
- exon/intron

5. AEGIS có thể so sánh annotation:

- so overlap giữa gene models
- so annotation versions
- dùng sequence homology, synteny, coordinate lift-over để infer orthology

Ví dụ minh họa ý tưởng chuẩn hóa annotation file:

Input GFF3 lộn xộn:

```text
chr1  GenBank  CDS   100  900  .  +  0  ID=cds1;Parent=geneA;product=ORF5 protein
chr1  GenBank  gene  100  900  .  +  .  ID=geneA;Name=GP5
chr1  GenBank  exon  100  900  .  +  .  Parent=cds1
```

Vấn đề:

- `CDS` xuất hiện trước `gene`
- quan hệ `Parent` chưa rõ
- annotation hierarchy có thể gây lỗi cho downstream tools

Về mặt ý tưởng, một annotation sạch hơn thường cần được tổ chức theo quan hệ rõ ràng hơn, ví dụ:

```text
chr1  AEGIS  gene        100  900  .  +  .  ID=geneA;Name=GP5
chr1  AEGIS  transcript  100  900  .  +  .  ID=geneA.t1;Parent=geneA
chr1  AEGIS  exon        100  900  .  +  .  ID=geneA.t1.exon1;Parent=geneA.t1
chr1  AEGIS  CDS         100  900  .  +  0  ID=geneA.t1.cds1;Parent=geneA.t1;product=ORF5 protein
```

Lưu ý quan trọng: ví dụ trên là **minh họa khái niệm**, không phải output được copy từ AEGIS. Output thực tế phụ thuộc command, input file, và rule xử lý của AEGIS.

## 7. Use case phù hợp với AEGIS

AEGIS phù hợp khi bài toán là:

- Có sẵn genome annotation ở dạng `GFF3/GTF`.
- Muốn clean hoặc validate annotation file.
- Muốn extract CDS/protein/promoter từ annotation.
- Muốn reformat/rename/merge annotation files.
- Muốn so sánh annotation giữa nhiều genome hoặc nhiều phiên bản annotation.
- Muốn phân tích gene model, transcript, exon, CDS, UTR.
- Muốn làm comparative genomics, orthology, synteny.

Ví dụ use case:

### Use case 1: Clean GFF3 trước khi dùng downstream tools

Input:

```text
genome.fasta
annotation.gff3
```

AEGIS tidy annotation để ID/Parent/format chuẩn hơn, giúp các tool khác đọc được.

### Use case 2: Extract protein từ annotation

Input:

```text
genome.fasta
annotation.gff3
```

AEGIS extract CDS/protein sequences từ các feature đã được annotation.

### Use case 3: So sánh annotation version

Input:

```text
Arabidopsis_v1.gff3
Arabidopsis_v2.gff3
genome.fasta
```

AEGIS so sánh gene models để xem gene nào thay đổi tọa độ, exon, transcript, hoặc cấu trúc.

## 8. ViraLift giải quyết bài toán gì?

ViraLift được thiết kế cho bài toán khác, cụ thể hơn:

> User có một reference virus chuẩn và nhiều query virus records. Query có thể có annotation lộn xộn hoặc thiếu annotation. Cần chuẩn hóa tên gene và nếu thiếu annotation thì lift annotation từ reference sang query.

Flow chính của ViraLift:

1. User upload reference GenBank và query GenBank/FASTA.
2. Tool detect virus và load alias config tương ứng.
3. Nếu query đã có annotation hữu ích:
   - extract feature trực tiếp
   - chuẩn hóa tên gene qua alias map
4. Nếu query thiếu annotation hoặc annotation không đủ:
   - lấy protein từ reference
   - chạy `tblastn` against query genome
   - suy ra tọa độ gene trên query
   - validate boundary/codon nếu phù hợp
5. Export kết quả chuẩn hóa.

ViraLift tập trung vào những vấn đề rất đặc thù của viral GenBank:

- tên gene không đồng nhất
- cùng một gene có nhiều tên khác nhau
- product/note có thể chứa thông tin hữu ích hoặc chỉ là mô tả chung
- một số virus dùng `mat_peptide`, một số dùng `CDS/ORF`
- query có thể thiếu annotation
- cần transfer annotation từ reference chuẩn

Ví dụ alias normalization trong ViraLift:

```text
GP5
ORF5 protein
major envelope glycoprotein
```

Tool cần quyết định:

```text
GP5 -> ORF5
ORF5 protein -> ORF5
major envelope glycoprotein -> có thể là alias hoặc ignored tùy virus/context
```

Đây không phải trọng tâm chính của AEGIS.

## 9. Nếu dùng AEGIS cho bài toán của ViraLift thì sao?

Giả sử có 100 GenBank viral records, như PRRSV/FMD/PED trong validation.

Nếu dùng AEGIS, flow có thể là:

1. Convert GenBank sang FASTA + GFF3.
2. Chạy AEGIS để tidy/validate GFF3.
3. Extract CDS/protein nếu annotation có sẵn.
4. Có thể so overlap hoặc compare nếu dữ liệu phù hợp.

Nhưng sẽ gặp các giới hạn:

- AEGIS không trực tiếp nhận GenBank viral workflow như ViraLift.
- AEGIS cần annotation có sẵn để clean/extract.
- Nếu query thiếu annotation, AEGIS không tự làm viral ref-guided tblastn annotation transfer theo flow của ViraLift.
- AEGIS không có alias map virus-specific để map `GP5`, `ORF5 protein`, `major envelope glycoprotein` về canonical `ORF5`.
- AEGIS không có UI để user approve alias suggestions.
- AEGIS không tập trung vào logic riêng cho `CDS`, `mat_peptide`, overlapping viral genes, ORF1a/ORF1b/ORF1ab convention.

Vì vậy AEGIS có thể hữu ích nếu mục tiêu là **clean GFF3 annotation format**, nhưng không thay thế trực tiếp ViraLift cho bài toán lab hiện tại.

## 10. Ví dụ: nếu dùng AEGIS để extract ORF5 thì sao?

Nếu dùng AEGIS để extract `ORF5`, nó sẽ hoạt động tốt khi annotation đã ghi rõ feature đó là `ORF5`.

Ví dụ GFF3 đã có tên rõ:

```text
chr1  GenBank  CDS  13788  14390  .  +  0  ID=cds5;Name=ORF5;product=major envelope glycoprotein
```

Trong trường hợp này, AEGIS có thể extract đoạn CDS/protein tại tọa độ:

```text
13788..14390
```

Tuy nhiên, trong viral GenBank thực tế, cùng một feature có thể được ghi bằng nhiều cách khác nhau:

```text
/gene="GP5"
/product="major envelope glycoprotein"
/note="ORF5 protein"
```

AEGIS có thể extract feature này nếu nó đã tồn tại trong annotation, nhưng nó không nhất thiết biết rằng:

```text
GP5 -> ORF5
ORF5 protein -> ORF5
major envelope glycoprotein -> alias hay ignored?
```

Đây là phần ViraLift xử lý bằng alias map theo từng virus.

Vì vậy:

- Nếu annotation ghi rõ `ORF5`, AEGIS có thể extract `ORF5`.
- Nếu annotation ghi `GP5`, `glycoprotein 5`, hoặc `major envelope glycoprotein`, AEGIS có thể extract feature, nhưng không chắc chuẩn hóa thành canonical `ORF5`.
- Nếu record không có annotation cho ORF5/CDS tương ứng, AEGIS không tự dùng reference để tìm ORF5 bằng `tblastn`.

ViraLift khác ở chỗ nó có thể:

```text
GP5 / ORF5 protein / glycoprotein 5
```

map về:

```text
ORF5
```

và nếu query thiếu annotation, ViraLift dùng protein `ORF5` từ reference chạy `tblastn` để lift tọa độ `ORF5` sang query.

## 11. So sánh nhanh AEGIS và ViraLift

| Tiêu chí | AEGIS | ViraLift |
|---|---|---|
| Mục tiêu chính | Manipulate, clean, extract, integrate genome annotations | Chuẩn hóa tên gene virus và transfer annotation từ reference sang query |
| Domain | General genome annotation, comparative genomics | Viral genomes |
| Input chính | GFF3/GTF + FASTA | GenBank/FASTA reference + query |
| Cần annotation có sẵn? | Có, chủ yếu xử lý annotation đã có | Không nhất thiết; nếu query thiếu annotation thì dùng tblastn lifting |
| Chuẩn hóa cái gì? | Chuẩn hóa cấu trúc/format annotation file | Chuẩn hóa tên gene về canonical names |
| Ví dụ chuẩn hóa | Fix ID/Parent, tidy GFF3, reformat GTF/GFF3 | `GP5`, `ORF5 protein` -> `ORF5` |
| Feature model | Gene/transcript/exon/CDS/UTR | CDS, mat_peptide, ORF, viral gene segments |
| Annotation transfer | Có comparative/liftover/orthology modules tổng quát | tblastn protein-based annotation lifting từ reference virus |
| Alias manager | Không phải trọng tâm | Có alias map, ignored names, ambiguous names, user review |
| Phù hợp nhất khi | Có annotation GFF3/GTF cần clean/compare | Có viral ref/query cần chuẩn hóa tên và annotate query thiếu feature |

## 12. Kết luận

AEGIS và ViraLift giống nhau ở tầng ý tưởng lớn: cả hai đều làm việc với genome annotation và đều muốn annotation trở nên dễ dùng hơn cho downstream analysis.

Nhưng hai tool giải quyết hai bài toán khác nhau:

- **AEGIS** là toolkit tổng quát để clean, extract, merge, reformat và compare genome annotation files như GFF3/GTF.
- **ViraLift** là workflow chuyên biệt cho virus, tập trung vào canonical gene-name standardization và reference-guided tblastn annotation transfer khi query records thiếu hoặc có annotation không đồng nhất.

Nói ngắn gọn:

> AEGIS chuẩn hóa và tích hợp annotation file/model. ViraLift chuẩn hóa tên gene virus và lift annotation từ reference chuẩn sang query.

Với bài toán hiện tại của lab, ViraLift phù hợp hơn vì input là viral GenBank/FASTA records, annotation thường không đồng nhất, và mục tiêu chính là chuẩn hóa tên gene cộng với annotate query dựa trên reference.

## 13. Mức độ kiểm chứng thông tin

Các thông tin sau đã được kiểm chứng từ paper/abstract Bioinformatics và GitHub README của AEGIS:

- AEGIS là **Annotation Extraction and Genomic Integration Suite**.
- AEGIS là Python toolkit cho manipulation, analysis, integration of genomic annotations.
- AEGIS làm việc chủ yếu với genome annotation như `GFF3/GTF`.
- AEGIS có chức năng parse, correct, standardise, validate genome annotations.
- AEGIS hỗ trợ feature extraction, ví dụ coding sequences và promoters.
- AEGIS có các command như extract, tidy, overlap, rename, summary, merge, subset, reformat, list.
- AEGIS có comparative/integrative module cho orthology, synteny, sequence homology, coordinate-based lift-over.
- AEGIS được demo trên Arabidopsis annotation versions và plant genomes.

Các phần sau là suy luận/so sánh dựa trên phạm vi công bố của AEGIS và thiết kế hiện tại của ViraLift:

- AEGIS không phải tool chuyên cho viral GenBank alias normalization.
- AEGIS không thay trực tiếp ViraLift cho bài toán map `GP5`, `ORF5 protein`, `major envelope glycoprotein` về canonical `ORF5`.
- AEGIS không được mô tả là workflow chuyên dụng cho reference-guided viral `tblastn` annotation transfer.
- Ví dụ GFF3 trong tài liệu này là ví dụ minh họa, không phải output thật của AEGIS.

Vì vậy, phần mô tả chức năng AEGIS là dựa trên nguồn thật; phần so sánh với ViraLift là phân tích use case, không phải claim từ tác giả AEGIS.
