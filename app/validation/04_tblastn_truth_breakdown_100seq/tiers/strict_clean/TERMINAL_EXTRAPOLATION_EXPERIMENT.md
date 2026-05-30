# FMD Terminal Extrapolation Experiment

Branch: `experiment/fmd-terminal-extrapolation`

Notebook:

- `fmd_terminal_extrapolation_eval.ipynb`

Output folder:

- `terminal_extrapolation_outputs/`

## Mục Tiêu

Thử cải thiện boundary cho FMD `mat_peptide` khi `tblastn` HSP không cover hết đầu/cuối protein.

Ý tưởng:

- HSP cho biết đoạn nào của protein ref được align qua `query_start/query_end`.
- Nếu HSP bắt đầu từ aa 5 thay vì aa 1, tức thiếu 4 aa ở N-terminal.
- Với nucleotide coordinate, thiếu `4 aa * 3 bp = 12 bp`.
- Tool thử extend boundary thêm phần missing terminal đó.

Điểm quan trọng:

- Đây không phải hardcode 12 bp.
- Extension được tính động từ HSP.
- Hiện chỉ áp dụng khi `validate_codons=False`, tức FMD `mat_peptide`.
- CDS/PRRS path chưa bị ảnh hưởng trong thử nghiệm này.

## Kết Quả Accuracy

Baseline:

- Total predictions: `1152`
- Exact match: `1088 / 1152 = 94.44%`
- Coordinate-correct: `1131 / 1152 = 98.18%`
- Failed coord/name: `21 / 1152 = 1.82%`

Sau terminal extrapolation:

- Total predictions: `1152`
- Exact match: `1114 / 1152 = 96.70%`
- Coordinate-correct: `1130 / 1152 = 98.09%`
- Failed coord/name: `22 / 1152 = 1.91%`

Diễn giải:

- Exact match tăng mạnh: `+26 cases`.
- Coordinate-correct giảm rất nhẹ: `-1 case`.
- Không có exact-match regression trong `fmd_regressed_cases.tsv`.
- Vì mục tiêu chính của experiment là sửa exact boundary, kết quả này rất tích cực.

## Kết Quả Theo Final Cause

Tool-side failures:

- Baseline: `36`
- Sau extrapolation: `9`
- Giảm: `27 cases`

Ref/truth artifacts:

- Baseline: `28`
- Sau extrapolation: `29`
- Tăng `1 case`, nghĩa là một case sau khi sửa boundary đã chuyển sang nhóm prediction match ref nhưng query truth khác convention.

## Những Nhóm Lỗi Được Fix

### `n_terminal_truncation_12bp`

- Baseline: `15`
- Sau extrapolation: `0`
- Giảm: `15`

Đây là finding quan trọng nhất.

Các case này chủ yếu thuộc `Lpro` và `3Cpro`:

- Trước: prediction start muộn 12 bp.
- Sau: extrapolation kéo start về đúng truth.
- Status mới: `ok_extrapolated`.

Kết luận:

- Giả thuyết HSP local alignment bị cụt N-terminal là đúng.
- Terminal extrapolation sửa đúng nhóm lỗi này.

### `short_peptide_boundary_offset`

- Baseline: `9`
- Sau extrapolation: `0`
- Giảm: `9`

Nhóm này chủ yếu là `2A`.

Vì `2A` rất ngắn, lệch 3 bp làm IoU bị phạt mạnh. Extrapolation sửa được các boundary offset nhỏ này.

### `minor_c_terminal_truncation_3bp`

- Baseline: `1`
- Sau extrapolation: `0`
- Giảm: `1`

Case này là `3B` thiếu 3 bp ở C-terminal. Extrapolation sửa được.

### `minor_vp1_boundary_offset`

- Baseline: `8`
- Sau extrapolation: `7`
- Giảm: `1`

VP1 chỉ cải thiện nhẹ. Các lỗi VP1 còn lại có thể không đơn giản là missing terminal từ HSP.

### `boundary_offset_matches_neither_ref_nor_truth`

- Baseline: `1`
- Sau extrapolation: `0`
- Giảm: `1`

Một mixed boundary case cũng được sửa.

## Những Nhóm Không Đổi

### `major_c_terminal_truncation`

- Baseline: `1`
- Sau extrapolation: `1`

Đây là case `3A` thiếu C-terminal `111 bp`.

Terminal extrapolation conservative không sửa case này, vì đây là lỗi lớn hơn nhiều so với missing terminal nhỏ.

### `no_hit`

- Baseline: `1`
- Sau extrapolation: `1`

Không có HSP thì không thể extrapolate.

### `truth_feature_absent`

- Baseline: `3`
- Sau extrapolation: `3`

Đây là ground-truth gap, không phải vấn đề tool boundary.

## Fixed Cases

File:

- `terminal_extrapolation_outputs/fmd_fixed_cases.tsv`

Có `26` cases chuyển từ non-exact sang exact.

Các pattern fixed:

- `Lpro`: nhiều case start được kéo upstream 12 bp và match truth.
- `3Cpro`: start được kéo upstream 12 bp và match truth.
- `2A`: boundary offset nhỏ được sửa.
- `3B`: C-terminal thiếu 3 bp được sửa.
- Một `VP1` case được sửa.

Không có regression trong:

- `terminal_extrapolation_outputs/fmd_regressed_cases.tsv`

## Kết Luận

Terminal extrapolation là một cải thiện đáng giữ để thử tiếp.

Nó giảm rõ rệt tool-side boundary failures:

- từ `36` xuống `9`
- exact match tăng từ `94.44%` lên `96.70%`
- không tạo exact-match regression trong FMD strict-clean run này

Nhóm được cải thiện đúng với hypothesis:

- HSP local alignment bỏ sót terminal amino acids.
- Tool cũ lấy HSP span làm boundary nên bị cụt đầu/cuối.
- Tool mới dùng query protein coordinates của HSP để extend terminal boundaries.

## Việc Nên Làm Tiếp

- Review `fmd_fixed_cases.tsv` để xác nhận các fixed cases hợp lý.
- Review `fmd_regressed_cases.tsv`; hiện file chỉ có header, tức chưa thấy exact regression.
- Investigate 9 tool-side cases còn lại:
  - `3A major_c_terminal_truncation`
  - `2A no_hit`
  - `VP1 residual boundary offsets`
- Cân nhắc giữ extrapolation chỉ cho `mat_peptide` trước.
- Sau khi ổn với FMD, mới test PRRS/CDS riêng vì PRRS có overlap và frameshift phức tạp hơn.
