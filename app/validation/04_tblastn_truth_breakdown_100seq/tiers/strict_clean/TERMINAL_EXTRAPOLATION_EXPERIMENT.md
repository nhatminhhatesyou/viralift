# FMD Terminal Extrapolation Experiment

Branch: `experiment/fmd-terminal-extrapolation-v2`

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

## Kết Quả Thử Nghiệm Trước Đó

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
- Không thấy exact-match regression trong FMD strict-clean run trước.
- Vì mục tiêu chính là sửa exact boundary, kết quả này tích cực.

## Nhóm Lỗi Được Fix Trong Run Trước

- `n_terminal_truncation_12bp`: `15 -> 0`
- `short_peptide_boundary_offset`: `9 -> 0`
- `minor_c_terminal_truncation_3bp`: `1 -> 0`
- `boundary_offset_matches_neither_ref_nor_truth`: `1 -> 0`
- `minor_vp1_boundary_offset`: `8 -> 7`

Những nhóm không đổi:

- `major_c_terminal_truncation`: `1 -> 1`
- `no_hit`: `1 -> 1`
- `truth_feature_absent`: `3 -> 3`

## Kết Luận

Terminal extrapolation là hướng đáng giữ để thử tiếp.

Nó giảm rõ rệt tool-side boundary failures trong FMD:

- từ `36` xuống `9`
- exact match tăng từ `94.44%` lên `96.70%`

Nhóm được cải thiện đúng với hypothesis:

- HSP local alignment bỏ sót terminal amino acids.
- Tool cũ lấy HSP span làm boundary nên bị cụt đầu/cuối.
- Tool mới dùng query protein coordinates của HSP để extend terminal boundaries.

## Việc Nên Làm Tiếp

- Chạy lại `fmd_terminal_extrapolation_eval.ipynb` trên branch này nếu muốn xác nhận với code/UI mới nhất.
- Review `fmd_fixed_cases.tsv` sau khi chạy lại.
- Review `fmd_regressed_cases.tsv`; nếu file chỉ có header là chưa thấy exact regression.
- Investigate các case còn lại:
  - `3A major_c_terminal_truncation`
  - `2A no_hit`
  - `VP1 residual boundary offsets`
- Giữ extrapolation cho `mat_peptide` trước.
- Sau khi ổn với FMD, mới test PRRS/CDS riêng vì PRRS có overlap và frameshift phức tạp hơn.
