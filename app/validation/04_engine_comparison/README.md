# 01 — Engine Benchmark: minimap2 vs tblastn

Question: **which alignment engine should ViraLift use for annotation transfer?**

This folder compares nucleotide-level minimap2 lifting with protein-guided tblastn lifting.

## Paper Use

Use this as the method-justification benchmark. The expected conclusion is that tblastn is more appropriate for viral gene annotation transfer when nucleotide-level similarity or boundary conservation is unreliable.

## Main Outputs

| File | Use |
|---|---|
| `outputs/engine_summary.tsv` | Overall comparison |
| `outputs/engine_summary_wide.tsv` | Paper/table-friendly summary |
| `outputs/engine_per_feature.tsv` | Per-feature comparison |
| `outputs/engine_accuracy.png` | Main figure |
| `outputs/engine_feature_heatmap.png` | Per-gene visual evidence |
| `outputs/engine_failure_modes.tsv` | Failure-mode breakdown |

## Interpretation Pattern

Report:

1. number of evaluated gene-record cases,
2. exact or IoU-based coordinate accuracy,
3. feature-level failure modes,
4. why tblastn was selected for the final pipeline.
