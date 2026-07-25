# 03 — End-To-End Pipeline Validation

Question: **does the complete ViraLift workflow behave correctly on real multi-record inputs?**

This is different from tblastn-only validation. It tests routing, alias resolution, direct extraction, tblastn lifting, status assignment, and export readiness together.

## Paper Use

Use this as a practical workflow validation, not as the main accuracy benchmark.

## Main Outputs

| File | Use |
|---|---|
| `outputs/e2e_summary.tsv` | Overall pipeline summary |
| `outputs/e2e_status_summary.tsv` | Status distribution |
| `outputs/e2e_per_prediction.tsv` | Detailed per-feature records |
| `outputs/e2e_accuracy.png` | Visual summary |
| `outputs/fmd/` | FMDV-specific outputs |
| `outputs/prrs/` | PRRSV-specific outputs |

## Interpretation Pattern

Report:

1. how many records/features were processed,
2. how many were direct vs tblastn,
3. how many needed review,
4. whether exported TSV/FASTA would be usable for downstream analyses.
