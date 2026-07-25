# 04 — tblastn Truth Validation And Failure Breakdown

Question: **how accurate is ViraLift's tblastn lifting when evaluated against annotated truth records?**

This is the main validation folder for FMDV and PRRSV.

## Paper Use

Use this folder for:

- overall exact and coordinate-correct accuracy,
- per-gene accuracy,
- non-exact failure analysis,
- failure attribution,
- biological interpretation,
- ablation evidence for rescue/extrapolation improvements.

## Recommended Reading Order

1. `tiers/strict_clean/summary_outputs/`  
   Overall FMDV + PRRSV accuracy.

2. `tiers/strict_clean/outputs_fmd/`  
   FMDV per-gene breakdown and terminal-boundary findings.

3. `tiers/strict_clean/outputs_prrsv/`  
   PRRSV per-gene breakdown, ORF1b start-boundary ambiguity, ORF7 rescue issue, ORF2b truth availability.

4. `tiers/strict_clean/terminal_extrapolation_outputs/`  
   FMDV before/after terminal extrapolation.

5. `tiers/strict_clean/prrsv_start_rescue_full_outputs/`  
   PRRSV before/after start-rescue improvement.

## Main Outputs

| File or folder | Use |
|---|---|
| `tiers/strict_clean/summary_outputs/fmd_prrsv_overall_accuracy.tsv` | Overall accuracy table |
| `tiers/strict_clean/summary_outputs/fmd_prrsv_overall_accuracy.png` | Overall accuracy figure |
| `tiers/strict_clean/summary_outputs/fmd_prrsv_per_gene_accuracy.tsv` | Per-gene accuracy table |
| `tiers/strict_clean/outputs_fmd/fmd_failure_final_summary.tsv` | FMDV failure attribution |
| `tiers/strict_clean/outputs_prrsv/prrsv_failure_final_summary.tsv` | PRRSV failure attribution |
| `tiers/strict_clean/prrsv_start_rescue_full_outputs/prrsv_start_rescue_per_gene_comparison.tsv` | PRRSV ablation evidence |

## Interpretation Pattern

For each virus:

1. show overall exact / coordinate-correct accuracy,
2. show per-gene accuracy on truth-available records only,
3. list main non-exact genes,
4. classify each failure as tool-side or reference/query annotation artifact,
5. show whether an implementation change improves the failure class.
