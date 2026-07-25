# 05 — External Baseline: GATU Export And Compare

Question: **how can representative ViraLift cases be prepared for comparison with GATU or manual annotation review?**

This folder contains exported FASTA inputs and manifests for external comparison.

## Paper Use

Use this as external-baseline support, especially for representative difficult cases rather than a fully automated benchmark.

## Main Outputs

| File | Use |
|---|---|
| `outputs/gatu_manifest.csv` | List of exported cases |
| `outputs/gatu_inputs/` | FASTA files for GATU input |
| `outputs/gatu_output/` | Place for manual/external GATU results |

## Interpretation Pattern

Report:

1. which cases were selected,
2. why they are representative,
3. whether GATU requires manual setup/annotation interpretation,
4. how ViraLift's reference-guided output differs.
