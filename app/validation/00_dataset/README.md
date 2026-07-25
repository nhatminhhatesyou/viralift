# 00 — Dataset Curation

Purpose: document the datasets used for validation before running alias or tblastn accuracy notebooks.

This folder is for accession lists, filtering rules, and strict-clean dataset definitions.

## Why This Exists

Raw GenBank records are not always suitable as ground truth. A record may be annotated, but still fail validation for reasons unrelated to ViraLift:

- missing canonical gene names,
- ORF1a/ORF1b vs ORF1ab granularity differences,
- reference and query using different boundary conventions,
- incomplete or partial genes,
- ambiguous product descriptions,
- absent truth annotation for a feature that the reference can still lift.

For paper-quality validation, keep dataset tiers explicit.

## Recommended Tiers

| Tier | Meaning | Use |
|---|---|---|
| Raw | All downloaded records | Descriptive statistics only |
| Usable-clean | Records with enough annotation to run pipeline checks | Pipeline-level validation |
| Strict-clean | Records with canonical truth for the evaluated reference genes | Accuracy / benchmark claims |

## Suggested Notebook

If this becomes a notebook later, name it:

```text
dataset_manifest_and_strict_clean_filtering.ipynb
```

It should output:

- accession manifest,
- reason for inclusion/exclusion,
- truth-available gene matrix,
- final record count per virus.
