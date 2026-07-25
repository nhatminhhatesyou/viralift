# 02 — Alias Coverage

Question: **how well do alias maps resolve heterogeneous GenBank gene names into ViraLift canonical names?**

This validation should be run before coordinate accuracy, because name-resolution noise can make correct coordinate predictions appear wrong.

## Paper Use

Use this section to show that ViraLift separates:

- useful aliases that map to canonical genes,
- excluded names that are too generic or unsafe,
- unresolved names that need user/LLM review.

## Main Outputs

| File | Use |
|---|---|
| `outputs/alias_summary_overall.tsv` | Overall canonical/excluded/unresolved coverage |
| `outputs/alias_summary.tsv` | Virus-level summary |
| `outputs/alias_summary_by_field.tsv` | Which qualifier fields are informative |
| `outputs/alias_top_raw_names.tsv` | Most frequent raw names |
| `outputs/alias_noncanonical_names.tsv` | Names worth reviewing |
| `outputs/alias_coverage.png` | Main coverage figure |
| `outputs/alias_coverage_heatmap.png` | Virus/field visual summary |

## Interpretation Pattern

Report:

1. canonical resolution rate,
2. actionable coverage after removing `excluded_names`,
3. top unresolved names,
4. whether remaining unresolved names affect tblastn validation.
