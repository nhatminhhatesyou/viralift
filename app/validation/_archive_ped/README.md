# 06 — PEDV Case Study

Question: **does the ViraLift workflow generalise to a newly configured virus?**

PEDV is used as a newer case study with two reference records and 100 query records.

## Main Notebooks

| Notebook | Purpose |
|---|---|
| `ped_alias_validation.ipynb` | Validate PEDV alias coverage and remaining names to review |
| `ped_tblastn_validation.ipynb` | Baseline PEDV tblastn-vs-truth validation |
| `ped_tblastn_validation_updated.ipynb` | Updated PEDV tblastn validation after alias/tool improvements |

## Main Outputs

| Folder | Use |
|---|---|
| `outputs/alias/` | Alias coverage tables and figures |
| `outputs/tblastn/` | Current tblastn validation outputs |
| `outputs/tblastn_preupdate/` | Baseline before selected updates |
| `outputs/tblastn_postupdate/` | After selected updates |

## Paper Use

Use PEDV as a generalisation case:

1. build or refine alias map,
2. validate alias coverage,
3. validate tblastn lifting against annotated records,
4. compare pre-update and post-update behavior,
5. discuss ORF1a/ORF1b/ORF1ab granularity if relevant.

## Figure/Table Candidates

- `outputs/alias/ped_alias_overall_resolution.png`
- `outputs/tblastn_postupdate/ped_tblastn_accuracy_by_ref_full100.png`
- `outputs/tblastn_postupdate/ped_tblastn_per_gene_accuracy_full100.tsv`
- `outputs/tblastn_postupdate/ped_tblastn_failure_summary_full100.tsv`
