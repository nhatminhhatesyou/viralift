# tblastn Truth Breakdown Tiers

This folder validates tblastn against annotated records using three input tiers.

- `strict_clean`: all truth features resolve to canonical names and canonical names are unique within the record.
- `usable_clean`: all truth features resolve to canonical names, but duplicate canonical names are allowed.
- `raw_realworld`: keeps annotated records as-is, for real-world noise comparison.

Build tier inputs with:

```bash
python app/validation/04_tblastn_truth_breakdown_100seq/scripts/build_validation_tiers.py
```

Each tier folder contains its own notebook, `inputs/`, and `outputs/`.
