# ViraLift Validation Workspace

This folder is split into small, purpose-specific notebooks. Each notebook has
its own `outputs/` folder so results do not mix across validation tasks.

## Notebooks

1. `01_engine_lifting_minimap_vs_tblastn/`
   - Compares coordinate lifting engines only: `tblastn` vs `minimap2`.
   - Uses fully annotated query records as ground truth.

2. `02_alias_coverage/`
   - Measures how well alias configs cover raw GenBank qualifier names.
   - Reports unresolved, ambiguous, and ignored names.

3. `03_end_to_end_pipeline_100seq/`
   - Runs the production pipeline behavior end-to-end.
   - Uses the real router: direct extraction when annotation is useful, `tblastn`
     when it is not.

4. `04_tblastn_truth_breakdown_100seq/`
   - Forces `tblastn` lifting and breaks down non-matching predictions.
   - Useful for diagnosing coordinate offsets, missing truth genes, and
     annotation disagreements.

5. `05_gatu_export_and_compare/`
   - Exports representative failed cases for GATU and compares GATU output
     against ViraLift predictions.

Shared validation-only helpers live in `_shared/validation_utils.py`. Pipeline
logic should stay in `app/src`; notebooks should call production APIs instead
of re-implementing the tool.

## Setup

Install validation-only notebook dependencies in the active environment:

```bash
pip install -r app/validation/requirements.txt
```
