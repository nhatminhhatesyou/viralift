# 07 — LLM Alias-Review Validation

Validates the LLM alias-assist feature (`app/src/llm/`) against ground truth,
following the same TSV-based house style as `06_ped_validation/`.

## Why this design

- PED alias-map correctness and tblastn coordinate accuracy are already
  validated separately (see `06_ped_validation/`, ~99% each). This
  validation isolates the **LLM naming/classification review layer** built
  on top of that, so it reuses the already-validated PED data as ground
  truth instead of re-running tblastn from scratch.
- **Track B** (`llm_suggestion_candidates_uncertain.tsv`): real raw
  qualifier values from the 100 PED query records, matched against the
  ref_1 gene model by coordinate overlap (a stand-in for tblastn's IoU
  match — justified since tblastn accuracy is already independently
  validated), scored with the real deterministic classifier
  (`classify_alias_candidate`), then filtered to the subset the real
  `needs_llm_review()` gate would actually send to the LLM.
- **Track A** (`track_a_ground_truth.py`): known excluded/unresolved raw names
  already characterized in `PED_VALIDATION_REPORT.md`, run through
  `review_unresolved_names` (the same function `ui/stages/resolve.py` calls).

## How to run

```bash
# 1. offline, no network needed — (re)builds the candidate dataset
python app/validation/07_llm_alias_validation/build_dataset.py

# 2. sanity-check the script mechanics without spending API calls
python app/validation/07_llm_alias_validation/run_llm_validation.py --mock

# 3. real run — needs OPENAI_API_KEY in .env and network access to
#    api.openai.com (this could NOT be executed from the sandbox that
#    prepared this validation — its outbound network is proxy-allowlisted
#    and api.openai.com is not on that allowlist)
python app/validation/07_llm_alias_validation/run_llm_validation.py
```

Outputs land in `outputs/`: `track_a_results.tsv`, `track_b_results.tsv`,
`summary.tsv`, plus the two candidate-dataset TSVs from step 1.

## Metrics computed

- `action_accuracy`: LLM recommendation vs expected action (`save_alias` /
  `ignore` / `skip`), per track.
- `canonical_accuracy_given_save`: when the LLM says `save_alias`, is the
  canonical name also correct?
- `dangerous_false_positive`: LLM said `save_alias` when it should not have
  — the costliest failure mode, since it would pre-tick "Save" in the UI
  for a name that shouldn't become a permanent alias.
- `operative_bucket`: what the row would actually become in
  `ui/stages/bootstrap_alias.py` (`save` / `exclude` / `skip`), applying the
  exact same confidence gating the UI uses
  (`llm_confidence in {"medium","high"}` for save/exclude; `skip` applies
  regardless of confidence).

## Current Action Semantics

The current alias model uses one exclusion bucket:

```text
excluded_names
```

So the live LLM/UI choices are:

| LLM action | UI meaning | Saved to config? |
|---|---|---|
| `save_alias` | Save raw name as an alias for a canonical gene | Yes, under `canonical_names` |
| `ignore` | Exclude this raw name from future alias matching | Yes, under `excluded_names` |
| `skip` | Do nothing for now | No |

Older cached LLM responses may still contain `move_to_ambiguous`; the current
UI treats that legacy action as `ignore` / exclude.
