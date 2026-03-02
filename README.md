# tocify — Weekly Journal ToC Digest (RSS -> OpenAI -> `digest.md`)
[![Weekly ToC Digest](https://github.com/ellchua/tocify/actions/workflows/weekly-digest.yml/badge.svg)](https://github.com/ellchua/tocify/actions/workflows/weekly-digest.yml)

This project is a personal fork of [voytek's original vibe-coded repository](https://github.com/voytek/tocify) (all credit for the core idea and initial implementation goes there).

This repo runs a GitHub Action once a week (or on-demand) that:

1. Pulls new items from a list of journal RSS feeds.
2. Uses OpenAI to triage which items match your research interests.
3. Writes a ranked digest to `digest.md` and commits digest artifacts back to the repo.

## Current behavior

- Two-stage model triage:
  - First-pass scoring with `gpt-5-nano`.
  - Re-ranking of top candidates with `gpt-5.2`.
- Persistent caching + incremental updates:
  - State saved at `output/tocify_state.json`.
  - Reuses previous model outputs when item content, model, and interests are unchanged.
- Category-based digest layout:
  - Sections are `Neuroblastoma`, `AI`, `Methods`, `Other`.
  - Single-best-match assignment (one paper appears in only one category).
  - Top 5 papers shown per category, sorted by score.
  - Empty categories are still shown.
- Snapshot outputs each run:
  - Root: `digest.md`, `latest.md`
  - Markdown history: `output/digest-YYYY-MM-DD_HHMMSS.md`
  - HTML history + latest: `html_outputs/digest-YYYY-MM-DD_HHMMSS.html`, `html_outputs/latest.html`

## Key environment variables

- Models:
  - `OPENAI_MODEL_CHEAP` (default: `gpt-5-nano`)
  - `OPENAI_MODEL_FINAL` (default: `gpt-5.2`)
- Throughput/quality:
  - `REFINE_TOP_N` (default: `40`)
  - `REFINE_BATCH_SIZE` (default: `20`)
  - `BATCH_SIZE` (default: `50`)
- Rendering/output:
  - `CATEGORY_MAX_ITEMS` (default: `5`)
  - `OUTPUT_DIR` (default: `output`)
  - `HTML_OUTPUT_DIR` (default: `html_outputs`)

## Release notes (2026-03-02)

- Added two-stage OpenAI triage (`nano` -> `gpt-5.2` re-rank).
- Added persistent cache/state to reduce repeated API calls.
- Added incremental processing for unchanged items.
- Added dated markdown and HTML snapshot outputs per run.
- Added categorized digest rendering with fixed section headers.
- Updated workflow env vars and artifact commit paths for new outputs.

## GitHub permissions note

If you push changes to `.github/workflows/weekly-digest.yml`, your GitHub token must include `workflow` scope.
