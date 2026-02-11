# tocify — Weekly Journal ToC Digest (RSS → OpenAI → `digest.md`)

This project is a personal fork of [voytek’s original vibe-coded repository](https://github.com/voytek/tocify) (all credit for the core idea and initial implementation goes there).

Essentially, this repo runs a GitHub Action once a week (or on-demand) that:

1. pulls new items from a list of journal RSS feeds  
2. uses OpenAI to triage which items match your research interests  
3. writes a ranked digest to `digest.md` and commits it back to the repo

The goal of this fork is to adapt the tool to my own research and interests by:
- Curating **feeds and sources** that match my work in cancer bioinformatics, multi‑omics, and AI-for-bio
- Tweaking the **interests and tagging logic** so I surface papers that are most relevant to my current projects
- Adjusting the **prompting/setup** so summaries and recommendations are aligned with how I like to read, triage, and revisit material over time

I’m planning to keep using this setup for about a month as an experiment and then reassess, especially with regard to **API usage and costs**. If it proves too expensive or unwieldy, I’ll iterate on the design rather than abandon the idea entirely.

## Ideas to implement

Currently, I have a few directions in mind to make the system more (cost-)efficient:

1. **Caching and memoization**
   - Cache responses for unchanged inputs (same URL / same paper / same query) to avoid repeated API calls
   - Store intermediate representations (embeddings, cleaned text) so re‑runs are cheaper and faster

2. **Incremental updates**
   - Only process **new** or **changed** content since the last run
   - Track per‑source history so I don’t re‑ingest entire feeds each time


Over time I’d like this repo to become a small, maintainable “reading copilot” tuned to my research life: opinionated about what to surface, conservative with API usage, and easy to adapt as my interests evolve!