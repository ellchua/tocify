import html
import json
import math
import os
import re
import time
import hashlib
from datetime import datetime, timezone, timedelta

import feedparser
import httpx
from dateutil import parser as dtparser
from openai import OpenAI, APITimeoutError, APIConnectionError, RateLimitError


# ---- config (env-tweakable) ----
CHEAP_MODEL = os.getenv("OPENAI_MODEL_CHEAP", "gpt-5-nano")
FINAL_MODEL = os.getenv("OPENAI_MODEL_FINAL", "gpt-5.2")
MAX_ITEMS_PER_FEED = int(os.getenv("MAX_ITEMS_PER_FEED", "50"))
MAX_TOTAL_ITEMS = int(os.getenv("MAX_TOTAL_ITEMS", "250"))
LOOKBACK_DAYS = int(os.getenv("LOOKBACK_DAYS", "7"))
INTERESTS_MAX_CHARS = int(os.getenv("INTERESTS_MAX_CHARS", "3000"))
SUMMARY_MAX_CHARS = int(os.getenv("SUMMARY_MAX_CHARS", "300"))
PREFILTER_KEEP_TOP = int(os.getenv("PREFILTER_KEEP_TOP", "100"))
BATCH_SIZE = int(os.getenv("BATCH_SIZE", "50"))
MIN_SCORE_READ = float(os.getenv("MIN_SCORE_READ", "0.7"))
MAX_RETURNED = int(os.getenv("MAX_RETURNED", "40"))
REFINE_TOP_N = int(os.getenv("REFINE_TOP_N", "40"))
REFINE_BATCH_SIZE = int(os.getenv("REFINE_BATCH_SIZE", "20"))
CATEGORY_MAX_ITEMS = int(os.getenv("CATEGORY_MAX_ITEMS", "5"))

OUTPUT_DIR = os.getenv("OUTPUT_DIR", "output")
HTML_OUTPUT_DIR = os.getenv("HTML_OUTPUT_DIR", "html_outputs")
STATE_PATH = os.getenv("STATE_PATH", os.path.join(OUTPUT_DIR, "tocify_state.json"))
STATE_RETENTION_DAYS = int(os.getenv("STATE_RETENTION_DAYS", "120"))

SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "week_of": {"type": "string"},
        "notes": {"type": "string"},
        "ranked": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "id": {"type": "string"},
                    "title": {"type": "string"},
                    "link": {"type": "string"},
                    "source": {"type": "string"},
                    "published_utc": {"type": ["string", "null"]},
                    "score": {"type": "number"},
                    "why": {"type": "string"},
                    "tags": {"type": "array", "items": {"type": "string"}},
                },
                "required": ["id", "title", "link", "source", "published_utc", "score", "why", "tags"],
            },
        },
    },
    "required": ["week_of", "notes", "ranked"],
}


# ---- tiny helpers ----
def load_feeds(path: str) -> list[dict]:
    """
    Supports:
    - blank lines
    - comments starting with #
    - optional naming via: Name | URL

    Returns list of:
    { "name": "...", "url": "..." }
    """
    feeds = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue

            # Named feed: "Name | URL"
            if "|" in s:
                name, url = [x.strip() for x in s.split("|", 1)]
            else:
                name, url = None, s

            feeds.append({"name": name, "url": url})

    return feeds


def read_text(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def load_prompt_template(path: str = "prompt.txt") -> str:
    if not os.path.exists(path):
        raise RuntimeError("prompt.txt not found in repo root")
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def sha1(s: str) -> str:
    return hashlib.sha1(s.encode("utf-8")).hexdigest()


def section(md: str, heading: str) -> str:
    m = re.search(rf"(?im)^\s*#{1,6}\s+{re.escape(heading)}\s*$", md)
    if not m:
        return ""
    rest = md[m.end() :]
    m2 = re.search(r"(?im)^\s*#{1,6}\s+\S", rest)
    return (rest[: m2.start()] if m2 else rest).strip()


def parse_interests_md(md: str) -> dict:
    keywords = []
    for line in section(md, "Keywords").splitlines():
        line = re.sub(r"^[\-\*\+]\s+", "", line.strip())
        if line:
            keywords.append(line)
    narrative = section(md, "Narrative").strip()
    if len(narrative) > INTERESTS_MAX_CHARS:
        narrative = narrative[:INTERESTS_MAX_CHARS] + "…"
    return {"keywords": keywords[:200], "narrative": narrative}


def item_fingerprint(it: dict) -> str:
    payload = {
        "source": it.get("source"),
        "title": it.get("title"),
        "link": it.get("link"),
        "published_utc": it.get("published_utc"),
        "summary": it.get("summary"),
    }
    return sha1(json.dumps(payload, sort_keys=True, ensure_ascii=False))


def upsert_best_ranked(rows: list[dict]) -> dict[str, dict]:
    best: dict[str, dict] = {}
    for r in rows:
        rid = r["id"]
        if rid not in best or r["score"] > best[rid]["score"]:
            best[rid] = r
    return best


# ---- state/cache ----
def load_state(path: str) -> dict:
    if not os.path.exists(path):
        return {"version": 1, "items": {}, "updated_utc": None}
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, dict):
            raise ValueError("invalid state format")
        data.setdefault("version", 1)
        data.setdefault("items", {})
        data.setdefault("updated_utc", None)
        return data
    except Exception:
        return {"version": 1, "items": {}, "updated_utc": None}


def save_state(path: str, state: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    state["updated_utc"] = datetime.now(timezone.utc).isoformat()
    with open(path, "w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=2, sort_keys=True)


def prune_state(state: dict, retention_days: int) -> None:
    cutoff = datetime.now(timezone.utc) - timedelta(days=retention_days)
    keep: dict[str, dict] = {}
    for item_id, entry in state.get("items", {}).items():
        last_seen = entry.get("last_seen_utc")
        if not last_seen:
            keep[item_id] = entry
            continue
        try:
            dt = dtparser.parse(last_seen)
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
            if dt >= cutoff:
                keep[item_id] = entry
        except Exception:
            keep[item_id] = entry
    state["items"] = keep


def read_cached_results(
    items: list[dict],
    state: dict,
    stage: str,
    model: str,
    interests_hash: str,
) -> tuple[list[dict], list[dict]]:
    cached, pending = [], []
    now_utc = datetime.now(timezone.utc).isoformat()

    for it in items:
        iid = it["id"]
        fp = item_fingerprint(it)
        entry = state.setdefault("items", {}).setdefault(iid, {})
        cached_fingerprint = entry.get("fingerprint")
        entry["last_seen_utc"] = now_utc

        cached_result = entry.get(f"{stage}_result")
        cached_model = entry.get(f"{stage}_model")
        cached_interests = entry.get(f"{stage}_interests_hash")
        if (
            cached_result
            and cached_model == model
            and cached_interests == interests_hash
            and cached_fingerprint == fp
        ):
            cached.append(cached_result)
        else:
            pending.append(it)
            entry["fingerprint"] = fp
    return cached, pending


def write_cached_results(
    state: dict,
    stage: str,
    model: str,
    interests_hash: str,
    items_by_id: dict[str, dict],
    ranked_rows: list[dict],
) -> None:
    now_utc = datetime.now(timezone.utc).isoformat()
    for r in ranked_rows:
        iid = r["id"]
        it = items_by_id.get(iid)
        if not it:
            continue
        entry = state.setdefault("items", {}).setdefault(iid, {})
        entry["last_seen_utc"] = now_utc
        entry["fingerprint"] = item_fingerprint(it)
        entry[f"{stage}_model"] = model
        entry[f"{stage}_interests_hash"] = interests_hash
        entry[f"{stage}_result"] = r


# ---- rss ----
def parse_date(entry) -> datetime | None:
    for attr in ("published_parsed", "updated_parsed"):
        t = getattr(entry, attr, None)
        if t:
            return datetime(*t[:6], tzinfo=timezone.utc)
    for key in ("published", "updated", "created"):
        val = entry.get(key)
        if val:
            try:
                dt = dtparser.parse(val)
                return dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
            except Exception:
                pass
    return None


def fetch_rss_items(feeds: list[dict]) -> list[dict]:
    cutoff = datetime.now(timezone.utc) - timedelta(days=LOOKBACK_DAYS)
    items = []
    for feed in feeds:
        url = feed["url"]
        d = feedparser.parse(url)

        # Priority: manual name > RSS title > URL
        source = (feed.get("name") or d.feed.get("title") or url).strip()
        for e in d.entries[:MAX_ITEMS_PER_FEED]:
            title = (e.get("title") or "").strip()
            link = (e.get("link") or "").strip()
            if not (title and link):
                continue
            dt = parse_date(e)
            if dt and dt < cutoff:
                continue
            summary = re.sub(r"\s+", " ", (e.get("summary") or e.get("description") or "").strip())
            if len(summary) > SUMMARY_MAX_CHARS:
                summary = summary[:SUMMARY_MAX_CHARS] + "…"
            items.append(
                {
                    "id": sha1(f"{source}|{title}|{link}"),
                    "source": source,
                    "title": title,
                    "link": link,
                    "published_utc": dt.isoformat() if dt else None,
                    "summary": summary,
                }
            )
    # dedupe + newest first
    items = list({it["id"]: it for it in items}.values())
    items.sort(key=lambda x: x["published_utc"] or "", reverse=True)
    return items[:MAX_TOTAL_ITEMS]


# ---- local prefilter ----
def keyword_prefilter(items: list[dict], keywords: list[str], keep_top: int) -> list[dict]:
    kws = [k.lower() for k in keywords if k.strip()]

    def hits(it):
        text = (it.get("title", "") + " " + it.get("summary", "")).lower()
        return sum(1 for k in kws if k in text)

    scored = [(hits(it), it) for it in items]
    matched = [it for s, it in scored if s > 0]
    if len(matched) < min(50, keep_top):
        return items[:keep_top]
    matched.sort(key=hits, reverse=True)
    return matched[:keep_top]


# ---- openai ----
def make_openai_client() -> OpenAI:
    key = os.environ.get("OPENAI_API_KEY", "").strip()
    if not key.startswith("sk-"):
        raise RuntimeError("OPENAI_API_KEY missing/invalid (expected to start with 'sk-').")
    http_client = httpx.Client(
        timeout=httpx.Timeout(connect=30.0, read=300.0, write=30.0, pool=30.0),
        http2=False,
        trust_env=False,
        headers={"Connection": "close", "Accept-Encoding": "gzip"},
    )
    return OpenAI(api_key=key, http_client=http_client)


def call_openai_triage(client: OpenAI, interests: dict, items: list[dict], model: str) -> dict:
    lean_items = [
        {
            "id": it["id"],
            "source": it["source"],
            "title": it["title"],
            "link": it["link"],
            "published_utc": it.get("published_utc"),
            "summary": (it.get("summary") or "")[:SUMMARY_MAX_CHARS],
        }
        for it in items
    ]

    template = load_prompt_template()

    prompt = (
        template.replace("{{KEYWORDS}}", json.dumps(interests["keywords"], ensure_ascii=False))
        .replace("{{NARRATIVE}}", interests["narrative"])
        .replace("{{ITEMS}}", json.dumps(lean_items, ensure_ascii=False))
    )

    last = None
    for attempt in range(6):
        try:
            resp = client.responses.create(
                model=model,
                input=prompt,
                text={"format": {"type": "json_schema", "name": "weekly_toc_digest", "schema": SCHEMA, "strict": True}},
            )
            return json.loads(resp.output_text)
        except (APITimeoutError, APIConnectionError, RateLimitError) as e:
            last = e
            time.sleep(min(60, 2**attempt))
    raise last


def triage_in_batches(client: OpenAI, interests: dict, items: list[dict], batch_size: int, model: str) -> dict:
    week_of = datetime.now(timezone.utc).date().isoformat()
    if not items:
        return {"week_of": week_of, "notes": "", "ranked": []}

    total = math.ceil(len(items) / batch_size)
    all_ranked, notes_parts = [], []

    for i in range(0, len(items), batch_size):
        batch = items[i : i + batch_size]
        print(f"Triage with {model}, batch {i // batch_size + 1}/{total} ({len(batch)} items)")
        res = call_openai_triage(client, interests, batch, model=model)
        if res.get("notes", "").strip():
            notes_parts.append(res["notes"].strip())
        all_ranked.extend(res.get("ranked", []))

    ranked = sorted(upsert_best_ranked(all_ranked).values(), key=lambda x: x["score"], reverse=True)
    return {"week_of": week_of, "notes": " ".join(dict.fromkeys(notes_parts))[:1000], "ranked": ranked}


def two_stage_triage(
    client: OpenAI,
    interests: dict,
    items: list[dict],
    state: dict,
) -> dict:
    week_of = datetime.now(timezone.utc).date().isoformat()
    interests_hash = sha1(json.dumps(interests, sort_keys=True, ensure_ascii=False))
    items_by_id = {it["id"]: it for it in items}

    cheap_cached, cheap_pending = read_cached_results(items, state, stage="cheap", model=CHEAP_MODEL, interests_hash=interests_hash)
    print(f"Cache hit (cheap model): {len(cheap_cached)} / {len(items)}")

    cheap_new = []
    cheap_notes = []
    if cheap_pending:
        cheap_res = triage_in_batches(client, interests, cheap_pending, batch_size=BATCH_SIZE, model=CHEAP_MODEL)
        cheap_new = cheap_res.get("ranked", [])
        if cheap_res.get("notes", "").strip():
            cheap_notes.append(cheap_res["notes"].strip())
        write_cached_results(
            state,
            stage="cheap",
            model=CHEAP_MODEL,
            interests_hash=interests_hash,
            items_by_id=items_by_id,
            ranked_rows=cheap_new,
        )

    cheap_map = upsert_best_ranked([*cheap_cached, *cheap_new])
    cheap_ranked = sorted(cheap_map.values(), key=lambda x: x["score"], reverse=True)

    refine_ids = [r["id"] for r in cheap_ranked[:REFINE_TOP_N]]
    refine_items = [items_by_id[rid] for rid in refine_ids if rid in items_by_id]

    final_cached, final_pending = read_cached_results(refine_items, state, stage="final", model=FINAL_MODEL, interests_hash=interests_hash)
    print(f"Cache hit (final model): {len(final_cached)} / {len(refine_items)}")

    final_new = []
    final_notes = []
    if final_pending:
        final_res = triage_in_batches(
            client,
            interests,
            final_pending,
            batch_size=REFINE_BATCH_SIZE,
            model=FINAL_MODEL,
        )
        final_new = final_res.get("ranked", [])
        if final_res.get("notes", "").strip():
            final_notes.append(final_res["notes"].strip())
        write_cached_results(
            state,
            stage="final",
            model=FINAL_MODEL,
            interests_hash=interests_hash,
            items_by_id=items_by_id,
            ranked_rows=final_new,
        )

    final_map = upsert_best_ranked([*final_cached, *final_new])

    merged = []
    for rid, cheap_row in cheap_map.items():
        merged.append(final_map.get(rid, cheap_row))

    notes_parts = [*cheap_notes, *final_notes]
    notes = " ".join(dict.fromkeys([x for x in notes_parts if x])).strip()[:1000]
    ranked = sorted(merged, key=lambda x: x["score"], reverse=True)
    return {"week_of": week_of, "notes": notes, "ranked": ranked}


# ---- render ----
def render_digest_md(result: dict, items_by_id: dict[str, dict]) -> str:
    week_of = result["week_of"]
    notes = result.get("notes", "").strip()
    ranked = result.get("ranked", [])
    kept = [r for r in ranked if r["score"] >= MIN_SCORE_READ]

    lines = [f"# Weekly ToC Digest (week of {week_of})", ""]
    if notes:
        lines += [notes, ""]
    lines += [
        f"**Included:** {len(kept)} (score ≥ {MIN_SCORE_READ:.2f})  ",
        f"**Scored:** {len(ranked)} total items",
        "",
        f"**Models:** `{CHEAP_MODEL}` first-pass, `{FINAL_MODEL}` re-rank top {REFINE_TOP_N}",
        "",
        "---",
        "",
    ]
    if not kept:
        return "\n".join(lines + ["_No items met the relevance threshold this week._", ""])

    categories = [
        {
            "name": "Neuroblastoma",
            "signals": [
                "neuroblastoma", "nb", "high-risk neuroblastoma", "mycn", "adrn", "mes",
                "alk", "atr", "relapse", "recurrence", "resistance", "refractory",
                "drug", "therapy", "biomarker", "pediatric oncology", "pediatric solid tumors",
            ],
            "must_any": ["neuroblastoma", "nb"],
            "boost": 6,
        },
        {
            "name": "AI",
            "signals": [
                "ai", "machine learning", "deep learning", "foundation model", "scfm",
                "representation learning", "perturbation modeling", "predictive model",
                "transfer learning", "single-cell foundation model",
            ],
            "boost": 3,
        },
        {
            "name": "Methods",
            "signals": [
                "method", "methods", "pipeline", "algorithm", "workflow", "benchmark", "reproducible",
                "multi-omics", "proteomics", "phosphoproteomics", "pathway", "gsva", "gsea",
                "viper", "msviper", "master regulator", "causal inference", "network",
                "single-cell", "scrna-seq", "snrna-seq",
            ],
            "boost": 2,
        },
        {
            "name": "Other",
            "signals": [],
            "boost": 0,
        },
    ]

    def normalized_text(row: dict) -> str:
        it = items_by_id.get(row["id"], {})
        parts = [
            row.get("title", ""),
            row.get("source", ""),
            row.get("why", ""),
            " ".join(row.get("tags", []) or []),
            it.get("summary", ""),
        ]
        return " ".join(parts).lower()

    def score_category(text: str, cfg: dict) -> int:
        score = 0
        signals = cfg.get("signals", [])
        excludes = cfg.get("exclude", [])
        must_any = cfg.get("must_any", [])
        must_all_group = cfg.get("must_all_group", [])

        for sig in signals:
            if sig and sig in text:
                score += 2

        for ex in excludes:
            if ex and ex in text:
                score -= 2

        if must_any and not any(x in text for x in must_any):
            return -999

        for group in must_all_group:
            if not any(x in text for x in group):
                return -999

        score += int(cfg.get("boost", 0))
        return score

    grouped: dict[str, list[dict]] = {c["name"]: [] for c in categories}
    for r in kept:
        text = normalized_text(r)
        scored = [(score_category(text, c), idx, c["name"]) for idx, c in enumerate(categories)]
        # best score wins; index preserves your priority order in ties
        best = sorted(scored, key=lambda x: (x[0], -x[1]), reverse=True)[0]
        target = best[2] if best[0] > -999 else "Other"
        grouped[target].append(r)

    for cname, rows in grouped.items():
        rows.sort(key=lambda x: x["score"], reverse=True)
        visible = rows[:CATEGORY_MAX_ITEMS]
        lines += [f"## {cname} ({len(visible)} shown / {len(rows)} total)", ""]

        for r in visible:
            it = items_by_id.get(r["id"], {})
            tags = ", ".join(r.get("tags", [])) if r.get("tags") else ""
            pub = r.get("published_utc")
            summary = (it.get("summary") or "").strip()

            lines += [
                f"### [{r['title']}]({r['link']})",
                f"*{r['source']}*  ",
                f"Score: **{r['score']:.2f}**" + (f"  \nPublished: {pub}" if pub else ""),
                (f"Tags: {tags}" if tags else ""),
                "",
                r["why"].strip(),
                "",
            ]
            if summary:
                lines += ["<details>", "<summary>RSS summary</summary>", "", summary, "", "</details>", ""]
            lines += ["---", ""]
    return "\n".join(lines)


def markdown_to_html(md: str, title: str) -> str:
    try:
        import markdown

        body = markdown.markdown(md, extensions=["extra", "sane_lists"])
    except Exception:
        body = f"<pre>{html.escape(md)}</pre>"

    return (
        "<!doctype html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        f"  <title>{html.escape(title)}</title>\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "  <style>"
        "body{max-width:900px;margin:40px auto;padding:0 16px;line-height:1.5;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;}"
        "h1,h2,h3{line-height:1.2;} code{background:#f5f5f5;padding:2px 4px;border-radius:4px;}"
        "details{border:1px solid #ddd;border-radius:6px;padding:8px 10px;background:#fafafa;}"
        "a{color:#0b5fff;text-decoration:none;}a:hover{text-decoration:underline;}"
        "</style>\n"
        "</head>\n"
        "<body>\n"
        f"{body}\n"
        "</body>\n"
        "</html>\n"
    )


def write_dated_outputs(md: str, week_of: str) -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(HTML_OUTPUT_DIR, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d_%H%M%S")

    md_snapshot = os.path.join(OUTPUT_DIR, f"digest-{stamp}.md")
    html_snapshot = os.path.join(HTML_OUTPUT_DIR, f"digest-{stamp}.html")

    with open(md_snapshot, "w", encoding="utf-8") as f:
        f.write(md)

    with open(html_snapshot, "w", encoding="utf-8") as f:
        f.write(markdown_to_html(md, title=f"Weekly ToC Digest ({week_of})"))

    # keep convenience pointers to latest render
    with open("latest.md", "w", encoding="utf-8") as f:
        f.write(md)
    with open(os.path.join(HTML_OUTPUT_DIR, "latest.html"), "w", encoding="utf-8") as f:
        f.write(markdown_to_html(md, title=f"Weekly ToC Digest ({week_of})"))

    print(f"Wrote snapshots: {md_snapshot}, {html_snapshot}")


def main():
    interests = parse_interests_md(read_text("interests.md"))
    feeds = load_feeds("feeds.txt")
    items = fetch_rss_items(feeds)
    print(f"Fetched {len(items)} RSS items (pre-filter)")

    today = datetime.now(timezone.utc).date().isoformat()
    if not items:
        md = f"# Weekly ToC Digest (week of {today})\n\n_No RSS items found in the last {LOOKBACK_DAYS} days._\n"
        with open("digest.md", "w", encoding="utf-8") as f:
            f.write(md)
        write_dated_outputs(md, week_of=today)
        print("No items; wrote digest.md")
        return

    items = keyword_prefilter(items, interests["keywords"], keep_top=PREFILTER_KEEP_TOP)
    print(f"Sending {len(items)} RSS items to model (post-filter)")

    items_by_id = {it["id"]: it for it in items}
    client = make_openai_client()

    state = load_state(STATE_PATH)
    result = two_stage_triage(client, interests, items, state)
    md = render_digest_md(result, items_by_id)

    with open("digest.md", "w", encoding="utf-8") as f:
        f.write(md)
    write_dated_outputs(md, week_of=result["week_of"])

    prune_state(state, retention_days=STATE_RETENTION_DAYS)
    save_state(STATE_PATH, state)

    print("Wrote digest.md")


if __name__ == "__main__":
    main()
