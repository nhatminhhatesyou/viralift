"""gbcrawler — Streamlit demo UI (query builder).

Run from the viralift/ folder:

    streamlit run gbcrawler/ui_app.py

The user fills structured fields (virus, length, date, region, …); the app
assembles the NCBI query, previews the match count, crawls GenBank records,
splits them by species against the ViraLift registry, and offers each species
file for download. Same engine as the CLI (gbcrawler/__main__.py).
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import streamlit as st

# Make the sibling modules importable when launched via `streamlit run`.
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import fetch as fetch_mod          # noqa: E402
import split as split_mod          # noqa: E402

_DEFAULT_REGISTRY = _HERE.parent / "app" / "config" / "virus_alias_registry.json"
_DEFAULT_OUT = _HERE.parent / "gbcrawler_out"

# Friendly virus picker -> proper [Organism] term (longest registry keyword).
_FALLBACK_VIRUSES = {
    "PRRSV": "Porcine reproductive and respiratory syndrome virus",
    "FMDV": "Foot-and-mouth disease virus",
    "PEDV": "Porcine epidemic diarrhea virus",
}


def _virus_options(registry_path: Path) -> dict:
    try:
        reg = split_mod.load_registry(registry_path)
        opts = {}
        for name, keywords in reg:
            organism = max(keywords, key=len) if keywords else name
            opts[name] = organism.strip()
        return opts or dict(_FALLBACK_VIRUSES)
    except Exception:  # noqa: BLE001
        return dict(_FALLBACK_VIRUSES)


def build_query(organism, complete_only, length, dates, region, extra) -> str:
    parts = []
    if organism:
        parts.append(f'"{organism}"[Organism]')
    if complete_only:
        parts.append("complete genome[Title]")
    if length:
        lo, hi = length
        parts.append(f'("{lo}"[SLEN] : "{hi}"[SLEN])')
    if dates:
        d0, d1 = dates
        parts.append(f'("{d0:%Y/%m/%d}"[PDAT] : "{d1:%Y/%m/%d}"[PDAT])')
    if region:
        parts.append(f'"{region}"[All Fields]')
    if extra and extra.strip():
        parts.append(extra.strip())
    return " AND ".join(parts)


# ══════════════════════════════════════════════════════════════════════
# Page + theme
# ══════════════════════════════════════════════════════════════════════
st.set_page_config(page_title="gbcrawler", page_icon="🧬", layout="wide")

st.markdown(
    """
    <style>
      :root { color-scheme: light dark; }
      .block-container { padding-top: 2.2rem; max-width: 1080px; }
      #MainMenu, footer { visibility: hidden; }

      .gc-hero {
        background:
          radial-gradient(circle at 88% 8%, rgba(255,255,255,.18), transparent 42%),
          linear-gradient(120deg, #174f49 0%, #1d6c63 52%, #2f8a73 100%);
        color: #f3faf7; border-radius: 18px;
        padding: 26px 30px; margin-bottom: 22px;
        box-shadow: 0 18px 48px rgba(23,79,73,.28);
      }
      .gc-hero h1 { color:#fff; font-size:1.9rem; margin:0 0 4px; font-weight:750; }
      .gc-hero p  { color:#d6ece6; margin:0; font-size:.97rem; }
      .gc-step { font-size:.78rem; letter-spacing:.08em; text-transform:uppercase;
                 color:#1d6c63; font-weight:700; margin-bottom:.3rem; }

      /* bordered containers -> cards */
      div[data-testid="stVerticalBlockBorderWrapper"] {
        border-radius: 14px; border:1px solid rgba(29,108,99,.18) !important;
        box-shadow: 0 8px 28px rgba(42,70,58,.07);
        background: rgba(255,255,255,.55);
      }
      .stButton>button {
        border-radius: 10px; font-weight: 650; border:1px solid rgba(29,108,99,.35);
      }
      .stButton>button[kind="primary"] { background:#1d6c63; border-color:#1d6c63; }
      div[data-testid="stMetricValue"] { color:#1d6c63; }
      code { color:#174f49; }
      .gc-q {
        background:#0f1a17; color:#a7e8d6; border-radius:10px;
        padding:12px 14px; font-family:ui-monospace,Menlo,monospace;
        font-size:.85rem; word-break:break-word; line-height:1.5;
      }
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
    <div class="gc-hero">
      <h1>🧬 gbcrawler</h1>
      <p>Crawl sequences from NCBI GenBank → split by species → input for ViraLift.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

ss = st.session_state
ss.setdefault("count_result", None)
ss.setdefault("crawl_summary", None)
ss.setdefault("out_dir", None)

# ── sidebar: NCBI identity + paths ────────────────────────────────────
with st.sidebar:
    st.subheader("⚙️ NCBI settings")
    email = st.text_input("Email (required)", placeholder="you@lab.org")
    api_key = st.text_input("API key (optional)", type="password",
                            help="With a key: 10 req/s instead of 3.")
    db = st.text_input("Database", value="nuccore")
    st.divider()
    registry_path = st.text_input("Registry (ViraLift)", value=str(_DEFAULT_REGISTRY))
    out_dir = st.text_input("Output folder", value=str(_DEFAULT_OUT))
    use_ledger = st.checkbox("Dedup ledger", value=True)
    ledger_path = st.text_input("Ledger file",
                                value=str(_HERE.parent / "gbcrawler_ledger.txt"),
                                disabled=not use_ledger)

viruses = _virus_options(Path(registry_path))

# ══════════════════════════════════════════════════════════════════════
# 1. Query builder
# ══════════════════════════════════════════════════════════════════════
st.markdown("<div class='gc-step'>Step 1 · Build query</div>", unsafe_allow_html=True)
with st.container(border=True):
    mode = st.radio("Input source", ["By criteria (query)", "By accession"],
                    horizontal=True, label_visibility="collapsed")

    query = ""
    accessions: list[str] = []

    if mode == "By criteria (query)":
        r1c1, r1c2 = st.columns([2, 1])
        with r1c1:
            vlabels = list(viruses.keys()) + ["Other (custom)…"]
            pick = st.selectbox("🦠 Virus", vlabels)
            if pick == "Other (custom)…":
                organism = st.text_input("Organism name",
                                         placeholder="e.g. Classical swine fever virus")
            else:
                organism = viruses[pick]
                st.caption(f"→ `\"{organism}\"[Organism]`")
        with r1c2:
            complete_only = st.toggle("Complete genomes only", value=True,
                                      help="Adds `complete genome[Title]`")

        st.markdown("**Optional filters**")
        f1, f2 = st.columns(2)
        with f1:
            use_len = st.checkbox("📏 Sequence length (bp)")
            length = None
            if use_len:
                lo, hi = st.slider("Length range", 0, 35000, (14000, 16000), step=500)
                length = (lo, hi)
        with f2:
            use_date = st.checkbox("📅 Publication date")
            dates = None
            if use_date:
                dates = st.date_input("Date range",
                                      value=(date(2020, 1, 1), date.today()),
                                      max_value=date.today())
                if not (isinstance(dates, (tuple, list)) and len(dates) == 2):
                    dates = None  # user mid-selection

        g1, g2 = st.columns(2)
        with g1:
            region_choice = st.selectbox(
                "🌍 Region / country",
                ["Global (no filter)", "Vietnam", "Custom…"],
            )
            region = ""
            if region_choice == "Vietnam":
                region = "Vietnam"
            elif region_choice == "Custom…":
                region = st.text_input("Country / region", placeholder="e.g. China")
            if region:
                st.caption("⚠️ Region filtering on nuccore is approximate — verify with Count.")
        with g2:
            extra = st.text_input("➕ Extra condition (NCBI syntax)",
                                  placeholder='e.g. NOT partial[Title]')

        query = build_query(organism, complete_only, length, dates, region, extra)

        st.markdown("**Query sent to NCBI:**")
        st.markdown(f"<div class='gc-q'>{query or '— not enough criteria yet —'}</div>",
                    unsafe_allow_html=True)
        with st.expander("✏️ Edit query manually (advanced)"):
            manual = st.text_area("Override query", value=query, height=80)
            if manual.strip() and manual.strip() != query:
                query = manual.strip()

    else:  # accession mode
        up = st.file_uploader("Accession file (.txt/.csv)", type=["txt", "csv"])
        text = st.text_area("…or paste accessions (one per line)", height=110,
                            placeholder="PP209408.1\nAF176348.2")
        raw = ""
        if up is not None:
            raw += up.getvalue().decode("utf-8", "ignore") + "\n"
        raw += text
        for line in raw.splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                accessions.extend(p.strip() for p in line.replace(",", " ").split())
        st.caption(f"{len(accessions)} accessions detected.")

    retmax = st.number_input("Max records to fetch (retmax)", 1, 100000, 1000, 100)

# ══════════════════════════════════════════════════════════════════════
# 2. Actions
# ══════════════════════════════════════════════════════════════════════
st.markdown("<div class='gc-step'>Step 2 · Preview count, then fetch</div>", unsafe_allow_html=True)
cc1, cc2 = st.columns(2)
do_count = cc1.button("① Count only (no download)", use_container_width=True,
                      disabled=(mode != "By criteria (query)"))
do_crawl = cc2.button("② Crawl + split by species", type="primary", use_container_width=True)


def _need_email() -> bool:
    if not email.strip():
        st.error("Enter your email in the sidebar — NCBI requires it.")
        return True
    return False


if do_count:
    if not _need_email():
        if not query.strip():
            st.error("Query has no criteria yet.")
        else:
            fetch_mod.configure(email, api_key or None)
            with st.spinner("Counting on NCBI…"):
                try:
                    ss.count_result = fetch_mod.count_query(query, db=db,
                                                            api_key=api_key or None)
                except Exception as e:  # noqa: BLE001
                    st.error(f"Query error: {e}")
                    ss.count_result = None

if ss.count_result is not None:
    n = ss.count_result
    st.metric("Records matching query", f"{n:,}")
    if n == 0:
        st.warning("0 records — query is wrong or too strict.")
    elif n > 50000:
        st.warning("Very many — the query may be too broad.")
    else:
        st.success("Looks reasonable. Press ② to download.")

if do_crawl and not _need_email():
    if mode == "By criteria (query)" and not query.strip():
        st.error("Query has no criteria yet.")
    elif mode == "By accession" and not accessions:
        st.error("No accessions provided.")
    else:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        raw_path = out / "raw_combined.gb"
        try:
            fetch_mod.configure(email, api_key or None)
            with st.spinner("Downloading from NCBI… (large genomes may take a few minutes)"):
                with open(raw_path, "w") as fh:
                    if mode == "By criteria (query)":
                        fetch_mod.fetch_from_query(query, fh, db=db, retmax=int(retmax),
                                                   api_key=api_key or None, quiet=True)
                    else:
                        fetch_mod.fetch_from_accessions(accessions, fh, db=db,
                                                        api_key=api_key or None, quiet=True)
            registry = split_mod.load_registry(Path(registry_path))
            seen = split_mod.load_ledger(Path(ledger_path)) if (use_ledger and ledger_path) else set()
            with st.spinner("Splitting by species…"):
                summary = split_mod.split_genbank(raw_path, registry, out, seen=seen)
            if use_ledger and ledger_path:
                split_mod.append_ledger(Path(ledger_path), summary["new_accessions"])
            ss.crawl_summary = summary
            ss.out_dir = str(out)
        except Exception as e:  # noqa: BLE001
            st.error(f"Crawl error: {e}")

# ══════════════════════════════════════════════════════════════════════
# 3. Results
# ══════════════════════════════════════════════════════════════════════
summary = ss.crawl_summary
if summary:
    st.markdown("<div class='gc-step'>Step 3 · Results</div>", unsafe_allow_html=True)
    with st.container(border=True):
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total fetched", summary["total"])
        m2.metric("New", summary["new"])
        m3.metric("Dup (batch)", summary["dup_in_batch"])
        m4.metric("Dup (ledger)", summary["dup_in_ledger"])

        st.markdown("**Species files — feed to ViraLift:**")
        for key, count in sorted(summary["counts"].items()):
            fpath = Path(summary["files"][key])
            label = key if key != split_mod.UNMATCHED else "(unmatched)"
            a, b = st.columns([3, 1])
            a.write(f"🧬 **{label}** — {count} records · `{fpath.name}`")
            if fpath.exists():
                b.download_button("⬇️ .gb", data=fpath.read_bytes(),
                                  file_name=fpath.name, key=f"dl_{key}",
                                  use_container_width=True)

        man = Path(summary["manifest"])
        if man.exists():
            st.markdown("**manifest.csv**")
            try:
                import pandas as pd
                st.dataframe(pd.read_csv(man), use_container_width=True, height=300)
            except Exception:  # noqa: BLE001
                st.text(man.read_text()[:4000])
            st.download_button("⬇️ manifest.csv", data=man.read_bytes(),
                               file_name="manifest.csv")

        st.info(
            "Next, inside `viralift/`:\n\n```\npython -m app.src.main "
            f"--reference <ref.gb> --query {ss.out_dir}/<virus>.gb "
            "--output output/run\n```"
        )
