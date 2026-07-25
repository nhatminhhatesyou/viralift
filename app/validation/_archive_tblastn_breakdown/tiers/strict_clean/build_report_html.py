#!/usr/bin/env python3
"""Build a self-contained, visually rich HTML report from the ViraLift validation report.

All referenced figures are embedded as base64 so the resulting HTML is a single
portable file that can be sent to a reviewer and opened anywhere.
"""
import base64
from pathlib import Path

BASE = Path(__file__).parent

IMAGES = {
    "overall": "summary_outputs/fmd_prrsv_overall_accuracy.png",
    "per_gene": "summary_outputs/fmd_prrsv_per_gene_accuracy.png",
    "fmd_extrap": "terminal_extrapolation_outputs/fmd_accuracy_comparison.png",
    "prrsv_rescue": "prrsv_start_rescue_full_outputs/prrsv_start_rescue_exact_comparison.png",
    "orf7": "orf7_start_rescue_outputs/orf7_rescue_accuracy_comparison.png",
    "fmd_blame": "outputs_fmd/fmd_final_blame_split.png",
    "prrsv_orf7_delta": "outputs_prrsv/prrsv_orf7_delta_patterns.png",
}


def img_data_uri(rel_path: str) -> str:
    data = (BASE / rel_path).read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    return f"data:image/png;base64,{b64}"


def main():
    imgs = {k: img_data_uri(v) for k, v in IMAGES.items()}

    html = TEMPLATE
    for key, uri in imgs.items():
        html = html.replace(f"{{{{IMG_{key.upper()}}}}}", uri)

    out = BASE / "VIRALIFT_TOOL_AND_VALIDATION_REPORT.html"
    out.write_text(html, encoding="utf-8")
    size_kb = out.stat().st_size / 1024
    print(f"Wrote {out} ({size_kb:.0f} KB)")


TEMPLATE = r"""<!DOCTYPE html>
<html lang="vi">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ViraLift — Báo Cáo Tool & Validation</title>
<style>
  :root {
    --bg: #0b0f1a;
    --bg-soft: #111726;
    --card: #151c2e;
    --card-hover: #1b2440;
    --border: #233049;
    --text: #e6ecf5;
    --text-dim: #98a6bd;
    --text-faint: #65728c;
    --accent: #4f9cff;
    --accent-2: #7c5cff;
    --green: #34d399;
    --green-soft: #6ee7b7;
    --red: #f87171;
    --amber: #fbbf24;
    --mono: "SF Mono", "JetBrains Mono", "Fira Code", ui-monospace, Menlo, Consolas, monospace;
    --sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  html { scroll-behavior: smooth; scroll-padding-top: 24px; }
  body {
    background: radial-gradient(1200px 800px at 80% -10%, #1a2240 0%, transparent 55%),
                radial-gradient(900px 700px at -10% 10%, #1d1538 0%, transparent 50%),
                var(--bg);
    color: var(--text);
    font-family: var(--sans);
    line-height: 1.65;
    font-size: 16px;
    -webkit-font-smoothing: antialiased;
  }
  .layout { display: flex; max-width: 1400px; margin: 0 auto; }

  /* Sidebar */
  .sidebar {
    width: 280px;
    flex-shrink: 0;
    position: sticky;
    top: 0;
    align-self: flex-start;
    height: 100vh;
    overflow-y: auto;
    padding: 32px 20px;
    border-right: 1px solid var(--border);
  }
  .sidebar .brand {
    display: flex; align-items: center; gap: 10px;
    font-weight: 800; font-size: 20px; letter-spacing: -0.02em;
    margin-bottom: 6px;
  }
  .brand .logo {
    width: 34px; height: 34px; border-radius: 9px;
    background: linear-gradient(135deg, var(--accent), var(--accent-2));
    display: grid; place-items: center; font-size: 18px;
    box-shadow: 0 6px 18px rgba(79,156,255,0.35);
  }
  .brand .grad { background: linear-gradient(90deg,#7fb6ff,#b79cff); -webkit-background-clip: text; background-clip: text; color: transparent; }
  .sidebar .tag { font-size: 12px; color: var(--text-faint); margin-bottom: 26px; padding-left: 2px; }
  .toc a {
    display: block;
    color: var(--text-dim);
    text-decoration: none;
    font-size: 14px;
    padding: 8px 12px;
    border-radius: 8px;
    border-left: 2px solid transparent;
    transition: all .15s ease;
  }
  .toc a:hover { background: var(--card); color: var(--text); }
  .toc a.active { background: var(--card-hover); color: #fff; border-left-color: var(--accent); }
  .toc .num { color: var(--text-faint); font-variant-numeric: tabular-nums; margin-right: 8px; }

  /* Main */
  main { flex: 1; min-width: 0; padding: 48px 56px 120px; }

  /* Hero */
  .hero {
    border: 1px solid var(--border);
    border-radius: 20px;
    padding: 44px 40px;
    background:
      radial-gradient(600px 300px at 90% -40%, rgba(124,92,255,0.20), transparent 60%),
      linear-gradient(180deg, rgba(79,156,255,0.07), transparent 70%),
      var(--bg-soft);
    margin-bottom: 14px;
    position: relative;
    overflow: hidden;
  }
  .hero .eyebrow { color: var(--accent); font-size: 13px; font-weight: 700; letter-spacing: .12em; text-transform: uppercase; }
  .hero h1 { font-size: 42px; line-height: 1.1; letter-spacing: -0.03em; margin: 10px 0 14px; font-weight: 800; }
  .hero h1 .grad { background: linear-gradient(90deg,#7fb6ff,#b79cff); -webkit-background-clip: text; background-clip: text; color: transparent; }
  .hero p { color: var(--text-dim); max-width: 680px; font-size: 17px; }
  .hero .pills { display: flex; gap: 10px; margin-top: 22px; flex-wrap: wrap; }
  .pill {
    font-size: 13px; padding: 7px 14px; border-radius: 999px;
    border: 1px solid var(--border); background: rgba(255,255,255,0.03);
    color: var(--text-dim); display: inline-flex; align-items: center; gap: 7px;
  }
  .pill b { color: var(--text); }
  .dot { width: 8px; height: 8px; border-radius: 50%; display: inline-block; }
  .dot.green { background: var(--green); box-shadow: 0 0 10px var(--green); }
  .dot.blue { background: var(--accent); box-shadow: 0 0 10px var(--accent); }
  .dot.purple { background: var(--accent-2); box-shadow: 0 0 10px var(--accent-2); }

  /* Headline stat band */
  .stat-band { display: grid; grid-template-columns: repeat(3, 1fr); gap: 16px; margin: 26px 0 8px; }
  .stat {
    border: 1px solid var(--border); border-radius: 16px; padding: 24px;
    background: linear-gradient(180deg, var(--card), var(--bg-soft));
    position: relative; overflow: hidden;
  }
  .stat::after {
    content: ""; position: absolute; inset: 0; pointer-events: none;
    background: radial-gradient(300px 120px at 100% 0%, rgba(52,211,153,0.10), transparent 70%);
  }
  .stat .label { font-size: 13px; color: var(--text-dim); text-transform: uppercase; letter-spacing: .08em; }
  .stat .value { font-size: 44px; font-weight: 800; letter-spacing: -0.03em; margin-top: 6px; font-variant-numeric: tabular-nums; }
  .stat .value.green { color: var(--green-soft); }
  .stat .value.blue { color: #8cc0ff; }
  .stat .value.purple { color: #b79cff; }
  .stat .sub { font-size: 13px; color: var(--text-faint); margin-top: 4px; }

  /* Sections */
  section { scroll-margin-top: 24px; padding-top: 56px; }
  section > h2 {
    font-size: 28px; letter-spacing: -0.02em; font-weight: 800;
    display: flex; align-items: baseline; gap: 14px; margin-bottom: 4px;
  }
  section > h2 .idx {
    font-size: 15px; color: var(--accent); font-weight: 700;
    border: 1px solid var(--border); border-radius: 8px; padding: 2px 10px;
    background: rgba(79,156,255,0.08); flex-shrink: 0;
  }
  .lead { color: var(--text-dim); font-size: 17px; margin: 10px 0 22px; max-width: 760px; }
  h3 { font-size: 19px; margin: 28px 0 12px; font-weight: 700; letter-spacing: -0.01em; }
  p { margin: 12px 0; color: #cdd6e6; }
  ul { margin: 12px 0 12px 4px; }
  li { margin: 7px 0; padding-left: 22px; position: relative; color: #cdd6e6; }
  li::before { content: ""; position: absolute; left: 4px; top: 11px; width: 6px; height: 6px; border-radius: 50%; background: var(--accent); }
  ol { margin: 12px 0 12px 22px; }
  ol li { padding-left: 8px; }
  ol li::before { display: none; }
  strong, b { color: #fff; font-weight: 700; }
  code {
    font-family: var(--mono); font-size: 0.88em;
    background: rgba(124,92,255,0.14); color: #c9bfff;
    padding: 2px 7px; border-radius: 6px; border: 1px solid rgba(124,92,255,0.18);
  }

  /* Cards */
  .card {
    border: 1px solid var(--border); border-radius: 16px;
    background: var(--card); padding: 26px 28px; margin: 18px 0;
  }
  .card.flow { background: linear-gradient(180deg, var(--card), var(--bg-soft)); }

  /* Pipeline steps */
  .steps { counter-reset: step; display: grid; gap: 12px; margin: 18px 0; }
  .step {
    display: flex; gap: 16px; align-items: flex-start;
    border: 1px solid var(--border); border-radius: 14px; padding: 16px 18px;
    background: var(--card); transition: background .15s, border-color .15s;
  }
  .step:hover { background: var(--card-hover); border-color: #33456a; }
  .step .n {
    counter-increment: step; flex-shrink: 0;
    width: 30px; height: 30px; border-radius: 9px;
    background: linear-gradient(135deg, var(--accent), var(--accent-2));
    color: #fff; font-weight: 800; display: grid; place-items: center; font-size: 14px;
  }
  .step .n::before { content: counter(step); }
  .step .body { flex: 1; }
  .step .body b { display: block; margin-bottom: 2px; }
  .step .body span { color: var(--text-dim); font-size: 14.5px; }

  /* Code blocks */
  pre {
    background: #0a0e18; border: 1px solid var(--border); border-radius: 12px;
    padding: 18px 20px; overflow-x: auto; margin: 16px 0;
    font-family: var(--mono); font-size: 13.5px; line-height: 1.6; color: #c7d2e6;
  }
  pre .c1 { color: #5fd3a3; }      /* good / highlight */
  pre .c2 { color: #f3a86b; }      /* annotation */
  pre .c3 { color: var(--text-faint); }
  pre .k { color: #7fb6ff; }

  /* Status legend chips */
  .legend { display: flex; gap: 10px; flex-wrap: wrap; margin: 16px 0; }
  .chip {
    font-size: 13px; padding: 6px 12px; border-radius: 8px;
    border: 1px solid var(--border); display: inline-flex; gap: 8px; align-items: center;
    font-family: var(--mono); background: rgba(255,255,255,0.02);
  }
  .chip .b { width: 9px; height: 9px; border-radius: 3px; }
  .chip.ok .b { background: var(--green); }
  .chip.rescue .b { background: var(--accent); }
  .chip.extrap .b { background: var(--accent-2); }
  .chip.nohit .b { background: var(--text-faint); }
  .chip.invalid .b { background: var(--red); }

  /* Tables */
  .table-wrap { overflow-x: auto; border: 1px solid var(--border); border-radius: 14px; margin: 18px 0; }
  table { width: 100%; border-collapse: collapse; font-size: 14.5px; }
  thead th {
    background: var(--bg-soft); color: var(--text-dim);
    text-align: left; padding: 13px 16px; font-weight: 600; font-size: 13px;
    text-transform: uppercase; letter-spacing: .05em; border-bottom: 1px solid var(--border);
    white-space: nowrap;
  }
  thead th.r, tbody td.r { text-align: right; font-variant-numeric: tabular-nums; }
  tbody td { padding: 12px 16px; border-bottom: 1px solid rgba(35,48,73,0.6); }
  tbody tr:last-child td { border-bottom: none; }
  tbody tr:hover { background: rgba(79,156,255,0.05); }
  tbody td.gene { font-family: var(--mono); color: #c9bfff; font-weight: 600; }
  .acc { font-weight: 700; font-variant-numeric: tabular-nums; }
  .acc.perfect { color: var(--green-soft); }
  .acc.high { color: #8cc0ff; }
  .acc.mid { color: var(--amber); }
  tr.total-row { font-weight: 700; }
  tr.total-row td { background: rgba(124,92,255,0.08); border-top: 1px solid var(--border); }
  /* mini accuracy bar */
  .bar { position: relative; height: 7px; border-radius: 4px; background: rgba(255,255,255,0.06); overflow: hidden; min-width: 90px; }
  .bar > span { position: absolute; left: 0; top: 0; bottom: 0; border-radius: 4px; }
  .bar .exact { background: var(--green); }
  .bar .coord { background: var(--accent); }
  .bar .fail { background: var(--red); }
  td .barcell { display: flex; align-items: center; gap: 10px; }

  /* Figure */
  figure {
    border: 1px solid var(--border); border-radius: 16px; overflow: hidden;
    margin: 22px 0; background: #fff;
  }
  figure img { display: block; width: 100%; height: auto; }
  figcaption {
    background: var(--bg-soft); color: var(--text-dim);
    padding: 12px 18px; font-size: 13.5px; border-top: 1px solid var(--border);
    display: flex; gap: 10px; align-items: center;
  }
  figcaption .fig-tag {
    font-size: 11px; font-weight: 700; letter-spacing: .06em; text-transform: uppercase;
    color: var(--accent); border: 1px solid var(--border); border-radius: 6px; padding: 2px 8px;
    flex-shrink: 0; background: rgba(79,156,255,0.08);
  }

  /* Callouts */
  .callout {
    border-radius: 14px; padding: 18px 22px; margin: 20px 0;
    border: 1px solid var(--border); border-left-width: 4px;
    background: var(--card);
  }
  .callout.find { border-left-color: var(--amber); background: linear-gradient(90deg, rgba(251,191,36,0.07), transparent 60%); }
  .callout.note { border-left-color: var(--accent); background: linear-gradient(90deg, rgba(79,156,255,0.07), transparent 60%); }
  .callout.win  { border-left-color: var(--green); background: linear-gradient(90deg, rgba(52,211,153,0.07), transparent 60%); }
  .callout .ttl { font-weight: 800; margin-bottom: 6px; display: flex; gap: 8px; align-items: center; }

  /* phase narrative (baseline -> finding -> fix -> after) */
  .phase-rail { display: flex; gap: 8px; margin: 22px 0 8px; flex-wrap: wrap; }
  .phase-rail .pstep {
    flex: 1; min-width: 130px; border: 1px solid var(--border); border-radius: 10px;
    padding: 10px 12px; background: var(--card); font-size: 12.5px; color: var(--text-dim);
    display: flex; align-items: center; gap: 8px;
  }
  .phase-rail .pstep .pn {
    width: 20px; height: 20px; border-radius: 6px; flex-shrink: 0;
    display: grid; place-items: center; font-size: 11px; font-weight: 800; color: #0b0f1a;
  }
  .phase-rail .pstep.p1 .pn { background: var(--text-faint); }
  .phase-rail .pstep.p2 .pn { background: var(--amber); }
  .phase-rail .pstep.p3 .pn { background: var(--accent); }
  .phase-rail .pstep.p4 .pn { background: var(--green); }
  .phase-rail .arrow { color: var(--text-faint); align-self: center; font-size: 14px; }

  .phase { display: flex; align-items: center; gap: 12px; margin: 34px 0 4px; }
  .phase .badge {
    flex-shrink: 0; width: 30px; height: 30px; border-radius: 9px;
    display: grid; place-items: center; font-weight: 800; font-size: 14px; color: #0b0f1a;
  }
  .phase.p1 .badge { background: var(--text-faint); }
  .phase.p2 .badge { background: var(--amber); }
  .phase.p3 .badge { background: var(--accent); }
  .phase.p4 .badge { background: var(--green); }
  .phase h3 { margin: 0; font-size: 19px; }
  .phase .kicker { font-size: 12px; text-transform: uppercase; letter-spacing: .1em; color: var(--text-faint); display: block; }

  /* phase block wrapper — makes each of the 4 steps a distinct boxed unit */
  .pblock { border: 1px solid var(--border); border-left-width: 4px; border-radius: 14px; padding: 4px 26px 22px; margin: 16px 0; background: var(--bg-soft); }
  .pblock.p1 { border-left-color: var(--text-faint); }
  .pblock.p2 { border-left-color: var(--amber); }
  .pblock.p3 { border-left-color: var(--accent); }
  .pblock.p4 { border-left-color: var(--green); }
  .pblock .phase { margin-top: 18px; }
  .pblock > p:first-of-type { margin-top: 6px; }

  /* gene status grid — per-gene picture, scannable instead of prose */
  .gene-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(165px, 1fr)); gap: 11px; margin: 16px 0; }
  .gcard { border: 1px solid var(--border); border-radius: 12px; padding: 13px 15px; background: var(--card); display: flex; flex-direction: column; gap: 6px; }
  .gcard.hot { border-color: rgba(248,113,113,0.55); box-shadow: 0 0 0 1px rgba(248,113,113,0.18); background: linear-gradient(180deg, rgba(248,113,113,0.06), var(--card)); }
  .gcard .gname { font-family: var(--mono); font-weight: 700; font-size: 15px; color: #fff; display: flex; align-items: center; justify-content: space-between; gap: 8px; }
  .gcard .gnum { font-family: var(--mono); font-size: 12.5px; color: var(--text-dim); font-variant-numeric: tabular-nums; font-weight: 600; }
  .gcard .gnote { font-size: 12.5px; color: var(--text-dim); line-height: 1.45; }
  .gstat { align-self: flex-start; font-size: 10.5px; font-weight: 800; padding: 3px 9px; border-radius: 999px; letter-spacing: .04em; text-transform: uppercase; }
  .gstat.clean  { background: rgba(52,211,153,0.15); color: var(--green-soft); }
  .gstat.conv   { background: rgba(79,156,255,0.15); color: #8cc0ff; }
  .gstat.minor  { background: rgba(251,191,36,0.16); color: var(--amber); }
  .gstat.invest { background: rgba(248,113,113,0.16); color: #fda4a4; }
  .gstat.fixed  { background: rgba(52,211,153,0.15); color: var(--green-soft); }
  .grid-legend { display: flex; gap: 14px; flex-wrap: wrap; margin: 4px 0 0; font-size: 12.5px; color: var(--text-dim); }
  .grid-legend span { display: inline-flex; align-items: center; gap: 6px; }
  .grid-legend i { width: 10px; height: 10px; border-radius: 3px; display: inline-block; }

  /* numbered findings */
  .findings { display: grid; gap: 14px; margin: 18px 0; }
  .finding {
    display: flex; gap: 18px; border: 1px solid var(--border); border-radius: 14px;
    padding: 18px 22px; background: var(--card);
  }
  .finding .fn {
    flex-shrink: 0; font-size: 26px; font-weight: 800; line-height: 1;
    background: linear-gradient(135deg,#7fb6ff,#b79cff); -webkit-background-clip: text; background-clip: text; color: transparent;
    width: 36px; font-variant-numeric: tabular-nums;
  }
  .finding .fbody b { color: #fff; }
  .finding .fbody p { margin: 4px 0 0; font-size: 14.5px; color: var(--text-dim); }

  /* before/after delta badges */
  .delta { display: inline-flex; align-items: center; gap: 6px; font-family: var(--mono); font-size: 14px; }
  .delta .from { color: var(--text-faint); text-decoration: line-through; }
  .delta .arrow { color: var(--text-faint); }
  .delta .to { color: var(--green-soft); font-weight: 700; }

  .quote {
    border-radius: 16px; padding: 26px 30px; margin: 24px 0;
    background: linear-gradient(135deg, rgba(79,156,255,0.10), rgba(124,92,255,0.10));
    border: 1px solid var(--border); font-size: 18px; line-height: 1.6; color: #e9eefb;
    position: relative; font-style: italic;
  }
  .quote::before { content: "“"; position: absolute; top: -6px; left: 16px; font-size: 64px; color: rgba(124,92,255,0.4); font-family: Georgia, serif; }
  .quote .body { padding-left: 30px; }

  footer { margin-top: 70px; padding-top: 26px; border-top: 1px solid var(--border); color: var(--text-faint); font-size: 13.5px; text-align: center; }

  /* responsive */
  @media (max-width: 1000px) {
    .sidebar { display: none; }
    main { padding: 32px 22px 90px; }
    .hero h1 { font-size: 32px; }
    .stat-band { grid-template-columns: 1fr; }
  }
</style>
</head>
<body>
<div class="layout">
  <aside class="sidebar">
    <div class="brand"><span class="logo">🧬</span><span>Vira<span class="grad">Lift</span></span></div>
    <div class="tag">Tool &amp; Validation Report</div>
    <nav class="toc">
      <a href="#s1"><span class="num">01</span>ViraLift Là Gì?</a>
      <a href="#s2"><span class="num">02</span>Tool Chạy Thế Nào?</a>
      <a href="#s3"><span class="num">03</span>Ví Dụ Sử Dụng</a>
      <a href="#s4"><span class="num">04</span>Validation Dataset</a>
      <a href="#s5"><span class="num">05</span>Cách Tính Accuracy</a>
      <a href="#s6"><span class="num">06</span>Kết Quả Tổng Quan</a>
      <a href="#s7"><span class="num">07</span>Breakdown: FMD</a>
      <a href="#s8"><span class="num">08</span>Breakdown: PRRSV</a>
      <a href="#s9"><span class="num">09</span>Finding Chính</a>
      <a href="#s10"><span class="num">10</span>Kết Luận</a>
    </nav>
  </aside>

  <main>
    <!-- HERO -->
    <div class="hero">
      <div class="eyebrow">Báo cáo Tool &amp; Validation</div>
      <h1>Vira<span class="grad">Lift</span><br>Chuẩn hóa &amp; lift annotation gene cho virus genome</h1>
      <p>Tool hỗ trợ chuẩn hóa tên gene từ GenBank và chuyển annotation sang genome chưa được annotate, bằng reference-guided <code>tblastn</code> lifting. Validate trên hai bộ dữ liệu strict-clean: <b>FMDV</b> và <b>PRRSV</b>.</p>
      <div class="pills">
        <span class="pill"><span class="dot green"></span>FMD <b>98.35%</b></span>
        <span class="pill"><span class="dot blue"></span>PRRSV <b>99.61%</b></span>
        <span class="pill"><span class="dot purple"></span>1910 gene-record cases</span>
        <span class="pill">tblastn lifting</span>
      </div>
    </div>

    <div class="stat-band">
      <div class="stat"><div class="label">FMD Accuracy</div><div class="value green">98.35%</div><div class="sub">1114 exact + 16 coord-only / 1149</div></div>
      <div class="stat"><div class="label">PRRSV Accuracy</div><div class="value blue">99.61%</div><div class="sub">674 exact + 84 coord-only / 761</div></div>
      <div class="stat"><div class="label">Tổng 2 tập</div><div class="value purple">98.85%</div><div class="sub">1788 exact + 100 coord-only / 1910</div></div>
    </div>

    <!-- 1 -->
    <section id="s1">
      <h2><span class="idx">01</span>ViraLift Là Gì?</h2>
      <p class="lead">ViraLift là tool hỗ trợ chuẩn hóa và chuyển annotation gene/peptide cho virus genome.</p>
      <p>Mục tiêu chính:</p>
      <ul>
        <li>Chuẩn hóa tên gene từ nhiều cách đặt tên khác nhau trong GenBank.</li>
        <li>Nếu query genome đã có annotation đủ tốt, tool trích xuất trực tiếp từ annotation.</li>
        <li>Nếu query genome chưa có annotation hoặc annotation thiếu, tool dùng reference chuẩn do user cung cấp để lift gene bằng <code>tblastn</code>.</li>
        <li>Xuất kết quả thành bảng TSV và FASTA để user kiểm tra hoặc dùng tiếp.</li>
      </ul>
      <p>Tool hiện được validate chính trên 2 nhóm virus:</p>
      <ul>
        <li><b>FMDV</b>: dùng <code>mat_peptide</code>, gene/peptide thường liên tiếp.</li>
        <li><b>PRRSV</b>: dùng <code>CDS</code>/ORF, có nhiều gene chồng lấp và convention annotation phức tạp hơn.</li>
      </ul>
    </section>

    <!-- 2 -->
    <section id="s2">
      <h2><span class="idx">02</span>Tool Chạy Như Thế Nào?</h2>
      <p class="lead">Pipeline chính từ reference GenBank đến bảng prediction + FASTA.</p>
      <div class="steps">
        <div class="step"><div class="n"></div><div class="body"><b>User cung cấp reference GenBank đã annotation chuẩn.</b></div></div>
        <div class="step"><div class="n"></div><div class="body"><b>Đọc reference &amp; xác định feature type hữu ích</b><span>ví dụ <code>CDS</code> hoặc <code>mat_peptide</code>.</span></div></div>
        <div class="step"><div class="n"></div><div class="body"><b>Chuẩn hóa tên gene bằng alias map.</b></div></div>
        <div class="step"><div class="n"></div><div class="body"><b>Với mỗi query genome</b><span>Nếu có annotation đủ hữu ích → direct extract. Nếu thiếu → dùng <code>tblastn</code> để map protein từ reference sang query.</span></div></div>
        <div class="step"><div class="n"></div><div class="body"><b>Validate boundary</b><span>Với <code>CDS</code>: kiểm tra start/stop codon &amp; frame. Với <code>mat_peptide</code>: không dùng start/stop codon vì cleavage product không nhất thiết có ATG/stop riêng.</span></div></div>
        <div class="step"><div class="n"></div><div class="body"><b>Xuất kết quả</b><span>bảng prediction · FASTA sequence · run summary · status cho từng gene.</span></div></div>
      </div>

      <div class="card flow">
        <h3 style="margin-top:0">Với tblastn lifting</h3>
        <ul>
          <li>Reference feature được dịch sang protein.</li>
          <li>Protein được search trên query genome bằng <code>tblastn</code>.</li>
          <li>HSPs được merge để suy ra tọa độ nucleotide.</li>
          <li>Tool chỉnh boundary nếu cần — ví dụ <b>terminal extrapolation</b> cho FMD hoặc <b>start-codon rescue</b> cho PRRSV.</li>
        </ul>
        <div class="callout note">
          <div class="ttl">💡 HSP là gì?</div>
          <span style="color:var(--text-dim)"><code>HSP</code> là một đoạn alignment tốt do BLAST trả về. Một gene có thể có một hoặc nhiều HSP. ViraLift dùng các HSP này để suy ra vùng nucleotide tương ứng trên query genome.</span>
        </div>
      </div>
    </section>

    <!-- 3 -->
    <section id="s3">
      <h2><span class="idx">03</span>Ví Dụ Sử Dụng</h2>
      <p class="lead">Chạy bằng CLI và đọc status output.</p>
      <pre><span class="k">python</span> -m app.src.main \
  --reference app/data/PRRS_ref_test.gb \
  --query app/data/PRRS_PP946131_noAnno.gb \
  --output output/prrs_example</pre>

      <h3>Ví dụ logic kết quả</h3>
      <pre><span class="c1">ORF5</span>  -> lifted bằng tblastn, tọa độ đúng, status <span class="c1">ok</span>
<span class="c1">ORF7</span>  -> lifted bằng tblastn, start được rescue, status <span class="k">ok_rescued</span>
<span class="c1">ORF2b</span> -> lifted đúng nhưng một số query truth có thể không annotate ORF2b riêng</pre>

      <h3>Ý nghĩa status thường gặp</h3>
      <div class="legend">
        <span class="chip ok"><span class="b"></span>ok — lift thành công, boundary hợp lệ</span>
        <span class="chip rescue"><span class="b"></span>ok_rescued — boundary đã được rescue lại</span>
        <span class="chip extrap"><span class="b"></span>ok_extrapolated — mở rộng bằng terminal extrapolation</span>
        <span class="chip nohit"><span class="b"></span>no_hit — không tìm được hit phù hợp</span>
        <span class="chip invalid"><span class="b"></span>invalid_boundaries — boundary/codon validation chưa đạt</span>
      </div>
      <div class="callout note">
        <span style="color:var(--text-dim)">Với <code>CDS</code>, <code>invalid_boundaries</code> có thể xảy ra nếu thiếu start codon, thiếu stop codon, hoặc CDS length không chia hết cho 3 (<code>in_frame = false</code>).</span>
      </div>
    </section>

    <!-- 4 -->
    <section id="s4">
      <h2><span class="idx">04</span>Validation Dataset</h2>
      <p class="lead">Validation dùng 2 bộ virus annotated records: <b>FMDV</b> và <b>PRRSV</b> strict-clean.</p>
      <h3>Nguyên tắc lọc strict-clean</h3>
      <ul>
        <li>Chỉ giữ record có annotation đủ để làm ground truth.</li>
        <li>Tên gene phải map được về canonical name qua alias map.</li>
        <li>Với PRRSV, không dùng nested-feature filter mù vì <code>ORF2b</code> có thể overlap/nest trong <code>ORF2a</code> nhưng vẫn là gene thật.</li>
        <li>Với gene không xuất hiện trong truth của một query record, không dùng record đó để chấm accuracy cho gene đó.</li>
      </ul>
      <div class="callout find">
        <div class="ttl">📌 Ví dụ truth-availability</div>
        <span style="color:var(--text-dim)">Reference có <code>ORF7</code>. Trong 95 PRRSV records, truth có <code>ORF7</code> ở 94 records. Vậy accuracy của <code>ORF7</code> chỉ tính trên 94 records đó — tránh đánh giá sai tool vì truth thiếu hoặc dùng convention khác.</span>
      </div>
    </section>

    <!-- 5 -->
    <section id="s5">
      <h2><span class="idx">05</span>Cách Tính Accuracy</h2>
      <p class="lead">Validation được tính theo từng gene.</p>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Metric</th><th>Ý nghĩa</th></tr></thead>
          <tbody>
            <tr><td class="gene">total</td><td>số query records thật sự có gene đó trong truth</td></tr>
            <tr><td class="gene">exact</td><td>prediction khớp hoàn toàn tên gene + start + end</td></tr>
            <tr><td class="gene">coord_only</td><td>tọa độ đúng theo IoU threshold nhưng không exact tuyệt đối</td></tr>
            <tr><td class="gene">failed</td><td>không exact và cũng không coordinate-correct</td></tr>
            <tr><td class="gene">accuracy_pct</td><td>(exact + coord_only) / total</td></tr>
            <tr><td class="gene">IoU</td><td>chỉ số đo mức overlap giữa tọa độ prediction và truth</td></tr>
          </tbody>
        </table>
      </div>

      <h3>IoU là gì?</h3>
      <p><code>IoU</code> = <b>Intersection over Union</b>: độ dài phần prediction &amp; truth chồng lên nhau, chia cho độ dài vùng bao phủ bởi cả hai.</p>
      <div class="stat-band" style="grid-template-columns:repeat(3,1fr); margin:18px 0;">
        <div class="card" style="margin:0">
          <b style="color:#fff">Khớp hoàn toàn</b>
          <pre style="margin:10px 0 0;font-size:12.5px">Truth:      100 - 199
Pred:       100 - 199
<span class="c1">IoU = 100/100 = 1.00</span></pre>
        </div>
        <div class="card" style="margin:0">
          <b style="color:#fff">Lệch một chút</b>
          <pre style="margin:10px 0 0;font-size:12.5px">Truth:      100 - 199
Pred:       103 - 199
<span class="c2">IoU = 97/100 = 0.97</span></pre>
        </div>
        <div class="card" style="margin:0">
          <b style="color:#fff">Lệch nhiều hơn</b>
          <pre style="margin:10px 0 0;font-size:12.5px">Truth:      100 - 199
Pred:       130 - 199
<span style="color:var(--red)">IoU = 70/100 = 0.70</span></pre>
        </div>
      </div>
      <div class="callout note">
        <div class="ttl">Ngưỡng coord_correct</div>
        <span style="color:var(--text-dim)">Một prediction là <code>coord_correct</code> nếu same-gene truth tồn tại <b>và</b> <code>IoU ≥ 0.90</code>. → IoU = 1.00: khớp hoàn toàn · IoU ≥ 0.90: đủ gần, coordinate-correct · IoU &lt; 0.90: failed về tọa độ.</span>
      </div>
      <div class="callout find">
        <div class="ttl">⚠️ Lưu ý về gene ngắn</div>
        <span style="color:var(--text-dim)">Với gene dài, lệch vài bp thường IoU vẫn cao. Nhưng với gene rất ngắn như FMD <code>2A</code>, chỉ lệch <b>6 bp</b> cũng có thể làm IoU tụt dưới 0.90. Một số failed cases ở gene ngắn không phải lift sai vùng lớn, mà là boundary precision bị strict scoring phạt mạnh.</span>
      </div>
      <p>Lý do tách <b>exact</b> và <b>coord_only</b>: <code>exact</code> đo boundary precision rất nghiêm ngặt; <code>coord_only</code> cho biết tool đã tìm đúng vùng gene, dù start/end lệch nhỏ hoặc annotation convention khác. Với virus annotation, nhiều case lệch vài bp có thể do convention, không nhất thiết là localization failure.</p>
    </section>

    <!-- 6 -->
    <section id="s6">
      <h2><span class="idx">06</span>Kết Quả Tổng Quan</h2>
      <p class="lead">Kết quả sau các cải thiện hiện tại.</p>
      <figure><img src="{{IMG_OVERALL}}" alt="Overall accuracy by virus"><figcaption><span class="fig-tag">Figure 1</span>Overall accuracy theo virus (FMD vs PRRSV).</figcaption></figure>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Virus</th><th class="r">Total cases</th><th class="r">Exact</th><th class="r">Coord only</th><th class="r">Failed</th><th class="r">Accuracy</th><th>Phân bố</th></tr></thead>
          <tbody>
            <tr><td class="gene">FMD</td><td class="r">1149</td><td class="r">1114</td><td class="r">16</td><td class="r">19</td><td class="r"><span class="acc high">98.35%</span></td><td><div class="barcell"><div class="bar"><span class="exact" style="width:96.95%"></span><span class="coord" style="left:96.95%;width:1.39%"></span><span class="fail" style="left:98.35%;width:1.65%"></span></div></div></td></tr>
            <tr><td class="gene">PRRSV</td><td class="r">761</td><td class="r">674</td><td class="r">84</td><td class="r">3</td><td class="r"><span class="acc perfect">99.61%</span></td><td><div class="barcell"><div class="bar"><span class="exact" style="width:88.57%"></span><span class="coord" style="left:88.57%;width:11.04%"></span><span class="fail" style="left:99.61%;width:0.39%"></span></div></div></td></tr>
            <tr class="total-row"><td class="gene">All</td><td class="r">1910</td><td class="r">1788</td><td class="r">100</td><td class="r">22</td><td class="r"><span class="acc perfect">98.85%</span></td><td><div class="barcell"><div class="bar"><span class="exact" style="width:93.61%"></span><span class="coord" style="left:93.61%;width:5.24%"></span><span class="fail" style="left:98.85%;width:1.15%"></span></div></div></td></tr>
          </tbody>
        </table>
      </div>
      <p style="color:var(--text-dim);font-size:13.5px"><span class="chip ok" style="font-size:12px"><span class="b"></span>exact</span> &nbsp; <span class="chip rescue" style="font-size:12px"><span class="b"></span>coord-only</span> &nbsp; <span class="chip invalid" style="font-size:12px"><span class="b"></span>failed</span></p>
      <ul>
        <li>FMD đạt <b>98.35%</b> khi tính exact + coordinate-correct.</li>
        <li>PRRSV đạt <b>99.61%</b> khi tính trên các gene thật sự có trong truth.</li>
        <li>Phần lớn lỗi còn lại là boundary issue nhỏ hoặc annotation/ref-truth mismatch.</li>
      </ul>
      <figure><img src="{{IMG_PER_GENE}}" alt="Per-gene accuracy"><figcaption><span class="fig-tag">Figure 2</span>Breakdown theo từng gene. Xanh đậm = exact, xanh nhạt = coord-only, đỏ = failed.</figcaption></figure>
    </section>

    <!-- 7 -->
    <section id="s7">
      <h2><span class="idx">07</span>Breakdown Theo Gene: FMD</h2>
      <p class="lead">FMD có 12 peptide (<code>mat_peptide</code>). Ta xem từng peptide ra sao ở baseline, khoanh vùng cái fail nhiều, rồi cải thiện.</p>
      <div class="phase-rail">
        <div class="pstep p1"><span class="pn">1</span>Từng gene ra sao</div><span class="arrow">→</span>
        <div class="pstep p2"><span class="pn">2</span>Tìm lỗi</div><span class="arrow">→</span>
        <div class="pstep p3"><span class="pn">3</span>Cải thiện</div><span class="arrow">→</span>
        <div class="pstep p4"><span class="pn">4</span>Kết quả sau</div>
      </div>

      <!-- FMD phase 1 -->
      <div class="pblock p1">
        <div class="phase p1"><span class="badge">1</span><div><span class="kicker">Bước 1 · Baseline</span><h3>Bức tranh từng peptide</h3></div></div>
        <p>Chạy tblastn với code gốc → exact <b>94.44%</b> (1088/1152). Đa số peptide đã gần như hoàn hảo; lỗi chỉ dồn vào vài chỗ:</p>
        <div class="gene-grid">
          <div class="gcard"><div class="gname">VP4</div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">VP2</div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">VP3</div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">2B</div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">2C</div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">3A</div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">3B</div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">3Cpro</div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">3Dpol</div><span class="gstat clean">Sạch</span></div>
          <div class="gcard hot"><div class="gname">Lpro</div><span class="gstat invest">Fail nhiều</span><div class="gnote">HSP cụt <b>12 bp</b> ở đầu N-terminal (14 ca lệch start)</div></div>
          <div class="gcard hot"><div class="gname">2A</div><span class="gstat invest">Fail nhiều</span><div class="gnote">peptide siêu ngắn → lệch vài bp đủ làm IoU &lt; 0.90</div></div>
          <div class="gcard"><div class="gname">VP1</div><span class="gstat minor">Lệch nhẹ</span><div class="gnote">coord đúng, boundary lệch theo convention</div></div>
        </div>
        <div class="grid-legend">
          <span><i style="background:var(--green)"></i>Sạch ~100%</span>
          <span><i style="background:var(--amber)"></i>Lệch boundary nhẹ</span>
          <span><i style="background:var(--red)"></i>Fail nhiều → cần đào sâu</span>
        </div>
        <p style="margin-top:14px"><b>Khoanh vùng:</b> 9/12 peptide đã sạch. Hai chỗ fail nhiều là <code>Lpro</code> (cụt N-terminal) và <code>2A</code> (peptide ngắn). <code>Lpro</code> là lỗi tool thật sự → đáng cải thiện.</p>
      </div>

      <!-- FMD phase 2 -->
      <div class="pblock p2">
        <div class="phase p2"><span class="badge">2</span><div><span class="kicker">Bước 2 · Finding</span><h3>tblastn align cụt vài amino acid ở terminal</h3></div></div>
        <p>Gom lỗi tool-side theo cơ chế (đã loại các case do ref/truth convention):</p>
        <div class="table-wrap">
          <table>
            <thead><tr><th>Cơ chế lỗi tool (baseline)</th><th class="r">Ca</th><th>Gene chính</th></tr></thead>
            <tbody>
              <tr><td class="gene">n_terminal_truncation_12bp</td><td class="r"><b>15</b></td><td>Lpro — thiếu 4 aa (12 bp) đầu N-terminal</td></tr>
              <tr><td class="gene">short_peptide_boundary_offset</td><td class="r"><b>9</b></td><td>peptide ngắn lệch boundary vài bp</td></tr>
              <tr><td class="gene">minor_vp1_boundary_offset</td><td class="r"><b>8</b></td><td>VP1 — lệch nhẹ (convention)</td></tr>
              <tr><td class="gene">minor_c_terminal_truncation_3bp</td><td class="r"><b>1</b></td><td>thiếu 1 aa ở C-terminal</td></tr>
            </tbody>
          </table>
        </div>
        <div class="callout find">
          <div class="ttl">💡 Vì sao không rescue bằng codon được?</div>
          <span style="color:var(--text-dim)"><code>mat_peptide</code> là cleavage product → <b>không có start/stop codon riêng</b>. tblastn align cụt vài residue đầu/cuối thì không thể tìm ATG/stop để rescue như CDS. Cần cơ chế khác.</span>
        </div>
      </div>

      <!-- FMD phase 3 -->
      <div class="pblock p3">
        <div class="phase p3"><span class="badge">3</span><div><span class="kicker">Bước 3 · Cải thiện</span><h3>Terminal extrapolation</h3></div></div>
        <p>Dùng <b>tọa độ query protein trong HSP</b> để biết phần amino acid bị cụt, rồi mở rộng boundary tới đúng đầu/cuối peptide. Chỉ áp dụng khi lượng thiếu nhỏ; case được mở rộng gắn status <code>ok_extrapolated</code>.</p>
      </div>

      <!-- FMD phase 4 -->
      <div class="pblock p4">
        <div class="phase p4"><span class="badge">4</span><div><span class="kicker">Bước 4 · Kết quả sau</span><h3>Lpro về 100%, exact 94.44% → 96.70%, 0 ca regress</h3></div></div>
        <div style="display:flex;gap:26px;flex-wrap:wrap;margin:6px 0 16px">
          <div><div style="font-size:12px;color:var(--text-faint);text-transform:uppercase;letter-spacing:.06em">Raw exact</div><div class="delta" style="font-size:18px"><span class="from">94.44%</span><span class="arrow">→</span><span class="to">96.70%</span></div></div>
          <div><div style="font-size:12px;color:var(--text-faint);text-transform:uppercase;letter-spacing:.06em">Fixed / regress</div><div class="delta" style="font-size:18px"><span class="to">+26</span> <span style="color:var(--text-faint)">/ 0</span></div></div>
          <div><div style="font-size:12px;color:var(--text-faint);text-transform:uppercase;letter-spacing:.06em">Cụt 12bp (Lpro)</div><div class="delta" style="font-size:18px"><span class="from">15</span><span class="arrow">→</span><span class="to">0</span></div></div>
          <div><div style="font-size:12px;color:var(--text-faint);text-transform:uppercase;letter-spacing:.06em">short_peptide</div><div class="delta" style="font-size:18px"><span class="from">9</span><span class="arrow">→</span><span class="to">0</span></div></div>
        </div>
        <figure><img src="{{IMG_FMD_EXTRAP}}" alt="FMD terminal extrapolation comparison"><figcaption><span class="fig-tag">Figure</span>FMD accuracy: baseline vs terminal extrapolation.</figcaption></figure>

        <h3>Bảng per-gene cuối</h3>
        <div class="table-wrap">
          <table>
            <thead><tr><th>Gene</th><th class="r">Total</th><th class="r">Exact</th><th class="r">Coord only</th><th class="r">Failed</th><th class="r">Accuracy</th></tr></thead>
            <tbody>
              <tr><td class="gene">Lpro</td><td class="r">96</td><td class="r">96</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">VP4</td><td class="r">95</td><td class="r">95</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">VP2</td><td class="r">96</td><td class="r">95</td><td class="r">0</td><td class="r">1</td><td class="r"><span class="acc high">98.96%</span></td></tr>
              <tr><td class="gene">VP3</td><td class="r">95</td><td class="r">95</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">VP1</td><td class="r">95</td><td class="r">79</td><td class="r">16</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">2A</td><td class="r">96</td><td class="r">79</td><td class="r">0</td><td class="r">17</td><td class="r"><span class="acc mid">82.29%</span></td></tr>
              <tr><td class="gene">2B</td><td class="r">96</td><td class="r">96</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">2C</td><td class="r">96</td><td class="r">96</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">3A</td><td class="r">96</td><td class="r">95</td><td class="r">0</td><td class="r">1</td><td class="r"><span class="acc high">98.96%</span></td></tr>
              <tr><td class="gene">3B</td><td class="r">96</td><td class="r">96</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">3Cpro</td><td class="r">96</td><td class="r">96</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">3Dpol</td><td class="r">96</td><td class="r">96</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
            </tbody>
          </table>
        </div>
        <div class="callout note">
          <div class="ttl">Lỗi còn lại — đều không phải truncation tool nữa</div>
          <ul style="margin-bottom:0">
            <li><code>VP1</code>: 16 coord-only — đúng vùng, boundary lệch theo convention.</li>
            <li><code>2A</code>: gene siêu ngắn, lệch vài bp = IoU &lt; 0.90 → bị strict scoring phạt, không phải lift sai vùng.</li>
            <li><code>VP2</code>, <code>3A</code>: mỗi gene 1 ca, cần manual review.</li>
          </ul>
        </div>
      </div>
    </section>

    <!-- 8 -->
    <section id="s8">
      <h2><span class="idx">08</span>Breakdown Theo Gene: PRRSV</h2>
      <p class="lead">PRRSV có 9 ORF (<code>CDS</code>), một số overlap nhau (<code>ORF2a/ORF2b</code>). Xem từng ORF ra sao, tìm ORF fail bất thường, rồi đào sâu.</p>
      <div class="phase-rail">
        <div class="pstep p1"><span class="pn">1</span>Từng gene ra sao</div><span class="arrow">→</span>
        <div class="pstep p2"><span class="pn">2</span>Tìm lỗi</div><span class="arrow">→</span>
        <div class="pstep p3"><span class="pn">3</span>Cải thiện</div><span class="arrow">→</span>
        <div class="pstep p4"><span class="pn">4</span>Kết quả sau</div>
      </div>

      <!-- PRRSV phase 1 -->
      <div class="pblock p1">
        <div class="phase p1"><span class="badge">1</span><div><span class="kicker">Bước 1 · Baseline</span><h3>Bức tranh từng ORF (exact / total)</h3></div></div>
        <p>Per-gene exact ở baseline (tính trên records có same-gene trong truth). Nhìn ngay ra <code>ORF7</code> lạc loài:</p>
        <div class="gene-grid">
          <div class="gcard"><div class="gname">ORF2a <span class="gnum">95/95</span></div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">ORF3 <span class="gnum">95/95</span></div><span class="gstat clean">Sạch</span></div>
          <div class="gcard"><div class="gname">ORF1a <span class="gnum">61/62</span></div><span class="gstat clean">Sạch</span><div class="gnote">1 ca do granularity ORF1ab</div></div>
          <div class="gcard"><div class="gname">ORF4 <span class="gnum">94/95</span></div><span class="gstat minor">Lệch nhẹ</span><div class="gnote">boundary offset vài bp</div></div>
          <div class="gcard"><div class="gname">ORF2b <span class="gnum">45/48</span></div><span class="gstat minor">Lệch nhẹ</span><div class="gnote">boundary offset</div></div>
          <div class="gcard"><div class="gname">ORF5 <span class="gnum">91/95</span></div><span class="gstat minor">Lệch nhẹ</span><div class="gnote">boundary offset</div></div>
          <div class="gcard"><div class="gname">ORF6 <span class="gnum">91/94</span></div><span class="gstat minor">Lệch nhẹ</span><div class="gnote">boundary offset</div></div>
          <div class="gcard"><div class="gname">ORF1b <span class="gnum">0/83</span></div><span class="gstat conv">Convention</span><div class="gnote">coord đúng 83/83; start lệch do frameshift → <b>không phải bug</b></div></div>
          <div class="gcard hot"><div class="gname">ORF7 <span class="gnum">66/94</span></div><span class="gstat invest">Fail nhiều</span><div class="gnote">tệ hơn hẳn các ORF khác → đào sâu</div></div>
        </div>
        <div class="grid-legend">
          <span><i style="background:var(--green)"></i>Sạch</span>
          <span><i style="background:var(--amber)"></i>Lệch boundary nhẹ</span>
          <span><i style="background:#4f9cff"></i>Convention (không phải bug)</span>
          <span><i style="background:var(--red)"></i>Fail nhiều → cần đào sâu</span>
        </div>
        <p style="margin-top:14px"><b>Khoanh vùng:</b> phần lớn ORF đã đúng hoặc chỉ lệch boundary nhẹ. <code>ORF1b</code> 0 exact nhưng là <b>convention</b> (coord đúng hết). Riêng <code>ORF7</code> fail tới 28/94 — bất thường → tách riêng investigate.</p>
      </div>

      <!-- PRRSV phase 2 -->
      <div class="pblock p2">
        <div class="phase p2"><span class="badge">2</span><div><span class="kicker">Bước 2 · Finding</span><h3>ORF7: start rescue chọn nhầm internal ATG</h3></div></div>
        <p>So delta tọa độ prediction vs truth cho riêng ORF7 → một pattern lặp lại rất rõ:</p>
        <div class="callout find">
          <div class="ttl">🔍 Pattern ORF7</div>
          <span style="color:var(--text-dim)"><code>delta_end = 0</code> nhưng <code>delta_start = +43 / +48</code> → tool tìm <b>đúng điểm kết thúc</b>, nhưng start lại nhảy vào một <b>ATG nội bộ</b> bên trong thay vì ATG upstream thật → CDS bị ngắn lại.</span>
        </div>
        <figure><img src="{{IMG_PRRSV_ORF7_DELTA}}" alt="ORF7 delta patterns"><figcaption><span class="fig-tag">Evidence</span>ORF7: end khớp (delta_end=0), start lệch dương đều đặn → rescue vào internal ATG.</figcaption></figure>
      </div>

      <!-- PRRSV phase 3 -->
      <div class="pblock p3">
        <div class="phase p3"><span class="badge">3</span><div><span class="kicker">Bước 3 · Cải thiện</span><h3>Frame &amp; ref-length-aware start rescue</h3></div></div>
        <p>Thay vì chọn ATG gần nhất một cách mù, start rescue giờ:</p>
        <ul style="margin-bottom:0">
          <li>kiểm tra frame: <code>len(CDS) % 3 == 0</code>;</li>
          <li>ưu tiên ATG tạo CDS <b>đúng frame</b> và có <b>độ dài gần reference</b> → bỏ qua các internal ATG tạo CDS quá ngắn.</li>
        </ul>
      </div>

      <!-- PRRSV phase 4 -->
      <div class="pblock p4">
        <div class="phase p4"><span class="badge">4</span><div><span class="kicker">Bước 4 · Kết quả sau</span><h3>ORF7 66 → 91; các ORF lệch nhẹ cũng về 100%</h3></div></div>
        <div style="display:flex;gap:22px;flex-wrap:wrap;margin:6px 0 16px">
          <div><div style="font-size:12px;color:var(--text-faint);text-transform:uppercase;letter-spacing:.06em">ORF7</div><div class="delta" style="font-size:18px"><span class="from">66</span><span class="arrow">→</span><span class="to">91</span></div></div>
          <div><div style="font-size:12px;color:var(--text-faint);text-transform:uppercase;letter-spacing:.06em">ORF5</div><div class="delta" style="font-size:18px"><span class="from">91</span><span class="arrow">→</span><span class="to">95</span></div></div>
          <div><div style="font-size:12px;color:var(--text-faint);text-transform:uppercase;letter-spacing:.06em">ORF6</div><div class="delta" style="font-size:18px"><span class="from">91</span><span class="arrow">→</span><span class="to">94</span></div></div>
          <div><div style="font-size:12px;color:var(--text-faint);text-transform:uppercase;letter-spacing:.06em">ORF2b</div><div class="delta" style="font-size:18px"><span class="from">45</span><span class="arrow">→</span><span class="to">48</span></div></div>
          <div><div style="font-size:12px;color:var(--text-faint);text-transform:uppercase;letter-spacing:.06em">ORF4</div><div class="delta" style="font-size:18px"><span class="from">94</span><span class="arrow">→</span><span class="to">95</span></div></div>
        </div>
        <figure><img src="{{IMG_PRRSV_RESCUE}}" alt="PRRSV start rescue exact comparison"><figcaption><span class="fig-tag">Figure</span>PRRSV exact per-gene: baseline vs start rescue.</figcaption></figure>

        <h3>Bảng per-gene cuối</h3>
        <div class="table-wrap">
          <table>
            <thead><tr><th>Gene</th><th class="r">Total</th><th class="r">Exact</th><th class="r">Coord only</th><th class="r">Failed</th><th class="r">Accuracy</th></tr></thead>
            <tbody>
              <tr><td class="gene">ORF1a</td><td class="r">62</td><td class="r">61</td><td class="r">1</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">ORF1b</td><td class="r">83</td><td class="r">0</td><td class="r">83</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">ORF2a</td><td class="r">95</td><td class="r">95</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">ORF2b</td><td class="r">48</td><td class="r">48</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">ORF3</td><td class="r">95</td><td class="r">95</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">ORF4</td><td class="r">95</td><td class="r">95</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">ORF5</td><td class="r">95</td><td class="r">95</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">ORF6</td><td class="r">94</td><td class="r">94</td><td class="r">0</td><td class="r">0</td><td class="r"><span class="acc perfect">100.00%</span></td></tr>
              <tr><td class="gene">ORF7</td><td class="r">94</td><td class="r">91</td><td class="r">0</td><td class="r">3</td><td class="r"><span class="acc high">96.81%</span></td></tr>
            </tbody>
          </table>
        </div>
        <div class="callout note">
          <div class="ttl">3 ca ORF7 còn lại — cố ý không auto-fix</div>
          <span style="color:var(--text-dim)">Cả 3 ca có <code>truth_len = 387 bp</code> trong khi reference ORF7 chỉ <code>372 bp</code>. Ref-length rescue cố ý <b>không kéo dài vượt reference</b>, nên để lại cho manual review / annotation convention — không nên auto-fix bằng độ dài ref.</span>
        </div>
        <figure><img src="{{IMG_ORF7}}" alt="ORF7 rescue accuracy comparison"><figcaption><span class="fig-tag">Evidence</span>ORF7-only experiment: fix đúng 25 ca ref/truth cùng 372 bp, còn lại 3 ca truth dài hơn.</figcaption></figure>
      </div>
    </section>

    <!-- 9 -->
    <section id="s9">
      <h2><span class="idx">09</span>Finding Chính</h2>
      <p class="lead">5 finding quan trọng rút ra từ validation.</p>
      <div class="findings">
        <div class="finding"><div class="fn">1</div><div class="fbody"><b>Accuracy phải tính trên gene thật sự có trong truth.</b><p>Nếu query truth không annotate gene đó, không thể dùng case đó để kết luận tool sai.</p></div></div>
        <div class="finding"><div class="fn">2</div><div class="fbody"><b>FMD và PRRSV cần xử lý khác nhau.</b><p>FMD dùng <code>mat_peptide</code> → phù hợp terminal extrapolation. PRRSV dùng <code>CDS</code> → phù hợp codon/frame-aware rescue.</p></div></div>
        <div class="finding"><div class="fn">3</div><div class="fbody"><b>tblastn thường tìm đúng vùng gene.</b><p>Lỗi chính không phải search sai vùng lớn — lỗi thường nằm ở boundary start/end.</p></div></div>
        <div class="finding"><div class="fn">4</div><div class="fbody"><b>ORF1b PRRSV cần interpretation riêng.</b><p>Coord đúng 83/83. Exact fail do frameshift/start-boundary convention.</p></div></div>
        <div class="finding"><div class="fn">5</div><div class="fbody"><b>Gene ngắn như FMD 2A dễ bị IoU phạt mạnh.</b><p>Lệch 6 bp ở gene rất ngắn có thể làm IoU dưới 0.90.</p></div></div>
      </div>
    </section>

    <!-- 10 -->
    <section id="s10">
      <h2><span class="idx">10</span>Kết Luận</h2>
      <p class="lead">ViraLift cho kết quả tốt trên cả FMD và PRRSV khi validation được tính đúng theo same-gene truth availability.</p>
      <div class="stat-band">
        <div class="stat"><div class="label">FMD</div><div class="value green">98.35%</div></div>
        <div class="stat"><div class="label">PRRSV</div><div class="value blue">99.61%</div></div>
        <div class="stat"><div class="label">Tổng 2 tập</div><div class="value purple">98.85%</div></div>
      </div>
      <h3>Kết luận kỹ thuật</h3>
      <ul>
        <li><code>tblastn</code> là hướng phù hợp cho annotation transfer khi query thiếu annotation.</li>
        <li>Các lỗi còn lại chủ yếu là boundary precision hoặc annotation convention mismatch.</li>
        <li>Terminal extrapolation cải thiện FMD.</li>
        <li>Frame/ref-length-aware start rescue cải thiện PRRSV, đặc biệt ORF7.</li>
      </ul>
      <div class="quote"><div class="body">ViraLift reliably transfers viral gene annotations using reference-guided tblastn lifting. On strict-clean FMDV and PRRSV validation datasets, the tool achieves high localization accuracy, with remaining failures mostly caused by short-feature boundary sensitivity or reference/query annotation convention differences.</div></div>
    </section>

    <footer>
      ViraLift — Tool &amp; Validation Report · Strict-clean FMDV &amp; PRRSV · Self-contained HTML
    </footer>
  </main>
</div>

<script>
  // Scrollspy: highlight active TOC item
  const links = Array.from(document.querySelectorAll('.toc a'));
  const map = new Map(links.map(a => [a.getAttribute('href').slice(1), a]));
  const obs = new IntersectionObserver((entries) => {
    entries.forEach(e => {
      if (e.isIntersecting) {
        links.forEach(l => l.classList.remove('active'));
        const a = map.get(e.target.id);
        if (a) a.classList.add('active');
      }
    });
  }, { rootMargin: '-20% 0px -70% 0px', threshold: 0 });
  document.querySelectorAll('section').forEach(s => obs.observe(s));
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
