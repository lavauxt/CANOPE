# ============================================================================
# Unified report assets: CSS + JavaScript + HTML assembly.
#
# Shared by the ECHO and CANOPE report generators so both packages emit
# byte-identical HTML structure, styling and client-side behaviour.  The only
# thing that differs between an ECHO and a CANOPE report is the JSON payload
# (the data blob embedded in <script id="unified-data">).
#
# Both echo_report.R and canope_report.R source() this file (or it is loaded
# once per package) and call report_html(payload, settings, out_file).
#
# Requires R >= 4.0 for raw strings (r"---(...)---").
# ============================================================================

# ----------------------------------------------------------------------------
# Stylesheet  (merged ECHO + CANOPE design system)
# ----------------------------------------------------------------------------

REPORT_CSS <- r"---(/* Unified interactive report -- design tokens, layout and components.
   Merges the ECHO and CANOPE stylesheets so both packages share one look.
   Categorical hues are applied by JS in fixed order; chrome (surface / ink /
   grid) swaps between light and dark.  Colour is never the only cue: every
   status colour is paired with an icon, symbol or text label. */
:root {
  color-scheme: light;
  --page: #f9f9f7;
  --surface: #fcfcfb;
  --text: #0b0b0b;
  --text-2: #52514e;
  --muted: #898781;
  --grid: #e1e0d9;
  --axis: #c3c2b7;
  --ring: rgba(11, 11, 11, 0.10);
  --accent: #2a78d6;
  --accent-ink: #ffffff;
  --accent-soft: #e7f0fc;
  --good: #006300;
  --warn: #8a5a00;
  --fail: #b3261e;
  --bad: #c42b2b;
  --ok-ink: #0b6b3a;
  --fail-soft: #fdeceb;
  --warn-soft: #fff3d6;
  --ok-soft: #e4f5ea;
  --info-soft: #e7f0fc;
  --bad-soft: #fdeaea;
  --c-del: #e34948; --c-dup: #2a78d6;
  --c-high: #1baf7a; --c-medium: #eda100; --c-low: #898781; --c-na: #c3c2b7;
  --radius: 10px;
  --font: system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
    --muted: #898781; --grid: #2c2c2a; --axis: #383835;
    --ring: rgba(255, 255, 255, 0.10); --accent: #3987e5; --accent-ink: #0b0b0b;
    --accent-soft: #1c2a3d; --good: #0ca30c; --warn: #fab219;
    --fail: #ff8a80; --bad: #f08080; --ok-ink: #5fd38d;
    --fail-soft: #3a1c1a; --warn-soft: #3a2d0b; --ok-soft: #12301f;
    --info-soft: #1c2a3d; --bad-soft: #3a1c1c;
    --c-del: #e66767; --c-dup: #3987e5; --c-high: #199e70; --c-medium: #c98500; --c-low: #898781; --c-na: #555552;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
  --muted: #898781; --grid: #2c2c2a; --axis: #383835;
  --ring: rgba(255, 255, 255, 0.10); --accent: #3987e5; --accent-ink: #0b0b0b;
  --accent-soft: #1c2a3d; --good: #0ca30c; --warn: #fab219;
  --fail: #ff8a80; --bad: #f08080; --ok-ink: #5fd38d;
  --fail-soft: #3a1c1a; --warn-soft: #3a2d0b; --ok-soft: #12301f;
  --info-soft: #1c2a3d; --bad-soft: #3a1c1c;
  --c-del: #e66767; --c-dup: #3987e5; --c-high: #199e70; --c-medium: #c98500; --c-low: #898781; --c-na: #555552;
}

* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body {
  margin: 0; background: var(--page); color: var(--text);
  font: 14px/1.5 var(--font); -webkit-font-smoothing: antialiased;
}
a { color: var(--accent); }
h1, h2, h3 { margin: 0; font-weight: 600; letter-spacing: -0.01em; }
h2 { font-size: 18px; } h3 { font-size: 14px; }
code, pre { font-family: var(--mono); font-size: 12px; }
pre { margin: 0; padding: 12px; overflow: auto; background: var(--page);
      border: 1px solid var(--ring); border-radius: 8px; max-height: 420px; }
.muted { color: var(--muted); } .small { font-size: 12px; }
.sr-only { position: absolute; width: 1px; height: 1px; overflow: hidden; clip: rect(0 0 0 0); }

/* ---- Header & tabs ---- */
.top { background: var(--surface); border-bottom: 1px solid var(--ring); }
.top-inner { max-width: 1440px; margin: 0 auto; padding: 14px 16px 0; }
.title-row { display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }
.title-row h1 { font-size: 20px; }
.badge { display: inline-flex; align-items: center; gap: 6px; padding: 2px 10px;
         border-radius: 999px; background: var(--accent-soft); color: var(--text-2);
         font-size: 12px; border: 1px solid var(--ring); }
.badge.demo { background: transparent; border-color: var(--warn); color: var(--warn); font-weight: 600; }
.badge.warn { background: var(--warn-soft); border-color: var(--warn); color: var(--warn); font-weight: 600; cursor: pointer; }
.badge.bad, .badge.fail { background: var(--fail-soft); border-color: var(--fail); color: var(--fail); font-weight: 600; cursor: pointer; }
.badge.ok { color: var(--good); }
.spacer { flex: 1; }
.tabs { display: flex; gap: 2px; overflow-x: auto; scrollbar-width: thin; margin-top: 10px; }
.tab { appearance: none; border: 0; background: none; color: var(--text-2); font: inherit;
       padding: 10px 14px; cursor: pointer; border-bottom: 2px solid transparent;
       white-space: nowrap; border-radius: 6px 6px 0 0; }
.tab:hover { background: var(--accent-soft); color: var(--text); }
.tab[aria-selected="true"] { color: var(--text); border-bottom-color: var(--accent); font-weight: 600; }
.tab:focus-visible, button:focus-visible, select:focus-visible, input:focus-visible {
  outline: 2px solid var(--accent); outline-offset: 2px; }
main { max-width: 1440px; margin: 0 auto; padding: 16px; }
.panel[hidden] { display: none; }

/* ---- Cards, grids, controls ---- */
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius);
        padding: 14px; min-width: 0; }
.card > header { display: flex; align-items: baseline; gap: 8px; margin-bottom: 8px; flex-wrap: wrap; }
.card > header .hint { margin-left: auto; }
.card.flagged { border-color: var(--warn); }
.grid { display: grid; gap: 12px; }
.g2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
.g3 { grid-template-columns: repeat(3, minmax(0, 1fr)); }
.g-main { grid-template-columns: minmax(0, 1fr) 360px; align-items: start; }
.g-viewer { grid-template-columns: 340px minmax(0, 1fr); align-items: start; }
.stack { display: flex; flex-direction: column; gap: 12px; }
.stack.tight { gap: 6px; }
.toolbar { display: flex; gap: 12px; flex-wrap: wrap; align-items: flex-end; margin-bottom: 12px; }
.toolbar:last-child { margin-bottom: 0; }
.field { display: flex; flex-direction: column; gap: 3px; font-size: 12px; color: var(--text-2); min-width: 0; max-width: 100%; }
.field.inline { flex-direction: row; align-items: center; gap: 6px; }
.field.inline input[type="checkbox"] { width: 16px; height: 16px; accent-color: var(--accent); }
select, input[type="text"], input[type="number"], input[type="search"] {
  font: inherit; color: var(--text); background: var(--surface); border: 1px solid var(--axis);
  border-radius: 6px; padding: 5px 8px; min-width: 120px; max-width: 280px; }
input[type="range"] { accent-color: var(--accent); }
button.btn { font: inherit; padding: 5px 12px; border-radius: 6px; cursor: pointer;
             border: 1px solid var(--axis); background: var(--surface); color: var(--text); }
button.btn:hover { background: var(--accent-soft); }
button.btn.primary { background: var(--accent); border-color: var(--accent); color: var(--accent-ink); }
button.btn[disabled] { opacity: 0.45; cursor: not-allowed; }
button.btn.small-btn { padding: 1px 9px; font-size: 12px; margin-left: 6px; }
.plot { width: 100%; min-height: 360px; }
.plot.tall { min-height: 560px; } .plot.short { min-height: 280px; }
.plot.strip { min-height: 70px; height: 70px; }
.plot.mid { min-height: 240px; }
.plot.viewer { height: 760px; }
.empty { padding: 28px 12px; text-align: center; color: var(--muted); }

/* ---- Stat tiles ---- */
.tiles { display: grid; gap: 12px; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); }
.tile { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius); padding: 12px 14px; }
.tile .v { font-size: 26px; font-weight: 600; letter-spacing: -0.02em; }
.tile .l { color: var(--text-2); font-size: 12px; }
.tile.ok { border-left: 4px solid var(--c-high); }
.tile.warn { border-left: 4px solid var(--warn); background: var(--warn-soft); }
.tile.fail, .tile.bad { border-left: 4px solid var(--fail); background: var(--fail-soft); }
.tile.warn .v { color: var(--warn); } .tile.fail .v, .tile.bad .v { color: var(--fail); } .tile.ok .v { color: var(--ok-ink); }
.scope { display: grid; gap: 12px; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); margin-bottom: 12px; }
.scope .tile .l { margin-top: 2px; }

/* ---- Legend (HTML so categories can be toggled) ---- */
.legend { display: flex; flex-wrap: wrap; gap: 4px 12px; margin-top: 8px; font-size: 12px; color: var(--text-2); }
.legend button { appearance: none; border: 0; background: none; font: inherit; color: inherit;
                 display: inline-flex; gap: 6px; align-items: center; cursor: pointer; padding: 2px 4px; border-radius: 4px; }
.legend button:hover { background: var(--accent-soft); }
.legend button[aria-pressed="false"] { opacity: 0.4; text-decoration: line-through; }
.swatch { width: 10px; height: 10px; border-radius: 3px; display: inline-block; }
.swatch.conf.HIGH { background: var(--c-high); } .swatch.conf.MEDIUM { background: var(--c-medium); }
.swatch.conf.LOW { background: var(--c-low); } .swatch.conf.NA { background: var(--c-na); }

/* ---- Tables ---- */
.tbl-wrap { overflow: auto; border: 1px solid var(--ring); border-radius: 8px; max-height: 520px; }
table.tbl { border-collapse: collapse; width: 100%; font-size: 12.5px; }
.tbl th, .tbl td { padding: 5px 10px; text-align: left; white-space: nowrap; border-bottom: 1px solid var(--grid); }
.tbl th { position: sticky; top: 0; background: var(--surface); cursor: pointer; user-select: none;
          font-weight: 600; color: var(--text-2); border-bottom: 1px solid var(--axis); z-index: 1; }
.tbl th[aria-sort="ascending"]::after { content: " \25B2"; font-size: 9px; }
.tbl th[aria-sort="descending"]::after { content: " \25BC"; font-size: 9px; }
.tbl td.num { text-align: right; font-variant-numeric: tabular-nums; }
.tbl td.wrap { white-space: normal; min-width: 260px; max-width: 520px; }
.tbl tbody tr:hover { background: var(--accent-soft); }
.tbl tbody tr.clickable, .tbl tbody tr.pick { cursor: pointer; }
.tbl tbody tr.sel { outline: 2px solid var(--accent); outline-offset: -2px; background: var(--accent-soft); }
.tbl td.warn { background: var(--warn-soft); color: var(--warn); font-weight: 600; }
.tbl td.fail { background: var(--fail-soft); color: var(--fail); font-weight: 600; }
.tbl-bar { display: flex; gap: 10px; align-items: center; margin: 8px 0; flex-wrap: wrap; }
.pager { display: inline-flex; gap: 6px; align-items: center; margin-left: auto; }

/* ---- Selection bars / kv / notes ---- */
.bars { display: flex; flex-direction: column; gap: 5px; }
.bar-row { display: grid; grid-template-columns: 110px 1fr 70px; gap: 8px; align-items: center; font-size: 12px; }
.bar-row .lab { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; color: var(--text-2); }
.bar-track { background: var(--grid); border-radius: 4px; height: 10px; overflow: hidden; }
.bar-fill { height: 100%; border-radius: 4px; }
.bar-fill.lv-ok { background: var(--c-high); }
.bar-fill.lv-warn { background: var(--c-medium); }
.bar-fill.lv-fail { background: var(--c-del); }
.note { font-size: 12px; color: var(--muted); }
.note-flag { font-size: 11.5px; }
.kv { display: grid; grid-template-columns: max-content 1fr; gap: 3px 14px; font-size: 12.5px; }
.kv dt { color: var(--text-2); } .kv dd { margin: 0; overflow-wrap: anywhere; }
.prose p { margin: 0 0 8px; max-width: 90ch; }

/* ---- Alerts (always an icon + text, never colour alone) ---- */
.alerts { display: flex; flex-direction: column; gap: 8px; margin-bottom: 12px; }
.alert { display: flex; gap: 10px; align-items: flex-start; padding: 9px 14px; border-radius: 8px;
         border: 1px solid var(--ring); border-left: 4px solid var(--muted); background: var(--surface); font-size: 13px; }
.alert .ico, .alert .ic { font-weight: 700; flex: none; width: 18px; text-align: center; line-height: 1.35; }
.alert .body { min-width: 0; overflow-wrap: anywhere; }
.alert .ttl { font-weight: 600; }
.alert.fail, .alert.error { border-left-color: var(--fail); background: var(--fail-soft); }
.alert.fail .ico, .alert.fail .ic, .alert.fail .ttl, .alert.error .ico, .alert.error .ic, .alert.error .ttl { color: var(--fail); }
.alert.warn, .alert.warning { border-left-color: var(--warn); background: var(--warn-soft); }
.alert.warn .ico, .alert.warn .ic, .alert.warn .ttl, .alert.warning .ico, .alert.warning .ic, .alert.warning .ttl { color: var(--warn); }
.alert.ok { border-left-color: var(--c-high); background: var(--ok-soft); }
.alert.ok .ico, .alert.ok .ic { color: var(--ok-ink); }
.alert.info { border-left-color: var(--accent); background: var(--info-soft); }
.alert.info .ico, .alert.info .ic, .alert.info .ttl { color: var(--accent); }
.alert a { cursor: pointer; }

/* ---- Collapsible alerts ---- */
details.alerts-details { margin-bottom: 12px; border-radius: var(--radius); }
details.alerts-details > summary {
  cursor: pointer; padding: 9px 14px; border-radius: var(--radius);
  background: var(--surface); border: 1px solid var(--ring);
  display: flex; gap: 8px; align-items: center; font-size: 13px; font-weight: 600;
  list-style: none; user-select: none;
}
details.alerts-details > summary::-webkit-details-marker { display: none; }
details.alerts-details > summary::before {
  content: "\25B8"; display: inline-block; transition: transform 0.15s;
  color: var(--muted); font-weight: 400;
}
details.alerts-details[open] > summary::before { transform: rotate(90deg); }
details.alerts-details > summary .ico { flex: none; width: 16px; text-align: center; }
details.alerts-details.fail   > summary { border-left: 4px solid var(--fail); }
details.alerts-details.warn   > summary { border-left: 4px solid var(--warn); }
details.alerts-details.error  > summary { border-left: 4px solid var(--fail); }
details.alerts-details.warning > summary { border-left: 4px solid var(--warn); }
details.alerts-details.ok     > summary { border-left: 4px solid var(--c-high); }
details.alerts-details.info   > summary { border-left: 4px solid var(--accent); }
details.alerts-details > .stack, details.alerts-details > .alerts { margin-top: 8px; }

/* ---- Chips ---- */
.chip { display: inline-flex; align-items: center; gap: 4px; padding: 0 8px; border-radius: 999px; font-size: 11.5px;
        font-weight: 600; line-height: 20px; border: 1px solid var(--ring); white-space: nowrap; background: var(--page); color: var(--text-2); }
.chip .ico { font-size: 11px; }
.chip.ok { color: var(--ok-ink); background: var(--ok-soft); }
.chip.warn, .chip.warning { color: var(--warn); background: var(--warn-soft); }
.chip.fail, .chip.error { color: var(--fail); background: var(--fail-soft); }
.chip.info { color: var(--accent); background: var(--info-soft); }
.chip.cat { font-weight: 500; }
.chip.type.DEL { color: var(--c-del); border-color: var(--c-del); background: transparent; }
.chip.type.DUP { color: var(--c-dup); border-color: var(--c-dup); background: transparent; }
.chip.DEL { color: var(--c-del); border-color: var(--c-del); }
.chip.DUP { color: var(--c-dup); border-color: var(--c-dup); }
.chip.conf.HIGH, .chip.HIGH { color: var(--ok-ink); background: var(--ok-soft); }
.chip.conf.MEDIUM, .chip.MEDIUM { color: var(--warn); background: var(--warn-soft); }
.chip.conf.LOW, .chip.LOW { color: var(--text-2); background: var(--page); }
.chip.conf.NA, .chip.NA { color: var(--text-2); background: var(--page); }
.chip-btn { appearance: none; font: inherit; font-size: 11.5px; font-weight: 600; display: inline-flex; gap: 4px; align-items: center;
            padding: 0 8px; line-height: 20px; border-radius: 999px; cursor: pointer; border: 1px solid var(--ring); background: var(--page); color: var(--text); }
.chip-btn:hover { background: var(--accent-soft); }
.chip-btn.warn { color: var(--warn); background: var(--warn-soft); }
.chip-btn.fail { color: var(--fail); background: var(--fail-soft); }
.chip-btn.link { color: var(--accent); }
.chips { display: flex; flex-wrap: wrap; gap: 6px; }
.count-badge { display: inline-block; margin-left: 6px; min-width: 18px; padding: 0 6px; border-radius: 999px; font-size: 11px; line-height: 18px;
               text-align: center; background: var(--warn); color: var(--page); font-weight: 700; }
.pill-list { display: flex; flex-wrap: wrap; gap: 6px; }

/* ---- Segmented control ---- */
.seg { display: inline-flex; border: 1px solid var(--axis); border-radius: 6px; overflow: hidden; }
.seg button { font: inherit; border: 0; background: var(--surface); color: var(--text-2); padding: 4px 10px; cursor: pointer;
              border-right: 1px solid var(--axis); }
.seg button:last-child { border-right: 0; }
.seg button[aria-pressed="true"] { background: var(--accent); color: var(--accent-ink); }

/* ---- Sub-tabs (inner switch within a tab) ---- */
.subtabs { display: flex; gap: 2px; border-bottom: 1px solid var(--ring); margin-bottom: 12px; overflow-x: auto; }
.subtab { appearance: none; border: 0; background: none; color: var(--text-2); font: inherit; font-size: 13px;
          padding: 8px 14px; cursor: pointer; border-bottom: 2px solid transparent; white-space: nowrap; border-radius: 6px 6px 0 0; }
.subtab:hover { background: var(--accent-soft); color: var(--text); }
.subtab[aria-selected="true"] { color: var(--text); border-bottom-color: var(--accent); font-weight: 600; }
.subpanel[hidden] { display: none; }

/* ---- Sample pickers ---- */
.sample-pick { display: flex; align-items: center; gap: 6px; font-size: 12px; color: var(--text-2); }
.sample-pick select { min-width: 200px; }
.sample-pick label { font-weight: 600; }
.sample-pick-inline { display: flex; align-items: center; gap: 8px; padding: 8px 12px; background: var(--surface);
                      border: 1px solid var(--ring); border-radius: 8px; margin-bottom: 12px; flex-wrap: wrap; }
.sample-pick-inline select { min-width: 200px; }
.sample-pick-inline label { font-size: 12px; color: var(--text-2); font-weight: 600; }

/* ---- Chromosome strip ---- */
.chrom-strip { display: grid; gap: 10px; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); }
.chrom-box { border: 1px solid var(--ring); border-left: 4px solid var(--muted); border-radius: 8px; padding: 8px 12px; display: flex; flex-direction: column; gap: 4px; align-items: flex-start; }
.chrom-box .nm { font-weight: 600; }
.chrom-box.ok { border-left-color: var(--c-high); }
.chrom-box.warn { border-left-color: var(--warn); background: var(--warn-soft); }
.chrom-box.fail { border-left-color: var(--fail); background: var(--fail-soft); }
.chrom-box.info { border-left-color: var(--accent); }

/* ---- Gallery (per-call plots) ---- */
.gallery { display: grid; gap: 12px; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); }
.shot { background: var(--surface); border: 1px solid var(--ring); border-radius: 8px; overflow: hidden; margin: 0; }
.shot button { all: unset; display: block; cursor: zoom-in; width: 100%; }
.shot img { display: block; width: 100%; height: 190px; object-fit: contain; background: #fff; }
.shot figcaption { padding: 6px 10px; font-size: 12px; color: var(--text-2); overflow-wrap: anywhere; }
.shot .nopreview { height: 190px; display: grid; place-items: center; color: var(--muted); background: var(--page); }
dialog.lightbox { border: 0; padding: 0; background: transparent; max-width: 96vw; max-height: 96vh; }
dialog.lightbox::backdrop { background: rgba(0, 0, 0, 0.7); }
dialog.lightbox img { max-width: 96vw; max-height: 92vh; background: #fff; border-radius: 8px; display: block; }

/* ---- Call viewer detail ---- */
.detail-head { display: flex; flex-wrap: wrap; gap: 8px 14px; align-items: baseline; margin-bottom: 8px; }
.detail-head h2 { font-size: 16px; }

@media (max-width: 1000px) { .g-main, .g-viewer, .g2, .g3 { grid-template-columns: minmax(0, 1fr); } }
@media (max-width: 720px) { main { padding: 12px; } .bar-row { grid-template-columns: 80px 1fr 56px; } }
@media print {
  .tabs, .toolbar, .btn, .subtabs, .sample-pick, .sample-pick-inline, .seg { display: none !important; }
  .panel[hidden] { display: block !important; }
  .subpanel[hidden] { display: block !important; }
  .card { break-inside: avoid; }
}
@media (forced-colors: active) { .bar-fill, .swatch { forced-color-adjust: none; } }
)---"


# ----------------------------------------------------------------------------
# Client-side application  (unified -- handles both ECHO and CANOPE payloads)
# ----------------------------------------------------------------------------

REPORT_JS <- r"---(/* Unified interactive report -- client-side application.
   Reads the JSON blob in #unified-data and builds 8 lazily-rendered tabs:
   Overview · QC & warnings · Samples · CNV calls · Cohort/PCA · Chromosomes ·
   Methods · Run.  The payload's meta.source field ("ECHO" or "CANOPE") selects
   which data-specific columns and panels are shown; all labels are identical
   between the two sources (combined labels such as "Q_SOME / BF" are used even
   when the source package only produces one of the two values).

   Architecture follows the SPARKS / ECHO / CANOPE report family: h() DOM
   helper, tableView, theme switch, hash navigation, collapsible alerts,
   sub-tabs, inline sample pickers.  No server and no network access are
   needed once plotly.js is loaded. */
(function () {
  'use strict';

  // ================================================================ data & utils
  var DATA = JSON.parse(document.getElementById('unified-data').textContent);
  var META = DATA.meta || {};
  var IS_ECHO = META.source === 'ECHO';
  var IS_CANOPE = META.source === 'CANOPE';
  var TH = DATA.thresholds || {};
  var SUM = DATA.summary || {};
  var S = DATA.samples || { cols: [], data: {}, n: 0 };
  var CALLS = DATA.calls || { cols: [], data: {}, n: 0 };
  var PLOTS = DATA.plots || DATA.windows || {};
  var WARNINGS = DATA.warnings || [];
  var ALERTS = DATA.alerts || [];
  var CHR = DATA.chromosomes || {};
  var SCOPE = DATA.scope || { cols: [], data: {}, n: 0 };
  var QC = DATA.qc || null;             // CANOPE QC metrics block
  var QC_FLAGS = DATA.qc_flags || null;  // ECHO QC flags table
  var CORR = DATA.corr || null;
  var PCA = DATA.pca || null;
  var DEPTH_HIST = DATA.depth_hist || null;
  var CONF_RULES = DATA.conf_rules || {};
  var PARAMS = DATA.params || {};
  var SCORING = DATA.scoring || null;
  var SETTINGS = DATA.settings || [];
  var LOG = DATA.log || {};

  var SAMPLES = (S.data.sample || S.data.Sample || []).slice();
  var SAMPLE_KEY = S.data.sample ? 'sample' : 'Sample';
  var SIDX = {}; SAMPLES.forEach(function (s, i) { SIDX[s] = i; });

  var $ = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };

  function h(tag, attrs) {
    var el = document.createElement(tag);
    if (attrs) Object.keys(attrs).forEach(function (k) {
      var v = attrs[k];
      if (v === null || v === undefined || v === false) return;
      if (k === 'class') el.className = v;
      else if (k === 'text') el.textContent = v;
      else if (k.slice(0, 2) === 'on' && typeof v === 'function') el.addEventListener(k.slice(2), v);
      else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
      else el.setAttribute(k, v === true ? '' : v);
    });
    for (var i = 2; i < arguments.length; i++) append(el, arguments[i]);
    return el;
  }
  function append(el, kid) {
    if (kid === null || kid === undefined || kid === false) return;
    if (Array.isArray(kid)) kid.forEach(function (k) { append(el, k); });
    else if (kid.nodeType) el.appendChild(kid);
    else el.appendChild(document.createTextNode(String(kid)));
  }
  function clear(el) { while (el.firstChild) el.removeChild(el.firstChild); return el; }
  function range(n) { var a = new Array(n); for (var i = 0; i < n; i++) a[i] = i; return a; }
  function uniq(a) { var seen = {}; return a.filter(function (v) { var k = String(v); if (seen[k]) return false; seen[k] = true; return true; }); }
  function asArray(x) { return x === null || x === undefined ? [] : Array.isArray(x) ? x : [x]; }
  function isNum(v) { return typeof v === 'number' && isFinite(v); }
  function fmtInt(n) { return Number(n).toLocaleString('en-US'); }
  function fmt(v) {
    if (v === null || v === undefined || (typeof v === 'number' && !isFinite(v))) return '';
    if (typeof v === 'boolean') return v ? 'yes' : 'no';
    if (typeof v !== 'number') return String(v);
    if (Number.isInteger(v)) return String(v);
    if (v !== 0 && Math.abs(v) < 1e-3) return v.toExponential(2);
    return String(Number(v.toPrecision(4)));
  }
  function orNa(v) { return (v === null || v === undefined || v === '' || (typeof v === 'number' && !isFinite(v))) ? 'n/a' : fmt(v); }
  function download(name, text, type) {
    var blob = new Blob([text], { type: type || 'text/csv;charset=utf-8' });
    var a = h('a', { href: URL.createObjectURL(blob), download: name });
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
  }
  function csvField(v) {
    var s = v === null || v === undefined ? '' : String(v);
    return /[",\n\r\t]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }
  function toCSV(cols, rowFn, n) {
    var lines = [cols.map(csvField).join(',')];
    for (var i = 0; i < n; i++) lines.push(cols.map(function (c) { return csvField(rowFn(c, i)); }).join(','));
    return lines.join('\n') + '\n';
  }
  function store(key, val) {
    try { if (val === undefined) return window.localStorage.getItem(key); window.localStorage.setItem(key, val); } catch (e) {}
    return null;
  }
  function plural(n, w, ws) { return n + ' ' + (n === 1 ? w : (ws || w + 's')); }
  function pretty(c) { return String(c).replace(/_/g, ' '); }

  // table helper: get a column array from a table object
  function tcol(T, c) { return (T && T.data && T.data[c]) || []; }
  function sval(s, col) { var i = SIDX[s]; return i === undefined ? null : tcol(S, col)[i]; }
  function hasCol(T, c) { return T && T.data && T.data[c] && T.data[c].length; }

  // ---- icons & levels ----
  var ICON = { fail: '\u2716', error: '\u2716', warn: '\u25b2', warning: '\u25b2', info: '\u2139', ok: '\u2714' };
  var LEVEL_LABEL = { fail: 'Fail', error: 'Fail', warn: 'Warning', warning: 'Warning', info: 'Note', ok: 'OK' };
  var SEV_ORDER = { error: 0, fail: 0, warning: 1, warn: 1, info: 2, ok: 3 };

  // ================================================================ theme & palette
  var CHROME = {
    light: { surface: '#fcfcfb', text: '#0b0b0b', secondary: '#52514e', muted: '#898781', grid: '#e1e0d9', axis: '#c3c2b7', mid: '#f0efec', band: '#f0efec' },
    dark: { surface: '#1a1a19', text: '#f0efec', secondary: '#c3c2b7', muted: '#898781', grid: '#2c2c2a', axis: '#383835', mid: '#383835', band: '#232321' }
  };
  var CAT = {
    light: ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'],
    dark: ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767']
  };
  var SEM = {
    light: { DEL: '#e34948', DUP: '#2a78d6', HIGH: '#006300', MEDIUM: '#eda100', LOW: '#898781', NA: '#c3c2b7',
             fail: '#b3261e', error: '#b3261e', warn: '#8a5a00', warning: '#8a5a00', ok: '#1baf7a', info: '#2a78d6',
             test: '#2a78d6', ref: '#a9a8a0', flank: '#898781', hit: '#e34948', del: '#e34948', dup: '#2a78d6',
             high: '#006300', medium: '#eda100', low: '#898781', bad: '#c42b2b' },
    dark: { DEL: '#e66767', DUP: '#3987e5', HIGH: '#0ca30c', MEDIUM: '#fab219', LOW: '#898781', NA: '#555552',
            fail: '#ff8a80', error: '#ff8a80', warn: '#fab219', warning: '#fab219', ok: '#199e70', info: '#3987e5',
            test: '#3987e5', ref: '#6f6e68', flank: '#898781', hit: '#e66767', del: '#e66767', dup: '#3987e5',
            high: '#0ca30c', medium: '#fab219', low: '#898781', bad: '#f08080' }
  };
  var SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'];
  var DIV = ['#e34948', '#f0a3a2', '#f0efec', '#9ec5f4', '#2a78d6'];
  var FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif";

  var themePref = store('unified-report-theme') || 'auto';
  function resolved() {
    if (themePref === 'light' || themePref === 'dark') return themePref;
    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  }
  function applyThemeAttr() {
    if (themePref === 'auto') document.documentElement.removeAttribute('data-theme');
    else document.documentElement.setAttribute('data-theme', themePref);
  }
  function tint(hex, amt, toward) {
    if (!amt) return hex;
    var n = function (s) { return [1, 3, 5].map(function (i) { return parseInt(s.slice(i, i + 2), 16); }); };
    var a = n(hex), b = n(toward);
    return '#' + a.map(function (v, i) { return ('0' + Math.round(v + (b[i] - v) * amt).toString(16)).slice(-2); }).join('');
  }
  function catColors(n, theme) {
    var base = CAT[theme], out = [], tints = [0, 0.38, 0.62];
    for (var i = 0; i < n; i++) out.push(tint(base[i % 8], tints[Math.min(Math.floor(i / 8), 2)], theme === 'light' ? '#ffffff' : '#000000'));
    return out;
  }
  function alpha(hex, a) {
    var v = [1, 3, 5].map(function (i) { return parseInt(hex.slice(i, i + 2), 16); });
    return 'rgba(' + v.join(',') + ',' + a + ')';
  }

  // ================================================================ plotly helpers
  var registry = [];
  var CONFIG = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['autoScale2d', 'lasso2d', 'select2d'] };

  function baseLayout(t, o) {
    o = o || {};
    var c = CHROME[t];
    var ax = function (label, extra) {
      return Object.assign({
        title: { text: label || '', font: { color: c.secondary } },
        gridcolor: c.grid, zerolinecolor: c.axis, linecolor: c.axis,
        tickfont: { color: c.muted }, automargin: true
      }, extra || {});
    };
    return {
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: FONT, size: 12, color: c.secondary },
      xaxis: ax(o.xlab, o.xaxis), yaxis: ax(o.ylab, o.yaxis),
      showlegend: !!o.legend,
      legend: Object.assign({ orientation: 'h', y: -0.2, font: { color: c.secondary } }, o.legendOpts || {}),
      margin: o.margin || { l: 56, r: 16, t: 16, b: 48 },
      hoverlabel: { font: { family: FONT } },
      hovermode: o.hovermode || 'closest',
      annotations: o.annotations, shapes: o.shapes, barmode: o.barmode, dragmode: o.dragmode, bargap: o.bargap
    };
  }
  function plot(div, build) {
    if (div._note) { clear(div); div._note = false; }
    div._render = function () {
      var s = build(resolved());
      var cfg = Object.assign({}, CONFIG, s.config || {});
      s.layout.autosize = true;
      window.Plotly.react(div, s.data, s.layout, cfg);
      if (s.onClick && div.on) {
        if (div.removeAllListeners) div.removeAllListeners('plotly_click');
        div.on('plotly_click', s.onClick);
      }
      div._dirty = false;
    };
    if (registry.indexOf(div) < 0) registry.push(div);
    div._render();
    return div;
  }
  function rerenderAll() {
    registry = registry.filter(function (d) { return document.body.contains(d) && d._render; });
    registry.forEach(function (d) { if (d.offsetParent !== null) d._render(); else d._dirty = true; });
  }
  function showNote(div, text) {
    if (window.Plotly && div.data) window.Plotly.purge(div);
    clear(div); div._render = null; div._note = true; div._dirty = false;
    div.appendChild(emptyNote(text));
  }
  function refreshVisible(root) {
    $$('.plot', root).forEach(function (d) {
      if (!d._render) return;
      if (d._dirty) d._render(); else if (window.Plotly && d.data) window.Plotly.Plots.resize(d);
    });
  }
  function plotBox(cls) { return h('div', { class: 'plot ' + (cls || '') }); }
  function emptyNote(text) { return h('div', { class: 'empty', text: text }); }

  // ================================================================ widgets
  function chip(level, text) {
    return h('span', { class: 'chip ' + level },
      h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[level] || '' }), text === undefined ? LEVEL_LABEL[level] : text);
  }
  function confChip(v) { return h('span', { class: 'chip conf ' + (v || 'NA') }, v === 'NA' || !v ? 'n/a' : v); }
  function typeChip(v) {
    if (v !== 'DEL' && v !== 'DUP') return fmt(v);
    return h('span', { class: 'chip type ' + v }, h('span', { class: 'ico', 'aria-hidden': 'true', text: v === 'DEL' ? '\u25BC' : '\u25B2' }), v);
  }
  function sevChip(v) { return v ? chip(v, v) : ''; }
  function card(title, hint, body, cls) {
    return h('section', { class: 'card' + (cls ? ' ' + cls : '') },
      h('header', null, h('h3', { text: title }), hint ? h('span', { class: 'muted small hint', text: hint }) : null), body);
  }
  function tile(value, label, level) {
    return h('div', { class: 'tile' + (level ? ' ' + level : '') }, h('div', { class: 'v', text: String(value) }), h('div', { class: 'l', text: label }));
  }
  function alertBox(level, title, text, extra) {
    return h('div', { class: 'alert ' + level, role: level === 'fail' || level === 'error' ? 'alert' : 'note' },
      h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[level] }),
      h('div', { class: 'body' },
        title ? h('span', { class: 'ttl', text: title }) : null, title && text ? ' ' : null, text || null, extra || null));
  }
  function kv(pairs) {
    var dl = h('dl', { class: 'kv' });
    pairs.forEach(function (p) {
      if (p[1] === null || p[1] === undefined || p[1] === '') return;
      dl.appendChild(h('dt', { text: p[0] })); dl.appendChild(h('dd', null, p[1]));
    });
    return dl;
  }
  function field(label, el) { return h('label', { class: 'field' }, label, el); }
  function selectOf(options, value, onchange, aria) {
    var sel = h('select', { 'aria-label': aria || 'Select' }, options.map(function (o) {
      return typeof o === 'string' ? h('option', { value: o, text: o }) : h('option', { value: o.value, text: o.text });
    }));
    if (value !== undefined) sel.value = value;
    if (onchange) sel.addEventListener('change', function () { onchange(sel.value); });
    return sel;
  }

  // ================================================================ tables
  function colType(T, c) {
    var arr = tcol(T, c);
    for (var i = 0; i < arr.length; i++) { if (arr[i] !== null && arr[i] !== undefined) return typeof arr[i]; }
    return 'string';
  }
  // T = {cols, data: {col: [..]}, n}. opts: filename, pageSize, render {col: fn(v,i)->node},
  // cellClass {col: fn(v,i)->class}, labels {col: header}, sort [col,'asc'|'desc'],
  // onRow fn(i), selected fn()->row index, wrap [cols], note
  function tableView(host, T, opts) {
    opts = opts || {};
    var n = T.n, pageSize = opts.pageSize || 25, page = 0, filter = '', pred = null;
    var cols = (T.cols || []).slice();
    var sortCol = opts.sort ? opts.sort[0] : null, sortDir = opts.sort ? opts.sort[1] : 'asc';
    var types = {}; cols.forEach(function (c) { types[c] = colType(T, c); });
    var idx = range(n);
    var search = h('input', { type: 'search', placeholder: 'Filter rows\u2026', 'aria-label': 'Filter rows' });
    var count = h('span', { class: 'muted small' });
    var pager = h('span', { class: 'pager' });
    var tbody = h('tbody'), thead = h('thead');
    var table = h('table', { class: 'tbl' }, thead, tbody);
    var wrap = h('div', { class: 'tbl-wrap' }, table);
    var dl = h('button', { class: 'btn', type: 'button', text: 'Download CSV', onclick: function () {
      download((opts.filename || 'table') + '.csv', toCSV(cols, function (c, k) { return tcol(T, c)[idx[k]]; }, idx.length));
    } });
    host.appendChild(h('div', { class: 'tbl-bar' }, search, count, dl, pager));
    host.appendChild(wrap);
    if (opts.note) host.appendChild(h('div', { class: 'note', text: opts.note }));

    function recompute() {
      var f = filter.toLowerCase();
      idx = [];
      for (var i = 0; i < n; i++) {
        if (pred && !pred(i)) continue;
        if (f) {
          var hit = false;
          for (var j = 0; j < cols.length && !hit; j++) {
            var v = tcol(T, cols[j])[i];
            if (v !== null && v !== undefined && String(v).toLowerCase().indexOf(f) >= 0) hit = true;
          }
          if (!hit) continue;
        }
        idx.push(i);
      }
      if (sortCol) {
        var arr = tcol(T, sortCol), dir = sortDir === 'asc' ? 1 : -1, num = types[sortCol] === 'number';
        idx.sort(function (a, b) {
          var x = arr[a], y = arr[b];
          if (x === null || x === undefined) return 1;
          if (y === null || y === undefined) return -1;
          return (num ? x - y : String(x).localeCompare(String(y), undefined, { numeric: true })) * dir;
        });
      }
      page = Math.min(page, Math.max(0, Math.ceil(idx.length / pageSize) - 1));
      draw();
    }
    function draw() {
      var tr = h('tr');
      cols.forEach(function (c) {
        tr.appendChild(h('th', {
          scope: 'col', tabindex: '0', text: (opts.labels && opts.labels[c]) || pretty(c),
          'aria-sort': sortCol === c ? (sortDir === 'asc' ? 'ascending' : 'descending') : 'none',
          onclick: function () { sortDir = sortCol === c && sortDir === 'asc' ? 'desc' : 'asc'; sortCol = c; recompute(); },
          onkeydown: function (e) { if (e.key === 'Enter') this.click(); }
        }));
      });
      clear(thead).appendChild(tr);
      clear(tbody);
      var sel = opts.selected ? opts.selected() : -1;
      var start = page * pageSize, end = Math.min(idx.length, start + pageSize);
      var wrapCols = opts.wrap || [];
      for (var k = start; k < end; k++) {
        var i = idx[k];
        var row = h('tr', { class: (opts.onRow ? 'clickable ' : '') + (i === sel ? 'sel' : '') });
        if (opts.onRow) (function (ii) {
          row.tabIndex = 0;
          row.addEventListener('click', function (e) { if (e.target.closest && e.target.closest('button,a')) return; opts.onRow(ii); });
          row.addEventListener('keydown', function (e) { if (e.key === 'Enter') opts.onRow(ii); });
        })(i);
        cols.forEach(function (c) {
          var v = tcol(T, c)[i], r = opts.render && opts.render[c];
          var cc = opts.cellClass && opts.cellClass[c] ? opts.cellClass[c](v, i) : '';
          var cls = ((types[c] === 'number' ? 'num ' : '') + (wrapCols.indexOf(c) >= 0 ? 'wrap ' : '') + (cc || '')).trim();
          row.appendChild(h('td', { class: cls || false }, r ? r(v, i) : fmt(v)));
        });
        tbody.appendChild(row);
      }
      count.textContent = fmtInt(idx.length) + (idx.length === n ? '' : ' of ' + fmtInt(n)) + ' rows' +
        (T.truncated ? ' (truncated from ' + fmtInt(T.n_total) + ')' : '');
      var pages = Math.max(1, Math.ceil(idx.length / pageSize));
      clear(pager);
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '\u2039', 'aria-label': 'Previous page', disabled: page === 0, onclick: function () { page--; draw(); } }));
      pager.appendChild(h('span', { class: 'small', text: (page + 1) + ' / ' + pages }));
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '\u203a', 'aria-label': 'Next page', disabled: page >= pages - 1, onclick: function () { page++; draw(); } }));
    }
    var timer = null;
    search.addEventListener('input', function () {
      clearTimeout(timer);
      timer = setTimeout(function () { filter = search.value; page = 0; recompute(); }, 120);
    });
    recompute();
    return { setPredicate: function (fn) { pred = fn; page = 0; recompute(); }, refresh: recompute, redraw: draw, count: function () { return idx.length; }, visible: function () { return idx.slice(); }, select: function (i) { draw(); } };
  }
  function simpleTable(T, opts) { var host = h('div'); tableView(host, T, opts); return host; }
  function subTable(T, cols) {
    var out = { cols: cols, data: {}, n: T.n };
    cols.forEach(function (c) { out.data[c] = tcol(T, c); });
    return out;
  }

  // ================================================================ global state (sample selection)
  var state = { sample: META.focus || '', call: null, listeners: [] };
  var headerPicker = null;
  var inlinePickers = [];
  function onState(cb) { state.listeners.push(cb); }
  function notify() {
    state.listeners.forEach(function (cb) { try { cb(); } catch (e) { if (window.console) console.error(e); } });
  }
  function setSample(s, quiet) {
    state.sample = s || '';
    if (headerPicker) headerPicker.value = state.sample;
    inlinePickers.forEach(function (p) { if (p.value !== state.sample) p.value = state.sample; });
    if (!quiet) notify();
  }
  // Build an inline sample picker (for QC and CNV calls tabs).
  // Renders after the collapsible alerts.  Synced with the header picker.
  function inlineSamplePicker(label) {
    if (!SAMPLES.length) return h('div', { class: 'sample-pick-inline', hidden: true });
    var sel = h('select', { 'aria-label': 'Sample' },
      h('option', { value: '', text: 'All samples (' + SAMPLES.length + ')' }),
      SAMPLES.map(function (s) { return h('option', { value: s, text: s }); }));
    sel.value = state.sample;
    sel.addEventListener('change', function () { setSample(sel.value); });
    inlinePickers.push(sel);
    return h('div', { class: 'sample-pick-inline' },
      h('label', { text: label || 'Sample' }), sel,
      h('span', { class: 'note', text: 'filters the tables and plots below to one sample' }));
  }
  function sampleSelect(label) {
    if (!headerPicker) return field(label || 'Sample', h('span', { class: 'muted small', text: 'Use the sample picker in the header or the inline picker above.' }));
    return field(label || 'Sample', h('span', { class: 'muted small', text: 'Use the sample picker above.' }));
  }

  // ================================================================ alerts (collapsible)
  // CANOPE-style warnings: {level, category, text, sample}
  function canopeAlerts(limit) {
    var fails = WARNINGS.filter(function (w) { return w.level === 'fail'; });
    var warns = WARNINGS.filter(function (w) { return w.level === 'warn'; });
    var infos = WARNINGS.filter(function (w) { return w.level === 'info'; });
    var worst = fails.length ? 'fail' : warns.length ? 'warn' : infos.length ? 'info' : 'ok';
    var summaryText = fails.length + ' fail \u00b7 ' + warns.length + ' warn \u00b7 ' + infos.length + ' note';
    if (worst === 'ok') {
      return h('div', { class: 'stack tight' }, alertBox('ok', 'No warnings.',
        'All sample metrics are within their thresholds and every chromosome class was analysed.'));
    }
    var box = h('div', { class: 'stack tight' });
    fails.slice(0, limit).forEach(function (w) { box.appendChild(alertBox('fail', w.category + ':', w.text, w.sample ? ' (' + w.sample + ')' : null)); });
    warns.slice(0, limit).forEach(function (w) { box.appendChild(alertBox('warn', w.category + ':', w.text, w.sample ? ' (' + w.sample + ')' : null)); });
    var shown = fails.length + warns.length;
    if (shown > limit * 2) box.appendChild(h('div', { class: 'note', text: (shown - limit * 2) + ' more warning(s) \u2014 open the sub-tab.' }));
    return h('details', { class: 'alerts-details ' + worst, open: worst === 'fail' ? '' : false },
      h('summary', null, h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[worst] }), h('span', { text: summaryText })),
      box);
  }
  // ECHO-style alerts: {kind, severity, title, text, tab, samples}
  function echoAlerts(kinds, minSev, withLinks) {
    var list = (ALERTS).filter(function (a) {
      if (kinds && kinds.indexOf(a.kind) < 0) return false;
      if (minSev && (SEV_ORDER[a.severity] || 3) > (SEV_ORDER[minSev] || 3)) return false;
      if (state.sample && a.samples && a.samples.length && a.samples.indexOf(state.sample) < 0) return false;
      return true;
    }).sort(function (a, b) { return (SEV_ORDER[a.severity] || 3) - (SEV_ORDER[b.severity] || 3); });
    if (!list.length) return null;
    var worst = list[0].severity;
    var items = list.map(function (a) {
      return h('div', { class: 'alert ' + a.severity, role: a.severity === 'info' ? 'note' : 'alert' },
        h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[a.severity] }),
        h('div', { class: 'body' },
          h('span', { class: 'ttl', text: a.title }), ' ', a.text,
          (withLinks !== false && a.tab && a.tab !== active)
            ? h('a', { href: '#tab=' + a.tab, text: ' Details \u2192', onclick: function (e) { e.preventDefault(); goto(a.tab); } })
            : null));
    });
    return h('details', { class: 'alerts-details ' + worst, open: worst === 'error' || worst === 'fail' ? '' : false },
      h('summary', null, h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[worst] }), h('span', { text: plural(list.length, 'warning') })),
      h('div', { class: 'alerts' }, items));
  }
  // Unified: show whichever alert source is present (both if both exist).
  function globalAlerts(limit) {
    var parts = [];
    var ea = echoAlerts(null, 'warning', true);
    if (ea) parts.push(ea);
    if (WARNINGS.length) parts.push(canopeAlerts(limit || 6));
    if (!parts.length) parts.push(h('div', { class: 'stack tight' }, alertBox('ok', 'No warnings.', 'All metrics within thresholds.')));
    return h('div', { class: 'stack tight' }, parts);
  }
  // Alerts filtered to a sample (for sample-specific tabs)
  function sampleAlerts() {
    if (IS_ECHO) return echoAlerts(['sample', 'cohort'], 'info', false);
    if (IS_CANOPE) {
      var s = state.sample, fl = s ? (DATA.flags && DATA.flags[s]) || [] : [];
      if (!fl.length) return h('div', { class: 'stack tight' }, alertBox('ok', null, 'No warnings for this sample.'));
      return h('div', { class: 'stack tight' }, fl.map(function (f) { return alertBox(f.level, null, f.text); }));
    }
    return null;
  }

  // ================================================================ sub-tabs
  // Builds an inner segmented control + panels.  items: [{id, label, build(root)}]
  function subTabs(items, defaultId) {
    var current = defaultId || items[0].id;
    var bar = h('div', { class: 'subtabs', role: 'tablist' });
    var panels = {};
    var container = h('div');
    items.forEach(function (it) {
      var btn = h('button', { class: 'subtab', role: 'tab', type: 'button', 'aria-selected': it.id === current ? 'true' : 'false',
        text: it.label, onclick: function () { switchTo(it.id); } });
      bar.appendChild(btn);
      var p = h('div', { class: 'subpanel', role: 'tabpanel', hidden: it.id === current ? false : true });
      panels[it.id] = { btn: btn, panel: p, built: false, build: it.build };
      container.appendChild(p);
    });
    function switchTo(id) {
      current = id;
      Object.keys(panels).forEach(function (k) {
        var on = k === id; panels[k].panel.hidden = !on;
        panels[k].btn.setAttribute('aria-selected', on ? 'true' : 'false');
      });
      if (!panels[id].built) { panels[id].built = true; try { panels[id].build(panels[id].panel); } catch (e) { panels[id].panel.appendChild(h('div', { class: 'empty', text: 'Error: ' + e.message })); if (window.console) console.error(e); } }
      else refreshVisible(panels[id].panel);
      setTimeout(function () { refreshVisible(panels[id].panel); }, 60);
    }
    container.insertBefore(bar, container.firstChild);
    // build the default panel immediately
    if (!panels[current].built) { panels[current].built = true; panels[current].build(panels[current].panel); }
    return { el: container, switchTo: switchTo, refresh: function () { if (panels[current].built) refreshVisible(panels[current].panel); } };
  }

  // ================================================================ call helpers
  var CALL_ROW = {};
  tcol(CALLS, 'id').forEach(function (id, i) { CALL_ROW[id] = i; });
  if (!hasCol(CALLS, 'id') && hasCol(CALLS, 'ID')) tcol(CALLS, 'ID').forEach(function (id, i) { CALL_ROW[id] = i; });

  function callsOf(s) {
    var out = [], sm = tcol(CALLS, 'sample');
    if (!sm.length) sm = tcol(CALLS, 'Sample');
    for (var i = 0; i < CALLS.n; i++) if (!s || sm[i] === s) out.push(i);
    return out;
  }
  var CONF_RANK = { HIGH: 3, MEDIUM: 2, LOW: 1 };
  function callOrder(rows) {
    var q = tcol(CALLS, 'q_some'), bf = tcol(CALLS, 'bf'), cf = tcol(CALLS, 'confidence');
    if (!cf.length) cf = tcol(CALLS, 'Confidence');
    return rows.slice().sort(function (a, b) {
      var d = (CONF_RANK[cf[b]] || 0) - (CONF_RANK[cf[a]] || 0);
      if (d) return d;
      var qa = isNum(q[a]) ? q[a] : isNum(bf[a]) ? bf[a] : -1;
      var qb = isNum(q[b]) ? q[b] : isNum(bf[b]) ? bf[b] : -1;
      return qb - qa;
    });
  }
  function hasPlot(i) {
    var hp = tcol(CALLS, 'has_plot'); if (hp.length) return !!hp[i];
    var v = tcol(CALLS, 'view'); return !!v[i];
  }
  function callId(i) { return tcol(CALLS, 'id')[i] !== undefined ? tcol(CALLS, 'id')[i] : tcol(CALLS, 'ID')[i]; }
  function callSample(i) { return tcol(CALLS, 'sample')[i] || tcol(CALLS, 'Sample')[i]; }
  function callType(i) { return tcol(CALLS, 'type')[i] || tcol(CALLS, 'Type')[i]; }
  function callConf(i) { return tcol(CALLS, 'confidence')[i] || tcol(CALLS, 'Confidence')[i]; }
  function callGene(i) { return tcol(CALLS, 'gene')[i] || tcol(CALLS, 'Gene')[i] || ''; }
  function callInterval(i) { return tcol(CALLS, 'interval')[i] || tcol(CALLS, 'Genomic')[i] || ''; }
  function goto(id) { activate(id, true); }
  function openCall(id) { state.call = id; var r = CALL_ROW[id]; if (r !== undefined) { var s = callSample(r); if (state.sample && state.sample !== s) setSample(s); } notify(); goto('calls'); }

  // ================================================================ TABS registry
  var TABS = [];
  function tab(id, label, available, build) {
    var t = { id: id, label: label, available: available, build: build, built: false };
    TABS.push(t); return t;
  }

  // ================================================================ OVERVIEW
  function levelCount(l) { return WARNINGS.filter(function (w) { return w.level === l; }).length; }
  function chromStrip() {
    var st = CHR.status, wrap = h('div', { class: 'chrom-strip' });
    if (!st || !st.n) return wrap;
    asArray(tcol(st, 'class')).forEach(function (k, i) {
      var inp = st.data.in_input ? st.data.in_input[i] : true, ana = st.data.analysed ? st.data.analysed[i] : true, reason = st.data.reason ? st.data.reason[i] : '';
      var level = ana ? 'ok' : inp ? (k === 'Y' ? 'warn' : 'fail') : 'info';
      wrap.appendChild(h('div', { class: 'chrom-box ' + level },
        h('div', { class: 'nm', text: st.data.label ? st.data.label[i] : k }),
        chip(level, ana ? 'analysed' : inp ? 'NOT analysed' : 'not in BED'),
        reason ? h('div', { class: 'small muted', text: reason }) : null));
    });
    return wrap;
  }
  function overviewBars() {
    var host = h('div'), box = plotBox('');
    var order = selectOf([{ value: 'input', text: 'Input order' }, { value: 'calls', text: 'Most calls first' }, { value: 'flags', text: 'Most flags first' }], 'calls', null, 'Sample order');
    host.appendChild(h('div', { class: 'toolbar' }, field('Order', order)));
    host.appendChild(box);
    if (!SAMPLES.length) { host.appendChild(emptyNote('No samples.')); return host; }
    function draw() {
      var nC = tcol(S, 'n_cnvs'), nC2 = hasCol(S, 'CNV calls') ? tcol(S, 'CNV calls') : nC;
      var nH = tcol(S, 'n_high'), nH2 = hasCol(S, 'HIGH') ? tcol(S, 'HIGH') : nH;
      var nM = tcol(S, 'n_medium'), nM2 = hasCol(S, 'MEDIUM') ? tcol(S, 'MEDIUM') : nM;
      var nL = tcol(S, 'n_low'), nL2 = hasCol(S, 'LOW') ? tcol(S, 'LOW') : nL;
      var ids = range(S.n);
      if (order.value === 'calls') ids.sort(function (a, b) { return (nC2[b] || 0) - (nC2[a] || 0); });
      plot(box, function (t) {
        var c = SEM[t], x = ids.map(function (i) { return SAMPLES[i]; });
        var mk = function (name, arr, col) { return { type: 'bar', name: name, x: x, y: ids.map(function (i) { return arr[i] || 0; }), marker: { color: col }, hovertemplate: '%{x}<br>' + name + ': %{y}<extra></extra>' }; };
        var data = [mk('HIGH', nH2.length ? nH2 : nH, c.HIGH || c.high), mk('MEDIUM', nM2.length ? nM2 : nM, c.MEDIUM || c.medium), mk('LOW', nL2.length ? nL2 : nL, c.LOW || c.low)];
        var L = baseLayout(t, { ylab: 'CNV calls', barmode: 'stack', legend: true, xaxis: { type: 'category', tickangle: -45 }, margin: { l: 56, r: 16, t: 16, b: 90 } });
        L.yaxis.rangemode = 'tozero';
        return { data: data, layout: L, onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) { setSample(p.x); goto('calls'); } } };
      });
    }
    order.addEventListener('change', draw);
    draw();
    return host;
  }
  function healthBars() {
    var host = h('div', { class: 'stack' }), bars = h('div', { class: 'bars' });
    var counts = { ok: 0, warn: 0, fail: 0 };
    var lvCol = tcol(S, 'level'), stCol = tcol(S, 'status');
    if (lvCol.length) lvCol.forEach(function (l) { if (counts[l] !== undefined) counts[l]++; });
    else if (stCol.length) stCol.forEach(function (l) { var m = { OK: 'ok', Warning: 'warn', Fail: 'fail', Excluded: 'warn' }[l] || 'ok'; counts[m] = (counts[m] || 0) + 1; });
    var total = Math.max(1, S.n);
    [['ok', 'No flags'], ['warn', 'Warnings'], ['fail', 'Fail']].forEach(function (p) {
      var col = SEM[resolved()][p[0]] || SEM[resolved()][p[0] === 'ok' ? 'ok' : p[0]];
      bars.appendChild(h('div', { class: 'bar-row' }, h('span', { class: 'lab', text: p[1] }),
        h('div', { class: 'bar-track' }, h('div', { class: 'bar-fill lv-' + p[0], style: { width: (100 * counts[p[0]] / total) + '%' } })),
        h('span', { class: 'small', text: counts[p[0]] + ' sample' + (counts[p[0]] === 1 ? '' : 's') })));
    });
    host.appendChild(bars);
    return host;
  }
  function runKv() {
    var p = PARAMS;
    return kv([
      ['Report generated', META.generated], ['Source', META.source], ['Version', META.version], ['Engine', META.engine_label || META.engine],
      ['Chromosome mode', META.mode_label || META.modechrom_label], ['Genome build', META.genome_build || META.genome_version],
      ['Coverage input', p.coverage_input || p.reads_file || null],
      ['Scope incomplete', META.scope_incomplete ? 'yes \u2013 some chromosome classes not analysed' : 'no']
    ]);
  }

  tab('overview', 'Overview', function () { return true; }, function (root) {
    var nF = levelCount('fail'), nW = levelCount('warn');
    var nErr = ALERTS.filter(function (a) { return a.severity === 'error' || a.severity === 'fail'; }).length;
    var nWarn = ALERTS.filter(function (a) { return a.severity === 'warning' || a.severity === 'warn'; }).length;
    var del = 0, dup = 0, high = 0;
    var typeCol = tcol(CALLS, 'type'); if (!typeCol.length) typeCol = tcol(CALLS, 'Type');
    typeCol.forEach(function (t) { if (t === 'DEL') del++; else if (t === 'DUP') dup++; });
    var confCol = tcol(CALLS, 'confidence'); if (!confCol.length) confCol = tcol(CALLS, 'Confidence');
    confCol.forEach(function (c) { if (c === 'HIGH') high++; });
    var both = SUM.both_chromosome_classes;
    var called = SUM.samples_called !== undefined ? SUM.samples_called : (META.n_analysed || 0);
    var requested = SUM.samples_requested !== undefined ? SUM.samples_requested : (META.n_samples || S.n);
    var tiles = h('div', { class: 'tiles' },
      tile(called + ' / ' + requested, 'samples analysed / requested', called < requested ? 'warn' : null),
      tile(fmtInt(CALLS.n), 'CNV calls (' + del + ' DEL, ' + dup + ' DUP)'),
      tile(fmtInt(high), 'HIGH-confidence calls'),
      tile(SUM.samples_with_calls || 0, 'samples with \u2265 1 call'),
      tile(fmtInt(META.n_targets || SUM.targets || 0), 'targets analysed'),
      tile(both === false ? 'Partial' : both === true ? 'Both' : (META.scope_incomplete ? 'Partial' : 'n/a'), 'autosomes + gonosomes', (both === false || META.scope_incomplete) ? 'fail' : both === true ? 'ok' : null),
      tile((nF + nW + nErr + nWarn), 'warnings (' + (nF + nErr) + ' fail, ' + (nW + nWarn) + ' warn)', (nF + nErr) ? 'fail' : (nW + nWarn) ? 'warn' : 'ok'));
    root.appendChild(h('div', { class: 'stack' },
      globalAlerts(6), tiles,
      card('Chromosome coverage of this run', META.mode_label || META.modechrom_label, h('div', null, chromStrip(),
        h('p', { class: 'note', text: 'One run analyses one chromosome mode (A = autosomes, XX = chrX with female references, XY = chrX/chrY with male references). Run the other modes to cover the whole genome.' }))),
      h('div', { class: 'grid g-main' }, card('Calls per sample', 'click a bar to list that sample\u2019s calls', overviewBars()), card('Sample health', null, healthBars())),
      card('Run', null, runKv())));
  });

  // ================================================================ QC & WARNINGS (sub-tabs)
  // QC_METRICS for the metric bars
  var QC_METRICS = [
    { col: 'median_depth', col2: 'Median depth', label: 'Median target depth', thr: function () { return TH.qc_min_cov || TH.min_cov; }, dir: 'below', hint: 'below the line = low coverage' },
    { col: 'total_reads', col2: 'Total reads', label: 'Total reads on the analysed targets', thr: function () { return TH.qc_min_total_reads || TH.min_total_reads; }, dir: 'below', hint: 'below the line = too few reads' },
    { col: 'max_corr', col2: 'Max correlation', label: 'Best Spearman correlation with another sample', thr: function () { return TH.qc_min_corr || TH.min_corr; }, dir: 'below', hint: 'below the line = no comparable reference' },
    { col: 'noise_z', col2: '', label: 'Sample noise (robust z of the exon-count dispersion)', thr: function () { return TH.qc_zscore; }, dir: 'above', hint: 'above the line = noisy sample' },
    { col: 'n_refs', col2: 'References', label: 'Reference samples used', thr: function () { return IS_CANOPE ? 5 : null; }, dir: 'below', hint: 'below the line = weak reference set' },
    { col: 'n_cnvs', col2: 'CNV calls', label: 'CNV calls', thr: function () { return null; }, dir: null, hint: 'outliers suggest a noisy sample / batch effect' }
  ];
  function qcTypeChip(v) {
    var lv = v === 'Missing Chromosome' ? 'fail' : v === 'Chromosome not analysed' ? 'warn' : v === 'Whole sample' ? 'warn' : 'info';
    return h('span', { class: 'chip ' + lv }, h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[lv] }), v);
  }
  function metricBars(m, orderSel) {
    var box = plotBox('short');
    var vals = tcol(S, m.col);
    if (!vals.length && m.col2) vals = tcol(S, m.col2);
    var ids = range(S.n).filter(function (i) { return isNum(vals[i]); });
    if (!ids.length) return null;
    function draw() {
      var o = orderSel.value, order = ids.slice();
      if (o === 'asc') order.sort(function (a, b) { return vals[a] - vals[b]; });
      else if (o === 'desc') order.sort(function (a, b) { return vals[b] - vals[a]; });
      plot(box, function (t) {
        var c = SEM[t], thr = m.thr();
        var bad = function (i) { return thr !== null && thr !== undefined && (m.dir === 'below' ? vals[i] < thr : vals[i] > thr); };
        var shapes = (thr !== null && thr !== undefined) ? [{ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: thr, y1: thr, line: { color: c.fail, dash: 'dash', width: 1.5 } }] : [];
        var d = [{ type: 'bar', x: order.map(function (i) { return SAMPLES[i]; }), y: order.map(function (i) { return vals[i]; }),
          marker: { color: order.map(function (i) { return bad(i) ? c.warn : CAT[t][0]; }) },
          hovertemplate: '%{x}<br>' + m.label + ': %{y}<extra></extra>' }];
        var lo = Math.min(0, Math.min.apply(null, ids.map(function (i) { return vals[i]; }))), hi = Math.max(Math.max.apply(null, ids.map(function (i) { return vals[i]; })), thr !== null && thr !== undefined ? thr : -Infinity);
        var L = baseLayout(t, { ylab: '', yaxis: { type: 'linear', range: [lo < 0 ? lo * 1.1 : 0, hi > 0 ? hi * 1.08 : 1] }, xaxis: { type: 'category', tickangle: -45, showticklabels: order.length <= 60 }, shapes: shapes, margin: { l: 60, r: 16, t: 10, b: order.length <= 60 ? 80 : 20 } });
        return { data: d, layout: L, onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) setSample(p.x); } };
      });
    }
    orderSel.addEventListener('change', draw);
    draw();
    var nbad = m.thr() === null || m.thr() === undefined ? 0 : ids.filter(function (i) { return m.dir === 'below' ? vals[i] < m.thr() : vals[i] > m.thr(); }).length;
    return card(m.label, (m.thr() !== null && m.thr() !== undefined ? 'threshold ' + fmt(m.thr()) + ' \u2013 ' : '') + m.hint + (nbad ? ' \u2013 ' + nbad + ' sample(s) flagged' : ''), box, nbad ? 'flagged' : '');
  }

  // ---- sub-tab: Alerts ----
  function buildAlertsSubtab(root) {
    var host = h('div', { class: 'stack tight' });
    if (IS_CANOPE) {
      var lvl = selectOf([{ value: '', text: 'All levels' }, { value: 'fail', text: 'Fail' }, { value: 'warn', text: 'Warning' }, { value: 'info', text: 'Note' }], '', null, 'Level');
      var cats = uniq(WARNINGS.map(function (w) { return w.category; }));
      var cat = selectOf([{ value: '', text: 'All categories' }].concat(cats.map(function (c) { return { value: c, text: c }; })), '', null, 'Category');
      var q = h('input', { type: 'search', placeholder: 'Search warnings\u2026', 'aria-label': 'Search warnings' });
      var count = h('span', { class: 'muted small' });
      var csvBtn = h('button', { class: 'btn', type: 'button', text: 'Download CSV', onclick: function () {
        download('warnings.csv', toCSV(['level', 'category', 'sample', 'text'], function (c, i) { return WARNINGS[i][c]; }, WARNINGS.length));
      } });
      root.appendChild(h('div', { class: 'toolbar' }, field('Level', lvl), field('Category', cat), field('Search', q), count, csvBtn));
      root.appendChild(host);
      function draw() {
        clear(host);
        var f = q.value.toLowerCase();
        var rows = WARNINGS.filter(function (w) {
          if (lvl.value && w.level !== lvl.value) return false;
          if (cat.value && w.category !== cat.value) return false;
          if (state.sample && w.sample !== state.sample) return false;
          if (f && (w.text + ' ' + w.category + ' ' + (w.sample || '')).toLowerCase().indexOf(f) < 0) return false;
          return true;
        });
        rows.sort(function (a, b) { return (SEV_ORDER[a.level] || 3) - (SEV_ORDER[b.level] || 3); });
        count.textContent = rows.length + ' of ' + WARNINGS.length;
        if (!rows.length) { host.appendChild(alertBox(WARNINGS.length ? 'info' : 'ok', WARNINGS.length ? 'No warning matches the filters.' : 'No warnings were raised.')); return; }
        [['fail', 'Fail', 'Fails'], ['warn', 'Warning', 'Warnings'], ['info', 'Note', 'Notes']].forEach(function (g) {
          var grp = rows.filter(function (w) { return w.level === g[0]; }).slice(0, 400);
          if (!grp.length) return;
          host.appendChild(h('details', { class: 'alerts-details ' + g[0], open: g[0] === rows[0].level ? '' : false },
            h('summary', null, h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[g[0]] }), h('span', { text: grp.length + ' ' + (grp.length === 1 ? g[1] : g[2]) })),
            h('div', { class: 'stack tight' }, grp.map(function (w) {
              return alertBox(w.level, null, null, [h('span', { class: 'chip cat', text: w.category }), ' ', w.sample ? h('button', { class: 'chip-btn link', type: 'button', onclick: function () { setSample(w.sample); goto('samples'); }, text: w.sample }) : null, w.sample ? ' ' : null, w.text]);
            }))));
        });
      }
      [lvl, cat].forEach(function (el) { el.addEventListener('change', draw); });
      q.addEventListener('input', draw);
      onState(draw);
      draw();
    } else {
      // ECHO-style alerts
      var sev = '', typ = '';
      var segBtns = ['', 'error', 'warning', 'info'].map(function (v) {
        return h('button', { type: 'button', 'aria-pressed': v === sev ? 'true' : 'false', text: v ? v[0].toUpperCase() + v.slice(1) : 'All', onclick: function () {
          sev = v; segBtns.forEach(function (b, k) { b.setAttribute('aria-pressed', ['', 'error', 'warning', 'info'][k] === v ? 'true' : 'false'); }); draw(); } });
      });
      var count = h('span', { class: 'muted small' });
      var csvBtn = h('button', { class: 'btn', type: 'button', text: 'Download CSV', onclick: function () {
        download('alerts.csv', toCSV(['severity', 'kind', 'title', 'text', 'samples'], function (c, i) {
          var a = ALERTS[i]; return c === 'samples' ? (a.samples || []).join(';') : a[c];
        }, ALERTS.length));
      } });
      root.appendChild(h('div', { class: 'toolbar' }, field('Severity', h('div', { class: 'seg', role: 'group' }, segBtns)), count, csvBtn));
      root.appendChild(host);
      function draw() {
        clear(host);
        var rows = ALERTS.filter(function (a) {
          if (sev && a.severity !== sev) return false;
          if (state.sample && a.samples && a.samples.length && a.samples.indexOf(state.sample) < 0) return false;
          return true;
        }).sort(function (a, b) { return (SEV_ORDER[a.severity] || 3) - (SEV_ORDER[b.severity] || 3); });
        count.textContent = rows.length + ' of ' + ALERTS.length;
        if (!rows.length) { host.appendChild(alertBox('ok', 'No alerts.')); return; }
        rows.forEach(function (a) {
          host.appendChild(alertBox(a.severity, a.title, a.text, a.tab ? h('a', { href: '#tab=' + a.tab, text: ' Details \u2192', onclick: function (e) { e.preventDefault(); goto(a.tab); } }) : null));
        });
      }
      onState(draw);
      draw();
    }
  }

  // ---- sub-tab: Metrics ----
  function buildMetricsSubtab(root) {
    var stack = h('div', { class: 'stack' });
    var orderSel = selectOf([{ value: 'input', text: 'Input order' }, { value: 'asc', text: 'Lowest first' }, { value: 'desc', text: 'Highest first' }], 'input', null, 'Order');
    stack.appendChild(h('div', { class: 'toolbar' }, field('Sort bars', orderSel),
      h('span', { class: 'note', text: 'Bars in amber breach the QC threshold (dashed line). Click a bar to select the sample.' })));
    var grid = h('div', { class: 'grid g2' });
    QC_METRICS.forEach(function (m) { var c = metricBars(m, orderSel); if (c) grid.appendChild(c); });
    stack.appendChild(grid);
    // scatter: depth vs correlation
    var depth = tcol(S, 'median_depth'); if (!depth.length) depth = tcol(S, 'Median depth');
    var corr = tcol(S, 'max_corr'); if (!corr.length) corr = tcol(S, 'Max correlation');
    var okIds = range(S.n).filter(function (i) { return isNum(depth[i]) && isNum(corr[i]); });
    if (okIds.length >= 2) {
      var sc = plotBox('short');
      stack.appendChild(card('Depth vs. best correlation', 'dashed lines = QC thresholds; click a point to select', sc));
      plot(sc, function (t) {
        var c = SEM[t], data = [];
        var lvCol = tcol(S, 'level'), stCol = tcol(S, 'status');
        ['ok', 'warn', 'fail'].forEach(function (lv) {
          var ids = okIds.filter(function (i) {
            if (lvCol.length) return lvCol[i] === lv;
            if (stCol.length) { var m = { OK: 'ok', Warning: 'warn', Fail: 'fail', Excluded: 'warn' }[stCol[i]] || 'ok'; return m === lv; }
            return lv === 'ok';
          });
          if (!ids.length) return;
          data.push({ type: 'scatter', mode: 'markers', name: LEVEL_LABEL[lv], x: ids.map(function (i) { return depth[i]; }), y: ids.map(function (i) { return corr[i]; }),
            text: ids.map(function (i) { return SAMPLES[i]; }), customdata: ids.map(function (i) { return SAMPLES[i]; }),
            marker: { color: c[lv], size: 8, symbol: lv === 'ok' ? 'circle' : lv === 'warn' ? 'diamond' : 'x', line: { width: 1, color: CHROME[t].surface } },
            hovertemplate: '<b>%{text}</b><br>depth %{x}<br>corr %{y}<extra>' + LEVEL_LABEL[lv] + '</extra>' });
        });
        var shapes = [];
        if (TH.qc_min_cov || TH.min_cov) shapes.push({ type: 'line', xref: 'x', yref: 'paper', x0: TH.qc_min_cov || TH.min_cov, x1: TH.qc_min_cov || TH.min_cov, y0: 0, y1: 1, line: { color: c.fail, dash: 'dash' } });
        if (TH.qc_min_corr || TH.min_corr) shapes.push({ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: TH.qc_min_corr || TH.min_corr, y1: TH.qc_min_corr || TH.min_corr, line: { color: c.fail, dash: 'dash' } });
        return { data: data, layout: baseLayout(t, { xlab: 'median depth', ylab: 'best correlation', legend: true, shapes: shapes }), onClick: function (ev) { var p = ev.points && ev.points[0]; if (p && p.customdata) setSample(p.customdata); } };
      });
    }
    root.appendChild(stack);
  }

  // ---- sub-tab: Flags ----
  function buildFlagsSubtab(root) {
    var stack = h('div', { class: 'stack' });
    // ECHO QC flags table
    if (QC_FLAGS && QC_FLAGS.n) {
      var host = h('div');
      var sev = '', typ = '';
      var tv = null;
      var bar = h('div', { class: 'toolbar' });
      var segBtns = ['', 'error', 'warning', 'info'].map(function (v) {
        return h('button', { type: 'button', 'aria-pressed': v === sev ? 'true' : 'false', text: v ? v[0].toUpperCase() + v.slice(1) : 'All', onclick: function () {
          sev = v; segBtns.forEach(function (b, k) { b.setAttribute('aria-pressed', ['', 'error', 'warning', 'info'][k] === v ? 'true' : 'false'); }); apply(); } });
      });
      var types = uniq(tcol(QC_FLAGS, 'Type')).sort();
      var typeSel = h('select', { 'aria-label': 'Flag type' }, h('option', { value: '', text: 'All types' }), types.map(function (x) { return h('option', { value: x, text: x }); }));
      typeSel.addEventListener('change', function () { typ = typeSel.value; apply(); });
      bar.appendChild(field('Severity', h('div', { class: 'seg', role: 'group' }, segBtns)));
      bar.appendChild(field('Type', typeSel));
      host.appendChild(bar);
      var tblHost = h('div'); host.appendChild(tblHost);
      tv = tableView(tblHost, QC_FLAGS, { filename: 'qc_flags', pageSize: 20, wrap: ['Details'], sort: ['Severity', 'asc'],
        render: { Severity: function (v) { return sevChip(v); }, Details: function (v) { return v; } } });
      function apply() {
        var bySample = function (i) { var s = state.sample; if (!s) return true; var v = tcol(QC_FLAGS, 'Sample')[i]; return v === s || v === 'All'; };
        tv.setPredicate(function (i) { return bySample(i) && (!sev || tcol(QC_FLAGS, 'Severity')[i] === sev) && (!typ || tcol(QC_FLAGS, 'Type')[i] === typ); });
      }
      apply();
      onState(apply);
      stack.appendChild(card('QC flags', 'Whole sample \u00b7 Whole exon \u00b7 Exon variability \u00b7 Cohort size \u00b7 Missing Chromosome', host));
    }
    // CANOPE QC metrics block
    if (QC && QC.rows && QC.rows.n) {
      var host2 = h('div'), typeBox = plotBox('short');
      var tcounts = tcol(QC.types, 'Freq'), tnames = tcol(QC.types, 'Type');
      stack.appendChild(h('div', { class: 'grid g-main' },
        card('QC flags by type', QC.n_total + ' flags in total', tnames.length ? typeBox : emptyNote('No QC flags.')),
        card('Genes with the most QC flags', 'top 60', QC.genes && QC.genes.n ? simpleTable(QC.genes, { filename: 'qc_genes', pageSize: 10 }) : emptyNote('No exon-level flags.'))));
      if (tnames.length) plot(typeBox, function (t) {
        return { data: [{ type: 'bar', orientation: 'h', y: tnames, x: tcounts, marker: { color: CAT[t][3] }, hovertemplate: '%{y}: %{x}<extra></extra>' }],
          layout: baseLayout(t, { xlab: 'flags', yaxis: { autorange: 'reversed' }, margin: { l: 170, r: 16, t: 10, b: 40 } }) };
      });
      var tv2 = tableView(host2, QC.rows, { filename: 'qc_metrics_flags', pageSize: 20, render: { Type: qcTypeChip },
        note: QC.capped ? 'Table capped at ' + fmtInt(QC.rows.n) + ' rows; full list in the QC TSV.' : null });
      var sm = tcol(QC.rows, 'Sample');
      var apply2 = function () { tv2.setPredicate(state.sample ? function (i) { return sm[i] === state.sample || sm[i] === 'All'; } : null); };
      onState(apply2); apply2();
      stack.appendChild(card('Exon-level QC flags', 'filtered by the selected sample', host2));
    }
    if (!stack.children.length) stack.appendChild(card('QC flags', null, emptyNote('No QC metrics file was provided.')));
    root.appendChild(stack);
  }

  // ---- sub-tab: Chromosome scope ----
  function buildScopeSubtab(root) {
    var stack = h('div', { class: 'stack' });
    // ECHO scope tiles
    if (SCOPE && SCOPE.n) {
      var scopeTiles = range(SCOPE.n).map(function (i) {
        var present = tcol(SCOPE, 'present')[i], analysed = tcol(SCOPE, 'analysed')[i];
        var status = !present ? 'not in BED' : analysed ? 'analysed' : 'NOT analysed';
        var cls = !present ? 'warn' : analysed ? 'ok' : 'fail';
        if (tcol(SCOPE, 'class')[i] === 'other') cls = '';
        return tile(tcol(SCOPE, 'label')[i], fmtInt(tcol(SCOPE, 'targets')[i]) + ' targets', cls, status);
      });
      stack.appendChild(card('Chromosome scope of this run', 'modechrom = ' + META.modechrom, h('div', null,
        h('p', { class: 'note', text: 'One run analyses one chromosome class: autosomes (A), chrX (X) or chrY (Y). Classes present in the BED but not analysed need their own run.' }),
        h('div', { class: 'scope' }, scopeTiles))));
    }
    // CANOPE chromosome status strip + inventory
    if (CHR.status && CHR.status.n) {
      stack.appendChild(card('Chromosome classes', null, chromStrip()));
      var inv = CHR.inventory;
      if (inv && inv.n) {
        var invHost = h('div');
        tableView(invHost, inv, { filename: 'chromosomes', pageSize: 30, labels: { targets_in_bed: 'targets in BED', targets_analysed: 'targets analysed' },
          render: { analysed: function (v) { return v ? chip('ok', 'analysed') : chip('warn', 'not analysed'); } },
          cellClass: { targets_analysed: function (v, i) { return tcol(inv, 'analysed')[i] ? '' : 'warn'; } } });
        var barBox = plotBox('short');
        stack.appendChild(h('div', { class: 'grid g2' }, card('Targets per chromosome', 'in BED vs. analysed', barBox), card('Chromosome inventory', null, invHost)));
        plot(barBox, function (t) {
          var c = tcol(inv, 'chromosome');
          return { data: [{ type: 'bar', name: 'in BED', x: c, y: tcol(inv, 'targets_in_bed'), marker: { color: SEM[t].NA } },
                          { type: 'bar', name: 'analysed', x: c, y: tcol(inv, 'targets_analysed'), marker: { color: CAT[t][0] } }],
            layout: baseLayout(t, { ylab: 'targets', barmode: 'overlay', legend: true, xaxis: { type: 'category', tickangle: -45 }, margin: { l: 56, r: 16, t: 10, b: 80 } }) };
        });
      }
    }
    if (!stack.children.length) stack.appendChild(card('Chromosome scope', null, emptyNote('No chromosome scope information available.')));
    root.appendChild(stack);
  }

  tab('qc', 'QC & warnings', function () { return true; }, function (root) {
    // Top: collapsible global alerts
    root.appendChild(globalAlerts(6));
    // Then: inline sample picker
    root.appendChild(inlineSamplePicker('Sample'));
    // Then: sub-tabs
    var st = subTabs([
      { id: 'alerts', label: 'Alerts', build: buildAlertsSubtab },
      { id: 'metrics', label: 'Metrics', build: buildMetricsSubtab },
      { id: 'flags', label: 'Flags', build: buildFlagsSubtab },
      { id: 'scope', label: 'Chromosome scope', build: buildScopeSubtab }
    ], 'alerts');
    root.appendChild(st.el);
    onState(function () { st.refresh(); });
  });

  // ================================================================ SAMPLES
  tab('samples', 'Samples', function () { return S.n > 0; }, function (root) {
    // Build a display table with both ECHO and CANOPE columns
    var dispCols = ['sample', 'status', 'level', 'n_cnvs', 'n_del', 'n_dup', 'n_high', 'n_medium', 'n_low',
                    'median_depth', 'total_reads', 'max_corr', 'n_refs', 'sex_inferred'];
    // Map ECHO column names
    var colMap = {};
    dispCols.forEach(function (c) {
      if (hasCol(S, c)) colMap[c] = c;
      else if (c === 'sample' && hasCol(S, 'Sample')) colMap[c] = 'Sample';
      else if (c === 'status' && hasCol(S, 'status')) colMap[c] = 'status';
      else if (c === 'level' && hasCol(S, 'Status')) colMap[c] = 'Status';
      else if (c === 'n_cnvs' && hasCol(S, 'CNV calls')) colMap[c] = 'CNV calls';
      else if (c === 'n_high' && hasCol(S, 'HIGH')) colMap[c] = 'HIGH';
      else if (c === 'n_medium' && hasCol(S, 'MEDIUM')) colMap[c] = 'MEDIUM';
      else if (c === 'n_low' && hasCol(S, 'LOW')) colMap[c] = 'LOW';
      else if (c === 'n_del' && hasCol(S, 'DEL')) colMap[c] = 'DEL';
      else if (c === 'n_dup' && hasCol(S, 'DUP')) colMap[c] = 'DUP';
      else if (c === 'median_depth' && hasCol(S, 'Median depth')) colMap[c] = 'Median depth';
      else if (c === 'total_reads' && hasCol(S, 'Total reads')) colMap[c] = 'Total reads';
      else if (c === 'max_corr' && hasCol(S, 'Max correlation')) colMap[c] = 'Max correlation';
      else if (c === 'n_refs' && hasCol(S, 'References')) colMap[c] = 'References';
    });
    var useCols = Object.keys(colMap);
    var T = { cols: useCols.map(function (c) { return colMap[c]; }), n: S.n, data: {} };
    useCols.forEach(function (c) { T.data[colMap[c]] = tcol(S, colMap[c]); });
    // add flags count column (CANOPE)
    if (IS_CANOPE) {
      var nf = SAMPLES.map(function (s) { var f = DATA.flags && DATA.flags[s]; return f ? f.filter(function (x) { return x.level === 'fail' || x.level === 'warn'; }).length : 0; });
      T.cols.push('flags'); T.data.flags = nf;
    }
    // Add sample-sheet attrs
    if (DATA.sample_attrs) Object.keys(DATA.sample_attrs).forEach(function (k) { T.cols.push(k); T.data[k] = asArray(DATA.sample_attrs[k]); });

    var tv = null, detail = h('div'), scatterBox = plotBox('');
    var labels = { n_cnvs: 'CNVs', n_del: 'DEL', n_dup: 'DUP', n_high: 'HIGH', n_medium: 'MEDIUM', n_low: 'LOW',
                   median_depth: 'median depth', total_reads: 'total reads', max_corr: 'best corr.', n_refs: 'refs', sex_inferred: 'sex (inferred)', flags: 'flags' };
    var cellClass = {};
    var render = {
      level: function (v) { return chip(v, LEVEL_LABEL[v] || v); },
      status: function (v) { return v === 'called' ? 'analysed' : v === 'Yes' ? 'analysed' : v === 'No' ? 'excluded' : v; },
      sample: function (v) { return h('strong', { text: v }); }
    };
    // numeric flag checks (CANOPE)
    if (IS_CANOPE) {
      var NF = {
        median_depth: function (v) { return isNum(v) && v < (TH.qc_min_cov || 0); },
        total_reads: function (v) { return isNum(v) && v < (TH.qc_min_total_reads || 0); },
        max_corr: function (v) { return isNum(v) && v < (TH.qc_min_corr || 0); },
        noise_z: function (v) { return isNum(v) && v > (TH.qc_zscore || 99); },
        n_refs: function (v) { return isNum(v) && v < 5; }
      };
      Object.keys(NF).forEach(function (k) { if (colMap[k]) cellClass[colMap[k]] = function (v) { return NF[k](v) ? 'warn' : ''; }; });
      cellClass.status = function (v) { return v === 'failed' || v === 'excluded' ? 'fail' : ''; };
    }
    var host = h('div');
    root.appendChild(h('div', { class: 'stack' },
      h('div', { class: 'toolbar' }, sampleSelect('Sample'),
        h('span', { class: 'note', text: 'Cells shaded amber/red are outside the QC thresholds.' })),
      card('Samples', 'click a row to select it', host),
      card('Selected sample', null, detail),
      card('Depth vs. best correlation', 'click a point to select the sample', scatterBox)));
    tv = tableView(host, T, { filename: 'samples', pageSize: 20, sort: ['flags', 'desc'], labels: labels, render: render, cellClass: cellClass,
      selected: function () { return state.sample ? SIDX[state.sample] : -1; },
      onRow: function (i) { setSample(SAMPLES[i] === state.sample ? '' : SAMPLES[i]); } });
    function drawDetail() {
      clear(detail);
      var s = state.sample;
      if (!s) { detail.appendChild(alertBox('ok', null, 'Select a sample to read its details and flags.')); return; }
      var lv = sval(s, colMap.level || 'level') || sval(s, 'Status');
      detail.appendChild(h('div', { class: 'toolbar' }, h('h3', { text: s }), lv ? chip(lv, LEVEL_LABEL[lv] || lv) : null,
        h('button', { class: 'btn', type: 'button', text: 'Show its CNV calls', onclick: function () { goto('calls'); } })));
      // flags
      if (IS_CANOPE) {
        var fl = (DATA.flags && DATA.flags[s]) || [];
        var list = h('div', { class: 'stack tight' });
        if (!fl.length) list.appendChild(alertBox('ok', null, 'No warnings for this sample.'));
        else fl.forEach(function (f) { list.appendChild(alertBox(f.level, null, f.text)); });
        detail.appendChild(list);
      }
      // key-value pairs
      var pairs = [];
      useCols.forEach(function (c) {
        if (c === 'sample' || c === 'level') return;
        var v = sval(s, colMap[c]); if (v === null || v === undefined || v === '') return;
        pairs.push([labels[c] || pretty(c), String(v)]);
      });
      if (pairs.length) detail.appendChild(h('div', { style: { marginTop: '10px' } }, kv(pairs)));
    }
    function drawScatter() {
      var depth = tcol(S, colMap.median_depth || 'median_depth'); if (!depth.length) depth = tcol(S, 'Median depth');
      var corr = tcol(S, colMap.max_corr || 'max_corr'); if (!corr.length) corr = tcol(S, 'Max correlation');
      var ok = range(S.n).filter(function (i) { return isNum(depth[i]) && isNum(corr[i]); });
      if (ok.length < 2) { showNote(scatterBox, 'Needs at least two analysed samples.'); return; }
      plot(scatterBox, function (t) {
        var c = SEM[t], data = [];
        ['ok', 'warn', 'fail'].forEach(function (lv) {
          var ids = ok.filter(function (i) {
            if (IS_CANOPE) return sval(SAMPLES[i], 'level') === lv;
            var st = tcol(S, 'Status')[i]; var m = { OK: 'ok', Warning: 'warn', Fail: 'fail', Excluded: 'warn' }[st] || 'ok'; return m === lv;
          });
          if (!ids.length) return;
          data.push({ type: 'scatter', mode: 'markers', name: LEVEL_LABEL[lv], x: ids.map(function (i) { return depth[i]; }), y: ids.map(function (i) { return corr[i]; }),
            text: ids.map(function (i) { return SAMPLES[i]; }), customdata: ids.map(function (i) { return SAMPLES[i]; }),
            marker: { color: c[lv], size: ids.map(function (i) { return 8 + Math.min(14, 2 * Math.sqrt((tcol(S, 'n_cnvs') || tcol(S, 'CNV calls'))[i] || 0)); }), symbol: lv === 'ok' ? 'circle' : lv === 'warn' ? 'diamond' : 'x', line: { width: 1, color: CHROME[t].surface } },
            hovertemplate: '<b>%{text}</b><br>depth %{x}<br>corr %{y}<extra>' + LEVEL_LABEL[lv] + '</extra>' });
        });
        return { data: data, layout: baseLayout(t, { xlab: 'median depth', ylab: 'best correlation', legend: true,
          shapes: [{ type: 'line', xref: 'x', yref: 'paper', x0: TH.qc_min_cov || TH.min_cov || 0, x1: TH.qc_min_cov || TH.min_cov || 0, y0: 0, y1: 1, line: { color: c.fail, dash: 'dash' } },
                   { type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: TH.qc_min_corr || TH.min_corr || 0, y1: TH.qc_min_corr || TH.min_corr || 0, line: { color: c.fail, dash: 'dash' } }] }),
          onClick: function (ev) { var p = ev.points && ev.points[0]; if (p && p.customdata) setSample(p.customdata); } };
      });
    }
    onState(function () { tv.redraw(); drawDetail(); drawScatter(); });
    drawDetail(); drawScatter();
  });

  // ================================================================ CNV CALLS (with folded 4-panel viewer)
  function chrNum(c) {
    var s = String(c).replace(/^chr/i, '').toUpperCase(), n = parseInt(s, 10);
    if (s === 'X') return 23; if (s === 'Y') return 24; if (s === 'M' || s === 'MT') return 25;
    return isFinite(n) ? n : 99;
  }
  function callFlags(i) {
    var out = [];
    if (hasCol(CALLS, 'homology_gene') && tcol(CALLS, 'homology_gene')[i]) out.push('homology gene');
    if (hasCol(CALLS, 'bg_flag') && tcol(CALLS, 'bg_flag')[i]) out.push('background');
    if (hasCol(CALLS, 'calibration') && tcol(CALLS, 'calibration')[i] === 'check') out.push('calibration');
    if (!hasPlot(i)) out.push('no plot');
    return out.join(', ');
  }
  // Combined-label cell accessors for the calls table
  function qsomeBfCell(i) {
    var q = tcol(CALLS, 'q_some')[i], bf = tcol(CALLS, 'bf')[i];
    if (IS_CANOPE) return orNa(q);
    if (IS_ECHO) return orNa(bf);
    return orNa(q !== undefined ? q : bf);
  }
  function mlcnCnCell(i) {
    var mlcn = tcol(CALLS, 'mlcn')[i], cn = tcol(CALLS, 'cn_label')[i];
    if (IS_CANOPE) return mlcn !== null && mlcn !== undefined ? (cn || '') + ' (MLCN ' + mlcn + ')' : orNa(mlcn);
    return '';
  }
  function corrRefsCell(i) {
    var corr = tcol(CALLS, 'correlation')[i], refs = tcol(CALLS, 'num_refs')[i], ncomp = tcol(CALLS, 'n_comp')[i];
    if (IS_CANOPE) return orNa(refs);
    if (IS_ECHO) return orNa(corr);
    return orNa(corr !== undefined ? corr : refs);
  }
  function exonsTargetsCell(i) {
    var nt = tcol(CALLS, 'num_targ')[i], ne = tcol(CALLS, 'nexons')[i];
    if (IS_CANOPE) return orNa(nt);
    if (IS_ECHO) return orNa(ne);
    return orNa(nt !== undefined ? nt : ne);
  }

  function callsTableDisplay() {
    // Build a unified display table with combined-label columns
    var cols = ['#', 'Sample', 'Chr', 'Start', 'End', 'Gene', 'Type', 'Exons / Targets', 'MLCN / CN', 'Q_SOME / BF', 'Corr / Refs', 'Confidence', 'Calibration', 'Plots'];
    var n = CALLS.n;
    var data = {};
    cols.forEach(function (c) { data[c] = new Array(n); });
    var idC = hasCol(CALLS, 'id') ? 'id' : 'ID';
    var smC = hasCol(CALLS, 'sample') ? 'sample' : 'Sample';
    var chrC = hasCol(CALLS, 'chr') ? 'chr' : 'Chr';
    var stC = hasCol(CALLS, 'start') ? 'start' : 'Start';
    var enC = hasCol(CALLS, 'end') ? 'end' : 'End';
    var geneC = hasCol(CALLS, 'gene') ? 'gene' : 'Gene';
    var typeC = hasCol(CALLS, 'type') ? 'type' : 'Type';
    var confC = hasCol(CALLS, 'confidence') ? 'confidence' : 'Confidence';
    for (var i = 0; i < n; i++) {
      data['#'][i] = tcol(CALLS, idC)[i];
      data['Sample'][i] = String(tcol(CALLS, smC)[i]);
      data['Chr'][i] = hasCol(CALLS, chrC) ? String(tcol(CALLS, chrC)[i]) : '';
      data['Start'][i] = hasCol(CALLS, stC) ? tcol(CALLS, stC)[i] : null;
      data['End'][i] = hasCol(CALLS, enC) ? tcol(CALLS, enC)[i] : null;
      data['Gene'][i] = String(callGene(i));
      data['Type'][i] = String(callType(i));
      data['Exons / Targets'][i] = exonsTargetsCell(i);
      data['MLCN / CN'][i] = mlcnCnCell(i);
      data['Q_SOME / BF'][i] = qsomeBfCell(i);
      data['Corr / Refs'][i] = corrRefsCell(i);
      data['Confidence'][i] = String(callConf(i));
      data['Calibration'][i] = hasCol(CALLS, 'calibration') ? tcol(CALLS, 'calibration')[i] : '';
      data['Plots'][i] = hasPlot(i);
    }
    return { cols: cols, n: n, data: data };
  }

  function genomeMap(box, rowsFn, yMode) {
    var ext = CHR.extent, chrs = ext ? asArray(ext.chr) : [], lens = ext ? asArray(ext.len) : [], off = {}, tot = 0, ticks = [], tt = [], shapes = [];
    // If no CANOPE extent, fall back to ECHO chr table
    if (!chrs.length && SCOPE && SCOPE.n) {
      // build extent from scope first/last
      var f = tcol(SCOPE, 'first'), l = tcol(SCOPE, 'last');
      chrs = range(SCOPE.n).map(function (i) { return tcol(SCOPE, 'chromosome')[i]; });
      lens = range(SCOPE.n).map(function (i) { return (l[i] || 0) - (f[i] || 0) + 1; });
    }
    chrs.forEach(function (c, i) {
      off[c] = tot;
      ticks.push(tot + lens[i] / 2); tt.push(String(c).replace(/^chr/, ''));
      if (i % 2 === 0) shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: tot, x1: tot + lens[i], y0: 0, y1: 1, line: { width: 0 }, fillcolor: 'rgba(128,128,128,0.08)', layer: 'below' });
      tot += lens[i] + 1;
    });
    plot(box, function (t) {
      var c = SEM[t], rows = rowsFn(), data = [];
      var chrCol = hasCol(CALLS, 'chr') ? 'chr' : 'Chr';
      var stCol = hasCol(CALLS, 'start') ? 'start' : 'Start';
      var enCol = hasCol(CALLS, 'end') ? 'end' : 'End';
      var gidxCol = hasCol(CALLS, 'gidx') ? 'gidx' : null;
      ['DEL', 'DUP'].forEach(function (type) {
        var ids = rows.filter(function (i) { return callType(i) === type; });
        if (!ids.length) return;
        var xs = ids.map(function (i) {
          if (gidxCol && isNum(tcol(CALLS, gidxCol)[i])) return tcol(CALLS, gidxCol)[i];
          var chr = String(tcol(CALLS, chrCol)[i]); if (off[chr] === undefined) return null;
          return off[chr] + (tcol(CALLS, stCol)[i] + tcol(CALLS, enCol)[i]) / 2;
        });
        var yv = ids.map(function (i) {
          if (yMode === 'mlcn') return tcol(CALLS, 'mlcn')[i];
          if (yMode === 'q_some') { var q = tcol(CALLS, 'q_some')[i]; return isNum(q) ? q : null; }
          return callSample(i);
        });
        data.push({ type: 'scatter', mode: 'markers', name: type, x: xs, y: yv,
          customdata: ids.map(function (i) { return callId(i); }),
          text: ids.map(function (i) { return '<b>' + callSample(i) + '</b> ' + type + ' ' + callInterval(i) + '<br>' + callGene(i) + '<br>Q_SOME/BF ' + qsomeBfCell(i) + ' \u00b7 ' + callConf(i); }),
          marker: { symbol: type === 'DEL' ? 'triangle-down' : 'triangle-up', color: c[type], size: 8, opacity: ids.map(function (i) { return { HIGH: 1, MEDIUM: 0.75, LOW: 0.45 }[callConf(i)] || 0.45; }), line: { width: 1, color: CHROME[t].surface } },
          hovertemplate: '%{text}<extra></extra>' });
      });
      var yaxis = yMode === 'sample' ? { type: 'category', autorange: 'reversed' } : { rangemode: 'tozero' };
      return { data: data, layout: baseLayout(t, { legend: true, ylab: yMode === 'sample' ? '' : yMode === 'mlcn' ? 'MLCN' : 'Q_SOME / BF', yaxis: yaxis, shapes: shapes,
        xaxis: { tickmode: 'array', tickvals: ticks, ticktext: tt, showgrid: false, title: { text: 'chromosome' } },
        margin: { l: yMode === 'sample' ? 90 : 56, r: 16, t: 10, b: 60 } }),
        onClick: function (ev) { var p = ev.points && ev.points[0]; if (p && p.customdata !== undefined) openCall(p.customdata); } };
    });
  }

  // ---- 4-panel call viewer ----
  function viewerFigure(P, row, t, showRefs) {
    var c = SEM[t], D = CALLS.data, type = callType(row), tc = c[type] || c.info || c.test;
    var px = asArray(P.px), aff = asArray(P.affected), n = px.length;
    var refLog = P.ref_log ? asArray(P.ref_log).map(asArray) : (P.ref_cov ? asArray(P.ref_cov).map(asArray) : []);
    var refZ = P.ref_z ? asArray(P.ref_z).map(asArray) : [];
    var testLog = asArray(P.test_log || P.test), ratio = asArray(P.ratio), lo = asArray(P.lo), hi = asArray(P.hi), z = asArray(P.z);
    var gene = asArray(P.gene), xlab = asArray(P.xlab), group = asArray(P.gene_group || P.group);
    var off = P.off ? asArray(P.off) : null;
    var txt = range(n).map(function (k) {
      var lab = xlab[k] !== undefined ? xlab[k] : (P.exon ? P.exon[k] : 'exon ' + (k + 1));
      return '<b>' + (gene[k] || '') + '</b> ' + (off && off[k] ? lab : 'exon ' + lab) + (aff[k] ? ' (called)' : '') + '<br>' + (P.pos ? P.pos[k] : '');
    });
    var markerCol = range(n).map(function (k) { return aff[k] ? tc : c.test; });
    var markerSym = range(n).map(function (k) { return aff[k] ? (type === 'DEL' ? 'triangle-down' : 'triangle-up') : 'circle'; });
    var markerSize = range(n).map(function (k) { return aff[k] ? 11 : 7; });
    var data = [];
    // Panel 0: gene strip (drawn separately, not here)
    // Panel 1: coverage
    if (showRefs && refLog.length) {
      var rlX = [], rlY = [];
      refLog.forEach(function (m) { for (var k = 0; k < n; k++) { rlX.push(px[k]); rlY.push(m[k]); } rlX.push(null); rlY.push(null); });
      data.push({ type: 'scatter', mode: 'markers', name: 'References (' + refLog.length + ')', x: rlX, y: rlY, xaxis: 'x', yaxis: 'y', marker: { color: c.ref, size: 4, opacity: 0.75 }, hoverinfo: 'skip', legendgroup: 'refs' });
    }
    data.push({ type: 'scatter', mode: 'markers', name: 'Test sample', x: px, y: testLog, xaxis: 'x', yaxis: 'y', text: txt,
      marker: { color: markerCol, symbol: markerSym, size: markerSize, line: { width: 1, color: CHROME[t].surface } },
      hovertemplate: '%{text}<br>log2(count) = %{y:.2f}<extra></extra>' });
    // Panel 2: observed/expected with interval
    var errUp = range(n).map(function (k) { return isNum(hi[k]) && isNum(ratio[k]) ? Math.max(0, hi[k] - ratio[k]) : 0; });
    var errDn = range(n).map(function (k) { return isNum(lo[k]) && isNum(ratio[k]) ? Math.max(0, ratio[k] - lo[k]) : 0; });
    data.push({ type: 'scatter', mode: 'markers', name: 'obs/exp', x: px, y: ratio, xaxis: 'x', yaxis: 'y2', showlegend: false, text: txt,
      error_y: { type: 'data', symmetric: false, array: errUp, arrayminus: errDn, color: c.flank, thickness: 1, width: 3 },
      marker: { color: markerCol, symbol: markerSym, size: markerSize, line: { width: 1, color: CHROME[t].surface } },
      hovertemplate: '%{text}<br>obs/exp = %{y:.2f}<extra></extra>' });
    // Panel 3: z-score
    if (showRefs && refZ.length) {
      var rzX = [], rzY = [];
      refZ.forEach(function (m) { for (var k = 0; k < n; k++) { rzX.push(px[k]); rzY.push(m[k]); } rzX.push(null); rzY.push(null); });
      data.push({ type: 'scatter', mode: 'markers', name: 'refs z', x: rzX, y: rzY, xaxis: 'x', yaxis: 'y3', marker: { color: c.ref, size: 4, opacity: 0.75 }, hoverinfo: 'skip', legendgroup: 'refs', showlegend: false });
    }
    data.push({ type: 'scatter', mode: 'markers', name: 'z', x: px, y: z, xaxis: 'x', yaxis: 'y3', showlegend: false, text: txt,
      marker: { color: markerCol, symbol: markerSym, size: markerSize, line: { width: 1, color: CHROME[t].surface } },
      hovertemplate: '%{text}<br>z = %{y:.2f}<extra></extra>' });
    data.push({ type: 'scatter', mode: 'markers', name: 'Called (' + type + ')', x: [null], y: [null], xaxis: 'x', yaxis: 'y',
      marker: { color: tc, symbol: type === 'DEL' ? 'triangle-down' : 'triangle-up', size: 11 } });

    var shapes = [];
    var affIdx = range(n).filter(function (k) { return aff[k]; });
    if (affIdx.length) {
      var x0 = Math.min.apply(null, affIdx.map(function (k) { return px[k]; })) - 0.5, x1 = Math.max.apply(null, affIdx.map(function (k) { return px[k]; })) + 0.5;
      shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: x0, x1: x1, y0: 0, y1: 1, fillcolor: alpha(tc, 0.1), line: { width: 0 }, layer: 'below' });
    }
    shapes.push({ type: 'line', xref: 'paper', yref: 'y2', x0: 0, x1: 1, y0: 0, y1: 0, line: { color: CHROME[t].axis, width: 1 } });
    [3, -3].forEach(function (v) { shapes.push({ type: 'line', xref: 'paper', yref: 'y3', x0: 0, x1: 1, y0: v, y1: v, line: { color: c.fail, dash: 'dot', width: 1 } }); });
    shapes.push({ type: 'line', xref: 'paper', yref: 'y3', x0: 0, x1: 1, y0: 0, y1: 0, line: { color: CHROME[t].axis, width: 1 } });

    var L = baseLayout(t, { legend: true, legendOpts: { y: 1.1, x: 0, xanchor: 'left', yanchor: 'bottom' }, shapes: shapes, margin: { l: 70, r: 16, t: 60, b: 70 },
      xaxis: { anchor: 'y3', domain: [0, 1], title: { text: 'exons (grouped by gene; gaps = gene boundaries)' } },
      yaxis: { domain: [0.72, 1], title: { text: 'log2(count + 0.5)' } } });
    L.yaxis = Object.assign({}, L.yaxis, { anchor: 'x', title: { text: 'log2(count + 0.5)', font: { color: CHROME[t].secondary } } });
    L.yaxis2 = Object.assign({}, L.yaxis, { domain: [0.38, 0.64], anchor: 'x', title: { text: 'obs / exp (95% interval)', font: { color: CHROME[t].secondary } } });
    L.yaxis3 = Object.assign({}, L.yaxis, { domain: [0, 0.30], anchor: 'x', title: { text: 'z-score', font: { color: CHROME[t].secondary } } });
    var pt = function (y, text) { return { xref: 'paper', yref: 'paper', x: 0, y: y, xanchor: 'left', yanchor: 'bottom', showarrow: false, text: text, font: { size: 12, color: CHROME[t].text } }; };
    L.annotations = [pt(1.0, '<b>1 \u00b7 Coverage</b> (test sample vs. references)'), pt(0.64, '<b>2 \u00b7 Observed / expected</b> (bars = 95% interval)'), pt(0.30, '<b>3 \u00b7 z-score</b> against references (dotted = \u00b13)')];
    return { data: data, layout: L };
  }
  function viewerGeneStrip(P, t) {
    // Panel 0: gene track with colored blocks per gene
    var px = asArray(P.px), gene = asArray(P.gene), group = asArray(P.gene_group || P.group), off = P.off ? asArray(P.off) : null;
    var n = px.length;
    return function () {
      var cc = catColors(8, t), c = CHROME[t], shapes = [], ann = [], k = 0, start = 0;
      for (var i = 1; i <= n; i++) {
        if (i === n || group[i] !== group[start]) {
          var g = gene[start] || '', isOff = !!(off && off[start]), cnt = i - start;
          var vert = cnt === 1 || (cnt === 2 && g.length > 5);
          shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: px[start] - 0.5, x1: px[i - 1] + 0.5, y0: 0.05, y1: 0.95, fillcolor: isOff ? c.axis : cc[k % 8], line: { width: 0 }, opacity: 0.9 });
          if (g) ann.push({ x: (px[start] + px[i - 1]) / 2, y: 0.5, xref: 'x', yref: 'paper', text: g, showarrow: false, textangle: vert ? -90 : 0, font: { color: isOff ? c.text : '#ffffff', size: vert ? 10 : 12 } });
          if (!isOff) k++;
          start = i;
        }
      }
      var lo = Math.min.apply(null, px) - 0.8, hi = Math.max.apply(null, px) + 0.8;
      var L = baseLayout(t, { shapes: shapes, annotations: ann, xaxis: { tickmode: 'array', tickvals: px, ticktext: asArray(P.xlab), tickangle: -45, range: [lo, hi], showgrid: false, showticklabels: false }, margin: { l: 56, r: 16, t: 4, b: 4 } });
      L.yaxis = { visible: false, range: [0, 1], fixedrange: true };
      return { data: [{ type: 'scatter', x: [px[0]], y: [0.5], mode: 'markers', marker: { opacity: 0 }, hoverinfo: 'skip' }], layout: L, config: { staticPlot: true } };
    };
  }

  tab('calls', 'CNV calls', function () { return true; }, function (root) {
    // Top: collapsible alerts
    root.appendChild(globalAlerts(6));
    // Then: inline sample picker
    root.appendChild(inlineSamplePicker('Sample'));

    var stack = h('div', { class: 'stack' });
    if (!CALLS.n) {
      stack.appendChild(alertBox('info', 'No CNV was called.', 'A run without calls is a valid result only for the chromosomes that were analysed (see the Chromosomes tab).'));
      root.appendChild(stack); return;
    }
    var f = { type: '', chr: '', minQ: 0, minT: 1, conf: { HIGH: true, MEDIUM: true, LOW: true, NA: true }, y: 'sample' };
    var typeSel = selectOf([{ value: '', text: 'DEL + DUP' }, { value: 'DEL', text: 'DEL' }, { value: 'DUP', text: 'DUP' }], '', null, 'Type');
    var chrCol = hasCol(CALLS, 'chr') ? 'chr' : 'Chr';
    var chrs = uniq(tcol(CALLS, chrCol)).sort(function (a, b) { return chrNum(a) - chrNum(b); });
    var chrSel = selectOf([{ value: '', text: 'All chromosomes' }].concat(chrs.map(function (c) { return { value: c, text: c }; })), '', null, 'Chromosome');
    var minQ = h('input', { type: 'number', min: 0, max: 99, value: 0, 'aria-label': 'Minimum Q_SOME / BF', style: { width: '80px', minWidth: '0' } });
    var minT = h('input', { type: 'number', min: 1, value: 1, 'aria-label': 'Minimum targets', style: { width: '80px', minWidth: '0' } });
    var ySel = selectOf([{ value: 'sample', text: 'Sample' }, { value: 'mlcn', text: 'MLCN' }, { value: 'q_some', text: 'Q_SOME / BF' }], 'sample', null, 'Map y axis');
    var confBox = h('div', { class: 'legend', role: 'group', 'aria-label': 'Confidence filter' });
    [['HIGH', 'HIGH'], ['MEDIUM', 'MEDIUM'], ['LOW', 'LOW'], ['NA', 'unscored']].forEach(function (p) {
      confBox.appendChild(h('button', { type: 'button', 'aria-pressed': 'true', onclick: function () {
        f.conf[p[0]] = !f.conf[p[0]]; this.setAttribute('aria-pressed', f.conf[p[0]] ? 'true' : 'false'); apply(); } },
        h('span', { class: 'swatch conf ' + p[0] }), p[1]));
    });
    var count = h('span', { class: 'muted small' });
    var tableHost = h('div'), mapBox = plotBox('');
    mapBox.style.height = Math.min(820, 320 + (SAMPLES.length > 8 ? 12 * Math.min(SAMPLES.length, 40) : 0)) + 'px';
    var Tdisp = callsTableDisplay();
    function pass(i) {
      if (state.sample && callSample(i) !== state.sample) return false;
      if (f.type && callType(i) !== f.type) return false;
      if (f.chr && String(tcol(CALLS, chrCol)[i]) !== f.chr) return false;
      var ck = callConf(i); ck = (ck === 'HIGH' || ck === 'MEDIUM' || ck === 'LOW') ? ck : 'NA';
      if (!f.conf[ck]) return false;
      if (f.minQ > 0) { var q = tcol(CALLS, 'q_some')[i]; var bf = tcol(CALLS, 'bf')[i]; var v = IS_CANOPE ? q : bf; if (!(isNum(v) && v >= f.minQ)) return false; }
      if (f.minT > 1) { var nt = tcol(CALLS, 'num_targ')[i]; var ne = tcol(CALLS, 'nexons')[i]; var v = IS_CANOPE ? nt : ne; if (!(isNum(v) && v >= f.minT)) return false; }
      return true;
    }
    var rowsNow = function () { var out = []; for (var i = 0; i < CALLS.n; i++) if (pass(i)) out.push(i); return out; };
    stack.appendChild(h('div', { class: 'toolbar' }, field('Type', typeSel), field('Chromosome', chrSel), field('Min Q_SOME / BF', minQ), field('Min targets', minT),
      h('div', { class: 'field' }, 'Confidence', confBox), count));
    var mapNote = '';
    stack.appendChild(card('Genome-wide map of calls', 'triangles: \u25BC deletion, \u25B2 duplication; click a call to open the viewer',
      h('div', null, h('div', { class: 'toolbar' }, field('Y axis', ySel)), mapBox, mapNote ? h('div', { class: 'note', text: mapNote }) : null)));
    // Calls table + call viewer side by side
    var detail = h('div');
    stack.appendChild(card('CNV calls', 'click a row to open the 4-panel call viewer', tableHost));
    stack.appendChild(detail);
    root.appendChild(stack);

    var tv = tableView(tableHost, Tdisp, { filename: 'cnv_calls', pageSize: 20, sort: ['Q_SOME / BF', 'desc'],
      render: { Type: typeChip, Confidence: confChip, Sample: function (v) { return h('strong', { text: v }); },
        'MLCN / CN': function (v) { return v || ''; },
        'Calibration': function (v) { return v === 'check' ? chip('warn', 'check') : v ? chip('ok', v) : ''; },
        'Plots': function (v, i) { return v ? h('button', { class: 'btn small-btn', type: 'button', text: 'View', onclick: function () { openCall(callId(i)); } }) : h('span', { class: 'muted', text: '\u2013' }); },
        Start: function (v) { return fmtInt(v); }, End: function (v) { return fmtInt(v); } },
      onRow: function (i) { if (hasPlot(i)) openCall(callId(i)); else { state.call = callId(i); notify(); } showCallDetail(detail, i); } });
    function apply() {
      tv.setPredicate(pass);
      count.textContent = tv.count() + ' of ' + CALLS.n + ' calls';
      if (mapBox._render) mapBox._render();
    }
    typeSel.addEventListener('change', function () { f.type = typeSel.value; apply(); });
    chrSel.addEventListener('change', function () { f.chr = chrSel.value; apply(); });
    minQ.addEventListener('input', function () { f.minQ = +minQ.value || 0; apply(); });
    minT.addEventListener('input', function () { f.minT = +minT.value || 1; apply(); });
    ySel.addEventListener('change', function () { f.y = ySel.value; genomeMap(mapBox, rowsNow, f.y); });
    genomeMap(mapBox, rowsNow, f.y);
    onState(function () { apply(); if (state.call !== null) { var r = CALL_ROW[state.call]; if (r !== undefined) showCallDetail(detail, r); } });
    apply();
    // Show first call detail if any
    var first = tv.visible()[0];
    if (first !== undefined) showCallDetail(detail, first);
  });

  function showCallDetail(host, i) {
    clear(host);
    var id = callId(i), P = PLOTS[String(id)];
    var head = h('div', { class: 'detail-head' },
      h('h2', { text: callSample(i) + ' \u00b7 ' + callGene(i) }),
      typeChip(callType(i)), confChip(callConf(i)),
      h('span', { class: 'muted', text: callInterval(i) }));
    var kvPairs = [
      ['Exons / Targets', exonsTargetsCell(i)],
      ['MLCN / CN', mlcnCnCell(i)],
      ['Q_SOME / BF', qsomeBfCell(i)],
      ['Corr / Refs', corrRefsCell(i)],
      ['Confidence', callConf(i)],
      ['Calibration', hasCol(CALLS, 'calibration') ? tcol(CALLS, 'calibration')[i] : ''],
      ['Flags', callFlags(i)]
    ];
    if (IS_ECHO) {
      kvPairs = kvPairs.concat([
        ['Observed', orNa(tcol(CALLS, 'observed')[i])],
        ['Expected', orNa(tcol(CALLS, 'expected')[i])],
        ['Ratio', orNa(tcol(CALLS, 'ratio')[i])],
        ['Bayes factor', orNa(tcol(CALLS, 'bf')[i])],
        ['Comparators', orNa(tcol(CALLS, 'comparators')[i])]
      ]);
    }
    if (IS_CANOPE) {
      kvPairs = kvPairs.concat([
        ['Interval', orNa(tcol(CALLS, 'interval')[i])],
        ['KB', orNa(tcol(CALLS, 'kb')[i])],
        ['Targets', orNa(tcol(CALLS, 'targets')[i])]
      ]);
    }
    var body = h('div', { class: 'stack' });
    if (!P) {
      body.appendChild(emptyNote(META.windows_limited ? 'Coverage profiles were only embedded for the first ' + META.windows_limit + ' calls.' : 'No coverage window for this call.'));
    } else {
      // background calibration alert
      if (P.bg && P.bg.flag) body.appendChild(alertBox('warn', 'Background calibration', fmt(P.bg.pct_outside) + '% of non-called exons lie outside the 95% interval.'));
      else if (P.subtitle) body.appendChild(alertBox('warn', 'Background calibration', P.subtitle));
      // 4-panel: strip + coverage + obs/exp + z
      var stripBox = plotBox('strip'), covBox = plotBox('mid'), ciBox = plotBox('mid'), zBox = plotBox('mid');
      body.appendChild(card('Gene track', null, stripBox));
      body.appendChild(card('Coverage (log2 count, test vs references)', null, covBox));
      body.appendChild(card('Observed / expected with 95% interval', null, ciBox));
      body.appendChild(card('Z-score vs references', null, zBox));
      var refsChk = h('input', { type: 'checkbox', checked: true, 'aria-label': 'Show references' });
      body.appendChild(h('label', { class: 'field inline' }, refsChk, 'show reference samples'));
      function drawPanels() {
        plot(stripBox, viewerGeneStrip(P, resolved()));
        plot(covBox, function (t) { return viewerFigure(P, i, t, refsChk.checked); });
        plot(ciBox, function (t) { return viewerFigure(P, i, t, refsChk.checked); });
        plot(zBox, function (t) { return viewerFigure(P, i, t, refsChk.checked); });
      }
      refsChk.addEventListener('change', drawPanels);
      drawPanels();
    }
    host.appendChild(card('Call details', '#' + id, h('div', null, head, h('div', { class: 'grid g-main' }, h('div', { class: 'stack' }, body), kv(kvPairs)))));
    host.scrollIntoView && host.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
  }

  // ================================================================ COHORT / PCA
  tab('cohort', 'Cohort / PCA', function () { return !!(PCA || CORR || DEPTH_HIST); }, function (root) {
    var stack = h('div', { class: 'stack' });
    var p1 = plotBox('tall'), p2 = plotBox('tall'), p3 = plotBox('short');
    var cards = [];
    if (PCA) cards.push(card('PCA of coverage profiles', (PCA.n_samples || PCA.names && PCA.names.length) + ' samples \u00b7 ' + fmtInt(PCA.n_targets || PCA.n_rows || 0) + ' targets (log2 counts)', p1));
    if (CORR) cards.push(card('Sample\u2013sample correlation', 'Spearman, hierarchically clustered', p2));
    if (cards.length) stack.appendChild(h('div', { class: 'grid g2' }, cards));
    if (!PCA) stack.appendChild(card('PCA of coverage profiles', null, emptyNote('PCA needs at least 3 samples with variable coverage.')));
    if (DEPTH_HIST && DEPTH_HIST.counts && DEPTH_HIST.counts.length)
      stack.appendChild(card('Target depth distribution', 'median depth of each target; dashed line = threshold', p3));
    root.appendChild(stack);
    setTimeout(function () {
      if (PCA) pcaPlot(p1);
      if (CORR) corrPlot(p2);
      if (DEPTH_HIST && DEPTH_HIST.counts && DEPTH_HIST.counts.length) depthHist(p3);
    }, 0);
  });
  function pcaPlot(div) {
    var P = PCA;
    plot(div, function (t) {
      var names = P.names || P.samples, pc1 = P.pc1 || P.pcs.PC1, pc2 = P.pc2 || P.pcs.PC2;
      var groups = P.group ? uniq(P.group) : ['All'], cc = catColors(Math.max(groups.length, 1), t), s = state.sample;
      var data = groups.map(function (g, k) {
        var ii = range(names.length).filter(function (i) { return (P.group ? P.group[i] : 'All') === g; });
        return { type: 'scatter', mode: 'markers+text', name: g, x: ii.map(function (i) { return pc1[i]; }), y: ii.map(function (i) { return pc2[i]; }),
          text: ii.map(function (i) { return names[i]; }), textposition: 'top right', textfont: { size: 10, color: CHROME[t].secondary },
          hovertemplate: '<b>%{text}</b><br>PC1 %{x:.2f}<br>PC2 %{y:.2f}<extra>' + g + '</extra>',
          marker: { color: cc[k], size: ii.map(function (i) { return s && names[i] === s ? 15 : 9; }), line: { width: ii.map(function (i) { return s && names[i] === s ? 3 : 0; }), color: CHROME[t].text } } };
      });
      var v1 = P.var1 || (P.var ? P.var[0] : 0), v2 = P.var2 || (P.var ? P.var[1] : 0);
      return { data: data, layout: baseLayout(t, { xlab: 'PC1 (' + fmt(v1) + '%)', ylab: 'PC2 (' + fmt(v2) + '%)', legend: groups.length > 1, margin: { l: 56, r: 24, t: 28, b: 48 } }) };
    });
  }
  function corrPlot(div) {
    var Cm = CORR;
    plot(div, function (t) {
      var z = Cm.z, names = Cm.names || Cm.samples;
      var lo = 1; z.forEach(function (r) { r.forEach(function (v) { if (isNum(v) && v < lo) lo = v; }); });
      return { data: [{ type: 'heatmap', x: names, y: names, z: z, zmin: Math.max(0.9, lo), zmax: 1,
        colorscale: SEQ.map(function (c, i) { return [i / (SEQ.length - 1), c]; }), colorbar: { thickness: 10, len: 0.6, title: { text: 'r' } },
        hovertemplate: '%{y} \u00d7 %{x}<br>r = %{z:.3f}<extra></extra>' }],
        layout: baseLayout(t, { xaxis: { type: 'category', tickangle: -45, showgrid: false }, yaxis: { type: 'category', autorange: 'reversed', showgrid: false }, margin: { l: 90, r: 16, t: 16, b: 90 } }) };
    });
  }
  function depthHist(div) {
    var Hh = DEPTH_HIST;
    plot(div, function (t) {
      var sem = SEM[t];
      var data = [{ type: 'bar', x: Hh.mid, y: Hh.counts, marker: { color: Hh.mid.map(function (m) { return (TH.qc_min_cov || TH.min_cov) && Math.pow(10, m) < (TH.qc_min_cov || TH.min_cov) ? sem.bad : sem.dup; }) },
        hovertemplate: 'median depth \u2248 %{customdata:.0f}<br>%{y} targets<extra></extra>', customdata: Hh.mid.map(function (m) { return Math.pow(10, m); }) }];
      var thr = Math.log10(Math.max(TH.qc_min_cov || TH.min_cov || 1, 1));
      return { data: data, layout: baseLayout(t, { xlab: 'Median depth per target (log)', ylab: 'Targets',
        shapes: (TH.qc_min_cov || TH.min_cov) ? [{ type: 'line', xref: 'x', yref: 'paper', x0: thr, x1: thr, y0: 0, y1: 1, line: { color: sem.warn, width: 2, dash: 'dash' } }] : [] }) };
    });
  }

  // ================================================================ CHROMOSOMES
  tab('chromosomes', 'Chromosomes', function () { return true; }, function (root) {
    var stack = h('div', { class: 'stack' });
    var both = SUM.both_chromosome_classes;
    if (both === false || META.scope_incomplete) stack.appendChild(alertBox('fail', 'Autosomes and gonosomes were NOT both analysed.', 'Results cover only part of the genome (mode: ' + (META.mode_label || META.modechrom_label) + ').'));
    else if (both === true) stack.appendChild(alertBox('ok', 'Autosomes and chrX were both analysed.', null));
    if (CHR.messages) asArray(CHR.messages).forEach(function (m) { stack.appendChild(alertBox('warn', null, m)); });
    if (CHR.status && CHR.status.n) stack.appendChild(card('Chromosome classes', null, chromStrip()));
    // Inventory
    var inv = CHR.inventory;
    if (inv && inv.n) {
      var invHost = h('div'), barBox = plotBox('short');
      tableView(invHost, inv, { filename: 'chromosomes', pageSize: 30, labels: { targets_in_bed: 'targets in BED', targets_analysed: 'targets analysed' },
        render: { analysed: function (v) { return v ? chip('ok', 'analysed') : chip('warn', 'not analysed'); } },
        cellClass: { targets_analysed: function (v, i) { return tcol(inv, 'analysed')[i] ? '' : 'warn'; } } });
      stack.appendChild(h('div', { class: 'grid g2' }, card('Targets per chromosome', 'in BED vs. analysed', barBox), card('Chromosome inventory', null, invHost)));
      plot(barBox, function (t) {
        var c = tcol(inv, 'chromosome');
        return { data: [{ type: 'bar', name: 'in BED', x: c, y: tcol(inv, 'targets_in_bed'), marker: { color: SEM[t].NA } },
                        { type: 'bar', name: 'analysed', x: c, y: tcol(inv, 'targets_analysed'), marker: { color: CAT[t][0] } }],
          layout: baseLayout(t, { ylab: 'targets', barmode: 'overlay', legend: true, xaxis: { type: 'category', tickangle: -45 }, margin: { l: 56, r: 16, t: 10, b: 80 } }) };
      });
    }
    // Depth heatmap
    var D = CHR.depth;
    if (D && D.chromosomes && asArray(D.chromosomes).length) {
      var hbox = plotBox('tall');
      stack.appendChild(card('Chromosome depth relative to autosomes', 'median target depth / sample autosomal median', hbox));
      plot(hbox, function (t) {
        var dz = asArray(D.z).map(asArray);
        return { data: [{ type: 'heatmap', x: asArray(D.samples), y: asArray(D.chromosomes), z: dz, zmin: 0, zmax: 2,
          colorscale: [[0, DIV[0]], [0.25, DIV[1]], [0.5, DIV[2]], [0.75, DIV[3]], [1, DIV[4]]],
          hovertemplate: '%{x}<br>%{y}: %{z:.2f}<extra></extra>', colorbar: { title: { text: 'ratio' } } }],
          layout: baseLayout(t, { xaxis: { type: 'category', tickangle: -45 }, yaxis: { type: 'category', autorange: 'reversed' }, margin: { l: 70, r: 16, t: 10, b: 90 } }),
          onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) setSample(p.x); } };
      });
    }
    // Sex check
    var X = CHR.sex;
    if (X && X.n) {
      var sbox = plotBox(''), sh = h('div');
      stack.appendChild(h('div', { class: 'grid g2' },
        card('Sex check', 'chrX vs chrY depth ratio', sbox), card('Sex check table', null, sh)));
      tableView(sh, X, { filename: 'sex_check', pageSize: 10, labels: { x_ratio: 'chrX ratio', y_ratio: 'chrY ratio' },
        cellClass: { mismatch: function (v) { return v ? 'warn' : ''; } }, onRow: function (i) { setSample(tcol(X, 'sample')[i]); } });
    }
    if (!stack.children.length) stack.appendChild(card('Chromosomes', null, emptyNote('No chromosome information available.')));
    root.appendChild(stack);
  });

  // ================================================================ METHODS
  tab('methods', 'Methods', function () { return true; }, function (root) {
    var stack = h('div', { class: 'stack' });
    var p = PARAMS, R = CONF_RULES, sc = SCORING;

    // 1. Confidence score (both ECHO scoring and CANOPE conf_rules)
    if (sc || R) {
      var confRulesList = h('ul', null,
        h('li', null, h('b', { text: 'HIGH' }), ' only if all HIGH criteria are met.'),
        h('li', null, h('b', { text: 'MEDIUM' }), ' if HIGH is not met but all MEDIUM criteria are.'),
        h('li', null, h('b', { text: 'LOW' }), ' otherwise.'));
      if (R && R.genes && asArray(R.genes).length) confRulesList.appendChild(h('li', null, h('b', { text: 'Always LOW' }), ' (homology): ' + asArray(R.genes).join(', ')));
      if (sc && sc.forced_low && sc.forced_low.length) confRulesList.appendChild(h('li', null, h('b', { text: 'Always LOW' }), ' (homology): ' + sc.forced_low.join(', ')));
      if (sc && sc.terminal) confRulesList.appendChild(h('li', { text: sc.terminal }));
      var confTable;
      if (sc && sc.rows) {
        confTable = h('table', { class: 'tbl' },
          h('thead', null, h('tr', null, ['Category', 'Correlation', 'References', 'Fold change'].map(function (c) { return h('th', { text: c }); }))),
          h('tbody', null, sc.rows.map(function (r) { return h('tr', null, r.map(function (v, j) { return h('td', null, j === 0 ? chip(v, v) : v); })); })));
      } else {
        confTable = h('table', { class: 'tbl' },
          h('thead', null, h('tr', null, ['Category', 'Q_SOME / BF', 'References', 'Targets'].map(function (c) { return h('th', { text: c }); }))),
          h('tbody', null,
            h('tr', null, h('td', null, chip('conf HIGH', 'HIGH')), h('td', { text: '>= ' + (R.high_q || '') }), h('td', { text: '>= ' + (R.high_refs || '') }), h('td', { text: '>= ' + (R.high_targ || '') })),
            h('tr', null, h('td', null, chip('conf MEDIUM', 'MEDIUM')), h('td', { text: '>= ' + (R.med_q || '') }), h('td', { text: '>= ' + (R.med_refs || '') }), h('td', { text: '>= ' + (R.med_targ || '') })),
            h('tr', null, h('td', null, chip('conf LOW', 'LOW')), h('td', { text: 'anything else' }), h('td', { text: 'any' }), h('td', { text: 'any' }))));
      }
      stack.appendChild(card('Confidence score', null, h('div', null, confRulesList, h('div', { class: 'tbl-wrap', style: { marginTop: '10px' } }, confTable))));
    }

    // 2. Reading the plots
    stack.appendChild(card('Reading the plots', null, h('div', { class: 'prose' },
      h('p', { text: 'Each CNV call shows the test sample against its selected reference samples over the affected exons \u00b1 5 neighbouring targets. The 4-panel viewer shows: (0) a gene track, (1) coverage on a log2 scale, (2) observed/expected with the 95% interval, (3) z-scores against the references.' }),
      h('p', { text: 'For ECHO the interval is the beta-binomial 95% interval (overdispersion phi fitted by ExomeDepth). For CANOPE it is the negative-binomial predictive interval used by the HMM. A call should sit clearly outside the interval while non-called neighbours cluster near zero.' }),
      h('p', { text: 'If a high fraction of neighbouring (non-called) exons fall outside the interval, the call is flagged "check" / "background calibration": confirm by another method or rerun with a better reference set.' }),
      h('p', { text: 'Points are never joined by lines: exons are discrete positions, and connecting them would imply a continuous signal (and draw through the inserted gene boundaries).' }))));

    // 3. Parameters
    if (Object.keys(p).length) stack.appendChild(card('Run parameters', null, kv(Object.keys(p).map(function (k) { return [k, Array.isArray(p[k]) ? p[k].join(', ') : String(p[k])]; }))));

    // 4. Settings (ECHO)
    if (SETTINGS && SETTINGS.length) {
      var skv = kv(SETTINGS.map(function (s) { return [s[0], s[1]]; }));
      stack.appendChild(card('Settings', 'from the config used for this run', skv));
    }

    // 5. QC thresholds
    var thrPairs = [
      ['Median depth \u2265', fmt(TH.qc_min_cov || TH.min_cov)],
      ['Total reads \u2265', fmtInt(TH.qc_min_total_reads || TH.min_total_reads)],
      ['Best correlation \u2265', fmt(TH.qc_min_corr || TH.min_corr)],
      ['Exon CV \u2264', fmt(TH.qc_max_exon_cv || TH.max_exon_cv)]
    ];
    if (TH.qc_zscore) thrPairs.push(['Noise robust z \u2264', fmt(TH.qc_zscore)]);
    thrPairs.push(['Reference samples \u2265', '5 (warning below)']);
    stack.appendChild(card('QC thresholds used for the warnings', null, kv(thrPairs)));

    // 6. How it works
    stack.appendChild(card(IS_ECHO ? 'How ECHO works' : 'How CANOPE works', null, h('div', { class: 'prose' },
      IS_ECHO
        ? [h('p', { text: 'ECHO is an ExomeDepth-based exome CNV caller. For every sample it selects the most correlated reference samples, combines them into an expected coverage per exon, and models the read count with a beta-binomial distribution whose overdispersion (phi) is fitted from the references. A Bayes factor scores each candidate call; confidence (HIGH / MEDIUM / LOW) is derived from the correlation, the number of references and the fold change.' }),
           h('p', { text: 'One run covers one chromosome mode: autosomes (A), chrX (X) or chrY (Y). Sex chromosomes are therefore never analysed together with the autosomes.' })]
        : [h('p', { text: 'CANOPE is a CANOES-derived exome CNV caller. For every sample it selects the most correlated reference samples, combines them with non-negative weights into an expected coverage per exon, and models the read count with a negative binomial whose variance is estimated from the references. A hidden Markov model with distance-decayed transitions (Viterbi) calls deletions and duplications; Q_SOME is the phred-scaled probability that at least one exon of the call is abnormal, and MLCN the most likely copy number.' }),
           h('p', { text: 'One run covers one chromosome mode: autosomes (A), chrX with female references (XX) or chrX/chrY with male references (XY).' })])));

    // 7. Warning glossary
    var glossary = [
      ['Chromosomes not analysed', 'The BED/count table contains targets that this run did not call. A missing call there is not a negative result.'],
      ['Sample failed / excluded', 'Failed samples raised an error while calling; excluded samples were dropped before calling. Neither has CNV calls, and excluded samples are not used as references.'],
      ['Low depth / low reads', 'Median target depth or total reads below the QC threshold: reduced sensitivity, especially for single-exon events.'],
      ['Low correlation / few references', 'No other sample resembles this one (batch, capture kit, tissue), so the expected coverage is poor and false calls are likely.'],
      ['Noisy sample', 'The exon-count dispersion is an outlier within the cohort (robust z-score).'],
      ['Genotyping skipped', 'Too many calls were produced; Q_SOME and MLCN were not computed and every call is LOW.'],
      ['Sex / mode mismatch', 'The annotated sex disagrees with chrX/chrY coverage, or the sample does not match the sex expected by the mode.'],
      ['Background calibration', 'A high share of non-called exons around a call fall outside their 95% interval, meaning the expected-coverage model is miscalibrated.']];
    stack.appendChild(card('Warning glossary', null, kv(glossary)));

    // 8. Sample table
    if (DATA.sample_table && DATA.sample_table.n) {
      var h1 = h('div');
      tableView(h1, DATA.sample_table, { filename: 'sample_table', pageSize: 15 });
      stack.appendChild(card('Sample table', null, h1));
    }
    // 9. Reference BAMs
    if (DATA.ref_bams && DATA.ref_bams.n) {
      var h2 = h('div');
      tableView(h2, DATA.ref_bams, { filename: 'reference_bams', pageSize: 15 });
      stack.appendChild(card('External reference BAMs', null, h2));
    }
    root.appendChild(stack);
  });

  // ================================================================ RUN
  tab('run', 'Run', function () { return true; }, function (root) {
    var L = LOG || {}, blocks = [];
    if (L.warnings && L.warnings.length) blocks.push(card('Warnings and errors from the log (' + L.warnings.length + ')', null, h('pre', { text: L.warnings.join('\n') })));
    if (DATA.config_yaml) blocks.push(card('Configuration used', null, h('pre', { text: DATA.config_yaml })));
    if (L.tail && L.tail.length) blocks.push(card('Pipeline log' + (L.truncated ? ' (tail)' : ''), null, h('pre', { text: L.tail.join('\n') })));
    if (L.session && L.session.length) blocks.push(card('Session info', null, h('pre', { text: L.session.join('\n') })));
    if (!blocks.length) blocks.push(emptyNote('No run information was recorded.'));
    root.appendChild(h('div', { class: 'stack' }, blocks));
  });

  // ================================================================ shell
  var panels = {}, tabButtons = {}, active = null;
  function activate(id, push) {
    var t = TABS.filter(function (x) { return x.id === id && x.available(); })[0] || TABS.filter(function (x) { return x.available(); })[0];
    if (!t) return;
    id = t.id;
    Object.keys(panels).forEach(function (k) {
      var on = k === id; panels[k].hidden = !on;
      tabButtons[k].setAttribute('aria-selected', on ? 'true' : 'false'); tabButtons[k].tabIndex = on ? 0 : -1;
    });
    if (!t.built) {
      t.built = true;
      try { t.build(panels[id]); }
      catch (e) { panels[id].appendChild(h('div', { class: 'empty', text: 'This tab could not be built: ' + e.message })); if (window.console) console.error(e); }
    } else refreshVisible(panels[id]);
    setTimeout(function () { refreshVisible(panels[id]); }, 60);
    active = id;
    if (push) { try { history.replaceState(null, '', '#tab=' + id); } catch (e) {} }
  }
  function init() {
    applyThemeAttr();
    document.title = META.title;
    var tabsEl = $('#tabs'), main = $('#panels');
    var avail = TABS.filter(function (t) { return t.available(); });
    avail.forEach(function (t, i) {
      var b = h('button', { class: 'tab', role: 'tab', id: 'tab-' + t.id, 'aria-controls': 'panel-' + t.id, 'aria-selected': 'false', tabindex: '-1', type: 'button', text: t.label,
        onclick: function () { activate(t.id, true); },
        onkeydown: function (e) {
          var k = e.key === 'ArrowRight' ? 1 : e.key === 'ArrowLeft' ? -1 : 0;
          if (!k) return; e.preventDefault();
          var nx = avail[(i + k + avail.length) % avail.length]; activate(nx.id, true); tabButtons[nx.id].focus();
        } });
      tabButtons[t.id] = b; tabsEl.appendChild(b);
      panels[t.id] = h('section', { class: 'panel', role: 'tabpanel', id: 'panel-' + t.id, 'aria-labelledby': 'tab-' + t.id, hidden: true });
      main.appendChild(panels[t.id]);
    });
    // header sample picker
    var pick = $('#sample-pick');
    if (pick && SAMPLES.length > 0) {
      headerPicker = h('select', { id: 'sample-sel', 'aria-label': 'Sample' },
        h('option', { value: '', text: 'All samples (' + SAMPLES.length + ')' }),
        SAMPLES.map(function (s) { return h('option', { value: s, text: s }); }));
      headerPicker.value = state.sample;
      headerPicker.addEventListener('change', function () { setSample(headerPicker.value); });
      pick.appendChild(h('label', { for: 'sample-sel', text: 'Sample ' }));
      pick.appendChild(headerPicker);
    } else if (pick) { pick.hidden = true; }

    var btn = $('#theme-btn');
    var labelFn = function () { btn.textContent = 'Theme: ' + themePref; };
    labelFn();
    btn.addEventListener('click', function () {
      themePref = themePref === 'auto' ? 'light' : themePref === 'light' ? 'dark' : 'auto';
      store('unified-report-theme', themePref); applyThemeAttr(); labelFn(); rerenderAll();
    });
    if (window.matchMedia) {
      var mq = window.matchMedia('(prefers-color-scheme: dark)');
      var onChange = function () { if (themePref === 'auto') rerenderAll(); };
      if (mq.addEventListener) mq.addEventListener('change', onChange);
    }
    var m = /tab=([\w-]+)/.exec(location.hash);
    activate(m ? m[1] : 'overview', false);
    window.addEventListener('hashchange', function () { var mm = /tab=([\w-]+)/.exec(location.hash); if (mm && mm[1] !== active) activate(mm[1], false); });
    document.documentElement.setAttribute('data-ready', '1');
  }

  function ensurePlotly(cb, timeoutMs) {
    timeoutMs = timeoutMs || 8000;
    var done = false;
    function has() { return (typeof window.Plotly !== 'undefined') && window.Plotly && (typeof window.Plotly.react === 'function'); }
    function finish(ok) {
      if (done) return; done = true; clearTimeout(timer);
      if (ok && has()) cb();
      else {
        var host = document.getElementById('panels');
        if (host) host.innerHTML = '<div class="empty">Plotly.js failed to load. Check your connection or reinstall the R package.</div>';
        if (window.console) console.error('[report] Plotly unavailable.');
      }
    }
    var timer = setTimeout(function () { finish(false); }, timeoutMs);
    if (has()) { finish(true); return; }
    var s = document.createElement('script');
    s.src = 'https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js';
    s.async = false; s.onload = function () { finish(true); }; s.onerror = function () { finish(false); };
    document.head.appendChild(s);
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', function () { ensurePlotly(init); });
  else ensurePlotly(init);
})();
)---"


# ----------------------------------------------------------------------------
# R helpers (shared by both packages)
# ----------------------------------------------------------------------------

#' CSS string for the unified report
#' @export
report_css <- function() REPORT_CSS

#' JS string for the unified report
#' @export
report_js <- function() REPORT_JS

#' Column-oriented table serialiser (data.frame -> list for JSON)
#' @export
rep_table <- function(df) {
  if (is.null(df)) return(list(cols = I(character()), n = 0L, data = list()))
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  data <- lapply(df, function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (is.double(v)) v <- round(v, 4)
    I(unname(v))
  })
  list(cols = I(colnames(df)), n = nrow(df), data = data)
}

#' Empty table placeholder
#' @export
rep_empty_table <- function(cols) {
  rep_table(stats::setNames(rep(list(character()), length(cols)), cols))
}

#' Inert JSON for a <script type="application/json"> element
#' @export
report_json <- function(x) {
  json <- as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null",
                                        digits = NA, force = TRUE))
  json <- gsub("</", "<\\/", json, fixed = TRUE)
  gsub("<!--", "<\\!--", json, fixed = TRUE)
}

#' Read the plotly.js bundle from the installed plotly R package (or NULL)
#' @export
report_plotly_js <- function() {
  if (!requireNamespace("plotly", quietly = TRUE)) return(NULL)
  deps <- tryCatch(suppressWarnings(suppressMessages(htmltools::findDependencies(plotly::plot_ly()))),
                   error = function(e) list())
  dep <- Filter(function(d) identical(d$name, "plotly-main"), deps)
  if (!length(dep)) return(NULL)
  dep <- dep[[1]]
  path <- file.path(dep$src$file, dep$script)
  if (!file.exists(path)) path <- system.file(dep$src$file, dep$script, package = dep$package %||% "plotly")
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  raw <- readBin(path, "raw", file.size(path))
  txt <- rawToChar(raw); Encoding(txt) <- "UTF-8"
  gsub("</script", "<\\/script", txt, fixed = TRUE)
}

#' HTML-escape a string
#' @export
report_esc <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' Assemble the full HTML document from a payload and write it to out_file.
#'
#' Both echo_report.R and canope_report.R call this with their payload
#' (a list with a \code{meta$source} field set to "ECHO" or "CANOPE").
#' The resulting HTML is byte-identical in structure and styling regardless
#' of source; only the embedded JSON blob differs.
#'
#' @param payload List (the report data, with meta$source set).
#' @param settings Optional list; \code{settings$report_plotly} can be "inline" (default) or "cdn".
#' @param out_file Path to write the HTML file to.
#' @return Invisibly, out_file.
#' @export
report_html <- function(payload, settings = list(), out_file) {
  css <- report_css()
  js  <- report_js()
  mode <- tolower(as.character(settings$report_plotly %||% "inline")[1])
  plotly_tag <- NULL
  if (mode != "cdn") {
    src <- report_plotly_js()
    if (is.null(src)) {
      message("[WARNING] plotly.js bundle not found; loading from CDN (needs internet).")
    } else {
      plotly_tag <- paste0("<script>\n", src, "\n</script>\n")
    }
  }
  if (is.null(plotly_tag))
    plotly_tag <- "<script src=\"https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js\"></script>\n"

  title <- report_esc(payload$meta$title)
  html <- paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>", title, "</title>\n<style>\n", css, "\n</style>\n</head>\n<body>\n",
    "<header class=\"top\"><div class=\"top-inner\">\n",
    "<div class=\"title-row\"><h1>", title, "</h1>",
    "<span class=\"badge\">", report_esc(payload$meta$source), "</span>",
    "<span class=\"badge\">", report_esc(payload$meta$engine_label %||% payload$meta$engine %||% ""), "</span>",
    "<span class=\"spacer\"></span>",
    "<span id=\"sample-pick\" class=\"sample-pick\"></span>",
    "<button id=\"theme-btn\" class=\"btn\" type=\"button\" aria-label=\"Colour theme\">Theme</button></div>\n",
    "<nav id=\"tabs\" class=\"tabs\" role=\"tablist\" aria-label=\"Report sections\"></nav>\n",
    "</div></header>\n<main id=\"panels\"></main>\n",
    "<noscript><p style=\"padding:16px\">This report needs JavaScript enabled.</p></noscript>\n",
    plotly_tag,
    "<script id=\"unified-data\" type=\"application/json\">", report_json(payload), "</script>\n",
    "<script>\n", js, "\n</script>\n</body>\n</html>\n")
  con <- file(out_file, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(html)), con)
  message("[INFO] ", payload$meta$source, " report written: ", out_file, " (",
          format(structure(file.size(out_file), class = "object_size"), units = "auto"), ")")
  invisible(out_file)
}

# null-coalesce helper
`%||%` <- function(a, b) if (is.null(a)) b else a
