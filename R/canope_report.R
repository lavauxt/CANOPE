#' Generate the interactive, multi-sample HTML report for CANOPE
#'
#' Builds \strong{one self-contained HTML file} (no server, no Shiny app, no
#' network needed) from a CANOPE workspace: tabs for an overview, a consolidated
#' warnings list, per-sample views, QC metrics, chromosome coverage /
#' autosome-gonosome status, PCA, a filterable table of all CNV calls with a
#' genome-wide map, and a call viewer with the three interactive panels
#' (coverage, log2 observed/expected with its negative-binomial interval, and
#' z-scores against the references). All data are embedded as JSON, plotly.js is
#' inlined once and every tab is drawn lazily when opened -- the same design as
#' the SPARKS report.
#'
#' The report handles any number of samples (calls, QC, references and plots
#' are all per sample), lists samples with \emph{no} call as well as samples
#' that failed or were excluded (so "no CNV" is never confused with "not
#' analysed"), and raises warnings for: per-sample metrics outside their
#' thresholds, noisy samples, genotyping skipped, near-zero-coverage runs that
#' could hide homozygous deletions, sex mismatches, calls in genes with known
#' homology problems, and -- explicitly -- autosomes or gonosomes that were not
#' analysed in this run.
#'
#' @param rdata_output   Path to the CANOPE workspace \code{.RData} written by
#'   \code{\link{run_canope}}.
#' @param qc_metrics_file Optional path to the QC-metrics TSV
#'   (\code{\link{run_canope_qc_metrics}}).
#' @param output_dir     Directory for the HTML (created if missing).
#' @param settings       Optional named list of thresholds.
#' @param sample_table   Optional path to a sample sheet.
#' @param log_file       Optional pipeline log.
#' @param prefix         Output name: \code{CANOPE_<prefix>_report.html}.
#' @param template_path  Deprecated and ignored.
#' @param max_call_plots Maximum number of calls whose interactive plots are
#'   embedded.
#' @param max_qc_rows    Row cap for the embedded exon-level QC flag table.
#' @param flank          Number of flanking exons shown each side of a call.
#' @param title          Report title.
#' @param plotly_source  \code{"inline"} (default) or \code{"cdn"}.
#'
#' @return Invisibly returns the path to the HTML file.
#' @export
generate_canope_report <- function(rdata_output,
                                   qc_metrics_file = NULL,
                                   output_dir,
                                   settings = NULL,
                                   sample_table = NULL,
                                   log_file = NULL,
                                   prefix = "CANOPE",
                                   template_path = NULL,
                                   max_call_plots = 300L,
                                   max_qc_rows = 3000L,
                                   flank = 5L,
                                   title = NULL,
                                   plotly_source = c("inline", "cdn")) {
  plotly_source <- match.arg(plotly_source)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for report generation.")
  if (!is.null(template_path))
    message("[INFO] 'template_path' is ignored: the report is no longer an R Markdown template.")
  if (!file.exists(rdata_output)) stop("[ERROR] Workspace not found: ", rdata_output)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("Failed to create output directory: ", output_dir)
  }

  env <- new.env()
  load(rdata_output, envir = env)
  required <- c("cnv_calls", "counts", "bed_file", "models", "refs")
  missing_vars <- setdiff(required, ls(env))
  if (length(missing_vars) > 0)
    stop("[ERROR] Missing variables in RData: ", paste(missing_vars, collapse = ", "))

  data <- .rp_build_data(env, qc_metrics_file, settings, sample_table, log_file,
                         prefix, title, max_call_plots, max_qc_rows, flank)
  json <- .rp_json(data)
  html <- .rp_html(data$meta, json, plotly_source)

  base <- sub("^CANOPE_?", "", prefix)
  out_file <- file.path(normalizePath(output_dir, mustWork = TRUE),
                        if (nzchar(base)) paste0("CANOPE_", base, "_report.html") else "CANOPE_report.html")
  con <- file(out_file, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(html)), con)
  message("[INFO] CANOPE report written: ", out_file, " (",
          format(structure(file.size(out_file), class = "object_size"), units = "auto"), ")")
  invisible(out_file)
}


# ============================================================================
# Embedded report assets (CSS + client-side JS)
#
# These used to live in inst/report/*.css and *.js and were read at runtime
# via .rp_asset().  They are embedded here as R raw strings so the rendered
# HTML is fully self-contained and the package works from a bare
# devtools::load_all() without any inst/report/ present.
#
# Raw-string delimiter: `r"---(` ... `)---"`. R only allows dash sequences as
# raw-string delimiters, so the terminator `)---"` must not occur inside the
# asset -- and it cannot: it is neither valid CSS nor valid JavaScript.
# ============================================================================

.CANOPE_REPORT_CSS <- r"---(
/* CANOPE interactive report -- design tokens, layout and components (same look as the SPARKS report).
   Colours follow the validated palette: categorical hues are applied by JS in
   fixed order; chrome (surface / ink / grid) swaps between light and dark. */
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
  --ok-ink: #0b6b3a;
  --fail-soft: #fdeceb;
  --warn-soft: #fff3d6;
  --ok-soft: #e4f5ea;
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
    --fail: #ff8a80; --ok-ink: #5fd38d; --fail-soft: #3a1c1a; --warn-soft: #3a2d0b; --ok-soft: #12301f;
    --c-del: #e66767; --c-dup: #3987e5; --c-high: #199e70; --c-medium: #c98500; --c-low: #898781; --c-na: #555552;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
  --muted: #898781; --grid: #2c2c2a; --axis: #383835;
  --ring: rgba(255, 255, 255, 0.10); --accent: #3987e5; --accent-ink: #0b0b0b;
  --accent-soft: #1c2a3d; --good: #0ca30c; --warn: #fab219;
    --fail: #ff8a80; --ok-ink: #5fd38d; --fail-soft: #3a1c1a; --warn-soft: #3a2d0b; --ok-soft: #12301f;
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
.grid { display: grid; gap: 12px; }
.g2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
.g3 { grid-template-columns: repeat(3, minmax(0, 1fr)); }
.g-main { grid-template-columns: minmax(0, 1fr) 360px; align-items: start; }
.g-viewer { grid-template-columns: 340px minmax(0, 1fr); align-items: start; }
.stack { display: flex; flex-direction: column; gap: 12px; }
.toolbar { display: flex; gap: 12px; flex-wrap: wrap; align-items: flex-end; margin-bottom: 12px; }
.field { display: flex; flex-direction: column; gap: 3px; font-size: 12px; color: var(--text-2); min-width: 0; max-width: 100%; }
.field.inline { flex-direction: row; align-items: center; gap: 6px; }
select, input[type="text"], input[type="number"], input[type="search"] {
  font: inherit; color: var(--text); background: var(--surface); border: 1px solid var(--axis);
  border-radius: 6px; padding: 5px 8px; min-width: 120px; max-width: 280px; }
input[type="range"] { accent-color: var(--accent); }
button.btn { font: inherit; padding: 5px 12px; border-radius: 6px; cursor: pointer;
             border: 1px solid var(--axis); background: var(--surface); color: var(--text); }
button.btn:hover { background: var(--accent-soft); }
button.btn.primary { background: var(--accent); border-color: var(--accent); color: var(--accent-ink); }
button.btn[disabled] { opacity: 0.45; cursor: not-allowed; }
.plot { width: 100%; min-height: 360px; }
.plot.tall { min-height: 560px; } .plot.short { min-height: 280px; }
.empty { padding: 28px 12px; text-align: center; color: var(--muted); }

/* ---- Stat tiles ---- */
.tiles { display: grid; gap: 12px; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); }
.tile { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius); padding: 12px 14px; }
.tile .v { font-size: 26px; font-weight: 600; letter-spacing: -0.02em; }
.tile .l { color: var(--text-2); font-size: 12px; }

/* ---- Legend (HTML so categories can be toggled) ---- */
.legend { display: flex; flex-wrap: wrap; gap: 4px 12px; margin-top: 8px; font-size: 12px; color: var(--text-2); }
.legend button { appearance: none; border: 0; background: none; font: inherit; color: inherit;
                 display: inline-flex; gap: 6px; align-items: center; cursor: pointer; padding: 2px 4px; border-radius: 4px; }
.legend button:hover { background: var(--accent-soft); }
.legend button[aria-pressed="false"] { opacity: 0.4; text-decoration: line-through; }
.swatch { width: 10px; height: 10px; border-radius: 3px; display: inline-block; }

/* ---- Tables ---- */
.tbl-wrap { overflow: auto; border: 1px solid var(--ring); border-radius: 8px; max-height: 520px; }
table.tbl { border-collapse: collapse; width: 100%; font-size: 12.5px; }
.tbl th, .tbl td { padding: 5px 10px; text-align: left; white-space: nowrap; border-bottom: 1px solid var(--grid); }
.tbl th { position: sticky; top: 0; background: var(--surface); cursor: pointer; user-select: none;
          font-weight: 600; color: var(--text-2); border-bottom: 1px solid var(--axis); z-index: 1; }
.tbl th[aria-sort="ascending"]::after { content: " \25B2"; font-size: 9px; }
.tbl th[aria-sort="descending"]::after { content: " \25BC"; font-size: 9px; }
.tbl td.num { text-align: right; font-variant-numeric: tabular-nums; }
.tbl tbody tr:hover { background: var(--accent-soft); }
.tbl-bar { display: flex; gap: 10px; align-items: center; margin: 8px 0; flex-wrap: wrap; }
.pager { display: inline-flex; gap: 6px; align-items: center; margin-left: auto; }

/* ---- Selection panel ---- */
.bars { display: flex; flex-direction: column; gap: 5px; }
.bar-row { display: grid; grid-template-columns: 110px 1fr 70px; gap: 8px; align-items: center; font-size: 12px; }
.bar-row .lab { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; color: var(--text-2); }
.bar-track { background: var(--grid); border-radius: 4px; height: 10px; overflow: hidden; }
.bar-fill { height: 100%; border-radius: 4px 4px 4px 4px; }
.note { font-size: 12px; color: var(--muted); }
.kv { display: grid; grid-template-columns: max-content 1fr; gap: 3px 14px; font-size: 12.5px; }
.kv dt { color: var(--text-2); } .kv dd { margin: 0; overflow-wrap: anywhere; }

/* ---- Gallery ---- */
.gallery { display: grid; gap: 12px; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); }
.shot { background: var(--surface); border: 1px solid var(--ring); border-radius: 8px; overflow: hidden; margin: 0; }
.shot button { all: unset; display: block; cursor: zoom-in; width: 100%; }
.shot img { display: block; width: 100%; height: 190px; object-fit: contain; background: #fff; }
.shot figcaption { padding: 6px 10px; font-size: 12px; color: var(--text-2); overflow-wrap: anywhere; }
.shot .nopreview { height: 190px; display: grid; place-items: center; color: var(--muted); background: var(--page); }
dialog.lightbox { border: 0; padding: 0; background: transparent; max-width: 96vw; max-height: 96vh; }
dialog.lightbox::backdrop { background: rgba(0, 0, 0, 0.7); }
dialog.lightbox img { max-width: 96vw; max-height: 92vh; background: #fff; border-radius: 8px; display: block; }

@media (max-width: 1000px) { .g-main, .g-viewer, .g2, .g3 { grid-template-columns: minmax(0, 1fr); } }
@media (max-width: 720px) { main { padding: 12px; } .bar-row { grid-template-columns: 80px 1fr 56px; } }
@media print {
  .tabs, .toolbar, .btn { display: none !important; }
  .panel[hidden] { display: block !important; }
  .card { break-inside: avoid; }
}
@media (forced-colors: active) { .bar-fill, .swatch { forced-color-adjust: none; } }


/* ---- Levels: alerts, chips, flagged cells (always an icon + text, never colour alone) ---- */
.alert { display: flex; gap: 10px; align-items: flex-start; padding: 9px 14px; border-radius: 8px;
         border: 1px solid var(--ring); border-left: 4px solid var(--muted); background: var(--surface); font-size: 13px; }
.alert .ico { font-size: 15px; line-height: 1.35; flex: none; }
.alert .body { min-width: 0; overflow-wrap: anywhere; }
.alert.fail { border-left-color: var(--fail); background: var(--fail-soft); }
.alert.fail .ico { color: var(--fail); }
.alert.warn { border-left-color: var(--warn); background: var(--warn-soft); }
.alert.warn .ico { color: var(--warn); }
.alert.ok { border-left-color: var(--c-high); background: var(--ok-soft); }
.alert.ok .ico { color: var(--ok-ink); }
.alert.info { border-left-color: var(--accent); background: var(--accent-soft); }
.alert.info .ico { color: var(--accent); }
.stack.tight { gap: 6px; }

.chip { display: inline-flex; align-items: center; gap: 4px; padding: 0 8px; border-radius: 999px; font-size: 11.5px;
        font-weight: 600; line-height: 20px; border: 1px solid var(--ring); white-space: nowrap; background: var(--page); color: var(--text-2); }
.chip.ok { color: var(--ok-ink); background: var(--ok-soft); }
.chip.warn { color: var(--warn); background: var(--warn-soft); }
.chip.fail { color: var(--fail); background: var(--fail-soft); }
.chip.info { color: var(--accent); background: var(--accent-soft); }
.chip.cat { font-weight: 500; }
.chip.type.DEL { color: var(--c-del); border-color: var(--c-del); background: transparent; }
.chip.type.DUP { color: var(--c-dup); border-color: var(--c-dup); background: transparent; }
.chip.conf.HIGH { color: var(--ok-ink); background: var(--ok-soft); }
.chip.conf.MEDIUM { color: var(--warn); background: var(--warn-soft); }
.chip.conf.LOW, .chip.conf.NA { color: var(--text-2); background: var(--page); }
.chip-btn { appearance: none; font: inherit; font-size: 11.5px; font-weight: 600; display: inline-flex; gap: 4px; align-items: center;
            padding: 0 8px; line-height: 20px; border-radius: 999px; cursor: pointer; border: 1px solid var(--ring); background: var(--page); color: var(--text); }
.chip-btn:hover { background: var(--accent-soft); }
.chip-btn.warn { color: var(--warn); background: var(--warn-soft); }
.chip-btn.fail { color: var(--fail); background: var(--fail-soft); }
.chip-btn.link { color: var(--accent); }
.chips { display: flex; flex-wrap: wrap; gap: 6px; }
.count-badge { display: inline-block; margin-left: 6px; min-width: 18px; padding: 0 6px; border-radius: 999px; font-size: 11px; line-height: 18px;
               text-align: center; background: var(--warn); color: var(--page); font-weight: 700; }
button.btn.small-btn { padding: 1px 9px; font-size: 12px; margin-left: 6px; }

.tile.ok { border-left: 4px solid var(--c-high); }
.tile.warn { border-left: 4px solid var(--warn); background: var(--warn-soft); }
.tile.fail { border-left: 4px solid var(--fail); background: var(--fail-soft); }
.card.flagged { border-color: var(--warn); }
.tbl td.warn { background: var(--warn-soft); color: var(--warn); font-weight: 600; }
.tbl td.fail { background: var(--fail-soft); color: var(--fail); font-weight: 600; }
.tbl tbody tr.clickable { cursor: pointer; }
.tbl tbody tr.sel { outline: 2px solid var(--accent); outline-offset: -2px; background: var(--accent-soft); }
.note-flag { font-size: 11.5px; }

.chrom-strip { display: grid; gap: 10px; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); }
.chrom-box { border: 1px solid var(--ring); border-left: 4px solid var(--muted); border-radius: 8px; padding: 8px 12px; display: flex; flex-direction: column; gap: 4px; align-items: flex-start; }
.chrom-box .nm { font-weight: 600; }
.chrom-box.ok { border-left-color: var(--c-high); }
.chrom-box.warn { border-left-color: var(--warn); background: var(--warn-soft); }
.chrom-box.fail { border-left-color: var(--fail); background: var(--fail-soft); }
.chrom-box.info { border-left-color: var(--accent); }

.bar-fill.lv-ok { background: var(--c-high); }
.bar-fill.lv-warn { background: var(--c-medium); }
.bar-fill.lv-fail { background: var(--c-del); }
.swatch.conf.HIGH { background: var(--c-high); } .swatch.conf.MEDIUM { background: var(--c-medium); }
.swatch.conf.LOW { background: var(--c-low); } .swatch.conf.NA { background: var(--c-na); }
.plot.viewer { height: 760px; }
.prose p { margin: 0 0 8px; max-width: 90ch; }
.field.inline input[type="checkbox"] { width: 16px; height: 16px; accent-color: var(--accent); }
.card > header h3 { font-size: 14px; }
.card .toolbar:last-child { margin-bottom: 0; }
)---"


.CANOPE_REPORT_JS <- r"---(
/* CANOPE interactive report -- client-side app.
   Reads the JSON blob in #canope-data, builds one lazily-rendered panel per tab and draws
   every figure with the inlined plotly.js. No server and no network access are needed.
   Design follows the SPARKS report (h() DOM helper, tableView, theme switch, hash navigation). */
(function () {
  'use strict';

  // ------------------------------------------------------------------ data & utils
  var DATA = JSON.parse(document.getElementById('canope-data').textContent);
  var META = DATA.meta, SUM = DATA.summary, THR = DATA.thresholds;
  var S = DATA.samples;
  var CALLS = DATA.calls || { cols: [], data: {}, n: 0 };
  var PLOTS = DATA.plots || {};
  var WARNINGS = DATA.warnings || [];
  var SAMPLES = (S.data.sample || []).slice();
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
    try { if (val === undefined) return window.localStorage.getItem(key); window.localStorage.setItem(key, val); } catch (e) { /* storage may be unavailable */ }
    return null;
  }
  function tcol(T, c) { return (T && T.data && T.data[c]) || []; }
  function sval(s, col) { var i = SIDX[s]; return i === undefined ? null : tcol(S, col)[i]; }
  function flagsOf(s) { var f = DATA.flags && DATA.flags[s]; return Array.isArray(f) ? f : []; }
  function hasFlag(s, re) { return flagsOf(s).some(function (f) { return re.test(f.text); }); }
  function pretty(c) { return String(c).replace(/_/g, ' '); }

  // ------------------------------------------------------------------ theme & palette
  var CHROME = {
    light: { surface: '#fcfcfb', text: '#0b0b0b', secondary: '#52514e', muted: '#898781', grid: '#e1e0d9', axis: '#c3c2b7', mid: '#f0efec' },
    dark: { surface: '#1a1a19', text: '#f0efec', secondary: '#c3c2b7', muted: '#898781', grid: '#2c2c2a', axis: '#383835', mid: '#383835' }
  };
  var CAT = {
    light: ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'],
    dark: ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767']
  };
  // meaning-bearing colours (always paired with a symbol / text, never colour alone)
  var SEM = {
    light: { DEL: '#e34948', DUP: '#2a78d6', HIGH: '#1baf7a', MEDIUM: '#eda100', LOW: '#898781', NA: '#c3c2b7',
             fail: '#d13b3a', warn: '#e09a00', ok: '#1baf7a', info: '#2a78d6', test: '#2a78d6', ref: '#c3c2b7', flank: '#898781' },
    dark: { DEL: '#e66767', DUP: '#3987e5', HIGH: '#199e70', MEDIUM: '#c98500', LOW: '#898781', NA: '#555552',
            fail: '#e66767', warn: '#c98500', ok: '#199e70', info: '#3987e5', test: '#3987e5', ref: '#4a4a46', flank: '#898781' }
  };
  var SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'];
  var DIV = ['#e34948', '#f0a3a2', '#f0efec', '#9ec5f4', '#2a78d6'];
  var FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif";

  var themePref = store('canope-theme') || 'auto';
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

  // ------------------------------------------------------------------ plotly helpers
  var registry = [];
  var themeListeners = [];
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
  // build(theme) -> {data, layout, config, onClick}
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
    themeListeners.forEach(function (cb) { cb(); });
  }
  function showNote(div, text) {
    if (window.Plotly) window.Plotly.purge(div);
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

  // ------------------------------------------------------------------ small widgets
  var ICON = { fail: '✖', warn: '⚠', info: 'ℹ', ok: '✓' };
  var LEVEL_LABEL = { fail: 'Fail', warn: 'Warning', info: 'Note', ok: 'OK' };
  function chip(level, text) {
    return h('span', { class: 'chip ' + level },
      h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[level] || '' }), text === undefined ? LEVEL_LABEL[level] : text);
  }
  function card(title, hint, body, cls) {
    return h('section', { class: 'card' + (cls ? ' ' + cls : '') },
      h('header', null, h('h3', { text: title }), hint ? h('span', { class: 'muted small hint', text: hint }) : null), body);
  }
  function tile(value, label, level) {
    return h('div', { class: 'tile' + (level ? ' ' + level : '') }, h('div', { class: 'v', text: value }), h('div', { class: 'l', text: label }));
  }
  function alertBox(level, title, text, extra) {
    return h('div', { class: 'alert ' + level, role: level === 'fail' ? 'alert' : 'note' },
      h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[level] }),
      h('div', { class: 'body' },
        title ? h('strong', { text: title }) : null, title && text ? ' ' : null, text || null, extra || null));
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

  // ------------------------------------------------------------------ global state (sample / call selection)
  var state = { sample: '', call: null, listeners: [] };
  var sampleSelects = [];
  function onState(cb) { state.listeners.push(cb); }
  function notify() {
    state.listeners.forEach(function (cb) { try { cb(); } catch (e) { if (window.console) console.error(e); } });
  }
  function setSample(s) {
    state.sample = s || '';
    sampleSelects.forEach(function (el) { el.value = state.sample; });
    notify();
  }
  function sampleSelect(label) {
    var sel = h('select', { 'aria-label': 'Sample' }, h('option', { value: '', text: 'All samples' }),
      SAMPLES.map(function (s) { return h('option', { value: s, text: s }); }));
    sel.value = state.sample;
    sel.addEventListener('change', function () { setSample(sel.value); });
    sampleSelects.push(sel);
    return field(label || 'Sample', sel);
  }

  // ------------------------------------------------------------------ tables
  function colType(T, c) {
    var arr = tcol(T, c);
    for (var i = 0; i < arr.length; i++) { if (arr[i] !== null && arr[i] !== undefined) return typeof arr[i]; }
    return 'string';
  }
  // T = {cols, data: {col: [..]}, n}.
  // opts: filename, pageSize, render {col: fn(value,i)->node}, cellClass {col: fn(value,i)->class},
  //       labels {col: header}, sort [col,'asc'|'desc'], onRow fn(i), selected fn()->row index, note
  function tableView(host, T, opts) {
    opts = opts || {};
    var n = T.n, pageSize = opts.pageSize || 25, page = 0, filter = '', pred = null;
    var cols = asArray(T.cols);
    var sortCol = opts.sort ? opts.sort[0] : null, sortDir = opts.sort ? opts.sort[1] : 'asc';
    var types = {}; cols.forEach(function (c) { types[c] = colType(T, c); });
    var idx = range(n);
    var search = h('input', { type: 'search', placeholder: 'Filter rows…', 'aria-label': 'Filter rows' });
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
          return (num ? x - y : String(x).localeCompare(String(y))) * dir;
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
          row.appendChild(h('td', { class: ((types[c] === 'number' ? 'num ' : '') + (cc || '')).trim() || false }, r ? r(v, i) : fmt(v)));
        });
        tbody.appendChild(row);
      }
      count.textContent = fmtInt(idx.length) + (idx.length === n ? '' : ' of ' + fmtInt(n)) + ' rows' +
        (T.truncated ? ' (truncated from ' + fmtInt(T.n_total) + ')' : '');
      var pages = Math.max(1, Math.ceil(idx.length / pageSize));
      clear(pager);
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '‹', 'aria-label': 'Previous page', disabled: page === 0, onclick: function () { page--; draw(); } }));
      pager.appendChild(h('span', { class: 'small', text: (page + 1) + ' / ' + pages }));
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '›', 'aria-label': 'Next page', disabled: page >= pages - 1, onclick: function () { page++; draw(); } }));
    }
    var timer = null;
    search.addEventListener('input', function () {
      clearTimeout(timer);
      timer = setTimeout(function () { filter = search.value; page = 0; recompute(); }, 120);
    });
    recompute();
    return { setPredicate: function (fn) { pred = fn; page = 0; recompute(); }, refresh: recompute, redraw: draw, count: function () { return idx.length; } };
  }
  function simpleTable(T, opts) { var host = h('div'); tableView(host, T, opts); return host; }
  function subTable(T, cols) {
    var out = { cols: cols, data: {}, n: T.n };
    cols.forEach(function (c) { out.data[c] = tcol(T, c); });
    return out;
  }

  // ------------------------------------------------------------------ derived sample data
  var NUMERIC_FLAGS = {
    median_depth: function (v) { return isNum(v) && v < THR.qc_min_cov; },
    total_reads: function (v) { return isNum(v) && v < THR.qc_min_total_reads; },
    max_corr: function (v) { return isNum(v) && v < THR.qc_min_corr; },
    noise_z: function (v) { return isNum(v) && v > THR.qc_zscore; },
    n_refs: function (v) { return isNum(v) && v < 5; }
  };
  function sampleTable() {
    var T = { cols: asArray(S.cols).slice(), data: {}, n: S.n };
    T.cols.forEach(function (c) { T.data[c] = tcol(S, c); });
    var nf = SAMPLES.map(function (s) { return flagsOf(s).filter(function (f) { return f.level === 'fail' || f.level === 'warn'; }).length; });
    T.cols.splice(T.cols.indexOf('level') + 1, 0, 'flags'); T.data.flags = nf;
    Object.keys(DATA.sample_attrs || {}).forEach(function (k) {
      if (T.data[k]) return;
      T.cols.push(k); T.data[k] = asArray(DATA.sample_attrs[k]);
    });
    return T;
  }
  function callsOf(s) {
    var out = [], sm = tcol(CALLS, 'sample');
    for (var i = 0; i < CALLS.n; i++) if (!s || sm[i] === s) out.push(i);
    return out;
  }
  var CALL_ROW = {};
  tcol(CALLS, 'id').forEach(function (id, i) { CALL_ROW[id] = i; });
  var CONF_RANK = { HIGH: 3, MEDIUM: 2, LOW: 1 };
  function callOrder(rows) {
    var q = tcol(CALLS, 'q_some'), cf = tcol(CALLS, 'confidence');
    return rows.slice().sort(function (a, b) {
      var d = (CONF_RANK[cf[b]] || 0) - (CONF_RANK[cf[a]] || 0);
      if (d) return d;
      return (isNum(q[b]) ? q[b] : -1) - (isNum(q[a]) ? q[a] : -1);
    });
  }
  function typeChip(v) {
    if (v !== 'DEL' && v !== 'DUP') return fmt(v);
    return h('span', { class: 'chip type ' + v }, h('span', { class: 'ico', 'aria-hidden': 'true', text: v === 'DEL' ? '▼' : '▲' }), v);
  }
  function confChip(v) {
    return h('span', { class: 'chip conf ' + (v || 'NA') }, v === 'NA' || !v ? 'n/a' : v);
  }
  function goto(id) { activate(id, true); }
  function openCall(id) { state.call = id; var r = CALL_ROW[id]; if (r !== undefined) { var s = tcol(CALLS, 'sample')[r]; if (state.sample && state.sample !== s) setSample(s); } notify(); goto('viewer'); }

  var TABS = [];
  function tab(id, label, available, build) {
    var t = { id: id, label: label, available: available, build: build, built: false };
    TABS.push(t); return t;
  }

  // ================================================================== OVERVIEW
  function levelCount(l) { return WARNINGS.filter(function (w) { return w.level === l; }).length; }
  var CHR = DATA.chromosomes;

  function chromStrip() {
    var st = CHR.status, wrap = h('div', { class: 'chrom-strip' });
    asArray(tcol(st, 'class')).forEach(function (k, i) {
      var inp = st.data.in_input[i], ana = st.data.analysed[i], reason = st.data.reason[i];
      var level = ana ? 'ok' : inp ? (k === 'Y' ? 'warn' : 'fail') : 'info';
      wrap.appendChild(h('div', { class: 'chrom-box ' + level },
        h('div', { class: 'nm', text: st.data.label[i] }),
        chip(level, ana ? 'analysed' : inp ? 'NOT analysed' : 'not in BED'),
        reason ? h('div', { class: 'small muted', text: reason }) : null));
    });
    return wrap;
  }

  function globalAlerts(limit) {
    var box = h('div', { class: 'stack tight' });
    var fails = WARNINGS.filter(function (w) { return w.level === 'fail'; });
    fails.slice(0, limit).forEach(function (w) { box.appendChild(alertBox('fail', w.category + ':', w.text)); });
    var rest = WARNINGS.filter(function (w) { return w.level === 'warn'; });
    if (fails.length > limit) rest = rest.concat(fails.slice(limit));
    if (rest.length) {
      box.appendChild(alertBox('warn', rest.length + ' more warning' + (rest.length > 1 ? 's' : ''), 'in sample metrics, chromosomes or calls. ',
        h('button', { class: 'btn small-btn', type: 'button', text: 'Review all warnings', onclick: function () { goto('warnings'); } })));
    }
    if (!fails.length && !rest.length) box.appendChild(alertBox('ok', 'No warnings.', 'All sample metrics are within their thresholds and every chromosome class was analysed.'));
    return box;
  }

  var tOverview = tab('overview', 'Overview', function () { return true; }, function (root) {
    var nF = levelCount('fail'), nW = levelCount('warn');
    var del = 0, dup = 0, high = 0;
    tcol(CALLS, 'type').forEach(function (t) { if (t === 'DEL') del++; else if (t === 'DUP') dup++; });
    tcol(CALLS, 'confidence').forEach(function (c) { if (c === 'HIGH') high++; });
    var both = SUM.both_chromosome_classes;
    var tiles = h('div', { class: 'tiles' },
      tile(SUM.samples_called + ' / ' + SUM.samples_requested, 'samples analysed / requested', SUM.samples_called < SUM.samples_requested ? 'warn' : null),
      tile(fmtInt(SUM.calls), 'CNV calls (' + del + ' DEL, ' + dup + ' DUP)'),
      tile(fmtInt(high), 'HIGH-confidence calls'),
      tile(SUM.samples_with_calls, 'samples with ≥ 1 call'),
      tile(fmtInt(SUM.targets), 'targets analysed'),
      tile(both === false ? 'Partial' : both === true ? 'Both' : 'n/a', 'autosomes + gonosomes', both === false ? 'fail' : both === true ? 'ok' : null),
      tile(nF + nW, 'warnings (' + nF + ' fail, ' + nW + ' warn)', nF ? 'fail' : nW ? 'warn' : 'ok'));
    root.appendChild(h('div', { class: 'stack' },
      globalAlerts(6), tiles,
      card('Chromosome coverage of this run', META.mode_label, h('div', null, chromStrip(),
        h('p', { class: 'note', text: 'One run analyses one chromosome mode (A = autosomes, XX = chrX with female references, XY = chrX/chrY with male references). Run the other modes to cover the whole genome.' }))),
      h('div', { class: 'grid g-main' }, card('Calls per sample', 'click a bar to list that sample’s calls', overviewBars()), card('Sample health', null, healthBars())),
      card('Run', null, runKv())));
  });
  tOverview.badge = null;

  function overviewBars() {
    var host = h('div'), box = plotBox('');
    var order = selectOf([{ value: 'input', text: 'Input order' }, { value: 'calls', text: 'Most calls first' }, { value: 'flags', text: 'Most flags first' }], 'calls', null, 'Sample order');
    host.appendChild(h('div', { class: 'toolbar' }, field('Order', order)));
    host.appendChild(box);
    if (!SAMPLES.length) { host.appendChild(emptyNote('No samples.')); return host; }
    function draw() {
      var ids = range(S.n), nf = SAMPLES.map(function (s) { return flagsOf(s).filter(function (f) { return f.level !== 'info'; }).length; });
      if (order.value === 'calls') ids.sort(function (a, b) { return S.data.n_cnvs[b] - S.data.n_cnvs[a]; });
      else if (order.value === 'flags') ids.sort(function (a, b) { return nf[b] - nf[a] || S.data.n_cnvs[b] - S.data.n_cnvs[a]; });
      plot(box, function (t) {
        var c = SEM[t], x = ids.map(function (i) { return SAMPLES[i]; });
        var na = ids.map(function (i) { return Math.max(0, S.data.n_cnvs[i] - S.data.n_high[i] - S.data.n_medium[i] - S.data.n_low[i]); });
        var mk = function (name, arr, col) { return { type: 'bar', name: name, x: x, y: arr, marker: { color: col }, hovertemplate: '%{x}<br>' + name + ': %{y}<extra></extra>' }; };
        var data = [mk('HIGH', ids.map(function (i) { return S.data.n_high[i]; }), c.HIGH), mk('MEDIUM', ids.map(function (i) { return S.data.n_medium[i]; }), c.MEDIUM),
                    mk('LOW', ids.map(function (i) { return S.data.n_low[i]; }), c.LOW)];
        if (na.some(function (v) { return v > 0; })) data.push(mk('unscored', na, c.NA));
        var L = baseLayout(t, { ylab: 'CNV calls', barmode: 'stack', legend: true, xaxis: { type: 'category', tickangle: -45 }, margin: { l: 56, r: 16, t: 16, b: 90 } });
        return { data: data, layout: L, onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) { setSample(p.x); goto('calls'); } } };
      });
    }
    order.addEventListener('change', draw);
    draw();
    return host;
  }

  function healthBars() {
    var counts = { ok: 0, warn: 0, fail: 0 };
    tcol(S, 'level').forEach(function (l) { counts[l] = (counts[l] || 0) + 1; });
    var total = Math.max(1, S.n);
    var host = h('div', { class: 'stack' }), bars = h('div', { class: 'bars' });
    [['ok', 'No flags'], ['warn', 'Warnings'], ['fail', 'Fail']].forEach(function (p) {
      var col = SEM[resolved()][p[0]];
      bars.appendChild(h('div', { class: 'bar-row' }, h('span', { class: 'lab', text: p[1] }),
        h('div', { class: 'bar-track' }, h('div', { class: 'bar-fill lv-' + p[0], style: { width: (100 * counts[p[0]] / total) + '%' } })),
        h('span', { class: 'small', text: counts[p[0]] + ' sample' + (counts[p[0]] === 1 ? '' : 's') })));
    });
    host.appendChild(bars);
    var st = {}; tcol(S, 'status').forEach(function (v) { st[v] = (st[v] || 0) + 1; });
    host.appendChild(kv(Object.keys(st).map(function (k) { return [k === 'called' ? 'analysed' : k, String(st[k])]; })));
    var flagged = SAMPLES.filter(function (s) { return sval(s, 'level') !== 'ok'; });
    if (flagged.length) {
      host.appendChild(h('div', { class: 'note', text: 'Samples with flags:' }));
      host.appendChild(h('div', { class: 'chips' }, flagged.slice(0, 40).map(function (s) {
        return h('button', { class: 'chip-btn ' + sval(s, 'level'), type: 'button', title: 'Open sample', onclick: function () { setSample(s); goto('samples'); } },
          h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[sval(s, 'level')] }), s);
      }), flagged.length > 40 ? h('span', { class: 'muted small', text: '+' + (flagged.length - 40) + ' more' }) : null));
    }
    return host;
  }

  function runKv() {
    var p = DATA.params || {};
    return kv([
      ['Report generated', META.generated], ['CANOPE version', META.version], ['HMM engine', META.engine_label],
      ['Chromosome mode', META.mode_label], ['Genome build', META.genome_build],
      ['Coverage input', p.coverage_input || p.reads_file || null],
      ['Reference selection', p.min_correlation !== undefined ? ('min correlation ' + p.min_correlation + ', up to ' + (p.num_refs || 30) + ' references') : null],
      ['Legacy workspace', META.legacy_workspace ? 'yes – this workspace has no run summary (older CANOPE); some warnings cannot be computed' : null]
    ]);
  }

  // ================================================================== WARNINGS
  tab('warnings', 'Warnings', function () { return true; }, function (root) {
    var host = h('div', { class: 'stack tight' });
    var lvl = selectOf([{ value: '', text: 'All levels' }, { value: 'fail', text: 'Fail' }, { value: 'warn', text: 'Warning' }, { value: 'info', text: 'Note' }], '', null, 'Level');
    var cats = uniq(WARNINGS.map(function (w) { return w.category; }));
    var cat = selectOf([{ value: '', text: 'All categories' }].concat(cats.map(function (c) { return { value: c, text: c }; })), '', null, 'Category');
    var q = h('input', { type: 'search', placeholder: 'Search warnings…', 'aria-label': 'Search warnings' });
    var count = h('span', { class: 'muted small' });
    var sel = sampleSelect('Sample');
    var csvBtn = h('button', { class: 'btn', type: 'button', text: 'Download CSV', onclick: function () {
      download('canope_warnings.csv', toCSV(['level', 'category', 'sample', 'text'], function (c, i) { return WARNINGS[i][c]; }, WARNINGS.length));
    } });
    root.appendChild(h('div', { class: 'toolbar' }, field('Level', lvl), field('Category', cat), sel, field('Search', q), count, csvBtn));
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
      rows.sort(function (a, b) { return (a.level === b.level ? 0 : { fail: 0, warn: 1, info: 2 }[a.level] - { fail: 0, warn: 1, info: 2 }[b.level]); });
      count.textContent = rows.length + ' of ' + WARNINGS.length;
      if (!rows.length) { host.appendChild(alertBox(WARNINGS.length ? 'info' : 'ok', WARNINGS.length ? 'No warning matches the filters.' : 'No warnings were raised.', null)); return; }
      rows.slice(0, 400).forEach(function (w) {
        host.appendChild(alertBox(w.level, null, null,
          [h('span', { class: 'chip cat', text: w.category }), ' ',
           w.sample ? h('button', { class: 'chip-btn link', type: 'button', title: 'Open sample', onclick: function () { setSample(w.sample); goto('samples'); }, text: w.sample }) : null,
           w.sample ? ' ' : null, w.text]));
      });
      if (rows.length > 400) host.appendChild(h('div', { class: 'note', text: 'Showing the first 400 warnings; use the filters or download the CSV.' }));
    }
    [lvl, cat].forEach(function (el) { el.addEventListener('change', draw); });
    q.addEventListener('input', draw);
    onState(draw);
    draw();
  });

  // ================================================================== SAMPLES
  var tSamples = tab('samples', 'Samples', function () { return true; }, function (root) {
    var T = sampleTable();
    var tv = null, detail = h('div'), scatterBox = plotBox('');
    var labels = { n_cnvs: 'CNVs', n_del: 'DEL', n_dup: 'DUP', n_high: 'HIGH', n_medium: 'MEDIUM', n_low: 'LOW', median_depth: 'median depth', total_reads: 'total reads',
                   max_corr: 'best corr.', noise_z: 'noise z', n_refs: 'refs', ref_cor: 'ref corr.', n_targets: 'targets', n_dropped: 'dropped', sex_inferred: 'sex (inferred)' };
    var cellClass = {};
    Object.keys(NUMERIC_FLAGS).forEach(function (k) { cellClass[k] = function (v) { return NUMERIC_FLAGS[k](v) ? 'warn' : ''; }; });
    cellClass.status = function (v) { return v === 'failed' || v === 'excluded' ? 'fail' : ''; };
    cellClass.n_cnvs = function (v, i) { return hasFlag(SAMPLES[i], /^Unusually many/) ? 'warn' : ''; };
    cellClass.n_dropped = function (v, i) { return hasFlag(SAMPLES[i], /removed before calling/) ? 'warn' : ''; };
    cellClass.sex_inferred = function (v, i) { return hasFlag(SAMPLES[i], /sex|modechrom/i) ? 'warn' : ''; };
    cellClass.flags = function (v, i) { return v > 0 ? (sval(SAMPLES[i], 'level') === 'fail' ? 'fail' : 'warn') : ''; };
    var render = {
      level: function (v) { return chip(v, LEVEL_LABEL[v]); },
      status: function (v) { return v === 'called' ? 'analysed' : v; },
      sample: function (v) { return h('strong', { text: v }); }
    };
    var host = h('div');
    root.appendChild(h('div', { class: 'stack' },
      h('div', { class: 'toolbar' }, sampleSelect('Sample'),
        h('span', { class: 'note', text: 'Cells shaded amber/red are outside the QC thresholds (' + THR.qc_min_cov + ' median depth, ' + fmtInt(THR.qc_min_total_reads) + ' reads, correlation ≥ ' + THR.qc_min_corr + ', noise z ≤ ' + THR.qc_zscore + ').' })),
      card('Samples', 'click a row to select it', host),
      card('Selected sample', null, detail),
      card('Depth vs. best correlation', 'dashed lines = QC thresholds; click a point to select the sample', scatterBox)));
    tv = tableView(host, T, { filename: 'canope_samples', pageSize: 20, sort: ['flags', 'desc'], labels: labels, render: render, cellClass: cellClass,
      selected: function () { return state.sample ? SIDX[state.sample] : -1; },
      onRow: function (i) { setSample(SAMPLES[i] === state.sample ? '' : SAMPLES[i]); } });
    function drawDetail() {
      clear(detail);
      var s = state.sample;
      if (!s) {
        var flagged = SAMPLES.filter(function (x) { return sval(x, 'level') !== 'ok'; });
        detail.appendChild(alertBox(flagged.length ? 'warn' : 'ok', flagged.length + ' of ' + SAMPLES.length + ' samples have warnings.', 'Select one to read its flags, or use the filter above.'));
        return;
      }
      var fl = flagsOf(s);
      var lv = sval(s, 'level');
      detail.appendChild(h('div', { class: 'toolbar' }, h('h3', { text: s }), chip(lv, LEVEL_LABEL[lv]),
        h('button', { class: 'btn', type: 'button', text: 'Show its CNV calls', onclick: function () { goto('calls'); } }),
        h('button', { class: 'btn', type: 'button', text: 'Open best call in viewer', disabled: !callOrder(callsOf(s)).filter(function (r) { return tcol(CALLS, 'has_plot')[r]; }).length,
          onclick: function () { var r = callOrder(callsOf(s)).filter(function (x) { return tcol(CALLS, 'has_plot')[x]; })[0]; if (r !== undefined) openCall(tcol(CALLS, 'id')[r]); } })));
      var list = h('div', { class: 'stack tight' });
      if (!fl.length) list.appendChild(alertBox('ok', null, 'No warnings for this sample: all metrics are within their thresholds.'));
      fl.forEach(function (f) { list.appendChild(alertBox(f.level, null, f.text)); });
      detail.appendChild(list);
      var pairs = [['Status', sval(s, 'status') === 'called' ? 'analysed' : sval(s, 'status')], ['CNV calls', sval(s, 'n_cnvs') + ' (' + sval(s, 'n_del') + ' DEL, ' + sval(s, 'n_dup') + ' DUP; ' +
        sval(s, 'n_high') + ' HIGH / ' + sval(s, 'n_medium') + ' MEDIUM / ' + sval(s, 'n_low') + ' LOW)'], ['Median target depth', fmt(sval(s, 'median_depth'))],
        ['Total reads', fmt(sval(s, 'total_reads'))], ['Best Spearman correlation', fmt(sval(s, 'max_corr'))], ['Median reference correlation', fmt(sval(s, 'ref_cor'))],
        ['References used', fmt(sval(s, 'n_refs'))], ['Targets used', fmt(sval(s, 'n_targets'))], ['Sex (inferred from chrX/chrY)', sval(s, 'sex_inferred')]];
      Object.keys(DATA.sample_attrs || {}).forEach(function (k) { pairs.push([pretty(k), asArray(DATA.sample_attrs[k])[SIDX[s]]]); });
      detail.appendChild(h('div', { style: { marginTop: '10px' } }, kv(pairs.map(function (p) { return [p[0], p[1] === null ? '' : String(p[1])]; }))));
    }
    function drawScatter() {
      var depth = tcol(S, 'median_depth'), corr = tcol(S, 'max_corr');
      var ok = range(S.n).filter(function (i) { return isNum(depth[i]) && isNum(corr[i]); });
      if (ok.length < 2) { showNote(scatterBox, 'Needs at least two analysed samples.'); return; }
      plot(scatterBox, function (t) {
        var c = SEM[t], data = [];
        ['ok', 'warn', 'fail'].forEach(function (lv) {
          var ids = ok.filter(function (i) { return S.data.level[i] === lv; });
          if (!ids.length) return;
          data.push({ type: 'scatter', mode: 'markers', name: LEVEL_LABEL[lv] === 'OK' ? 'No flags' : LEVEL_LABEL[lv], x: ids.map(function (i) { return depth[i]; }), y: ids.map(function (i) { return corr[i]; }),
            text: ids.map(function (i) { return SAMPLES[i]; }), customdata: ids.map(function (i) { return SAMPLES[i]; }),
            marker: { color: c[lv], size: ids.map(function (i) { return 8 + Math.min(14, 2 * Math.sqrt(S.data.n_cnvs[i])); }), symbol: lv === 'ok' ? 'circle' : lv === 'warn' ? 'diamond' : 'x', line: { width: 1, color: CHROME[t].surface } },
            hovertemplate: '<b>%{text}</b><br>depth %{x}<br>best corr. %{y}<extra>' + LEVEL_LABEL[lv] + '</extra>' });
        });
        if (state.sample && SIDX[state.sample] !== undefined && ok.indexOf(SIDX[state.sample]) >= 0) {
          var k = SIDX[state.sample];
          data.push({ type: 'scatter', mode: 'markers', name: 'selected', x: [depth[k]], y: [corr[k]], hoverinfo: 'skip',
            marker: { color: 'rgba(0,0,0,0)', size: 24, line: { color: CHROME[t].text, width: 2 } } });
        }
        var L = baseLayout(t, { xlab: 'median target depth (log scale)', ylab: 'best Spearman correlation with another sample', legend: true, xaxis: { type: 'log' },
          shapes: [{ type: 'line', xref: 'x', yref: 'paper', x0: THR.qc_min_cov, x1: THR.qc_min_cov, y0: 0, y1: 1, line: { color: c.fail, dash: 'dash', width: 1 } },
                   { type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: THR.qc_min_corr, y1: THR.qc_min_corr, line: { color: c.fail, dash: 'dash', width: 1 } }] });
        return { data: data, layout: L, onClick: function (ev) { var p = ev.points && ev.points[0]; if (p && p.customdata) setSample(p.customdata); } };
      });
    }
    onState(function () { tv.redraw(); drawDetail(); drawScatter(); });
    drawDetail(); drawScatter();
  });
  tSamples.badge = function () {
    var n = tcol(S, 'level').filter(function (l) { return l !== 'ok'; }).length;
    return n ? h('span', { class: 'count-badge', text: n, title: n + ' sample(s) with warnings' }) : null;
  };

  // ================================================================== QC METRICS
  var QC_METRICS = [
    { col: 'median_depth', label: 'Median target depth', thr: function () { return THR.qc_min_cov; }, dir: 'below', hint: 'below the line = low coverage' },
    { col: 'total_reads', label: 'Total reads on the analysed targets', thr: function () { return THR.qc_min_total_reads; }, dir: 'below', hint: 'below the line = too few reads' },
    { col: 'max_corr', label: 'Best Spearman correlation with another sample', thr: function () { return THR.qc_min_corr; }, dir: 'below', hint: 'below the line = no comparable reference' },
    { col: 'noise_z', label: 'Sample noise (robust z of the exon-count dispersion)', thr: function () { return THR.qc_zscore; }, dir: 'above', hint: 'above the line = noisy sample' },
    { col: 'n_refs', label: 'Reference samples used', thr: function () { return 5; }, dir: 'below', hint: 'below the line = weak reference set' },
    { col: 'n_cnvs', label: 'CNV calls', thr: function () { return null; }, dir: null, hint: 'outliers suggest a noisy sample / batch effect' }
  ];

  function qcTypeChip(v) {
    var lv = v === 'Missing Chromosome' ? 'fail' : v === 'Chromosome not analysed' ? 'warn' : v === 'Whole sample' ? 'warn' : 'info';
    return h('span', { class: 'chip ' + lv }, h('span', { class: 'ico', 'aria-hidden': 'true', text: ICON[lv] }), v);
  }

  function metricBars(m, orderSel) {
    var box = plotBox('short');
    var vals = tcol(S, m.col);
    var ids = range(S.n).filter(function (i) { return isNum(vals[i]); });
    if (!ids.length) return null;
    function draw() {
      var o = orderSel.value, order = ids.slice();
      if (o === 'asc') order.sort(function (a, b) { return vals[a] - vals[b]; });
      else if (o === 'desc') order.sort(function (a, b) { return vals[b] - vals[a]; });
      plot(box, function (t) {
        var c = SEM[t], thr = m.thr();
        var bad = function (i) { return thr !== null && (m.dir === 'below' ? vals[i] < thr : vals[i] > thr); };
        var shapes = thr !== null ? [{ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: thr, y1: thr, line: { color: c.fail, dash: 'dash', width: 1.5 } }] : [];
        var d = [{ type: 'bar', x: order.map(function (i) { return SAMPLES[i]; }), y: order.map(function (i) { return vals[i]; }),
          marker: { color: order.map(function (i) { return bad(i) ? c.warn : CAT[t][0]; }), pattern: { shape: order.map(function (i) { return bad(i) ? '/' : ''; }) } },
          customdata: order.map(function (i) { return bad(i) ? ' ⚠ outside threshold' : ''; }),
          hovertemplate: '%{x}<br>' + m.label + ': %{y}%{customdata}<extra></extra>' }];
        // linear axis anchored at zero and always containing the threshold line
        var lo = Math.min(0, Math.min.apply(null, ids.map(function (i) { return vals[i]; }))), hi = Math.max(Math.max.apply(null, ids.map(function (i) { return vals[i]; })), thr !== null ? thr : -Infinity);
        var L = baseLayout(t, { ylab: '', yaxis: { type: 'linear', range: [lo < 0 ? lo * 1.1 : 0, hi > 0 ? hi * 1.08 : 1] }, xaxis: { type: 'category', tickangle: -45, showticklabels: order.length <= 60 },
          shapes: shapes, margin: { l: 60, r: 16, t: 10, b: order.length <= 60 ? 80 : 20 } });
        return { data: d, layout: L, onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) setSample(p.x); } };
      });
    }
    orderSel.addEventListener('change', draw);
    draw();
    var nbad = m.thr() === null ? 0 : ids.filter(function (i) { return m.dir === 'below' ? vals[i] < m.thr() : vals[i] > m.thr(); }).length;
    return card(m.label, (m.thr() !== null ? 'threshold ' + fmt(m.thr()) + ' – ' : '') + m.hint + (nbad ? ' – ' + nbad + ' sample(s) flagged' : ''), box, nbad ? 'flagged' : '');
  }

  tab('qc', 'QC metrics', function () { return true; }, function (root) {
    var stack = h('div', { class: 'stack' });
    var orderSel = selectOf([{ value: 'input', text: 'Input order' }, { value: 'asc', text: 'Lowest first' }, { value: 'desc', text: 'Highest first' }], 'input', null, 'Order');
    stack.appendChild(h('div', { class: 'toolbar' }, field('Sort bars', orderSel), sampleSelect('Highlight in table'),
      h('span', { class: 'note', text: 'Bars hatched in amber breach the QC threshold (dashed red line). Click a bar to select the sample.' })));
    var grid = h('div', { class: 'grid g2' });
    QC_METRICS.forEach(function (m) { var c = metricBars(m, orderSel); if (c) grid.appendChild(c); });
    stack.appendChild(grid);

    var Q = DATA.qc;
    if (!Q) {
      stack.appendChild(card('Exon-level QC flags', null, emptyNote('No QC-metrics file was provided to the report (run with qc_metrics enabled to get exon-level flags).')));
    } else {
      var host = h('div'), typeBox = plotBox('short');
      var tcounts = tcol(Q.types, 'Freq'), tnames = tcol(Q.types, 'Type');
      stack.appendChild(h('div', { class: 'grid g-main' },
        card('QC flags by type', Q.n_total + ' flags in total', tnames.length ? typeBox : emptyNote('No QC flags: every exon and sample passed.')),
        card('Genes with the most QC flags', 'top 60', Q.genes && Q.genes.n ? simpleTable(Q.genes, { filename: 'canope_qc_genes', pageSize: 10 }) : emptyNote('No exon-level flags.'))));
      if (tnames.length) plot(typeBox, function (t) {
        return { data: [{ type: 'bar', orientation: 'h', y: tnames, x: tcounts, marker: { color: CAT[t][3] }, hovertemplate: '%{y}: %{x}<extra></extra>' }],
          layout: baseLayout(t, { xlab: 'flags', yaxis: { autorange: 'reversed' }, margin: { l: 170, r: 16, t: 10, b: 40 } }) };
      });
      var tv = tableView(host, Q.rows, { filename: 'canope_qc_flags', pageSize: 20, render: { Type: qcTypeChip },
        note: Q.capped ? 'The table is capped at ' + fmtInt(Q.rows.n) + ' rows in this report; the full list is in the QC TSV file written by the pipeline.' : null });
      var sm = Q.rows.data.Sample;
      var apply = function () { tv.setPredicate(state.sample ? function (i) { return sm[i] === state.sample || sm[i] === 'All'; } : null); };
      onState(apply); apply();
      stack.appendChild(card('QC flags (samples, chromosomes, exons)', 'filtered by the selected sample', host));
    }
    root.appendChild(stack);
  });

  // ================================================================== CHROMOSOMES
  tab('chromosomes', 'Chromosomes', function () { return true; }, function (root) {
    var stack = h('div', { class: 'stack' });
    var both = SUM.both_chromosome_classes;
    if (both === false) stack.appendChild(alertBox('fail', 'Autosomes and gonosomes were NOT both analysed in this run.',
      'Results cover only part of the genome (mode: ' + META.mode_label + '). “No call” on the missing chromosomes does not mean “no CNV”.'));
    else if (both === true) stack.appendChild(alertBox('ok', 'Autosomes and chrX were both analysed.', null));
    else stack.appendChild(alertBox('info', 'Chromosome coverage is unknown.', 'This workspace comes from an older CANOPE version without a run summary.'));
    asArray(CHR.messages).forEach(function (m) { stack.appendChild(alertBox('warn', null, m)); });

    stack.appendChild(card('Chromosome classes', null, chromStrip()));

    var inv = CHR.inventory;
    var invHost = h('div');
    tableView(invHost, inv, { filename: 'canope_chromosomes', pageSize: 30, labels: { targets_in_bed: 'targets in BED', targets_analysed: 'targets analysed' },
      render: { analysed: function (v) { return v ? chip('ok', 'analysed') : chip('warn', 'not analysed'); } },
      cellClass: { targets_analysed: function (v, i) { return inv.data.analysed[i] ? '' : 'warn'; } } });
    var barBox = plotBox('short');
    stack.appendChild(h('div', { class: 'grid g2' }, card('Targets per chromosome', 'in the input BED vs. analysed in this run', barBox), card('Chromosome inventory', null, invHost)));
    plot(barBox, function (t) {
      var c = tcol(inv, 'chromosome');
      return { data: [{ type: 'bar', name: 'in BED', x: c, y: tcol(inv, 'targets_in_bed'), marker: { color: SEM[t].NA } },
                      { type: 'bar', name: 'analysed', x: c, y: tcol(inv, 'targets_analysed'), marker: { color: CAT[t][0] } }],
        layout: baseLayout(t, { ylab: 'targets', barmode: 'overlay', legend: true, xaxis: { type: 'category', tickangle: -45 }, margin: { l: 56, r: 16, t: 10, b: 80 } }) };
    });

    var D = CHR.depth;
    if (D && D.chromosomes && asArray(D.chromosomes).length) {
      var hbox = plotBox('tall');
      stack.appendChild(card('Chromosome depth relative to the autosomes', 'median target depth / sample autosomal median (1 ≈ diploid autosomes, ~0.5 chrX in males, ~0 chrY in females)', hbox));
      plot(hbox, function (t) {
        var dz = asArray(D.z).map(asArray);
        return { data: [{ type: 'heatmap', x: asArray(D.samples), y: asArray(D.chromosomes), z: dz, zmin: 0, zmax: 2, colorscale: [[0, DIV[0]], [0.25, DIV[1]], [0.5, DIV[2]], [0.75, DIV[3]], [1, DIV[4]]],
          hovertemplate: '%{x}<br>%{y}: %{z:.2f}<extra></extra>', colorbar: { title: { text: 'ratio' }, tickfont: { color: CHROME[t].muted } } }],
          layout: baseLayout(t, { xaxis: { type: 'category', tickangle: -45, showticklabels: asArray(D.samples).length <= 80 }, yaxis: { type: 'category', autorange: 'reversed' }, margin: { l: 70, r: 16, t: 10, b: 90 } }),
          onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) setSample(p.x); } };
      });
    } else {
      stack.appendChild(card('Chromosome depth relative to the autosomes', null, emptyNote(
        both === null || both === undefined ? 'Not available for this workspace.' : 'Depth per chromosome is only computed on the chromosomes that were analysed.')));
    }

    var X = CHR.sex;
    if (X && X.n) {
      var sbox = plotBox(''), sh = h('div');
      stack.appendChild(h('div', { class: 'grid g2' },
        card('Sex check', 'chrX vs chrY depth ratio; dashed lines = decision thresholds', sbox), card('Sex check table', null, sh)));
      tableView(sh, X, { filename: 'canope_sex_check', pageSize: 10, labels: { x_ratio: 'chrX ratio', y_ratio: 'chrY ratio' },
        cellClass: { mismatch: function (v) { return v ? 'warn' : ''; }, mode_mismatch: function (v) { return v ? 'warn' : ''; } }, onRow: function (i) { setSample(X.data.sample[i]); } });
      plot(sbox, function (t) {
        var c = SEM[t], data = [], cs = { F: CAT[t][4], M: CAT[t][0], unclear: c.NA };
        ['F', 'M', 'unclear'].forEach(function (g) {
          var ids = range(X.n).filter(function (i) { return X.data.inferred_sex[i] === g; });
          if (!ids.length) return;
          data.push({ type: 'scatter', mode: 'markers', name: g === 'unclear' ? 'unclear' : g === 'F' ? 'looks female' : 'looks male',
            x: ids.map(function (i) { return X.data.x_ratio[i]; }), y: ids.map(function (i) { return isNum(X.data.y_ratio[i]) ? X.data.y_ratio[i] : 0; }),
            customdata: ids.map(function (i) { return X.data.sample[i]; }), text: ids.map(function (i) { return X.data.sample[i] + (X.data.mismatch[i] ? ' (annotation mismatch)' : ''); }),
            marker: { color: cs[g], size: 10, symbol: g === 'F' ? 'circle' : g === 'M' ? 'square' : 'diamond',
              line: { width: ids.map(function (i) { return X.data.mismatch[i] ? 3 : 1; }), color: ids.map(function (i) { return X.data.mismatch[i] ? c.fail : CHROME[t].surface; }) } },
            hovertemplate: '%{text}<br>chrX %{x:.2f}, chrY %{y:.2f}<extra></extra>' });
        });
        return { data: data, layout: baseLayout(t, { xlab: 'chrX depth ratio', ylab: 'chrY depth ratio', legend: true,
          shapes: [{ type: 'line', xref: 'x', yref: 'paper', x0: 0.75, x1: 0.75, y0: 0, y1: 1, line: { color: c.fail, dash: 'dash', width: 1 } },
                   { type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: 0.1, y1: 0.1, line: { color: c.fail, dash: 'dash', width: 1 } }] }),
          onClick: function (ev) { var p = ev.points && ev.points[0]; if (p && p.customdata) setSample(p.customdata); } };
      });
    }
    root.appendChild(stack);
  });

  // ================================================================== PCA & CORRELATION
  tab('pca', 'PCA & similarity', function () { return !!(DATA.pca || DATA.corr); }, function (root) {
    var stack = h('div', { class: 'stack' }), P = DATA.pca;
    if (P) {
      var pcs = Object.keys(P.pcs), vars = asArray(P.var);
      var xs = selectOf(pcs.map(function (p, i) { return { value: p, text: p + ' (' + fmt(vars[i]) + '%)' }; }), pcs[0], null, 'X axis');
      var ys = selectOf(pcs.map(function (p, i) { return { value: p, text: p + ' (' + fmt(vars[i]) + '%)' }; }), pcs[Math.min(1, pcs.length - 1)], null, 'Y axis');
      var opts = [{ value: 'level', text: 'Warning level' }, { value: 'status', text: 'Status' }, { value: 'n_cnvs', text: 'Number of CNV calls' }, { value: 'sex_inferred', text: 'Sex (inferred)' }]
        .concat(Object.keys(DATA.sample_attrs || {}).map(function (k) { return { value: 'attr:' + k, text: pretty(k) }; }));
      var cb = selectOf(opts, 'level', null, 'Colour by');
      var box = plotBox('tall');
      stack.appendChild(h('div', { class: 'toolbar' }, field('X', xs), field('Y', ys), field('Colour by', cb), sampleSelect('Highlight')));
      stack.appendChild(card('PCA of exon counts', 'log2(count + 1), ' + fmtInt(P.n_targets) + ' targets; samples far from the cloud are poor references / likely batch outliers', box));
      var psamples = asArray(P.samples);
      var valueOf = function (s) {
        var v = cb.value, i = SIDX[s];
        if (i === undefined) return 'unknown';
        if (v === 'level') return LEVEL_LABEL[S.data.level[i]] === 'OK' ? 'No flags' : LEVEL_LABEL[S.data.level[i]];
        if (v.indexOf('attr:') === 0) return String(asArray(DATA.sample_attrs[v.slice(5)])[i]);
        var x = S.data[v][i]; return x === null || x === undefined ? 'NA' : String(x);
      };
      var draw = function () {
        plot(box, function (t) {
          var X = asArray(P.pcs[xs.value]), Y = asArray(P.pcs[ys.value]), cont = cb.value === 'n_cnvs', data = [];
          var xi = pcs.indexOf(xs.value), yi = pcs.indexOf(ys.value);
          if (cont) {
            data.push({ type: 'scatter', mode: 'markers', x: X, y: Y, text: psamples, customdata: psamples, hovertemplate: '<b>%{text}</b><br>%{marker.color} calls<extra></extra>',
              marker: { size: 11, color: psamples.map(function (s) { return sval(s, 'n_cnvs'); }), colorscale: SEQ.map(function (c2, i) { return [i / (SEQ.length - 1), c2]; }), showscale: true,
                colorbar: { title: { text: 'CNV calls' }, tickfont: { color: CHROME[t].muted } }, line: { width: 1, color: CHROME[t].surface } } });
          } else {
            var groups = uniq(psamples.map(valueOf)), cols = catColors(groups.length, t);
            if (cb.value === 'level') cols = groups.map(function (g) { return g === 'No flags' ? SEM[t].ok : g === 'Warning' ? SEM[t].warn : SEM[t].fail; });
            groups.forEach(function (g, gi) {
              var ids = range(psamples.length).filter(function (i) { return valueOf(psamples[i]) === g; });
              data.push({ type: 'scatter', mode: 'markers', name: g, x: ids.map(function (i) { return X[i]; }), y: ids.map(function (i) { return Y[i]; }),
                text: ids.map(function (i) { return psamples[i]; }), customdata: ids.map(function (i) { return psamples[i]; }),
                marker: { size: 11, color: cols[gi], symbol: g === 'Fail' ? 'x' : g === 'Warning' ? 'diamond' : 'circle', line: { width: 1, color: CHROME[t].surface } },
                hovertemplate: '<b>%{text}</b><br>' + g + '<extra></extra>' });
            });
          }
          var k = psamples.indexOf(state.sample);
          if (k >= 0) data.push({ type: 'scatter', mode: 'markers', name: 'selected', x: [X[k]], y: [Y[k]], hoverinfo: 'skip', showlegend: false,
            marker: { size: 24, color: 'rgba(0,0,0,0)', line: { color: CHROME[t].text, width: 2 } } });
          return { data: data, layout: baseLayout(t, { xlab: xs.value + ' (' + fmt(vars[xi]) + '%)', ylab: ys.value + ' (' + fmt(vars[yi]) + '%)', legend: !cont }),
            onClick: function (ev) { var p = ev.points && ev.points[0]; if (p && p.customdata) setSample(p.customdata); } };
        });
      };
      [xs, ys, cb].forEach(function (el) { el.addEventListener('change', draw); });
      onState(function () { if (box._render) box._render(); });
      draw();
    } else {
      stack.appendChild(card('PCA of exon counts', null, emptyNote('PCA needs at least 3 samples.')));
    }
    var C = DATA.corr;
    if (C) {
      var cbox = plotBox('tall');
      stack.appendChild(card('Sample-to-sample Spearman correlation', 'exon counts, hierarchically clustered; low blocks = batches that cannot serve as references for each other', cbox));
      plot(cbox, function (t) {
        var z = asArray(C.z).map(asArray), lo = 1;
        z.forEach(function (r, i) { r.forEach(function (v, j) { if (i !== j && isNum(v) && v < lo) lo = v; }); });
        return { data: [{ type: 'heatmap', x: asArray(C.samples), y: asArray(C.samples), z: z, zmin: Math.max(0, Math.floor(lo * 50) / 50), zmax: 1,
          colorscale: SEQ.map(function (c2, i) { return [i / (SEQ.length - 1), c2]; }), hovertemplate: '%{x} vs %{y}<br>rho = %{z:.3f}<extra></extra>',
          colorbar: { title: { text: 'rho' }, tickfont: { color: CHROME[t].muted } } }],
          layout: baseLayout(t, { xaxis: { type: 'category', tickangle: -45, showticklabels: asArray(C.samples).length <= 60 }, yaxis: { type: 'category', autorange: 'reversed', showticklabels: asArray(C.samples).length <= 60 },
            margin: { l: 80, r: 16, t: 10, b: 90 } }),
          onClick: function (ev) { var p = ev.points && ev.points[0]; if (p) setSample(p.x); } };
      });
    }
    root.appendChild(stack);
  });

  // ================================================================== CNV CALLS
  function chrNum(c) {
    var s = String(c).replace(/^chr/i, '').toUpperCase(), n = parseInt(s, 10);
    if (s === 'X') return 23; if (s === 'Y') return 24; if (s === 'M' || s === 'MT') return 25;
    return isFinite(n) ? n : 99;
  }
  function callFlags(i) {
    var out = [];
    if (tcol(CALLS, 'homology_gene')[i]) out.push('homology gene');
    if (tcol(CALLS, 'bg_flag')[i]) out.push('background');
    if (!tcol(CALLS, 'has_plot')[i]) out.push('no plot');
    var q = tcol(CALLS, 'q_some')[i];
    if (q === null || q === undefined) out.push('no Q_SOME');
    return out.join(', ');
  }
  function confReasons(i) {
    var R = DATA.conf_rules, c = CALLS.data, out = [];
    if (c.confidence[i] === 'HIGH') return out;
    if (c.homology_gene[i]) { out.push('the gene is on the pseudogene / homology list (always LOW): ' + asArray(R.genes).join(', ')); return out; }
    var q = c.q_some[i], refs = c.num_refs[i], nt = c.num_targ[i], ml = c.mlcn[i];
    if (!isNum(q)) out.push('Q_SOME is missing (genotyping was skipped for this sample)');
    else if (q < R.med_q) out.push('Q_SOME ' + q + ' is below the MEDIUM threshold (' + R.med_q + ')');
    else if (q < R.high_q) out.push('Q_SOME ' + q + ' is below the HIGH threshold (' + R.high_q + ')');
    if (isNum(refs) && refs < R.med_refs) out.push(refs + ' reference samples (< ' + R.med_refs + ' needed for MEDIUM)');
    else if (isNum(refs) && refs < R.high_refs) out.push(refs + ' reference samples (< ' + R.high_refs + ' needed for HIGH)');
    if (isNum(nt) && nt < R.med_targ) out.push(nt + ' target(s) (< ' + R.med_targ + ' needed)');
    if (c.confidence[i] !== 'HIGH' && isNum(q) && q >= R.high_q && isNum(ml) && !((c.type[i] === 'DEL' && ml <= 1) || (c.type[i] === 'DUP' && ml >= 3)))
      out.push('the most likely copy number (' + ml + ') is not consistent with a ' + c.type[i]);
    return out;
  }

  function callsTable() {
    var C = CALLS, cols = ['id', 'sample', 'interval', 'gene', 'type', 'kb', 'num_targ', 'mlcn', 'cn_label', 'q_some', 'confidence', 'num_refs', 'flags', 'view'];
    var T = { cols: cols, n: C.n, data: {} };
    cols.forEach(function (c) { T.data[c] = c === 'flags' ? range(C.n).map(callFlags) : c === 'view' ? range(C.n).map(function (i) { return C.data.has_plot[i]; }) : tcol(C, c); });
    return T;
  }

  function genomeMap(box, rowsFn, yMode, chrSet) {
    var ext = CHR.extent, chrs = asArray(ext.chr), lens = asArray(ext.len), off = {}, tot = 0, ticks = [], tt = [], shapes = [];
    chrs.forEach(function (c, i) {
      off[c] = tot;
      ticks.push(tot + lens[i] / 2); tt.push(String(c).replace(/^chr/, ''));
      if (i % 2 === 0) shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: tot, x1: tot + lens[i], y0: 0, y1: 1, line: { width: 0 }, fillcolor: 'rgba(128,128,128,0.08)', layer: 'below' });
      tot += lens[i] + 1;
    });
    plot(box, function (t) {
      var c = SEM[t], rows = rowsFn(), data = [];
      var ysamples = SAMPLES;
      ['DEL', 'DUP'].forEach(function (type) {
        var ids = rows.filter(function (i) { return CALLS.data.type[i] === type && off[CALLS.data.chr[i]] !== undefined; });
        if (!ids.length) return;
        var yv = ids.map(function (i) {
          return yMode === 'mlcn' ? CALLS.data.mlcn[i] : yMode === 'q_some' ? CALLS.data.q_some[i] : CALLS.data.sample[i];
        });
        data.push({ type: 'scatter', mode: 'markers', name: type, x: ids.map(function (i) { return off[CALLS.data.chr[i]] + (CALLS.data.start[i] + CALLS.data.end[i]) / 2; }), y: yv,
          customdata: ids.map(function (i) { return CALLS.data.id[i]; }),
          text: ids.map(function (i) { return '<b>' + CALLS.data.sample[i] + '</b> ' + type + ' ' + CALLS.data.interval[i] + '<br>' + (CALLS.data.gene[i] || '') + '<br>' + (CALLS.data.cn_label[i] || '') +
            ' · Q_SOME ' + (isNum(CALLS.data.q_some[i]) ? CALLS.data.q_some[i] : 'n/a') + ' · ' + CALLS.data.confidence[i] + ' · ' + CALLS.data.num_targ[i] + ' target(s)'; }),
          marker: { symbol: type === 'DEL' ? 'triangle-down' : 'triangle-up', color: c[type], size: ids.map(function (i) { return 8 + Math.min(12, 2 * Math.sqrt(CALLS.data.num_targ[i])); }),
            opacity: ids.map(function (i) { return { HIGH: 1, MEDIUM: 0.75, LOW: 0.45 }[CALLS.data.confidence[i]] || 0.45; }), line: { width: 1, color: CHROME[t].surface } },
          hovertemplate: '%{text}<extra></extra>' });
      });
      var yaxis = yMode === 'sample' ? { type: 'category', categoryorder: 'array', categoryarray: ysamples, autorange: 'reversed' } : { rangemode: 'tozero' };
      var L = baseLayout(t, { legend: true, ylab: yMode === 'sample' ? '' : yMode === 'mlcn' ? 'most likely copy number' : 'Q_SOME (phred)', yaxis: yaxis, shapes: shapes,
        xaxis: { tickmode: 'array', tickvals: ticks, ticktext: tt, showgrid: false, title: { text: 'chromosome (analysed chromosomes only, scaled to their target extent)' } },
        margin: { l: yMode === 'sample' ? 90 : 56, r: 16, t: 10, b: 60 } });
      return { data: data, layout: L, onClick: function (ev) { var p = ev.points && ev.points[0]; if (p && p.customdata !== undefined) openCall(p.customdata); } };
    });
  }

  tab('calls', 'CNV calls', function () { return true; }, function (root) {
    var stack = h('div', { class: 'stack' });
    if (!CALLS.n) {
      stack.appendChild(alertBox(SUM.samples_called ? 'info' : 'warn', SUM.samples_called ? 'No CNV was called.' : 'No sample was analysed.',
        SUM.samples_called ? 'A run without calls is a valid result only for the chromosomes that were analysed (see the Chromosomes tab).' : 'Check the Warnings tab.'));
      root.appendChild(stack); return;
    }
    if (SUM.both_chromosome_classes === false)
      stack.appendChild(alertBox('fail', 'Partial genome:', 'autosomes and gonosomes were not both analysed in this run (' + META.mode_label + '). Calls cover only that part.'));
    var f = { type: '', chr: '', minQ: 0, minT: 1, conf: { HIGH: true, MEDIUM: true, LOW: true, NA: true }, y: 'sample' };
    var typeSel = selectOf([{ value: '', text: 'DEL + DUP' }, { value: 'DEL', text: 'DEL' }, { value: 'DUP', text: 'DUP' }], '', null, 'Type');
    var chrs = uniq(tcol(CALLS, 'chr')).sort(function (a, b) { return chrNum(a) - chrNum(b); });
    var chrSel = selectOf([{ value: '', text: 'All chromosomes' }].concat(chrs.map(function (c) { return { value: c, text: c }; })), '', null, 'Chromosome');
    var minQ = h('input', { type: 'number', min: 0, max: 99, value: 0, 'aria-label': 'Minimum Q_SOME', style: { width: '80px', minWidth: '0' } });
    var minT = h('input', { type: 'number', min: 1, value: 1, 'aria-label': 'Minimum targets', style: { width: '80px', minWidth: '0' } });
    var ySel = selectOf([{ value: 'sample', text: 'Sample' }, { value: 'mlcn', text: 'Copy number (MLCN)' }, { value: 'q_some', text: 'Q_SOME' }], 'sample', null, 'Map y axis');
    var confBox = h('div', { class: 'legend', role: 'group', 'aria-label': 'Confidence filter' });
    [['HIGH', 'HIGH'], ['MEDIUM', 'MEDIUM'], ['LOW', 'LOW'], ['NA', 'unscored']].forEach(function (p) {
      confBox.appendChild(h('button', { type: 'button', 'aria-pressed': 'true', onclick: function () {
        f.conf[p[0]] = !f.conf[p[0]]; this.setAttribute('aria-pressed', f.conf[p[0]] ? 'true' : 'false'); apply(); } },
        h('span', { class: 'swatch conf ' + p[0] }), p[1]));
    });
    var count = h('span', { class: 'muted small' });
    var tableHost = h('div'), mapBox = plotBox('');
    mapBox.style.height = Math.min(820, 320 + (SAMPLES.length > 8 ? 12 * Math.min(SAMPLES.length, 40) : 0)) + 'px';
    var T = callsTable();
    var rowsNow = function () {
      var out = [], d = CALLS.data;
      for (var i = 0; i < CALLS.n; i++) if (pass(i)) out.push(i);
      return out;
    };
    function pass(i) {
      var d = CALLS.data;
      if (state.sample && d.sample[i] !== state.sample) return false;
      if (f.type && d.type[i] !== f.type) return false;
      if (f.chr && d.chr[i] !== f.chr) return false;
      var ck = d.confidence[i] === 'HIGH' || d.confidence[i] === 'MEDIUM' || d.confidence[i] === 'LOW' ? d.confidence[i] : 'NA';
      if (!f.conf[ck]) return false;
      if (f.minQ > 0 && !(isNum(d.q_some[i]) && d.q_some[i] >= f.minQ)) return false;
      if (f.minT > 1 && !(d.num_targ[i] >= f.minT)) return false;
      return true;
    }
    stack.appendChild(h('div', { class: 'toolbar' }, sampleSelect('Sample'), field('Type', typeSel), field('Chromosome', chrSel), field('Min Q_SOME', minQ), field('Min targets', minT),
      h('div', { class: 'field' }, 'Confidence', confBox), count));
    var mapNote = '';
    var notAn = asArray(tcol(CHR.inventory, 'chromosome')).filter(function (c, i) { return !CHR.inventory.data.analysed[i]; });
    if (notAn.length) mapNote = 'Not analysed in this run (absent from the map): ' + notAn.join(', ') + '.';
    stack.appendChild(card('Genome-wide map of calls', 'triangles: ▼ deletion, ▲ duplication; size = number of targets, opacity = confidence; click a call to open it in the viewer',
      h('div', null, h('div', { class: 'toolbar' }, field('Y axis', ySel)), mapBox, mapNote ? h('div', { class: 'note', text: mapNote }) : null)));
    stack.appendChild(card('CNV calls', 'click a row to open the call viewer', tableHost));
    root.appendChild(stack);

    var tv = tableView(tableHost, T, { filename: 'canope_cnv_calls', pageSize: 20, sort: ['q_some', 'desc'],
      labels: { id: '#', num_targ: 'targets', mlcn: 'MLCN', cn_label: 'copy number', q_some: 'Q_SOME', num_refs: 'refs', kb: 'kb', view: 'plots' },
      render: { type: typeChip, confidence: confChip, sample: function (v) { return h('strong', { text: v }); },
        view: function (v, i) { return v ? h('button', { class: 'btn small-btn', type: 'button', text: 'View', onclick: function () { openCall(CALLS.data.id[i]); } }) : h('span', { class: 'muted', text: '–' }); },
        flags: function (v) { return v ? h('span', { class: 'note-flag', text: v }) : ''; } },
      cellClass: { confidence: function (v) { return v === 'LOW' ? 'warn' : ''; }, flags: function (v) { return v ? 'warn' : ''; } },
      onRow: function (i) { if (tcol(CALLS, 'has_plot')[i]) openCall(CALLS.data.id[i]); else { state.call = CALLS.data.id[i]; } } });
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
    onState(apply);
    apply();
  });

  // ================================================================== CALL VIEWER
  function hasPlot(i) { return !!tcol(CALLS, 'has_plot')[i]; }

  function viewerFigure(P, row, t, showRefs) {
    var c = SEM[t], D = CALLS.data, type = D.type[row], tc = c[type] || c.info;
    var px = asArray(P.px), aff = asArray(P.affected), n = px.length, i;
    var refLog = asArray(P.ref_log).map(asArray), refZ = asArray(P.ref_z).map(asArray);
    var txt = range(n).map(function (k) {
      return '<b>' + (asArray(P.gene)[k] || '') + '</b> exon ' + asArray(P.exon)[k] + (aff[k] ? ' (called)' : '') + '<br>' + asArray(P.pos)[k];
    });
    // Concatenate one x/y pair per reference sample into a single trace (null
    // separators are skipped by plotly in markers mode).
    var series = function (mats) {
      var X = [], Y = [];
      mats.forEach(function (m) { for (var k = 0; k < n; k++) { X.push(px[k]); Y.push(m[k]); } X.push(null); Y.push(null); });
      return { x: X, y: Y };
    };
    var testLog = asArray(P.test_log), ratio = asArray(P.ratio), lo = asArray(P.lo), hi = asArray(P.hi), z = asArray(P.z);
    // Points are never joined by lines, in every panel -- same convention as
    // the static PDF plots (generate_plots): exons are ordered genomic
    // positions, and a connecting line would suggest a continuous signal
    // between them (and run straight across the inserted gene gaps).
    var markerCol = range(n).map(function (k) { return aff[k] ? tc : c.test; });
    var markerSym = range(n).map(function (k) { return aff[k] ? (type === 'DEL' ? 'triangle-down' : 'triangle-up') : 'circle'; });
    var markerSize = range(n).map(function (k) { return aff[k] ? 11 : 7; });
    var data = [];
    var nref = refLog.length;
    if (showRefs && nref) {
      var rl = series(refLog), rz = series(refZ);
      data.push({ type: 'scatter', mode: 'markers', name: 'Reference samples (' + nref + ')', x: rl.x, y: rl.y, xaxis: 'x', yaxis: 'y', marker: { color: c.ref, size: 4, opacity: 0.75 }, hoverinfo: 'skip', legendgroup: 'refs' });
      data.push({ type: 'scatter', mode: 'markers', name: 'refs z', x: rz.x, y: rz.y, xaxis: 'x', yaxis: 'y3', marker: { color: c.ref, size: 4, opacity: 0.75 }, hoverinfo: 'skip', legendgroup: 'refs', showlegend: false });
    }
    data.push({ type: 'scatter', mode: 'markers', name: 'Test sample ' + P.sample, x: px, y: testLog, xaxis: 'x', yaxis: 'y', text: txt,
      marker: { color: markerCol, symbol: markerSym, size: markerSize, line: { width: 1, color: CHROME[t].surface } },
      hovertemplate: '%{text}<br>log2(count + 0.5) = %{y:.2f}<extra></extra>' });
    var errUp = range(n).map(function (k) { return isNum(hi[k]) && isNum(ratio[k]) ? Math.max(0, hi[k] - ratio[k]) : 0; });
    var errDn = range(n).map(function (k) { return isNum(lo[k]) && isNum(ratio[k]) ? Math.max(0, ratio[k] - lo[k]) : 0; });
    var outside = range(n).map(function (k) { return isNum(ratio[k]) && isNum(lo[k]) && isNum(hi[k]) && (ratio[k] < lo[k] || ratio[k] > hi[k]); });
    data.push({ type: 'scatter', mode: 'markers', name: 'log2 obs/exp', x: px, y: ratio, xaxis: 'x', yaxis: 'y2', showlegend: false, text: txt,
      customdata: outside.map(function (o, k) { return (o ? 'outside' : 'inside') + ' the 95% interval [' + fmt(lo[k]) + ', ' + fmt(hi[k]) + ']'; }),
      error_y: { type: 'data', symmetric: false, array: errUp, arrayminus: errDn, color: c.flank, thickness: 1, width: 3 },
      marker: { color: markerCol, symbol: markerSym, size: markerSize, line: { width: 1, color: CHROME[t].surface } },
      hovertemplate: '%{text}<br>log2(obs/exp) = %{y:.2f}<br>%{customdata}<extra></extra>' });
    data.push({ type: 'scatter', mode: 'markers', name: 'z', x: px, y: z, xaxis: 'x', yaxis: 'y3', showlegend: false, text: txt,
      marker: { color: markerCol, symbol: markerSym, size: markerSize, line: { width: 1, color: CHROME[t].surface } },
      hovertemplate: '%{text}<br>z = %{y:.2f}<extra></extra>' });
    // legend entry for the called exons
    data.push({ type: 'scatter', mode: 'markers', name: 'Called exons (' + type + ')', x: [null], y: [null], xaxis: 'x', yaxis: 'y',
      marker: { color: tc, symbol: type === 'DEL' ? 'triangle-down' : 'triangle-up', size: 11 } });

    var affIdx = range(n).filter(function (k) { return aff[k]; });
    var shapes = [];
    if (affIdx.length) {
      var x0 = Math.min.apply(null, affIdx.map(function (k) { return px[k]; })) - 0.5, x1 = Math.max.apply(null, affIdx.map(function (k) { return px[k]; })) + 0.5;
      shapes.push({ type: 'rect', xref: 'x', yref: 'paper', x0: x0, x1: x1, y0: 0, y1: 1, fillcolor: alpha(tc, 0.1), line: { width: 0 }, layer: 'below' });
    }
    shapes.push({ type: 'line', xref: 'paper', yref: 'y2', x0: 0, x1: 1, y0: 0, y1: 0, line: { color: CHROME[t].axis, width: 1 } });
    [3, -3].forEach(function (v) { shapes.push({ type: 'line', xref: 'paper', yref: 'y3', x0: 0, x1: 1, y0: v, y1: v, line: { color: c.fail, dash: 'dot', width: 1 } }); });
    shapes.push({ type: 'line', xref: 'paper', yref: 'y3', x0: 0, x1: 1, y0: 0, y1: 0, line: { color: CHROME[t].axis, width: 1 } });

    // x ticks: one label per gene block
    var groups = asArray(P.group), gene = asArray(P.gene), seen = {}, tv = [], tx = [];
    for (i = 0; i < n; i++) {
      var g = groups[i];
      if (!seen[g]) { seen[g] = { s: 0, k: 0, gene: gene[i] }; }
      seen[g].s += px[i]; seen[g].k++;
    }
    Object.keys(seen).forEach(function (g) { tv.push(seen[g].s / seen[g].k); tx.push(seen[g].gene || ''); });

    var L = baseLayout(t, { legend: true, legendOpts: { y: 1.1, x: 0, xanchor: 'left', yanchor: 'bottom' }, shapes: shapes, margin: { l: 70, r: 16, t: 60, b: 70 },
      xaxis: { anchor: 'y3', domain: [0, 1], tickmode: 'array', tickvals: tv, ticktext: tx, tickangle: -30, title: { text: 'exons (grouped by gene; gaps = gene boundaries)' } },
      yaxis: { domain: [0.72, 1], title: { text: 'log2(count + 0.5)' } } });
    var sub = function (domain, label) { return Object.assign({}, L.yaxis, { domain: domain, anchor: 'x', title: { text: label, font: { color: CHROME[t].secondary } } }); };
    L.yaxis = Object.assign({}, L.yaxis, { anchor: 'x', title: { text: 'log2(count + 0.5)', font: { color: CHROME[t].secondary } } });
    L.yaxis2 = sub([0.38, 0.64], 'log2 (obs / exp)');
    L.yaxis3 = sub([0, 0.30], 'z-score');
    var pt = function (y, text) { return { xref: 'paper', yref: 'paper', x: 0, y: y, xanchor: 'left', yanchor: 'bottom', showarrow: false, text: text, font: { size: 12, color: CHROME[t].text } }; };
    L.annotations = [pt(1.0, '<b>1 · Coverage</b> (test sample vs. its reference samples)'), pt(0.64, '<b>2 · Observed / expected</b> (bars = 95% negative-binomial interval)'),
      pt(0.30, '<b>3 · z-score</b> against the references (dotted = ±3)')];
    return { data: data, layout: L };
  }

  tab('viewer', 'Call viewer', function () { return CALLS.n > 0; }, function (root) {
    var stack = h('div', { class: 'stack' }), infoBox = h('div'), alertsBox = h('div', { class: 'stack tight' }), plotHost = h('div', { class: 'viewer-plot' }), box = plotBox('viewer');
    var callSel = h('select', { 'aria-label': 'Call', style: { maxWidth: 'min(640px, 100%)' } });
    var refs = h('input', { type: 'checkbox', checked: true, 'aria-label': 'Show reference samples' });
    var prev = h('button', { class: 'btn', type: 'button', text: '‹ Previous' }), next = h('button', { class: 'btn', type: 'button', text: 'Next ›' });
    var current = [];
    stack.appendChild(h('div', { class: 'toolbar' }, sampleSelect('Sample'), field('Call', callSel), prev, next, h('label', { class: 'field inline' }, refs, 'show references')));
    stack.appendChild(h('div', { class: 'grid g-viewer' }, card('Call', null, h('div', null, infoBox, alertsBox)), card('Interactive plots', 'zoom/pan is linked across the three panels; the shaded band is the called region', box)));
    root.appendChild(stack);
    function fillList() {
      var rows = callOrder(callsOf(state.sample));
      current = rows;
      clear(callSel);
      rows.forEach(function (r) {
        var d = CALLS.data;
        callSel.appendChild(h('option', { value: d.id[r], text: '#' + d.id[r] + ' ' + (state.sample ? '' : d.sample[r] + ' · ') + d.type[r] + ' ' + d.interval[r] + ' ' + (d.gene[r] || '') + ' · ' + d.confidence[r] +
          (isNum(d.q_some[r]) ? ' · Q' + d.q_some[r] : '') + (hasPlot(r) ? '' : ' (no plot)') }));
      });
      if (!rows.length) return;
      if (state.call === null || current.map(function (r) { return CALLS.data.id[r]; }).indexOf(state.call) < 0) {
        var first = rows.filter(hasPlot)[0]; state.call = CALLS.data.id[first !== undefined ? first : rows[0]];
      }
      callSel.value = state.call;
    }
    function show() {
      clear(infoBox); clear(alertsBox);
      var row = CALL_ROW[state.call];
      if (row === undefined) { infoBox.appendChild(emptyNote('This sample has no CNV call.')); showNote(box, 'No call selected.'); return; }
      var d = CALLS.data;
      infoBox.appendChild(h('div', { class: 'toolbar' }, typeChip(d.type[row]), confChip(d.confidence[row]), h('strong', { text: d.sample[row] })));
      infoBox.appendChild(kv([['Interval', d.interval[row]], ['Gene(s)', d.gene[row]], ['Copy number', (d.cn_label[row] || '') + (isNum(d.mlcn[row]) ? ' (MLCN ' + d.mlcn[row] + ')' : '')],
        ['Q_SOME', isNum(d.q_some[row]) ? String(d.q_some[row]) : 'n/a'], ['Targets / size', d.num_targ[row] + ' targets, ' + d.kb[row] + ' kb'], ['Reference samples', fmt(d.num_refs[row])], ['Target ids', d.targets[row]]]));
      var reasons = confReasons(row);
      if (reasons.length) alertsBox.appendChild(alertBox('warn', 'Why ' + d.confidence[row] + '?', reasons.join('; ') + '.'));
      var P = PLOTS[String(state.call)];
      if (!P) {
        alertsBox.appendChild(alertBox('info', null, 'No interactive plot is embedded for this call (plot limit reached or no model data). The PDF plots written by the pipeline are unaffected.'));
        showNote(box, 'No plot data for this call.'); return;
      }
      if (P.bg && P.bg.flag) alertsBox.appendChild(alertBox('warn', 'Background calibration:',
        fmt(P.bg.pct_outside) + '% of the non-called exons in the window (' + P.bg.n_outside + ' of ' + P.bg.n_background + ') lie outside their 95% interval. The sample is noisy here; confirm the call by another method.'));
      else if (P.bg && isNum(P.bg.pct_outside)) alertsBox.appendChild(alertBox('ok', null, fmt(P.bg.pct_outside) + '% of the neighbouring exons are outside their 95% interval (' + P.bg.n_outside + ' of ' + P.bg.n_background + '): background looks calibrated.'));
      plot(box, function (t) { return viewerFigure(P, row, t, refs.checked); });
      setTimeout(function () { if (window.Plotly && box.data) window.Plotly.Plots.resize(box); }, 30);
    }
    function step(k) {
      var ids = current.map(function (r) { return CALLS.data.id[r]; }), at = ids.indexOf(state.call);
      if (!ids.length) return;
      state.call = ids[(at + k + ids.length) % ids.length]; callSel.value = state.call; show();
    }
    callSel.addEventListener('change', function () { state.call = +callSel.value; show(); });
    refs.addEventListener('change', show);
    prev.addEventListener('click', function () { step(-1); }); next.addEventListener('click', function () { step(1); });
    onState(function () { fillList(); show(); });
    fillList(); show();
  });

  // ================================================================== METHOD
  tab('method', 'Method & thresholds', function () { return true; }, function (root) {
    var p = DATA.params || {}, R = DATA.conf_rules, stack = h('div', { class: 'stack' });
    stack.appendChild(card('How CANOPE works', null, h('div', { class: 'prose' },
      h('p', { text: 'CANOPE is a CANOES-derived exome CNV caller. For every sample it selects the most correlated reference samples, combines them with non-negative weights into an expected coverage per exon, and models the read count with a negative binomial whose variance is estimated from the references. A hidden Markov model with distance-decayed transitions (Viterbi) calls deletions and duplications; Q_SOME is the phred-scaled probability that at least one exon of the call is abnormal, and MLCN the most likely copy number.' }),
      h('p', { text: 'One run covers one chromosome mode: autosomes (A), chrX with female references (XX) or chrX/chrY with male references (XY). Sex chromosomes are therefore never analysed together with the autosomes, which is why this report states which classes were covered.' }))));
    stack.appendChild(h('div', { class: 'grid g2' },
      card('Run parameters', null, kv(Object.keys(p).map(function (k) { return [k, Array.isArray(p[k]) ? p[k].join(', ') : String(p[k])]; }))),
      card('QC thresholds used for the warnings', null, kv([
        ['Median depth ≥', fmt(THR.qc_min_cov)], ['Total reads ≥', fmtInt(THR.qc_min_total_reads)], ['Best correlation ≥', fmt(THR.qc_min_corr)],
        ['Exon CV ≤', fmt(THR.qc_max_exon_cv)], ['Noise robust z ≤', fmt(THR.qc_zscore)], ['Reference samples ≥', '5 (warning below)']]))));
    stack.appendChild(card('Confidence tiers', null, kv([
      ['HIGH', 'Q_SOME ≥ ' + R.high_q + ', ≥ ' + R.high_refs + ' references, ≥ ' + R.high_targ + ' target(s), copy number consistent with the call, gene not on the homology list'],
      ['MEDIUM', 'Q_SOME ≥ ' + R.med_q + ', ≥ ' + R.med_refs + ' references, ≥ ' + R.med_targ + ' target(s), gene not on the homology list'],
      ['LOW', 'everything else, including genes with known pseudogene / homology issues: ' + asArray(R.genes).join(', ')]])));
    var glossary = [
      ['Chromosomes not analysed', 'The BED/count table contains autosomal or gonosomal targets that this run did not call (modechrom / removeY / QC filters). A missing call there is not a negative result.'],
      ['Sample failed / excluded', 'Failed samples raised an error while calling; excluded samples were dropped before calling (e.g. too few reads). Neither has CNV calls, and excluded samples are not used as references.'],
      ['Low depth / low reads', 'Median target depth or total reads below the QC threshold: reduced sensitivity, especially for single-exon events.'],
      ['Low correlation / few references', 'No other sample resembles this one (batch, capture kit, tissue), so the expected coverage is poor and false calls are likely.'],
      ['Noisy sample', 'The exon-count dispersion is an outlier within the cohort (robust z-score).'],
      ['Genotyping skipped', 'Too many calls were produced; Q_SOME and MLCN were not computed and every call is LOW.'],
      ['Low-coverage run removed', 'Runs of consecutive targets with ~no reads were removed before calling: a homozygous deletion or a capture dropout there cannot be called.'],
      ['Sex / mode mismatch', 'The annotated sex disagrees with chrX/chrY coverage, or the sample does not match the sex expected by modechrom: chrX calls are unreliable.'],
      ['Background calibration', 'A high share of non-called exons around a call fall outside their 95% interval, meaning the expected-coverage model is miscalibrated for this sample/locus.']];
    stack.appendChild(card('Warning glossary', null, kv(glossary)));
    root.appendChild(stack);
  });

  // ================================================================== RUN & LOG
  tab('run', 'Run & log', function () { return true; }, function (root) {
    var L = DATA.log || {}, blocks = [];
    if (asArray(L.warnings).length) blocks.push(card('Pipeline warnings (' + asArray(L.warnings).length + ')', null, h('pre', { text: asArray(L.warnings).join('\n') })));
    if (asArray(L.tail).length) blocks.push(card('Pipeline log (tail)', null, h('pre', { text: asArray(L.tail).join('\n') })));
    if (asArray(L.session).length) blocks.push(card('Session info', null, h('pre', { text: asArray(L.session).join('\n') })));
    if (!blocks.length) blocks.push(emptyNote('No pipeline log was supplied to the report.'));
    root.appendChild(h('div', { class: 'stack' }, blocks));
  });

  // ------------------------------------------------------------------ shell
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
    if (push) { try { history.replaceState(null, '', '#tab=' + id); } catch (e) { /* file:// may refuse */ } }
  }
  function init() {
    applyThemeAttr();
    document.title = META.title;
    var tabsEl = $('#tabs'), main = $('#panels');
    var avail = TABS.filter(function (t) { return t.available(); });
    avail.forEach(function (t, i) {
      var label = t.label;
      var b = h('button', { class: 'tab', role: 'tab', id: 'tab-' + t.id, 'aria-controls': 'panel-' + t.id, 'aria-selected': 'false', tabindex: '-1', type: 'button',
        onclick: function () { activate(t.id, true); },
        onkeydown: function (e) {
          var k = e.key === 'ArrowRight' ? 1 : e.key === 'ArrowLeft' ? -1 : 0;
          if (!k) return; e.preventDefault();
          var nx = avail[(i + k + avail.length) % avail.length]; activate(nx.id, true); tabButtons[nx.id].focus();
        } }, label, t.badge ? t.badge() : null);
      tabButtons[t.id] = b; tabsEl.appendChild(b);
      panels[t.id] = h('section', { class: 'panel', role: 'tabpanel', id: 'panel-' + t.id, 'aria-labelledby': 'tab-' + t.id, hidden: true });
      main.appendChild(panels[t.id]);
    });
    var btn = $('#theme-btn');
    var labelFn = function () { btn.textContent = 'Theme: ' + themePref; };
    labelFn();
    btn.addEventListener('click', function () {
      themePref = themePref === 'auto' ? 'light' : themePref === 'light' ? 'dark' : 'auto';
      store('canope-theme', themePref); applyThemeAttr(); labelFn(); rerenderAll();
    });
    if (window.matchMedia) {
      var mq = window.matchMedia('(prefers-color-scheme: dark)');
      var onChange = function () { if (themePref === 'auto') rerenderAll(); };
      if (mq.addEventListener) mq.addEventListener('change', onChange);
    }
    var m = /tab=([\w-]+)/.exec(location.hash);
    var cm = /call=(\d+)/.exec(location.hash);
    if (cm && CALL_ROW[+cm[1]] !== undefined) state.call = +cm[1];
    activate(m ? m[1] : 'overview', false);
    window.addEventListener('hashchange', function () { var mm = /tab=([\w-]+)/.exec(location.hash); if (mm && mm[1] !== active) activate(mm[1], false); });
    document.documentElement.setAttribute('data-ready', '1');
  }
  // Robustly make sure Plotly is available before init() runs.
  // 1) prefer whatever was already loaded (inline or CDN <script>)
  // 2) otherwise lazily load the CDN build
  // 3) otherwise show a clear message instead of a cryptic "Plotly.react" error
  function ensurePlotly(cb, timeoutMs) {
    timeoutMs = timeoutMs || 8000;
    var done = false;
    function has() {
      return (typeof window.Plotly !== 'undefined') &&
             window.Plotly && (typeof window.Plotly.react === 'function');
    }
    function finish(ok) {
      if (done) return;
      done = true;
      clearTimeout(timer);
      if (ok && has()) cb();
      else {
        var host = document.getElementById('panels');
        if (host) host.innerHTML = '<div class="empty">' +
          'Plotly.js failed to load. Please check your internet connection ' +
          'or reinstall the R package with <code>install.packages("plotly")</code>.</div>';
        if (window.console) console.error('[CANOPE] Plotly unavailable.');
      }
    }
    var timer = setTimeout(function () { finish(false); }, timeoutMs);
    if (has()) { finish(true); return; }
    if (window.console) console.warn('[CANOPE] Plotly not found inline; loading from jsDelivr CDN...');
    var s = document.createElement('script');
    s.src = 'https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js';
    s.async = false;
    s.onload = function () { finish(true); };
    s.onerror = function () { finish(false); };
    document.head.appendChild(s);
  }

  if (document.readyState === 'loading')
    document.addEventListener('DOMContentLoaded', function () { ensurePlotly(init); });
  else
    ensurePlotly(init);
})();
)---"


# ---------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------

.rp_meta_cols <- c("target", "chromosome", "start", "end", "gc", "GENE", "mean", "exon_number")

.rp_arr <- function(x) I(as.vector(x))

.rp_num <- function(x, digits = 3) .rp_arr(round(as.numeric(x), digits))

.rp_table <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  data <- lapply(df, function(col) {
    if (is.factor(col)) col <- as.character(col)
    if (is.numeric(col)) .rp_arr(col) else if (is.logical(col)) .rp_arr(col) else .rp_arr(as.character(col))
  })
  list(cols = .rp_arr(names(df)), data = data, n = nrow(df))
}

.rp_esc <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

.rp_chr_order <- function(chr) {
  s <- toupper(sub("^chr", "", as.character(chr), ignore.case = TRUE))
  n <- suppressWarnings(as.numeric(s))
  n[s == "X"] <- 23; n[s == "Y"] <- 24; n[s %in% c("M", "MT")] <- 25
  n[is.na(n)] <- 99
  n
}

.rp_conf_rank <- function(x) {
  r <- c(HIGH = 3L, MEDIUM = 2L, LOW = 1L)
  out <- unname(r[toupper(as.character(x))]); out[is.na(out)] <- 0L; out
}

# ---------------------------------------------------------------------------------
# per-call plot window
# ---------------------------------------------------------------------------------

#' Plot data for one CNV call (consumed by the report's call viewer)
#'
#' Uses the very objects the HMM used (\code{models[[sample]]}: expected mean,
#' variance, the sample's counts and the scaled reference matrix), so the three
#' panels agree with each other and with the static PDF plots.
#' @noRd
.rp_call_window <- function(call_row, bed, bed_chr, models, flank = 5L, gene_gap = 1, exon_index = NULL) {
  s <- as.character(call_row$SAMPLE)
  m <- models[[s]]
  if (is.null(m) || is.null(m$target) || is.null(m$ref_matrix)) return(NULL)
  tb <- parse_canope_targets(call_row$TARGETS)
  if (length(tb) < 2L) return(NULL)
  chr_call <- format_chr_label(call_row$CHR)

  in_model  <- bed$target %in% m$target
  call_rows <- which(bed$target >= tb[1] & bed$target <= tb[2] & bed_chr == chr_call & in_model)
  if (length(call_rows) == 0L) return(NULL)
  cand <- which(in_model & bed_chr == chr_call)
  lo_i <- match(min(call_rows), cand); hi_i <- match(max(call_rows), cand)
  win  <- cand[seq(max(1L, lo_i - flank), min(length(cand), hi_i + flank))]

  mi <- match(bed$target[win], m$target)
  mean_e <- m$mean[mi]; var_e <- m$var_estimate[mi]; test <- m$test_counts[mi]
  refmat <- m$ref_matrix[mi, , drop = FALSE]
  pc <- 0.5
  expected_safe <- pmax(mean_e, 1) + pc
  ratio <- log2((test + pc) / expected_safe)
  disp  <- var_e - mean_e
  size  <- ifelse(disp > 0, mean_e^2 / disp, Inf)
  lo <- log2((stats::qnbinom(0.025, mu = mean_e, size = size) + pc) / expected_safe)
  hi <- log2((stats::qnbinom(0.975, mu = mean_e, size = size) + pc) / expected_safe)

  z_mean <- rowMeans(refmat, na.rm = TRUE)
  z_sd   <- sqrt(pmax(var_e, 1)); z_sd[!is.finite(z_sd) | z_sd < 1e-6] <- 1e-6
  z      <- (test - z_mean) / z_sd
  ref_z  <- (refmat - z_mean) / z_sd

  gp <- compute_gene_gap_positions(bed, win, gap = gene_gap)
  affected <- win %in% call_rows
  bg <- check_background_calibration(ratio, lo, hi, affected)

  rows_by_ref <- function(mat, digits) lapply(seq_len(ncol(mat)), function(j) .rp_num(mat[, j], digits))
  list(
    sample = s,
    exon  = .rp_arr(if (!is.null(exon_index)) exon_index[win] else seq_along(win)),
    gene  = .rp_arr(as.character(bed$GENE[win])),
    pos   = .rp_arr(paste0(bed_chr[win], ":", bed$start[win], "-", bed$end[win])),
    px    = .rp_num(gp$px, 2),
    group = .rp_arr(gp$gene_group),
    affected = .rp_arr(as.integer(affected)),
    test_log = .rp_num(log2(test + pc), 3),
    ratio = .rp_num(ratio, 3), lo = .rp_num(lo, 3), hi = .rp_num(hi, 3),
    z = .rp_num(z, 2),
    ref_names = .rp_arr(colnames(refmat)),
    ref_log = rows_by_ref(log2(refmat + pc), 2),
    ref_z = rows_by_ref(ref_z, 2),
    bg = list(n_background = bg$n_background, n_outside = bg$n_outside,
              pct_outside = if (is.na(bg$pct_outside)) NULL else round(bg$pct_outside, 1),
              flag = isTRUE(bg$flag))
  )
}

# ---------------------------------------------------------------------------------
# main data builder
# ---------------------------------------------------------------------------------

.rp_build_data <- function(env, qc_metrics_file, settings, sample_table, log_file,
                           prefix, title, max_call_plots, max_qc_rows, flank) {
  cnv    <- env$cnv_calls
  counts <- env$counts
  bed    <- env$bed_file
  models <- env$models
  refs   <- env$refs
  run    <- if ("run_info" %in% ls(env)) env$run_info else list()
  engine <- if ("hmm_engine_used" %in% ls(env)) env$hmm_engine_used else (run$engine %||% "new")
  settings <- settings %||% list()
  set <- function(name, default) settings[[name]] %||% run$thresholds[[name]] %||% default

  if (is.null(cnv) || ncol(cnv) == 0L) cnv <- .empty_cnv_calls(scored = TRUE)
  if (!"Confidence" %in% colnames(cnv)) cnv$Confidence <- rep(NA_character_, nrow(cnv))
  if (!"CN_label"   %in% colnames(cnv)) cnv$CN_label   <- rep(NA_character_, nrow(cnv))

  thr <- list(qc_min_cov = as.numeric(set("qc_min_cov", 100)),
              qc_min_total_reads = as.numeric(set("qc_min_total_reads", 3e5)),
              qc_min_corr = as.numeric(set("qc_min_corr", 0.98)),
              qc_max_exon_cv = as.numeric(set("qc_max_exon_cv", 0.5)),
              qc_zscore = as.numeric(run$params$qc_zscore %||% 3))
  gene_gap <- as.numeric(settings$plot_gene_gap %||% 1)

  # ---- samples ------------------------------------------------------------------
  sample_cols <- setdiff(colnames(counts), .rp_meta_cols)
  requested <- unique(c(run$samples_requested, sample_cols, as.character(cnv$SAMPLE)))
  requested <- requested[!is.na(requested) & nzchar(requested)]

  depth <- reads <- setNames(rep(NA_real_, length(requested)), requested)
  if (length(sample_cols) > 0) {
    depth[sample_cols] <- vapply(sample_cols, function(s) stats::median(counts[[s]], na.rm = TRUE), numeric(1))
    reads[sample_cols] <- vapply(sample_cols, function(s) sum(counts[[s]], na.rm = TRUE), numeric(1))
  }

  # sample-vs-sample Spearman (also gives each sample's best match)
  corr <- NULL; max_corr <- setNames(rep(NA_real_, length(requested)), requested)
  if (length(sample_cols) >= 2L && length(sample_cols) <= 400L) {
    mat <- as.matrix(counts[, sample_cols, drop = FALSE]); storage.mode(mat) <- "numeric"
    cm <- tryCatch(suppressWarnings(stats::cor(mat, method = "spearman",
                                               use = if (anyNA(mat)) "pairwise.complete.obs" else "everything")),
                   error = function(e) NULL)
    if (!is.null(cm)) {
      cm2 <- cm; diag(cm2) <- NA_real_
      max_corr[sample_cols] <- apply(cm2, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
      ord <- if (length(sample_cols) >= 3L && !anyNA(cm))
        tryCatch(stats::hclust(stats::as.dist(1 - cm), method = "average")$order, error = function(e) seq_along(sample_cols))
      else seq_along(sample_cols)
      if (length(sample_cols) <= 150L)
        corr <- list(samples = .rp_arr(sample_cols[ord]),
                     z = lapply(ord, function(i) .rp_num(cm[i, ord], 3)))
    }
  }

  noise <- run$sample_noise
  if (is.null(noise) && length(sample_cols) >= 2L)
    noise <- tryCatch(detect_outlier_samples(counts[, sample_cols, drop = FALSE], z_threshold = thr$qc_zscore),
                      error = function(e) NULL)

  status_df <- run$sample_status
  get_status <- function(s, col, default = NA) {
    if (is.null(status_df) || nrow(status_df) == 0 || !s %in% status_df$sample) return(default)
    v <- status_df[[col]][match(s, status_df$sample)]
    if (is.null(v)) default else v
  }

  n_by <- function(s, f) sum(cnv$SAMPLE == s & f)
  conf <- toupper(as.character(cnv$Confidence))
  sex_check <- run$sex_check

  # ---- per-sample flags ------------------------------------------------------------
  flags <- setNames(vector("list", length(requested)), requested)
  add_flag <- function(s, level, text) flags[[s]][[length(flags[[s]]) + 1L]] <<- list(level = level, text = text)

  call_counts <- vapply(requested, function(s) sum(cnv$SAMPLE == s), numeric(1))
  analysed_calls <- call_counts[requested %in% (run$samples_analysed %||% sample_cols)]
  mad_calls <- if (length(analysed_calls) >= 5) stats::mad(analysed_calls) else NA_real_
  med_calls <- if (length(analysed_calls) >= 5) stats::median(analysed_calls) else NA_real_
  many_calls_cut <- if (!is.na(mad_calls)) max(10, med_calls + 3 * max(mad_calls, 1)) else NA_real_

  low_genes <- as.character(unlist(settings$low_confidence_genes %||%
                                     eval(formals(score_canope_confidence)$low_confidence_genes)))
  lcr <- run$low_coverage_runs

  for (s in requested) {
    st <- get_status(s, "status", if (s %in% sample_cols) "called" else "unknown")
    msg <- get_status(s, "message", "")
    if (identical(st, "failed"))
      add_flag(s, "fail", paste0("CNV calling FAILED for this sample (", msg, "). 'No calls' here does not mean 'no CNV'."))
    if (identical(st, "excluded"))
      add_flag(s, "fail", paste0("Excluded from the analysis before calling: ", msg, ". It has no CNV calls and was not used as a reference."))
    if (identical(st, "called") && nzchar(msg) && !is.na(msg))
      add_flag(s, "info", paste0("Warnings raised while calling: ", msg))
    if (!is.na(depth[[s]]) && depth[[s]] < thr$qc_min_cov)
      add_flag(s, "warn", sprintf("Low median target depth: %.1f (threshold %g).", depth[[s]], thr$qc_min_cov))
    if (!is.na(reads[[s]]) && reads[[s]] < thr$qc_min_total_reads)
      add_flag(s, "warn", sprintf("Low total reads on the analysed targets: %s (threshold %s).",
                                  format(round(reads[[s]]), big.mark = ",", scientific = FALSE), format(thr$qc_min_total_reads, big.mark = ",", scientific = FALSE)))
    if (!is.na(max_corr[[s]]) && max_corr[[s]] < thr$qc_min_corr)
      add_flag(s, "warn", sprintf("Best Spearman correlation with any other sample is only %.3f (threshold %g): few comparable references.",
                                  max_corr[[s]], thr$qc_min_corr))
    if (!is.null(noise) && s %in% noise$sample) {
      z <- noise$robust_z[match(s, noise$sample)]
      if (!is.na(z) && z > thr$qc_zscore && !identical(st, "excluded"))
        add_flag(s, "warn", sprintf("Noisy sample (robust z = %.2f > %g): calls are less reliable.", z, thr$qc_zscore))
    }
    nr <- get_status(s, "n_refs", NA_integer_)
    if (!is.na(nr) && nr < 5)
      add_flag(s, "warn", sprintf("Only %d reference sample(s) were used.", as.integer(nr)))
    if (isFALSE(get_status(s, "genotyped", NA)))
      add_flag(s, "fail", "Genotyping was skipped (too many calls): Q_SOME is missing and every call is LOW. The sample is probably too noisy to call.")
    if (!is.na(many_calls_cut) && call_counts[[s]] > many_calls_cut)
      add_flag(s, "warn", sprintf("Unusually many CNV calls (%d vs a cohort median of %g): possible noisy sample / batch effect.",
                                  as.integer(call_counts[[s]]), med_calls))
    if (!is.null(lcr) && s %in% lcr$SAMPLE) {
      r <- lcr[lcr$SAMPLE == s, , drop = FALSE]
      add_flag(s, "warn", sprintf(
        "%d run(s) of >= 3 consecutive targets with ~no reads were removed before calling (%s). A homozygous deletion or capture dropout there cannot be called (min_sample_count = %s).",
        nrow(r), paste(utils::head(paste0(r$CHR, ":", r$START, "-", r$END, if_else_gene(r$GENE)), 3), collapse = "; "),
        as.character(run$params$min_sample_count %||% 5)))
    }
    if (!is.null(sex_check) && s %in% sex_check$sample) {
      k <- match(s, sex_check$sample)
      if (isTRUE(sex_check$mismatch[k]))
        add_flag(s, "warn", sprintf("Annotated sex (%s) disagrees with the chrX/chrY coverage (looks %s; chrX ratio %.2f).",
                                    sex_check$annotated_sex[k], sex_check$inferred_sex[k], sex_check$x_ratio[k]))
      if (isTRUE(sex_check$mode_mismatch[k]) && !identical(st, "excluded"))
        add_flag(s, "warn", sprintf("modechrom = '%s' expects %s samples, but this sample looks %s: its chrX calls are not reliable.",
                                    run$modechrom, if (identical(run$modechrom, "XX")) "female" else "male",
                                    if (!is.na(sex_check$annotated_sex[k])) sex_check$annotated_sex[k] else sex_check$inferred_sex[k]))
    }
    own <- cnv[cnv$SAMPLE == s, , drop = FALSE]
    if (nrow(own) > 0) {
      g <- unique(trimws(unlist(strsplit(as.character(own$GENE), "[;, ]"))))
      hit <- intersect(g, low_genes)
      if (length(hit) > 0)
        add_flag(s, "info", paste0("CNV call(s) involve genes with known pseudogene / homology issues (always scored LOW): ",
                                   paste(hit, collapse = ", "), ". Consider orthogonal validation."))
    }
  }

  level_of <- function(f) if (any(vapply(f, function(x) x$level == "fail", logical(1)))) "fail" else
    if (any(vapply(f, function(x) x$level == "warn", logical(1)))) "warn" else "ok"
  sample_level <- vapply(flags, function(f) if (length(f)) level_of(f) else "ok", character(1))

  st_chr <- vapply(requested, function(s) as.character(get_status(s, "status", if (s %in% sample_cols) "called" else "unknown")), character(1))
  samples_df <- data.frame(
    sample = requested, status = st_chr, level = unname(sample_level),
    n_cnvs = as.integer(call_counts),
    n_del = vapply(requested, function(s) n_by(s, cnv$CNV == "DEL"), numeric(1)),
    n_dup = vapply(requested, function(s) n_by(s, cnv$CNV == "DUP"), numeric(1)),
    n_high = vapply(requested, function(s) n_by(s, conf == "HIGH"), numeric(1)),
    n_medium = vapply(requested, function(s) n_by(s, conf == "MEDIUM"), numeric(1)),
    n_low = vapply(requested, function(s) n_by(s, conf == "LOW"), numeric(1)),
    median_depth = round(as.numeric(depth), 1), total_reads = round(as.numeric(reads)),
    max_corr = round(as.numeric(max_corr), 4),
    noise_z = if (!is.null(noise)) round(noise$robust_z[match(requested, noise$sample)], 2) else NA_real_,
    n_refs = vapply(requested, function(s) as.numeric(get_status(s, "n_refs", length(refs[[s]]) %||% NA)), numeric(1)),
    ref_cor = round(vapply(requested, function(s) as.numeric(get_status(s, "ref_cor_median", NA_real_)), numeric(1)), 3),
    n_targets = vapply(requested, function(s) as.numeric(get_status(s, "n_targets", NA_integer_)), numeric(1)),
    n_dropped = vapply(requested, function(s) as.numeric(get_status(s, "n_low_count_dropped", NA_integer_)), numeric(1)),
    sex_inferred = if (!is.null(sex_check)) sex_check$inferred_sex[match(requested, sex_check$sample)] else NA_character_,
    stringsAsFactors = FALSE, row.names = NULL
  )
  if (is.null(sex_check)) samples_df$sex_inferred <- NA_character_

  # sample sheet columns (used for PCA colouring)
  st_tbl <- run$sample_table
  if (!is.null(sample_table) && is.character(sample_table) && file.exists(sample_table)) {
    tt <- tryCatch(utils::read.table(sample_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                                     check.names = FALSE, comment.char = ""), error = function(e) NULL)
    if (!is.null(tt) && "sample_name" %in% colnames(tt)) st_tbl <- tt
  } else if (is.data.frame(sample_table) && "sample_name" %in% colnames(sample_table)) {
    st_tbl <- sample_table
  }
  attrs <- list()
  if (!is.null(st_tbl)) {
    for (cn in setdiff(colnames(st_tbl), "sample_name")) {
      v <- as.character(st_tbl[[cn]][match(requested, st_tbl$sample_name)])
      nu <- length(unique(stats::na.omit(v)))
      if (nu >= 2L && nu <= 20L) attrs[[cn]] <- .rp_arr(ifelse(is.na(v), "NA", v))
    }
  }

  # ---- chromosomes ----------------------------------------------------------------
  if (!is.null(run$chromosome_status)) {
    chr_status <- run$chromosome_status
    chr_messages <- as.character(run$chromosome_messages %||% character(0))
  } else {
    an <- as.character(bed$chromosome)
    cs <- canope_chromosome_status(an, an, modechrom = "A", removeY = TRUE)
    chr_status <- cs$status; chr_messages <- cs$messages
    chr_messages <- c(chr_messages, "This workspace has no run summary (older CANOPE version): chromosomes missing from the analysed table are reported as absent, not as 'not analysed'.")
  }
  chr_depth <- run$chrom_depth
  chr_depth_json <- NULL
  if (!is.null(chr_depth) && !is.null(dim(chr_depth)) && ncol(chr_depth) > 0) {
    o <- order(.rp_chr_order(rownames(chr_depth)))
    chr_depth <- chr_depth[o, , drop = FALSE]
    chr_depth_json <- list(chromosomes = .rp_arr(rownames(chr_depth)), samples = .rp_arr(colnames(chr_depth)),
                           z = lapply(seq_len(nrow(chr_depth)), function(i) .rp_num(chr_depth[i, ], 3)))
  }
  inv <- run$chromosome_inventory
  analysed_chr <- unique(as.character(bed$chromosome))
  inv_df <- if (!is.null(inv)) {
    data.frame(chromosome = inv$chromosome, class = inv$class, targets_in_bed = inv$n_targets_input,
               targets_analysed = vapply(inv$chromosome, function(ch) sum(as.character(bed$chromosome) == ch), numeric(1)),
               analysed = inv$chromosome %in% analysed_chr, stringsAsFactors = FALSE)
  } else {
    tb <- table(as.character(bed$chromosome))
    data.frame(chromosome = names(tb), class = canope_chrom_class(names(tb)), targets_in_bed = as.integer(tb),
               targets_analysed = as.integer(tb), analysed = TRUE, stringsAsFactors = FALSE)
  }
  inv_df <- inv_df[order(.rp_chr_order(inv_df$chromosome)), , drop = FALSE]
  both_ok <- if (!is.null(run$chromosome_status)) {
    all(chr_status$analysed[chr_status$class %in% c("autosome", "X")] | !chr_status$in_input[chr_status$class %in% c("autosome", "X")]) &&
      any(chr_status$analysed[chr_status$class == "autosome"]) && any(chr_status$analysed[chr_status$class == "X"])
  } else NA
  ext <- tapply(as.numeric(bed$end), format_chr_label(bed$chromosome), max)
  ext <- ext[order(.rp_chr_order(names(ext)))]

  # ---- QC metrics file --------------------------------------------------------------
  qc <- NULL
  if (!is.null(qc_metrics_file) && file.exists(qc_metrics_file))
    qc <- tryCatch(data.table::fread(qc_metrics_file, header = TRUE, fill = TRUE, data.table = FALSE),
                   error = function(e) NULL)
  qc_json <- NULL
  if (!is.null(qc) && nrow(qc) > 0) {
    type_n <- as.data.frame(table(Type = qc$Type), stringsAsFactors = FALSE)
    keep_first <- qc$Type %in% c("Whole sample", "Missing Chromosome", "Chromosome not analysed")
    exon_rows <- which(!keep_first)
    capped <- length(exon_rows) > max_qc_rows
    rows <- c(which(keep_first), utils::head(exon_rows, max_qc_rows))
    gene_n <- NULL
    ex <- qc[!keep_first, , drop = FALSE]
    if (nrow(ex) > 0) {
      g <- as.data.frame(table(Gene = ex$Gene, Type = ex$Type), stringsAsFactors = FALSE)
      g <- g[g$Freq > 0, , drop = FALSE]
      g <- g[order(-g$Freq), , drop = FALSE]
      gene_n <- .rp_table(utils::head(g, 60))
    }
    qc_json <- list(types = .rp_table(type_n), rows = .rp_table(qc[rows, , drop = FALSE]),
                    n_total = nrow(qc), capped = capped, genes = gene_n)
  } else if (!is.null(qc)) {
    qc_json <- list(types = .rp_table(data.frame(Type = character(0), Freq = integer(0))),
                    rows = .rp_table(qc), n_total = 0L, capped = FALSE, genes = NULL)
  }

  # ---- PCA -----------------------------------------------------------------------------
  pca_json <- NULL
  if (length(sample_cols) >= 3L) {
    pca_json <- tryCatch({
      lm <- log2(as.matrix(counts[, sample_cols, drop = FALSE]) + 1)
      rv <- apply(lm, 1, stats::var, na.rm = TRUE)
      lm <- lm[!is.na(rv) & rv > 0, , drop = FALSE]
      if (nrow(lm) < 2L) NULL else {
        pc <- stats::prcomp(t(lm), scale. = TRUE, center = TRUE)
        k <- min(5L, ncol(pc$x))
        ve <- summary(pc)$importance[2, seq_len(k)] * 100
        list(samples = .rp_arr(sample_cols),
             pcs = setNames(lapply(seq_len(k), function(i) .rp_num(pc$x[, i], 3)), paste0("PC", seq_len(k))),
             var = .rp_num(ve, 1), n_targets = nrow(lm))
      }
    }, error = function(e) { message("[WARNING] report PCA skipped: ", conditionMessage(e)); NULL })
  }

  # ---- calls ------------------------------------------------------------------------------
  calls_json <- NULL; plots <- list(); n_flagged_bg <- 0L; n_plot_missing <- 0L
  if (nrow(cnv) > 0) {
    iv <- lapply(as.character(cnv$INTERVAL), parse_canope_interval)
    gene_hit <- vapply(as.character(cnv$GENE), function(g)
      any(trimws(unlist(strsplit(g, "[;, ]"))) %in% low_genes), logical(1))
    cdf <- data.frame(
      id = seq_len(nrow(cnv)), sample = as.character(cnv$SAMPLE),
      chr = format_chr_label(cnv$CHR),
      start = vapply(iv, function(p) as.numeric(p$start), numeric(1)),
      end = vapply(iv, function(p) as.numeric(p$end), numeric(1)),
      interval = as.character(cnv$INTERVAL), gene = as.character(cnv$GENE), type = as.character(cnv$CNV),
      kb = round(as.numeric(cnv$KB), 2), num_targ = as.integer(cnv$NUM_TARG),
      mlcn = as.integer(cnv$MLCN), cn_label = as.character(cnv$CN_label),
      q_some = as.integer(cnv$Q_SOME), num_refs = as.integer(cnv$NUM_REFS),
      confidence = ifelse(is.na(cnv$Confidence), "NA", toupper(as.character(cnv$Confidence))),
      homology_gene = gene_hit, targets = as.character(cnv$TARGETS),
      stringsAsFactors = FALSE)
    ord <- order(-.rp_conf_rank(cdf$confidence), -ifelse(is.na(cdf$q_some), -1L, cdf$q_some))
    embed <- utils::head(ord, max_call_plots)
    bed_chr <- format_chr_label(bed$chromosome)
    exon_index <- compute_exon_index(bed)
    cdf$bg_flag <- FALSE; cdf$has_plot <- FALSE
    for (i in embed) {
      w <- tryCatch(.rp_call_window(cnv[i, , drop = FALSE], bed, bed_chr, models, flank = flank,
                                    gene_gap = gene_gap, exon_index = exon_index),
                    error = function(e) NULL)
      if (!is.null(w)) {
        plots[[as.character(cdf$id[i])]] <- w
        cdf$has_plot[i] <- TRUE; cdf$bg_flag[i] <- isTRUE(w$bg$flag)
      }
    }
    n_flagged_bg <- sum(cdf$bg_flag)
    n_plot_missing <- sum(!cdf$has_plot)
    calls_json <- .rp_table(cdf)
  }

  # ---- global warnings -------------------------------------------------------------------------
  warns <- list()
  add_w <- function(level, category, text, sample = NULL)
    warns[[length(warns) + 1L]] <<- list(level = level, category = category, text = text, sample = sample)

  for (m in chr_messages) add_w("warn", "Chromosomes", m)
  if (isTRUE(both_ok == FALSE) || (!is.na(both_ok) && !both_ok)) {
    # message already included in chr_messages; keep a dedicated top-level entry
    add_w("fail", "Chromosomes",
          "Autosomes and gonosomes were NOT both analysed in this run: results cover only part of the genome. Do not read 'no call' on the missing chromosomes as 'no CNV'.")
  }
  fl <- samples_df$sample[samples_df$status == "failed"]
  if (length(fl)) add_w("fail", "Samples", sprintf("%d sample(s) FAILED calling and have no result: %s.", length(fl), paste(fl, collapse = ", ")))
  ex <- samples_df$sample[samples_df$status == "excluded"]
  if (length(ex)) add_w("warn", "Samples", sprintf("%d sample(s) were excluded before calling: %s.", length(ex), paste(ex, collapse = ", ")))
  n_an <- sum(samples_df$status == "called")
  if (n_an > 0 && n_an < 4) add_w("warn", "Samples", sprintf("Only %d sample(s) analysed: CANOPE needs >= 3 reference samples per sample, so sensitivity and specificity are limited.", n_an))
  for (s in requested) for (f in flags[[s]]) if (f$level %in% c("fail", "warn")) add_w(f$level, "Sample metrics", f$text, s)
  if (n_flagged_bg > 0)
    add_w("warn", "Calls", sprintf("%d call(s) have a high share of non-called exons outside the 95%% interval (background calibration flag): check them manually.", n_flagged_bg))
  if (n_plot_missing > 0 && !is.null(calls_json))
    add_w("info", "Calls", sprintf("%d call(s) have no embedded interactive plot (limit of %d plots, or no model data). The PDF plots are unaffected.", n_plot_missing, max_call_plots))
  if (!is.null(qc_json)) {
    nex <- sum(qc$Type %in% c("Whole exon", "Exon variability"))
    if (nex > 0) add_w("info", "QC metrics", sprintf("%d exon-level QC flag(s) (low depth / high variability); see the QC metrics tab.", nex))
  }
  if (identical(engine, "legacy_canoes")) add_w("info", "Method", "The legacy CANOES HMM engine was used.")

  # ---- log ------------------------------------------------------------------------------------------
  log_lines <- character(0)
  if (!is.null(log_file) && is.character(log_file) && file.exists(log_file))
    log_lines <- readLines(log_file, warn = FALSE)
  log_warn <- log_lines[startsWith(log_lines, "[WARNING]")]
  if (length(log_warn) > 0) add_w("info", "Pipeline log", sprintf("%d warning line(s) in the pipeline log (see the 'Run & log' tab).", length(log_warn)))
  si <- tryCatch(utils::capture.output(print(utils::sessionInfo())), error = function(e) character(0))

  prm <- run$params %||% list()
  prm <- lapply(prm, function(x) if (length(x) > 1) .rp_arr(x) else x)
  conf_rules <- list(
    high_q = settings$high_q_score %||% 50, med_q = settings$med_q_score %||% 20,
    high_refs = settings$high_num_refs %||% 10, med_refs = settings$med_num_refs %||% 5,
    high_targ = settings$high_num_targ %||% 1, med_targ = settings$med_num_targ %||% 1,
    genes = .rp_arr(low_genes))

  mode_label <- switch(as.character(run$modechrom %||% "?"),
                       A = "Autosomes (A)", XX = "chrX, female references (XX)", XY = "chrX/chrY, male references (XY)", "unknown")
  list(
    meta = list(
      title = title %||% paste("CANOPE CNV report -", sub("^CANOPE_?", "", prefix) %||% prefix),
      generated = format(Sys.time(), "%Y-%m-%d %H:%M"),
      version = run$canope_version %||% tryCatch(as.character(utils::packageVersion("CANOPE")), error = function(e) "dev"),
      engine = engine, engine_label = if (identical(engine, "legacy_canoes")) "Legacy CANOES (original HMM core)" else "CANOPE (modern HMM engine)",
      genome_build = run$genome_build %||% "unknown", modechrom = run$modechrom %||% NA_character_, mode_label = mode_label,
      legacy_workspace = is.null(run$sample_status)),
    thresholds = thr,
    conf_rules = conf_rules,
    params = prm,
    summary = list(
      samples_requested = length(requested), samples_called = n_an,
      samples_with_calls = length(unique(cnv$SAMPLE)), calls = nrow(cnv),
      targets = nrow(bed), warnings = length(warns),
      n_fail = sum(vapply(warns, function(w) w$level == "fail", logical(1))),
      n_warn = sum(vapply(warns, function(w) w$level == "warn", logical(1))),
      both_chromosome_classes = if (is.na(both_ok)) NULL else both_ok),
    samples = .rp_table(samples_df),
    flags = lapply(flags, function(f) unname(f)),
    sample_attrs = attrs,
    warnings = unname(warns),
    chromosomes = list(status = .rp_table(chr_status), messages = .rp_arr(chr_messages),
                       inventory = .rp_table(inv_df), depth = chr_depth_json,
                       extent = list(chr = .rp_arr(names(ext)), len = .rp_arr(as.numeric(ext))),
                       sex = if (!is.null(sex_check)) .rp_table(sex_check) else NULL),
    qc = qc_json,
    corr = corr,
    pca = pca_json,
    calls = calls_json,
    plots = plots,
    log = list(warnings = .rp_arr(utils::head(log_warn, 300L)),
               tail = .rp_arr(utils::tail(log_lines[!startsWith(log_lines, "[WARNING]")], 300L)),
               session = .rp_arr(utils::head(si, 120L)))
  )
}

# "chr:start-end (GENE)" helper for the low-coverage flag text
if_else_gene <- function(g) ifelse(is.na(g) | g == "", "", paste0(" (", g, ")"))


# ---------------------------------------------------------------------------------
# serialisation & page
# ---------------------------------------------------------------------------------

.rp_json <- function(data) {
  json <- as.character(jsonlite::toJSON(data, auto_unbox = TRUE, null = "null", na = "null",
                                        digits = NA, force = TRUE))
  # keep the JSON inert inside a <script> element
  json <- gsub("</", "<\\/", json, fixed = TRUE)
  gsub("<!--", "<\\!--", json, fixed = TRUE)
}

.rp_plotly_js <- function() {
  if (!requireNamespace("plotly", quietly = TRUE)) return(NULL)
  deps <- tryCatch(suppressWarnings(suppressMessages(htmltools::findDependencies(plotly::plot_ly()))),
                   error = function(e) NULL)
  dep <- Filter(function(d) identical(d$name, "plotly-main"), deps)
  if (!length(dep)) return(NULL)
  dep <- dep[[1L]]
  path <- file.path(dep$src$file, dep$script)
  if (!file.exists(path)) path <- system.file(dep$src$file, dep$script, package = dep$package %||% "plotly")
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  js <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  gsub("</script", "<\\/script", js, fixed = TRUE)
}

.rp_html <- function(meta, json, plotly_source = "inline") {
  css <- .CANOPE_REPORT_CSS
  js  <- .CANOPE_REPORT_JS
  plotly_tag <- NULL
  if (identical(plotly_source, "inline")) {
    pj <- .rp_plotly_js()
    if (!is.null(pj)) plotly_tag <- paste0("<script>\n", pj, "\n</script>\n")
    else message("[WARNING] plotly.js not found locally (install the 'plotly' package); using the jsDelivr CDN instead.")
  }
  if (is.null(plotly_tag))
    plotly_tag <- "<script src=\"https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js\"></script>\n"

  paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>", .rp_esc(meta$title), "</title>\n",
    "<style>\n", css, "\n</style>\n</head>\n<body>\n",
    "<header class=\"top\"><div class=\"top-inner\">\n",
    "<div class=\"title-row\"><h1>", .rp_esc(meta$title), "</h1>",
    "<span class=\"badge\">", .rp_esc(meta$engine_label), "</span>",
    "<span class=\"spacer\"></span>",
    "<button id=\"theme-btn\" class=\"btn\" type=\"button\" aria-label=\"Toggle colour theme\">Theme</button></div>\n",
    "<nav id=\"tabs\" class=\"tabs\" role=\"tablist\" aria-label=\"Report sections\"></nav>\n",
    "</div></header>\n<main id=\"panels\"></main>\n",
    "<noscript><p style=\"padding:16px\">This report needs JavaScript enabled.</p></noscript>\n",
    "<script id=\"canope-data\" type=\"application/json\">", json, "</script>\n",
    plotly_tag,
    "<script>\n", js, "\n</script>\n</body>\n</html>\n"
  )
}