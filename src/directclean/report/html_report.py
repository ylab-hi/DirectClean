"""
HTML report generator for DirectClean pipeline results.

Produces a self-contained interactive HTML file with embedded
Chart.js visualizations.  All data is serialized as JSON and
embedded directly in the page — no server required.

Usage::

    from directclean.report import HtmlReportGenerator

    generator = HtmlReportGenerator(
        report=pipeline_report,
        config=pipeline_config,
        input_fastq=Path("raw.fastq"),
        output_dir=Path("results/"),
    )
    generator.write(Path("results/directclean.report.html"))
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path

from directclean.pipeline import PipelineConfig, PipelineReport

logger = logging.getLogger(__name__)


class HtmlReportGenerator:
    """Generate a self-contained HTML report from pipeline results.

    The report includes interactive Chart.js charts for each pipeline
    stage plus an overall read-flow waterfall chart.

    Args:
        report:      PipelineReport with all stage statistics.
        config:      PipelineConfig with pipeline parameters.
        input_fastq: Path to the original input FASTQ.
        output_dir:  Output directory used by the pipeline.
        prefix:      Filename prefix used by the pipeline.
    """

    def __init__(
        self,
        report: PipelineReport,
        config: PipelineConfig,
        input_fastq: Path,
        output_dir: Path,
        prefix: str = "directclean",
    ) -> None:
        self.report = report
        self.config = config
        self.input_fastq = Path(input_fastq)
        self.output_dir = Path(output_dir)
        self.prefix = prefix

    # ------------------------------------------------------------------
    # Data serialization
    # ------------------------------------------------------------------

    def _build_report_data(self) -> dict:
        """Serialize all report data into a JSON-safe dictionary."""
        rpt = self.report
        cfg = self.config

        # -- Meta --
        try:
            from directclean import __version__
        except ImportError:
            __version__ = "unknown"

        data: dict = {
            "meta": {
                "version": __version__,
                "generated_at": datetime.now().isoformat(timespec="seconds"),
                "input_fastq": str(self.input_fastq),
                "output_dir": str(self.output_dir),
                "prefix": self.prefix,
                "elapsed_seconds": rpt.elapsed_seconds,
            },
            "config": {
                "threads": cfg.threads,
                "min_confidence": cfg.min_confidence,
                "context_window": cfg.context_window,
                "min_mapq": cfg.min_mapq,
                "max_edit_distance": cfg.adapter_config.max_edit_distance,
                "min_segment_length": cfg.adapter_config.min_segment_length,
                "scan_window": cfg.homopolymer_config.scan_window,
                "density_threshold": cfg.homopolymer_config.density_threshold,
                "min_run": cfg.homopolymer_config.min_run,
            },
        }

        # -- Stage 1: Breakinator --
        br = rpt.break_report
        if br is not None:
            data["breakinator"] = {
                "total_breakpoints": br.total_breakpoints,
                "foldback_count": br.foldback_count,
                "chimeric_count": br.chimeric_count,
                "pass_count": br.pass_count,
                "input_reads": br.input_reads,
                "removed_reads": br.removed_reads,
                "kept_reads": br.kept_reads,
            }

        # -- Stage 2: Restrander --
        rs = rpt.restrander_report
        if rs is not None:
            data["restrander"] = {
                "total_input": rs.total_input,
                "forward": rs.forward,
                "reverse": rs.reverse,
                "unknown": rs.unknown,
                "rtp_rtp": rs.rtp_rtp,
                "tso_tso": rs.tso_tso,
                "total_artefacts": rs.total_artefacts,
                "output_reads": rs.output_reads,
            }

        # -- Stage 2.5: Unknowns Rescue --
        ur = rpt.unknowns_rescue_report
        if ur is not None:
            data["unknowns_rescue"] = {
                "total_unknowns": ur.total_unknowns,
                "reads_with_adapter": ur.reads_with_adapter,
                "reads_without": ur.reads_without,
                "segments_produced": ur.segments_produced,
                "segments_discarded_short": ur.segments_discarded_short,
                "oriented_forward": ur.oriented_forward,
                "oriented_reverse": ur.oriented_reverse,
                "oriented_unknown": ur.oriented_unknown,
                "output_reads": ur.output_reads,
            }

        # -- Stage 3: Rescuer --
        rc = rpt.rescue_report
        if rc is not None:
            data["rescue"] = {
                "total_reads": rc.total_reads,
                "reads_with_internal": rc.reads_with_internal,
                "reads_without": rc.reads_without,
                "segments_rescued": rc.segments_rescued,
                "segments_discarded": rc.segments_discarded,
                "total_segments": rc.total_segments,
            }

        # -- Stage 5: Homopolymer Rescue --
        fr = rpt.filter_report
        if fr is not None:
            data["filter"] = {
                "total_chimeric_reads": fr.total_chimeric_reads,
                "artifact_reads": fr.artifact_reads,
                "clean_chimeric_reads": fr.clean_chimeric_reads,
                "total_junctions": fr.total_junctions,
                "artifact_junctions": fr.artifact_junctions,
                "skipped_junctions": fr.skipped_junctions,
                "total_reads_in_fastq": fr.total_reads_in_fastq,
                "output_reads": fr.output_reads,
                "segments_rescued": fr.segments_rescued,
                "segments_discarded": fr.segments_discarded,
            }

        # -- Waterfall data for the read-flow chart --
        waterfall = []
        if br is not None:
            waterfall.append(
                {
                    "label": "Raw Input",
                    "value": br.input_reads,
                    "type": "total",
                }
            )
            waterfall.append(
                {
                    "label": "Breakinator Removed",
                    "value": -br.removed_reads,
                    "type": "loss",
                }
            )

        if rs is not None:
            removed_by_restrander = rs.unknown + rs.total_artefacts
            waterfall.append(
                {
                    "label": "Restrander Removed",
                    "value": -removed_by_restrander,
                    "type": "loss",
                }
            )

        if ur is not None and ur.output_reads > 0:
            waterfall.append(
                {
                    "label": "Unknowns Rescued",
                    "value": ur.output_reads,
                    "type": "gain",
                }
            )

        if rc is not None:
            net_rescue = rc.total_segments - rc.total_reads
            if net_rescue != 0:
                waterfall.append(
                    {
                        "label": "Adapter Rescue (net)",
                        "value": net_rescue,
                        "type": "gain" if net_rescue > 0 else "loss",
                    }
                )

        if fr is not None:
            net_filter = fr.output_reads - fr.total_reads_in_fastq
            if net_filter != 0:
                waterfall.append(
                    {
                        "label": "Homopolymer Rescue (net)",
                        "value": net_filter,
                        "type": "gain" if net_filter > 0 else "loss",
                    }
                )
            waterfall.append(
                {
                    "label": "Final Output",
                    "value": fr.output_reads,
                    "type": "total",
                }
            )

        data["waterfall"] = waterfall

        return data

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def write(self, output_path: Path) -> None:
        """Render the HTML report and write to disk.

        Args:
            output_path: Destination file path for the HTML report.
        """
        data = self._build_report_data()
        data_json = json.dumps(data, indent=2)

        html = _HTML_TEMPLATE.replace("__REPORT_DATA__", data_json)

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(html, encoding="utf-8")

        logger.info(f"HTML report written: {output_path}")


# ======================================================================
# HTML Template
# ======================================================================

_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>DirectClean Report</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4"></script>
<style>
/* ---- CSS Reset & Variables ---- */
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
:root {
  --bg: #f4f6f9;
  --card-bg: #ffffff;
  --text: #1e293b;
  --text-muted: #64748b;
  --border: #e2e8f0;
  --primary: #2563eb;
  --primary-light: #dbeafe;
  --green: #16a34a;
  --green-light: #dcfce7;
  --red: #dc2626;
  --red-light: #fee2e2;
  --orange: #ea580c;
  --orange-light: #fff7ed;
  --purple: #7c3aed;
  --purple-light: #f3e8ff;
  --teal: #0d9488;
  --teal-light: #ccfbf1;
  --radius: 10px;
  --shadow: 0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.06);
  --shadow-md: 0 4px 6px rgba(0,0,0,0.07), 0 2px 4px rgba(0,0,0,0.06);
}

body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
  background: var(--bg);
  color: var(--text);
  line-height: 1.6;
  padding: 0;
}

/* ---- Header ---- */
.header {
  background: linear-gradient(135deg, #1e3a5f 0%, #2563eb 100%);
  color: #fff;
  padding: 2rem 2rem 1.5rem;
}
.header h1 { font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }
.header .subtitle { opacity: 0.85; font-size: 0.95rem; margin-top: 0.3rem; }
.header-meta {
  display: flex; flex-wrap: wrap; gap: 1.5rem;
  margin-top: 1rem; font-size: 0.85rem; opacity: 0.9;
}
.header-meta span { display: inline-flex; align-items: center; gap: 0.3rem; }

/* ---- Layout ---- */
.container { max-width: 1200px; margin: 0 auto; padding: 1.5rem; }

/* ---- Summary cards ---- */
.summary-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
  gap: 1rem;
  margin-bottom: 1.5rem;
}
.summary-card {
  background: var(--card-bg);
  border-radius: var(--radius);
  padding: 1.2rem 1.4rem;
  box-shadow: var(--shadow);
  border-left: 4px solid var(--primary);
}
.summary-card.green  { border-left-color: var(--green); }
.summary-card.red    { border-left-color: var(--red); }
.summary-card.orange { border-left-color: var(--orange); }
.summary-card.purple { border-left-color: var(--purple); }
.summary-card .label { font-size: 0.8rem; color: var(--text-muted); text-transform: uppercase; letter-spacing: 0.5px; font-weight: 600; }
.summary-card .value { font-size: 1.7rem; font-weight: 700; margin-top: 0.25rem; }
.summary-card .sub   { font-size: 0.8rem; color: var(--text-muted); margin-top: 0.15rem; }

/* ---- Cards ---- */
.card {
  background: var(--card-bg);
  border-radius: var(--radius);
  box-shadow: var(--shadow);
  margin-bottom: 1.25rem;
  overflow: hidden;
}
.card-header {
  display: flex; align-items: center; justify-content: space-between;
  padding: 0.9rem 1.4rem;
  cursor: pointer;
  user-select: none;
  border-bottom: 1px solid var(--border);
  background: #fafbfc;
  transition: background 0.15s;
}
.card-header:hover { background: #f1f5f9; }
.card-header h2 { font-size: 1rem; font-weight: 600; }
.card-header .badge {
  display: inline-block; font-size: 0.7rem; font-weight: 600;
  padding: 0.15rem 0.5rem; border-radius: 9999px;
  margin-left: 0.6rem;
}
.badge-blue   { background: var(--primary-light); color: var(--primary); }
.badge-green  { background: var(--green-light); color: var(--green); }
.badge-red    { background: var(--red-light); color: var(--red); }
.badge-orange { background: var(--orange-light); color: var(--orange); }
.badge-purple { background: var(--purple-light); color: var(--purple); }
.badge-teal   { background: var(--teal-light); color: var(--teal); }

.card-header .chevron {
  font-size: 0.8rem; color: var(--text-muted);
  transition: transform 0.2s;
}
.card-header.collapsed .chevron { transform: rotate(-90deg); }
.card-body { padding: 1.4rem; }
.card-body.collapsed { display: none; }

/* ---- Stage content layout ---- */
.stage-content {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 1.5rem;
  align-items: start;
}
@media (max-width: 768px) { .stage-content { grid-template-columns: 1fr; } }

/* ---- Tables ---- */
.stats-table { width: 100%; border-collapse: collapse; font-size: 0.88rem; }
.stats-table td { padding: 0.45rem 0.6rem; border-bottom: 1px solid var(--border); }
.stats-table td:first-child { color: var(--text-muted); font-weight: 500; white-space: nowrap; }
.stats-table td:last-child  { text-align: right; font-weight: 600; font-variant-numeric: tabular-nums; }
.stats-table tr.separator td { border-bottom: 2px solid var(--border); padding: 0.2rem; }
.stats-table tr:last-child td { border-bottom: none; }

/* ---- Chart containers ---- */
.chart-wrap {
  position: relative;
  width: 100%;
  max-width: 320px;
  margin: 0 auto;
}
.chart-wrap-wide {
  position: relative;
  width: 100%;
  max-height: 400px;
}
.chart-wrap-bar {
  position: relative;
  width: 100%;
  max-height: 280px;
}
.full-width { grid-column: 1 / -1; }

/* ---- Config table ---- */
.config-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
  gap: 0.5rem 2rem;
}
.config-grid .cfg-item {
  display: flex; justify-content: space-between;
  padding: 0.4rem 0; border-bottom: 1px solid var(--border);
  font-size: 0.88rem;
}
.config-grid .cfg-item .cfg-key { color: var(--text-muted); }
.config-grid .cfg-item .cfg-val { font-weight: 600; font-variant-numeric: tabular-nums; }

/* ---- Footer ---- */
.footer {
  text-align: center; padding: 1.5rem;
  font-size: 0.8rem; color: var(--text-muted);
  border-top: 1px solid var(--border);
  margin-top: 1rem;
}

/* ---- Print ---- */
@media print {
  body { background: #fff; }
  .header { background: #1e3a5f !important; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
  .card-body.collapsed { display: block !important; }
  .card-header .chevron { display: none; }
}
</style>
</head>
<body>

<!-- Header -->
<div class="header">
  <h1>DirectClean Pipeline Report</h1>
  <div class="subtitle" id="headerSubtitle"></div>
  <div class="header-meta" id="headerMeta"></div>
</div>

<div class="container">

  <!-- Executive Summary -->
  <div class="summary-grid" id="summaryGrid"></div>

  <!-- Waterfall Chart -->
  <div class="card">
    <div class="card-header" onclick="toggleCard(this)">
      <h2>Read Flow Through Pipeline <span class="badge badge-blue">Waterfall</span></h2>
      <span class="chevron">&#9660;</span>
    </div>
    <div class="card-body">
      <div class="chart-wrap-wide">
        <canvas id="waterfallChart"></canvas>
      </div>
    </div>
  </div>

  <!-- Stage 1: Breakinator -->
  <div class="card" id="breakCard" style="display:none">
    <div class="card-header" onclick="toggleCard(this)">
      <h2>Stage 1: Breakinator <span class="badge badge-red">Foldback Removal</span></h2>
      <span class="chevron">&#9660;</span>
    </div>
    <div class="card-body">
      <div class="stage-content">
        <div><table class="stats-table" id="breakTable"></table></div>
        <div class="chart-wrap"><canvas id="breakChart"></canvas></div>
      </div>
      <div style="margin-top:1.2rem">
        <div class="chart-wrap-bar"><canvas id="breakBarChart"></canvas></div>
      </div>
    </div>
  </div>

  <!-- Stage 2: Restrander -->
  <div class="card" id="restranderCard" style="display:none">
    <div class="card-header" onclick="toggleCard(this)">
      <h2>Stage 2: Restrander <span class="badge badge-purple">Strand Correction</span></h2>
      <span class="chevron">&#9660;</span>
    </div>
    <div class="card-body">
      <div class="stage-content">
        <div><table class="stats-table" id="restranderTable"></table></div>
        <div class="chart-wrap"><canvas id="restranderChart"></canvas></div>
      </div>
      <div style="margin-top:1.2rem">
        <div class="chart-wrap-bar"><canvas id="restranderBarChart"></canvas></div>
      </div>
    </div>
  </div>

  <!-- Stage 2.5: Unknowns Rescue -->
  <div class="card" id="unknownsCard" style="display:none">
    <div class="card-header" onclick="toggleCard(this)">
      <h2>Stage 2.5: Unknowns Rescue <span class="badge badge-teal">Recovery</span></h2>
      <span class="chevron">&#9660;</span>
    </div>
    <div class="card-body">
      <div class="stage-content">
        <div><table class="stats-table" id="unknownsTable"></table></div>
        <div class="chart-wrap"><canvas id="unknownsChart"></canvas></div>
      </div>
      <div style="margin-top:1.2rem">
        <div class="chart-wrap-bar"><canvas id="unknownsBarChart"></canvas></div>
      </div>
    </div>
  </div>

  <!-- Stage 3: Rescuer -->
  <div class="card" id="rescueCard" style="display:none">
    <div class="card-header" onclick="toggleCard(this)">
      <h2>Stage 3: Rescuer <span class="badge badge-orange">Adapter Detection</span></h2>
      <span class="chevron">&#9660;</span>
    </div>
    <div class="card-body">
      <div class="stage-content">
        <div><table class="stats-table" id="rescueTable"></table></div>
        <div class="chart-wrap"><canvas id="rescueChart"></canvas></div>
      </div>
      <div style="margin-top:1.2rem">
        <div class="chart-wrap-bar"><canvas id="rescueBarChart"></canvas></div>
      </div>
    </div>
  </div>

  <!-- Stage 5: Homopolymer Rescue -->
  <div class="card" id="filterCard" style="display:none">
    <div class="card-header" onclick="toggleCard(this)">
      <h2>Stage 5: Homopolymer Rescue <span class="badge badge-green">RT Artifact Chopping</span></h2>
      <span class="chevron">&#9660;</span>
    </div>
    <div class="card-body">
      <div class="stage-content">
        <div><table class="stats-table" id="filterTable"></table></div>
        <div class="chart-wrap"><canvas id="filterChart"></canvas></div>
      </div>
      <div style="margin-top:1.2rem">
        <div class="chart-wrap-bar"><canvas id="filterBarChart"></canvas></div>
      </div>
    </div>
  </div>

  <!-- Configuration -->
  <div class="card">
    <div class="card-header collapsed" onclick="toggleCard(this)">
      <h2>Configuration <span class="badge badge-blue">Parameters</span></h2>
      <span class="chevron">&#9660;</span>
    </div>
    <div class="card-body collapsed" id="configBody"></div>
  </div>

</div>

<div class="footer" id="footer"></div>

<script>
// ---- Embedded Report Data ----
const DATA = __REPORT_DATA__;

// ---- Helpers ----
function fmt(n) {
  if (n === undefined || n === null) return 'N/A';
  return n.toLocaleString();
}
function pct(n, total) {
  if (!total || total === 0) return 'N/A';
  return (n / total * 100).toFixed(1) + '%';
}
function fmtTime(s) {
  if (!s) return 'N/A';
  const h = Math.floor(s / 3600);
  const m = Math.floor((s % 3600) / 60);
  const sec = Math.floor(s % 60);
  if (h > 0) return h + 'h ' + m + 'm ' + sec + 's';
  if (m > 0) return m + 'm ' + sec + 's';
  return sec + 's';
}

function toggleCard(headerEl) {
  headerEl.classList.toggle('collapsed');
  const body = headerEl.nextElementSibling;
  if (body) body.classList.toggle('collapsed');
}

function statsRow(label, value) {
  return '<tr><td>' + label + '</td><td>' + value + '</td></tr>';
}
function separatorRow() {
  return '<tr class="separator"><td colspan="2"></td></tr>';
}

// Chart.js defaults
Chart.defaults.font.family = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif";
Chart.defaults.font.size = 12;
Chart.defaults.plugins.legend.labels.usePointStyle = true;
Chart.defaults.plugins.legend.labels.padding = 14;

// Color palette
const C = {
  blue: '#2563eb', blueA: 'rgba(37,99,235,0.7)',
  green: '#16a34a', greenA: 'rgba(22,163,74,0.7)',
  red: '#dc2626', redA: 'rgba(220,38,38,0.7)',
  orange: '#ea580c', orangeA: 'rgba(234,88,12,0.7)',
  purple: '#7c3aed', purpleA: 'rgba(124,58,237,0.7)',
  teal: '#0d9488', tealA: 'rgba(13,148,136,0.7)',
  gray: '#94a3b8', grayA: 'rgba(148,163,184,0.7)',
  yellow: '#ca8a04', yellowA: 'rgba(202,138,4,0.7)',
};

// ---- Header ----
document.getElementById('headerSubtitle').textContent =
  'Version ' + DATA.meta.version + '  |  Generated ' + DATA.meta.generated_at;
document.getElementById('headerMeta').innerHTML =
  '<span>Input: ' + DATA.meta.input_fastq + '</span>' +
  '<span>Output: ' + DATA.meta.output_dir + '</span>' +
  '<span>Prefix: ' + DATA.meta.prefix + '</span>';

// ---- Executive Summary ----
(function() {
  const grid = document.getElementById('summaryGrid');
  const inputReads = DATA.breakinator ? DATA.breakinator.input_reads : 0;
  const outputReads = DATA.filter ? DATA.filter.output_reads : 0;
  const totalRemoved = inputReads - outputReads;
  const retentionPct = inputReads > 0 ? (outputReads / inputReads * 100).toFixed(1) + '%' : 'N/A';

  grid.innerHTML =
    '<div class="summary-card">' +
      '<div class="label">Input Reads</div>' +
      '<div class="value">' + fmt(inputReads) + '</div>' +
      '<div class="sub">Raw FASTQ</div>' +
    '</div>' +
    '<div class="summary-card green">' +
      '<div class="label">Final Output</div>' +
      '<div class="value">' + fmt(outputReads) + '</div>' +
      '<div class="sub">Cleaned + rescued reads</div>' +
    '</div>' +
    '<div class="summary-card red">' +
      '<div class="label">Net Change</div>' +
      '<div class="value">' + fmt(totalRemoved > 0 ? -totalRemoved : '+' + Math.abs(totalRemoved)) + '</div>' +
      '<div class="sub">Retention: ' + retentionPct + '</div>' +
    '</div>' +
    '<div class="summary-card purple">' +
      '<div class="label">Processing Time</div>' +
      '<div class="value">' + fmtTime(DATA.meta.elapsed_seconds) + '</div>' +
      '<div class="sub">' + (DATA.config.threads || '?') + ' threads</div>' +
    '</div>';
})();

// ---- Waterfall Chart ----
(function() {
  const wf = DATA.waterfall || [];
  if (wf.length === 0) return;

  const labels = wf.map(w => w.label);
  const bgColors = [];
  const values = [];
  const bases = [];    // invisible base for floating bars
  let running = 0;

  for (let i = 0; i < wf.length; i++) {
    const w = wf[i];
    if (w.type === 'total') {
      bases.push(0);
      values.push(w.value);
      running = w.value;
      bgColors.push(C.blue);
    } else if (w.type === 'loss') {
      const absVal = Math.abs(w.value);
      running += w.value;  // subtract
      bases.push(running);
      values.push(absVal);
      bgColors.push(C.red);
    } else {
      // gain
      bases.push(running);
      values.push(w.value);
      running += w.value;
      bgColors.push(C.green);
    }
  }

  new Chart(document.getElementById('waterfallChart'), {
    type: 'bar',
    data: {
      labels: labels,
      datasets: [
        {
          label: 'Base',
          data: bases,
          backgroundColor: 'transparent',
          borderWidth: 0,
          barPercentage: 0.6,
        },
        {
          label: 'Value',
          data: values,
          backgroundColor: bgColors,
          borderRadius: 4,
          barPercentage: 0.6,
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        x: { stacked: true, grid: { display: false } },
        y: {
          stacked: true,
          beginAtZero: true,
          ticks: { callback: v => fmt(v) },
          title: { display: true, text: 'Read Count' },
        },
      },
      plugins: {
        legend: { display: false },
        tooltip: {
          callbacks: {
            label: function(ctx) {
              if (ctx.datasetIndex === 0) return null;
              const w = wf[ctx.dataIndex];
              const prefix = w.type === 'loss' ? '-' : (w.type === 'gain' ? '+' : '');
              return prefix + fmt(Math.abs(w.value)) + ' reads';
            },
          },
        },
      },
    },
  });
})();

// ---- Stage 1: Breakinator ----
(function() {
  const d = DATA.breakinator;
  if (!d) return;
  document.getElementById('breakCard').style.display = '';

  // Table
  const t = document.getElementById('breakTable');
  t.innerHTML =
    statsRow('Breakpoints total', fmt(d.total_breakpoints)) +
    statsRow('  Foldback', fmt(d.foldback_count)) +
    statsRow('  Chimeric', fmt(d.chimeric_count)) +
    statsRow('  Pass', fmt(d.pass_count)) +
    separatorRow() +
    statsRow('Input reads', fmt(d.input_reads)) +
    statsRow('Foldback reads removed', fmt(d.removed_reads) + ' (' + pct(d.removed_reads, d.input_reads) + ')') +
    statsRow('Reads kept', fmt(d.kept_reads));

  // Doughnut: breakpoint classification
  new Chart(document.getElementById('breakChart'), {
    type: 'doughnut',
    data: {
      labels: ['Foldback', 'Chimeric', 'Pass'],
      datasets: [{
        data: [d.foldback_count, d.chimeric_count, d.pass_count],
        backgroundColor: [C.red, C.orange, C.green],
        borderWidth: 2, borderColor: '#fff',
      }],
    },
    options: {
      responsive: true,
      plugins: {
        title: { display: true, text: 'Breakpoint Classification' },
        tooltip: {
          callbacks: { label: ctx => ctx.label + ': ' + fmt(ctx.raw) + ' (' + pct(ctx.raw, d.total_breakpoints) + ')' }
        },
      },
    },
  });

  // Bar: kept vs removed
  new Chart(document.getElementById('breakBarChart'), {
    type: 'bar',
    data: {
      labels: ['Reads'],
      datasets: [
        { label: 'Kept', data: [d.kept_reads], backgroundColor: C.greenA, borderRadius: 4 },
        { label: 'Removed', data: [d.removed_reads], backgroundColor: C.redA, borderRadius: 4 },
      ],
    },
    options: {
      indexAxis: 'y', responsive: true, maintainAspectRatio: false,
      scales: {
        x: { stacked: true, ticks: { callback: v => fmt(v) } },
        y: { stacked: true, display: false },
      },
      plugins: {
        title: { display: true, text: 'Reads Kept vs Removed' },
        tooltip: { callbacks: { label: ctx => ctx.dataset.label + ': ' + fmt(ctx.raw) } },
      },
    },
  });
})();

// ---- Stage 2: Restrander ----
(function() {
  const d = DATA.restrander;
  if (!d) return;
  document.getElementById('restranderCard').style.display = '';

  const t = document.getElementById('restranderTable');
  t.innerHTML =
    statsRow('Total input reads', fmt(d.total_input)) +
    statsRow('Forward (+)', fmt(d.forward)) +
    statsRow('Reverse (-)', fmt(d.reverse)) +
    statsRow('Unknown (?)', fmt(d.unknown) + ' (' + pct(d.unknown, d.total_input) + ')') +
    statsRow('Artefacts', fmt(d.total_artefacts) + ' (' + pct(d.total_artefacts, d.total_input) + ')') +
    statsRow('  RTP-RTP', fmt(d.rtp_rtp)) +
    statsRow('  TSO-TSO', fmt(d.tso_tso)) +
    separatorRow() +
    statsRow('Output reads', fmt(d.output_reads) + ' (' + pct(d.output_reads, d.total_input) + ')');

  // Pie: classification
  new Chart(document.getElementById('restranderChart'), {
    type: 'doughnut',
    data: {
      labels: ['Forward (+)', 'Reverse (-)', 'Unknown (?)', 'RTP-RTP', 'TSO-TSO'],
      datasets: [{
        data: [d.forward, d.reverse, d.unknown, d.rtp_rtp, d.tso_tso],
        backgroundColor: [C.green, C.blue, C.gray, C.red, C.orange],
        borderWidth: 2, borderColor: '#fff',
      }],
    },
    options: {
      responsive: true,
      plugins: {
        title: { display: true, text: 'Read Classification' },
        tooltip: {
          callbacks: { label: ctx => ctx.label + ': ' + fmt(ctx.raw) + ' (' + pct(ctx.raw, d.total_input) + ')' }
        },
      },
    },
  });

  // Stacked bar: output breakdown
  new Chart(document.getElementById('restranderBarChart'), {
    type: 'bar',
    data: {
      labels: ['Reads'],
      datasets: [
        { label: 'Forward', data: [d.forward], backgroundColor: C.greenA, borderRadius: 4 },
        { label: 'Reverse', data: [d.reverse], backgroundColor: C.blueA, borderRadius: 4 },
        { label: 'Unknown', data: [d.unknown], backgroundColor: C.grayA, borderRadius: 4 },
        { label: 'Artefacts', data: [d.total_artefacts], backgroundColor: C.redA, borderRadius: 4 },
      ],
    },
    options: {
      indexAxis: 'y', responsive: true, maintainAspectRatio: false,
      scales: {
        x: { stacked: true, ticks: { callback: v => fmt(v) } },
        y: { stacked: true, display: false },
      },
      plugins: {
        title: { display: true, text: 'Read Composition' },
        tooltip: { callbacks: { label: ctx => ctx.dataset.label + ': ' + fmt(ctx.raw) } },
      },
    },
  });
})();

// ---- Stage 2.5: Unknowns Rescue ----
(function() {
  const d = DATA.unknowns_rescue;
  if (!d) return;
  document.getElementById('unknownsCard').style.display = '';

  const t = document.getElementById('unknownsTable');
  t.innerHTML =
    statsRow('Total unknowns scanned', fmt(d.total_unknowns)) +
    statsRow('Reads with internal adapter', fmt(d.reads_with_adapter) + ' (' + pct(d.reads_with_adapter, d.total_unknowns) + ')') +
    statsRow('Reads without (skipped)', fmt(d.reads_without)) +
    separatorRow() +
    statsRow('Sub-reads produced', fmt(d.segments_produced)) +
    statsRow('Discarded (too short)', fmt(d.segments_discarded_short)) +
    statsRow('Oriented forward', fmt(d.oriented_forward)) +
    statsRow('Oriented reverse (RC)', fmt(d.oriented_reverse)) +
    statsRow('Orientation unknown', fmt(d.oriented_unknown)) +
    separatorRow() +
    statsRow('Output reads', fmt(d.output_reads));

  // Doughnut: orientation
  new Chart(document.getElementById('unknownsChart'), {
    type: 'doughnut',
    data: {
      labels: ['Forward', 'Reverse (RC)', 'Unknown'],
      datasets: [{
        data: [d.oriented_forward, d.oriented_reverse, d.oriented_unknown],
        backgroundColor: [C.green, C.blue, C.gray],
        borderWidth: 2, borderColor: '#fff',
      }],
    },
    options: {
      responsive: true,
      plugins: {
        title: { display: true, text: 'Sub-read Orientation' },
        tooltip: {
          callbacks: {
            label: ctx => {
              const total = d.oriented_forward + d.oriented_reverse + d.oriented_unknown;
              return ctx.label + ': ' + fmt(ctx.raw) + ' (' + pct(ctx.raw, total) + ')';
            }
          }
        },
      },
    },
  });

  // Bar: segments
  new Chart(document.getElementById('unknownsBarChart'), {
    type: 'bar',
    data: {
      labels: ['Sub-reads'],
      datasets: [
        { label: 'Output reads', data: [d.output_reads], backgroundColor: C.greenA, borderRadius: 4 },
        { label: 'Discarded (short)', data: [d.segments_discarded_short], backgroundColor: C.orangeA, borderRadius: 4 },
        { label: 'Orientation unknown', data: [d.oriented_unknown], backgroundColor: C.grayA, borderRadius: 4 },
      ],
    },
    options: {
      indexAxis: 'y', responsive: true, maintainAspectRatio: false,
      scales: {
        x: { stacked: true, ticks: { callback: v => fmt(v) } },
        y: { stacked: true, display: false },
      },
      plugins: {
        title: { display: true, text: 'Sub-read Outcome' },
        tooltip: { callbacks: { label: ctx => ctx.dataset.label + ': ' + fmt(ctx.raw) } },
      },
    },
  });
})();

// ---- Stage 3: Rescuer ----
(function() {
  const d = DATA.rescue;
  if (!d) return;
  document.getElementById('rescueCard').style.display = '';

  const t = document.getElementById('rescueTable');
  t.innerHTML =
    statsRow('Total reads processed', fmt(d.total_reads)) +
    statsRow('Reads with internal adapter', fmt(d.reads_with_internal) + ' (' + pct(d.reads_with_internal, d.total_reads) + ')') +
    statsRow('Reads passed unchanged', fmt(d.reads_without)) +
    separatorRow() +
    statsRow('Segments rescued', fmt(d.segments_rescued)) +
    statsRow('Segments discarded', fmt(d.segments_discarded)) +
    statsRow('Total output reads', fmt(d.total_segments));

  // Doughnut: adapter detection
  new Chart(document.getElementById('rescueChart'), {
    type: 'doughnut',
    data: {
      labels: ['With internal adapter', 'No adapter (unchanged)'],
      datasets: [{
        data: [d.reads_with_internal, d.reads_without],
        backgroundColor: [C.orange, C.green],
        borderWidth: 2, borderColor: '#fff',
      }],
    },
    options: {
      responsive: true,
      plugins: {
        title: { display: true, text: 'Internal Adapter Detection' },
        tooltip: {
          callbacks: { label: ctx => ctx.label + ': ' + fmt(ctx.raw) + ' (' + pct(ctx.raw, d.total_reads) + ')' }
        },
      },
    },
  });

  // Bar: segments
  new Chart(document.getElementById('rescueBarChart'), {
    type: 'bar',
    data: {
      labels: ['Segments'],
      datasets: [
        { label: 'Rescued', data: [d.segments_rescued], backgroundColor: C.greenA, borderRadius: 4 },
        { label: 'Discarded', data: [d.segments_discarded], backgroundColor: C.redA, borderRadius: 4 },
      ],
    },
    options: {
      indexAxis: 'y', responsive: true, maintainAspectRatio: false,
      scales: {
        x: { stacked: true, ticks: { callback: v => fmt(v) } },
        y: { stacked: true, display: false },
      },
      plugins: {
        title: { display: true, text: 'Chopped Segment Outcome' },
        tooltip: { callbacks: { label: ctx => ctx.dataset.label + ': ' + fmt(ctx.raw) } },
      },
    },
  });
})();

// ---- Stage 5: Homopolymer Rescue ----
(function() {
  const d = DATA.filter;
  if (!d) return;
  document.getElementById('filterCard').style.display = '';

  const t = document.getElementById('filterTable');
  t.innerHTML =
    statsRow('Chimeric reads examined', fmt(d.total_chimeric_reads)) +
    statsRow('Artifact reads (chopped)', fmt(d.artifact_reads) + ' (' + pct(d.artifact_reads, d.total_chimeric_reads) + ')') +
    statsRow('Clean chimeric reads', fmt(d.clean_chimeric_reads)) +
    separatorRow() +
    statsRow('Total junctions', fmt(d.total_junctions)) +
    statsRow('Artifact junctions', fmt(d.artifact_junctions)) +
    statsRow('Skipped junctions (filtered)', fmt(d.skipped_junctions)) +
    separatorRow() +
    statsRow('Input FASTQ reads', fmt(d.total_reads_in_fastq)) +
    statsRow('Segments rescued', fmt(d.segments_rescued)) +
    statsRow('Segments discarded (<100bp)', fmt(d.segments_discarded)) +
    statsRow('Output reads', fmt(d.output_reads));

  // Doughnut: chimeric classification
  new Chart(document.getElementById('filterChart'), {
    type: 'doughnut',
    data: {
      labels: ['Artifact (chopped)', 'Clean chimeric'],
      datasets: [{
        data: [d.artifact_reads, d.clean_chimeric_reads],
        backgroundColor: [C.red, C.green],
        borderWidth: 2, borderColor: '#fff',
      }],
    },
    options: {
      responsive: true,
      plugins: {
        title: { display: true, text: 'Chimeric Read Classification' },
        tooltip: {
          callbacks: { label: ctx => ctx.label + ': ' + fmt(ctx.raw) + ' (' + pct(ctx.raw, d.total_chimeric_reads) + ')' }
        },
      },
    },
  });

  // Grouped bar: junctions and segments
  new Chart(document.getElementById('filterBarChart'), {
    type: 'bar',
    data: {
      labels: ['Junctions', 'Segments'],
      datasets: [
        {
          label: 'Artifact / Rescued',
          data: [d.artifact_junctions, d.segments_rescued],
          backgroundColor: C.redA,
          borderRadius: 4,
        },
        {
          label: 'Skipped / Discarded',
          data: [d.skipped_junctions, d.segments_discarded],
          backgroundColor: C.orangeA,
          borderRadius: 4,
        },
        {
          label: 'Clean / N/A',
          data: [d.total_junctions - d.artifact_junctions - d.skipped_junctions, 0],
          backgroundColor: C.greenA,
          borderRadius: 4,
        },
      ],
    },
    options: {
      responsive: true, maintainAspectRatio: false,
      scales: {
        x: { grid: { display: false } },
        y: { ticks: { callback: v => fmt(v) }, title: { display: true, text: 'Count' } },
      },
      plugins: {
        title: { display: true, text: 'Junctions & Segments Breakdown' },
        tooltip: { callbacks: { label: ctx => ctx.dataset.label + ': ' + fmt(ctx.raw) } },
      },
    },
  });
})();

// ---- Configuration ----
(function() {
  const c = DATA.config;
  if (!c) return;
  const items = [
    ['Threads', c.threads],
    ['Min confidence', c.min_confidence],
    ['Max edit distance', c.max_edit_distance],
    ['Min segment length', c.min_segment_length + ' bp'],
    ['Context window', c.context_window + ' bp'],
    ['Min MAPQ', c.min_mapq],
    ['Scan window', c.scan_window + ' bp'],
    ['Density threshold', c.density_threshold],
    ['Min homopolymer run', c.min_run + ' bp'],
  ];
  let html = '<div class="config-grid">';
  items.forEach(([k, v]) => {
    html += '<div class="cfg-item"><span class="cfg-key">' + k + '</span><span class="cfg-val">' + v + '</span></div>';
  });
  html += '</div>';
  document.getElementById('configBody').innerHTML = html;
})();

// ---- Footer ----
document.getElementById('footer').innerHTML =
  'DirectClean v' + DATA.meta.version +
  ' &mdash; Generated ' + DATA.meta.generated_at +
  ' &mdash; <a href="https://github.com/your-org/directclean" style="color:var(--primary)">Documentation</a>';
</script>
</body>
</html>
"""
