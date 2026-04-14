#!/usr/bin/env python3
"""
lncRNA-seq Pipeline Report Generator
=====================================
Generates HTML (single-file, offline-ready with inline base64 images)
and Markdown reports from the pipeline output directory.

Usage
-----
python generate_report.py \\
    --config  <path/to/config.yaml> \\
    --outdir  <pipeline_output_dir> \\
    --report  <report_output_dir>   \\
    [--no-inline-images]
"""

import argparse
import base64
import csv
import datetime
import json
import os
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Optional YAML import (standard library fallback via simple parser)
# ---------------------------------------------------------------------------
try:
    import yaml  # PyYAML
    def load_yaml(path):
        with open(path) as fh:
            return yaml.safe_load(fh)
except ImportError:
    def load_yaml(path):
        """Minimal YAML loader sufficient for config.yaml (no anchors)."""
        import ast
        result = {}
        with open(path) as fh:
            lines = fh.readlines()
        # Use Python's json for anything that survived; otherwise parse line by line.
        # This is a best-effort fallback – PyYAML is always preferred.
        stack = [result]
        indent_stack = [-1]
        for raw in lines:
            stripped = raw.rstrip()
            if not stripped or stripped.lstrip().startswith("#"):
                continue
            indent = len(raw) - len(raw.lstrip())
            key_val = stripped.lstrip()
            if ":" in key_val:
                key, _, val = key_val.partition(":")
                key = key.strip()
                val = val.strip()
                while indent <= indent_stack[-1]:
                    stack.pop()
                    indent_stack.pop()
                parent = stack[-1]
                if isinstance(parent, dict):
                    if val == "" or val.startswith("#"):
                        new = {}
                        parent[key] = new
                        stack.append(new)
                        indent_stack.append(indent)
                    else:
                        val = val.split("#")[0].strip()
                        try:
                            parent[key] = ast.literal_eval(val)
                        except Exception:
                            parent[key] = val
        return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def safe_read(path):
    """Return file contents or None if the file does not exist."""
    p = Path(path)
    if p.is_file():
        try:
            return p.read_text(errors="replace")
        except Exception:
            return None
    return None


def image_to_base64(path, fmt=None):
    """Return a data-URI string for an image, or None if unavailable."""
    p = Path(path)
    if not p.is_file():
        return None
    ext = (fmt or p.suffix.lower().lstrip("."))
    mime_map = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg",
                "svg": "image/svg+xml", "gif": "image/gif", "pdf": "application/pdf"}
    mime = mime_map.get(ext, "image/png")
    data = base64.b64encode(p.read_bytes()).decode()
    return f"data:{mime};base64,{data}"


def html_img(src_uri, alt="", style="max-width:100%;border-radius:4px;"):
    if src_uri is None:
        return f'<p class="placeholder">⚠ 图片未找到 / Image not available: {alt}</p>'
    return f'<img src="{src_uri}" alt="{alt}" style="{style}">'


def parse_tsv(path):
    """Parse a TSV/TXT file and return (headers, rows)."""
    content = safe_read(path)
    if content is None:
        return [], []
    lines = [l for l in content.splitlines() if l.strip()]
    if not lines:
        return [], []
    reader = csv.reader(lines, delimiter="\t")
    rows = list(reader)
    return rows[0], rows[1:]


def parse_csv_file(path):
    """Parse a CSV file; return (headers, rows)."""
    content = safe_read(path)
    if content is None:
        return [], []
    lines = [l for l in content.splitlines() if l.strip()]
    if not lines:
        return [], []
    # Strip surrounding quotes from each cell
    reader = csv.reader(lines)
    rows = list(reader)
    return rows[0], rows[1:]


# ---------------------------------------------------------------------------
# MultiQC data parsers
# ---------------------------------------------------------------------------

def load_multiqc_table(path):
    """Return a dict: sample_name -> {column: value}."""
    headers, rows = parse_tsv(path)
    result = {}
    for row in rows:
        if not row:
            continue
        sample = row[0]
        result[sample] = dict(zip(headers[1:], row[1:]))
    return result


def parse_quality_scores(path):
    """
    Parse fastqc_per_sequence_quality_scores_plot.txt.
    Returns dict: sample -> {q20_pct, q30_pct}
    """
    content = safe_read(path)
    if not content:
        return {}
    result = {}
    for line in content.splitlines():
        if not line.strip() or line.startswith("Sample"):
            continue
        parts = line.split("\t")
        sample = parts[0]
        counts_by_q = {}
        for token in parts[1:]:
            token = token.strip()
            if not token:
                continue
            m = re.match(r"\((\d+\.?\d*),\s*(\d+\.?\d*)\)", token)
            if m:
                q = float(m.group(1))
                c = float(m.group(2))
                counts_by_q[q] = c
        total = sum(counts_by_q.values())
        if total == 0:
            continue
        q20 = sum(v for k, v in counts_by_q.items() if k >= 20)
        q30 = sum(v for k, v in counts_by_q.items() if k >= 30)
        result[sample] = {
            "q20_pct": round(q20 / total * 100, 2),
            "q30_pct": round(q30 / total * 100, 2),
        }
    return result


# ---------------------------------------------------------------------------
# QC table builder
# ---------------------------------------------------------------------------

def build_qc_table(config, output_dir):
    """
    Build a list of dicts (one per sample) for the three-line QC table.
    Columns: sample, condition, raw_reads, clean_reads, gc_pct, q20_pct,
             q30_pct, mapping_rate, unique_mapping_rate
    Falls back to "N/A" when data is missing.
    """
    multiqc_dir = Path(output_dir) / "results" / "multiqc_report_data"

    gen_stats = load_multiqc_table(multiqc_dir / "multiqc_general_stats.txt")
    star_stats = load_multiqc_table(multiqc_dir / "multiqc_star.txt")
    cutadapt_stats = load_multiqc_table(multiqc_dir / "multiqc_cutadapt.txt")
    fastqc_stats = load_multiqc_table(multiqc_dir / "multiqc_fastqc.txt")
    quality_scores = parse_quality_scores(
        multiqc_dir / "fastqc_per_sequence_quality_scores_plot.txt")

    samples = config["samples"]
    rows = []
    for sample_name, sample_info in samples.items():
        condition = sample_info.get("condition", "N/A")
        prefix = sample_info.get("prefix", sample_name)

        # --- raw reads ---
        raw_reads = "N/A"
        for key in [sample_name, prefix, f"{prefix}_R1"]:
            d = cutadapt_stats.get(key, {})
            if d.get("r_processed"):
                try:
                    # The value in multiqc is in millions
                    raw_reads = f"{float(d['r_processed']):.3f}M"
                except Exception:
                    pass
                break
            d = gen_stats.get(key, {})
            if d.get("fastqc-total_sequences"):
                try:
                    raw_reads = f"{float(d['fastqc-total_sequences']):.3f}M"
                except Exception:
                    pass
                break

        # --- clean reads ---
        clean_reads = "N/A"
        for key in [sample_name, prefix]:
            d = cutadapt_stats.get(key, {})
            if d.get("r_written"):
                try:
                    clean_reads = f"{float(d['r_written']):.3f}M"
                except Exception:
                    pass
                break

        # --- GC% ---
        gc_pct = "N/A"
        for key in [sample_name, f"{prefix}_R1", prefix]:
            d = fastqc_stats.get(key) or gen_stats.get(key, {})
            v = (d.get("%GC") or d.get("fastqc-percent_gc"))
            if v:
                try:
                    gc_pct = f"{float(v):.1f}%"
                except Exception:
                    pass
                break

        # --- Q20 / Q30 ---
        q20, q30 = "N/A", "N/A"
        for key in [sample_name, f"{prefix}_R1", prefix]:
            qs = quality_scores.get(key)
            if qs:
                q20 = f"{qs['q20_pct']:.2f}%"
                q30 = f"{qs['q30_pct']:.2f}%"
                break

        # --- mapping rate ---
        map_rate = "N/A"
        unique_rate = "N/A"
        for key in [sample_name, prefix]:
            d = star_stats.get(key, {})
            if d.get("uniquely_mapped_percent"):
                try:
                    unique_rate = f"{float(d['uniquely_mapped_percent']):.2f}%"
                    mapped_pct = d.get("mapped_percent") or d.get("uniquely_mapped_percent")
                    map_rate = f"{float(mapped_pct):.2f}%"
                except Exception:
                    pass
                break
            d = gen_stats.get(key, {})
            if d.get("star-uniquely_mapped_percent"):
                try:
                    unique_rate = f"{float(d['star-uniquely_mapped_percent']):.2f}%"
                    map_rate = unique_rate
                except Exception:
                    pass
                break

        rows.append({
            "sample": sample_name,
            "condition": condition,
            "raw_reads": raw_reads,
            "clean_reads": clean_reads,
            "gc_pct": gc_pct,
            "q20_pct": q20,
            "q30_pct": q30,
            "mapping_rate": map_rate,
            "unique_mapping_rate": unique_rate,
        })
    return rows


# ---------------------------------------------------------------------------
# Collect lncRNA / DESeq2 plots
# ---------------------------------------------------------------------------

def collect_lncrna_plots(output_dir):
    """Return list of (filename, title) for lncRNA downstream analysis plots."""
    plot_dir = Path(output_dir) / "lncRNA_analysis_output" / "plots"
    mapping = {
        "01_log2FC_distribution.png": "log₂FC Distribution",
        "02_neglog10pvalue_distribution.png": "-log₁₀(P-value) Distribution",
        "03_neglog10padj_distribution.png": "-log₁₀(Padj) Distribution",
        "04_baseMean_violin_by_class.png": "Base Mean by Gene Class",
        "05_ExonCount_distribution.png": "Exon Count Distribution",
        "06_GeneLength_distribution.png": "Gene Length Distribution",
        "07_MA_plot.png": "MA Plot",
        "08_Volcano_lncRNA.png": "lncRNA Volcano Plot",
        "12_Expr_vs_absLFC.png": "Expression vs |log₂FC|",
        "13_Expr_vs_neglog10padj.png": "Expression vs -log₁₀(Padj)",
        "15_GSEA_top_enrichment_curve.png": "GSEA Top Enrichment Curves",
    }
    results = []
    for fname, title in mapping.items():
        p = plot_dir / fname
        if p.is_file():
            results.append((str(p), title))
    return results


# ---------------------------------------------------------------------------
# Collect ERCC results
# ---------------------------------------------------------------------------

def collect_ercc_results(config, output_dir):
    """Return list of (sample, curve_pdf_path, mpc_table_path)."""
    samples = list(config["samples"].keys())
    results = []
    for s in samples:
        pdf = Path(output_dir) / "results" / f"{s}_ERCC_standard_curve.pdf"
        mpc = Path(output_dir) / "results" / f"{s}_mpc.txt"
        results.append({
            "sample": s,
            "curve_pdf": str(pdf) if pdf.is_file() else None,
            "mpc_table": str(mpc) if mpc.is_file() else None,
        })
    return results


def read_mpc_table(path):
    """Read a _mpc.txt file; return (headers, rows)."""
    if path is None:
        return [], []
    return parse_tsv(path)


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

CSS = """
/* ===== Reset & base ===== */
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: 'Segoe UI', Arial, 'PingFang SC', 'Microsoft YaHei', sans-serif;
    font-size: 14px;
    line-height: 1.65;
    color: #2c3e50;
    background: #f7f9fc;
}
a { color: #1a73e8; text-decoration: none; }
a:hover { text-decoration: underline; }

/* ===== Layout ===== */
.report-header {
    background: linear-gradient(135deg, #1a237e 0%, #283593 50%, #1565c0 100%);
    color: #fff;
    padding: 36px 48px 32px;
}
.report-header h1 { font-size: 1.9em; font-weight: 700; letter-spacing: .5px; margin-bottom: 6px; }
.report-header .meta { font-size: .88em; opacity: .85; }
.report-header .meta span { margin-right: 18px; }

.nav-bar {
    background: #fff;
    border-bottom: 1px solid #e0e6ed;
    padding: 0 48px;
    position: sticky;
    top: 0;
    z-index: 100;
    display: flex;
    gap: 2px;
    overflow-x: auto;
    box-shadow: 0 1px 4px rgba(0,0,0,.06);
}
.nav-bar a {
    display: inline-block;
    padding: 12px 16px;
    font-size: .82em;
    font-weight: 600;
    color: #455a64;
    border-bottom: 3px solid transparent;
    white-space: nowrap;
}
.nav-bar a:hover { color: #1565c0; border-bottom-color: #1565c0; text-decoration: none; }

.content { max-width: 1100px; margin: 0 auto; padding: 32px 24px 80px; }

/* ===== Section cards ===== */
.section {
    background: #fff;
    border-radius: 8px;
    box-shadow: 0 1px 6px rgba(0,0,0,.07);
    margin-bottom: 28px;
    overflow: hidden;
}
.section-header {
    background: #e8eaf6;
    padding: 14px 24px 12px;
    border-left: 5px solid #3949ab;
    display: flex;
    align-items: center;
    gap: 10px;
}
.section-header h2 {
    font-size: 1.05em;
    font-weight: 700;
    color: #1a237e;
    margin: 0;
}
.section-number {
    background: #3949ab;
    color: #fff;
    border-radius: 50%;
    width: 24px;
    height: 24px;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: .75em;
    font-weight: 700;
    flex-shrink: 0;
}
.section-body { padding: 22px 26px; }

/* ===== Sub-sections ===== */
.subsection { margin-bottom: 22px; }
.subsection h3 {
    font-size: .95em;
    color: #37474f;
    font-weight: 700;
    border-bottom: 1px solid #eceff1;
    padding-bottom: 5px;
    margin-bottom: 12px;
}

/* ===== Three-line table (三线表) ===== */
.three-line-table {
    width: 100%;
    border-collapse: collapse;
    font-size: .84em;
    margin: 8px 0;
}
.three-line-table thead tr th {
    border-top: 2px solid #2c3e50;
    border-bottom: 1px solid #2c3e50;
    padding: 8px 10px;
    text-align: center;
    font-weight: 700;
    background: transparent;
    color: #1a237e;
    white-space: nowrap;
}
.three-line-table tbody tr td {
    padding: 7px 10px;
    text-align: center;
    border: none;
    color: #37474f;
}
.three-line-table tbody tr:nth-child(even) td { background: #f8f9fa; }
.three-line-table tbody tr:last-child td {
    border-bottom: 2px solid #2c3e50;
}
.three-line-table caption {
    caption-side: bottom;
    font-size: .8em;
    color: #78909c;
    margin-top: 5px;
    text-align: left;
}

/* ===== Generic table ===== */
.data-table {
    width: 100%;
    border-collapse: collapse;
    font-size: .82em;
    margin: 8px 0;
    overflow-x: auto;
}
.data-table th {
    background: #e8eaf6;
    color: #283593;
    padding: 7px 10px;
    text-align: left;
    font-weight: 700;
    border: 1px solid #c5cae9;
}
.data-table td {
    padding: 6px 10px;
    border: 1px solid #eceff1;
    color: #455a64;
    vertical-align: top;
}
.data-table tbody tr:hover td { background: #f3f4f9; }

/* ===== Params grid ===== */
.params-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(240px, 1fr));
    gap: 12px;
    margin: 4px 0 8px;
}
.param-item {
    background: #f5f7ff;
    border: 1px solid #e3e8f5;
    border-radius: 6px;
    padding: 10px 14px;
}
.param-label { font-size: .77em; color: #78909c; font-weight: 600; text-transform: uppercase; }
.param-value { font-size: .9em; color: #1a237e; font-weight: 700; margin-top: 2px; word-break: break-all; }

/* ===== Image gallery ===== */
.figure-row {
    display: flex;
    flex-wrap: wrap;
    gap: 18px;
    margin-top: 8px;
}
.figure-box {
    flex: 1 1 44%;
    min-width: 260px;
    border: 1px solid #e0e6ed;
    border-radius: 6px;
    overflow: hidden;
    background: #fafbff;
}
.figure-box img { width: 100%; display: block; }
.figure-caption {
    font-size: .78em;
    color: #607d8b;
    padding: 5px 10px;
    text-align: center;
    font-style: italic;
    border-top: 1px solid #e8eaf6;
}

/* ===== Badges ===== */
.badge {
    display: inline-block;
    padding: 2px 10px;
    border-radius: 12px;
    font-size: .76em;
    font-weight: 700;
    letter-spacing: .3px;
}
.badge-on  { background: #e8f5e9; color: #2e7d32; border: 1px solid #a5d6a7; }
.badge-off { background: #fce4ec; color: #c62828; border: 1px solid #ef9a9a; }

/* ===== Placeholder ===== */
.placeholder {
    color: #90a4ae;
    font-style: italic;
    padding: 10px 0;
    font-size: .85em;
}

/* ===== File index ===== */
.file-index { list-style: none; padding: 0; }
.file-index li {
    padding: 6px 0;
    border-bottom: 1px dashed #eceff1;
    font-size: .84em;
}
.file-index li:last-child { border: none; }
.file-path { font-family: monospace; color: #37474f; }
.file-desc { color: #78909c; margin-left: 8px; }

/* ===== Footer ===== */
.report-footer {
    text-align: center;
    color: #90a4ae;
    font-size: .78em;
    padding: 20px;
    margin-top: 8px;
    border-top: 1px solid #eceff1;
}

/* ===== Print ===== */
@media print {
    .nav-bar { display: none; }
    .section { box-shadow: none; border: 1px solid #ccc; page-break-inside: avoid; }
    body { background: #fff; font-size: 11px; }
}
"""

NAV_ITEMS = [
    ("overview",     "项目概览"),
    ("qc",           "数据质控"),
    ("alignment",    "比对与定量"),
    ("lncrna",       "lncRNA 分析"),
    ("ercc",         "ERCC 分析"),
    ("deseq2",       "差异表达"),
    ("file-index",   "结果文件"),
]


def section_open(sec_id, title, num, hidden=False):
    hidden_attr = ' style="display:none"' if hidden else ""
    return (
        f'<div class="section" id="{sec_id}"{hidden_attr}>\n'
        f'  <div class="section-header">'
        f'<span class="section-number">{num}</span>'
        f'<h2>{title}</h2>'
        f'</div>\n'
        f'  <div class="section-body">\n'
    )


def section_close():
    return "  </div>\n</div>\n"


def three_line_table(headers, rows, caption=""):
    head_cells = "".join(f"<th>{h}</th>" for h in headers)
    body_rows = ""
    for row in rows:
        cells = "".join(f"<td>{c}</td>" for c in row)
        body_rows += f"<tr>{cells}</tr>\n"
    cap = f"<caption>{caption}</caption>" if caption else ""
    return (
        f'<div style="overflow-x:auto">'
        f'<table class="three-line-table">{cap}'
        f'<thead><tr>{head_cells}</tr></thead>'
        f'<tbody>{body_rows}</tbody>'
        f'</table></div>'
    )


def data_table_html(headers, rows, max_rows=20):
    if not headers:
        return '<p class="placeholder">暂无数据 / No data available.</p>'
    head_cells = "".join(f"<th>{h}</th>" for h in headers)
    body_rows = ""
    for row in rows[:max_rows]:
        cells = "".join(f"<td>{c}</td>" for c in row)
        body_rows += f"<tr>{cells}</tr>\n"
    note = ""
    if len(rows) > max_rows:
        note = (f'<p style="font-size:.8em;color:#78909c;margin-top:4px">'
                f'仅展示前 {max_rows} 行，共 {len(rows)} 行数据。'
                f'</p>')
    return (
        f'<div style="overflow-x:auto">'
        f'<table class="data-table">'
        f'<thead><tr>{head_cells}</tr></thead>'
        f'<tbody>{body_rows}</tbody>'
        f'</table></div>{note}'
    )


def param_grid(items):
    """items: list of (label, value)"""
    cells = "".join(
        f'<div class="param-item">'
        f'<div class="param-label">{label}</div>'
        f'<div class="param-value">{value}</div>'
        f'</div>'
        for label, value in items
    )
    return f'<div class="params-grid">{cells}</div>'


def figure_row(images, inline=True):
    """images: list of (path, caption)"""
    boxes = ""
    for path, cap in images:
        src = image_to_base64(path) if inline else path
        boxes += (
            f'<div class="figure-box">'
            f'{html_img(src, cap)}'
            f'<div class="figure-caption">{cap}</div>'
            f'</div>'
        )
    return f'<div class="figure-row">{boxes}</div>'


# ---------------------------------------------------------------------------
# Build HTML report
# ---------------------------------------------------------------------------

def build_html(config, output_dir, inline_images=True):
    """Return the complete HTML string."""
    report_cfg = config.get("report", {})
    title = report_cfg.get("title", "lncRNA-seq Analysis Report")
    author = report_cfg.get("author", "—")
    project = report_cfg.get("project", "—")
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    samples = config["samples"]
    contrast = config["analysis"].get("contrast", ["Treated", "Control"])
    ercc_enabled = bool(config.get("ercc", {}).get("use_ercc", False))
    deseq2_enabled = bool(config.get("analysis", {}).get("deseq2", True))
    use_rsem = bool(config.get("analysis", {}).get("use_rsem", False))

    # Collect data
    qc_rows = build_qc_table(config, output_dir)
    lncrna_plots = collect_lncrna_plots(output_dir) if deseq2_enabled else []
    ercc_results = collect_ercc_results(config, output_dir) if ercc_enabled else []

    # DESeq2 files
    contrast_str = f"{contrast[0]}_vs_{contrast[1]}"
    deseq2_csv = Path(output_dir) / "results" / f"{contrast_str}.deseq2_results.csv"
    volcano_png = Path(output_dir) / "results" / f"{contrast_str}.volcano_plot.png"
    go_png = Path(output_dir) / "results" / f"{contrast_str}.GO_dotplot.png"

    # lncRNA tables
    qc_summary_path = Path(output_dir) / "lncRNA_analysis_output" / "tables" / "QC_summary.csv"
    lnc_deg_path = Path(output_dir) / "lncRNA_analysis_output" / "tables" / "DEG_annotated_lncRNA_only.csv"
    lnc_sig_path = Path(output_dir) / "lncRNA_analysis_output" / "tables" / "SIG_lncRNA_padj0.05_lfc1.csv"

    # ----- Navigation -----
    nav_links = ""
    for sec_id, label in NAV_ITEMS:
        visible = True
        if sec_id == "ercc" and not ercc_enabled:
            visible = False
        if sec_id in ("deseq2", "lncrna") and not deseq2_enabled:
            visible = False
        if visible:
            nav_links += f'<a href="#{sec_id}">{label}</a>'

    # ----- SECTION 1: Project Overview -----
    sample_rows = []
    for s, info in samples.items():
        pe = "Paired-end" if info.get("R2") else "Single-end"
        sample_rows.append([s, info.get("condition", "—"), pe])

    ref_cfg = config.get("ref", {})
    params_cfg = config.get("params", {})

    params_items = [
        ("参考基因组 (FASTA)", Path(ref_cfg.get("fasta", "N/A")).name),
        ("注释文件 (GTF)", Path(ref_cfg.get("gtf", "N/A")).name),
        ("比较对照", f"{contrast[0]} vs {contrast[1]}"),
        ("FDR 阈值", config["analysis"].get("fdr_threshold", "N/A")),
        ("|log₂FC| 阈值", config["analysis"].get("lfc_threshold", "N/A")),
        ("链特异性", params_cfg.get("featurecounts_strandedness", "N/A")),
        ("CPU 线程数", config.get("threads", "N/A")),
        ("RSEM 定量", '<span class="badge badge-on">启用</span>' if use_rsem else '<span class="badge badge-off">未启用（featureCounts）</span>'),
        ("ERCC 质控", '<span class="badge badge-on">启用</span>' if ercc_enabled else '<span class="badge badge-off">未启用</span>'),
        ("DESeq2 分析", '<span class="badge badge-on">启用</span>' if deseq2_enabled else '<span class="badge badge-off">未启用</span>'),
        ("报告生成时间", now),
        ("物种注释库", config["analysis"].get("r_annotation_db", "N/A")),
    ]
    if use_rsem and ref_cfg.get("rsem_grp"):
        params_items.insert(2, ("RSEM 参考", Path(ref_cfg["rsem_grp"]).name))

    sec1 = section_open("overview", "项目与参数概览", 1)
    sec1 += '<div class="subsection"><h3>样本信息</h3>'
    sec1 += three_line_table(
        ["样本名称", "实验组", "测序类型"],
        sample_rows,
        caption="表1  样本基本信息"
    )
    sec1 += "</div>"
    sec1 += '<div class="subsection"><h3>关键参数</h3>'
    sec1 += param_grid(params_items)
    sec1 += "</div>"
    sec1 += section_close()

    # ----- SECTION 2: QC -----
    qc_headers = [
        "样本名称", "实验组", "原始 Reads (M)", "过滤后 Reads (M)",
        "GC%", "Q20%", "Q30%", "比对率", "唯一比对率"
    ]
    qc_data_rows = [
        [r["sample"], r["condition"], r["raw_reads"], r["clean_reads"],
         r["gc_pct"], r["q20_pct"], r["q30_pct"],
         r["mapping_rate"], r["unique_mapping_rate"]]
        for r in qc_rows
    ]

    sec2 = section_open("qc", "数据质控", 2)
    sec2 += '<div class="subsection"><h3>质控指标汇总</h3>'
    sec2 += three_line_table(
        qc_headers, qc_data_rows,
        caption="表2  质控关键指标（三线表）。原始 Reads 单位为百万（M）；Q20/Q30 由 FastQC 序列质量分布估算；比对率来自 STAR。"
    )
    sec2 += "</div>"
    # MultiQC plot mention
    multiqc_report = Path(output_dir) / "results" / "multiqc_report.html"
    if multiqc_report.is_file():
        sec2 += (
            '<p style="margin-top:6px;font-size:.84em;color:#455a64">'
            '💡 详细质控报告可在 MultiQC 报告中查看：'
            f'<code>results/multiqc_report.html</code></p>'
        )
    sec2 += section_close()

    # ----- SECTION 3: Alignment & Quantification -----
    sec3 = section_open("alignment", "比对与定量结果", 3)
    # STAR stats summary
    multiqc_dir = Path(output_dir) / "results" / "multiqc_report_data"
    star_headers, star_rows = parse_tsv(multiqc_dir / "multiqc_star.txt")
    if star_headers:
        display_cols = ["Sample", "total_reads", "uniquely_mapped",
                        "uniquely_mapped_percent", "multimapped",
                        "multimapped_percent", "unmapped_tooshort_percent"]
        idx = [i for i, h in enumerate(star_headers) if h in display_cols
               or star_headers[i] == "Sample"]
        filtered_h = [star_headers[i] for i in idx]
        filtered_r = [[row[i] if i < len(row) else "" for i in idx]
                      for row in star_rows
                      if row and row[0] in samples]
        sec3 += '<div class="subsection"><h3>STAR 比对统计</h3>'
        sec3 += three_line_table(
            ["样本", "总 Reads", "唯一比对 Reads", "唯一比对率 %",
             "多重比对", "多重比对率 %", "未比对（过短）%"],
            filtered_r,
            caption="表3  STAR 比对统计摘要"
        )
        sec3 += "</div>"

    # featureCounts
    fc_headers, fc_rows = parse_tsv(multiqc_dir / "multiqc_featurecounts.txt")
    if fc_headers:
        display_fc = ["Sample", "Total", "Assigned", "percent_assigned"]
        idx_fc = [i for i, h in enumerate(fc_headers) if h in display_fc]
        filtered_fch = [fc_headers[i] for i in idx_fc]
        filtered_fcr = [[row[i] if i < len(row) else "" for i in idx_fc]
                        for row in fc_rows]
        sec3 += '<div class="subsection"><h3>featureCounts 定量统计</h3>'
        sec3 += three_line_table(
            ["样本", "总 Reads", "已分配 Reads", "分配率 %"],
            filtered_fcr,
            caption="表4  featureCounts 基因定量统计"
        )
        sec3 += "</div>"

    # TPM matrix info
    tpm_path = Path(output_dir) / "results" / "RSEM.gene_tpm.symbol.tsv"
    counts_path = Path(output_dir) / "results" / "featureCounts.gene_counts.symbol.tsv"
    quant_grid_items = [
        ("featureCounts 计数矩阵", "results/featureCounts.gene_counts.symbol.tsv" if counts_path.is_file() else "未产出"),
    ]
    if use_rsem:
        quant_grid_items.insert(0, ("RSEM TPM 矩阵", "results/RSEM.gene_tpm.symbol.tsv" if tpm_path.is_file() else "未产出"))
    sec3 += '<div class="subsection"><h3>定量结果文件</h3>'
    sec3 += param_grid(quant_grid_items)
    sec3 += "</div>"
    sec3 += section_close()

    # ----- SECTION 4: lncRNA Analysis (conditional) -----
    sec4_hidden = not deseq2_enabled
    sec4 = section_open("lncrna", "lncRNA 相关结果汇总", 4, hidden=sec4_hidden)
    if not deseq2_enabled:
        sec4 += '<p class="placeholder">⚠ lncRNA 下游分析未启用（analysis.deseq2: false）。</p>'
    else:
        # QC summary table
        qs_h, qs_r = parse_csv_file(qc_summary_path)
        if qs_h:
            col_labels = {
                "total_genes": "总基因数",
                "genes_with_padj": "有 Padj 基因",
                "sig_all_padj005": "显著差异基因 (Padj<0.05)",
                "total_lnc": "lncRNA 总数",
                "lnc_with_padj": "有 Padj lncRNA",
                "sig_lnc_padj005": "显著 lncRNA (Padj<0.05)",
                "sig_lnc_strict_n": "显著 lncRNA (严格)",
                "sig_lnc_moderate_n": "显著 lncRNA (宽松)",
                "up_lnc_strict": "上调 lncRNA (严格)",
                "down_lnc_strict": "下调 lncRNA (严格)",
            }
            display_h = [col_labels.get(h, h) for h in qs_h]
            sec4 += '<div class="subsection"><h3>lncRNA 质控摘要</h3>'
            sec4 += three_line_table(display_h, qs_r, caption="表5  lncRNA 分析质控摘要")
            sec4 += "</div>"

        # lncRNA DEG table preview
        lnc_h, lnc_r = parse_csv_file(lnc_deg_path)
        if lnc_h:
            sec4 += '<div class="subsection"><h3>lncRNA 差异分析结果（前 20 行）</h3>'
            sec4 += data_table_html(lnc_h, lnc_r, max_rows=20)
            sec4 += "</div>"

        # Plots
        if lncrna_plots:
            sec4 += '<div class="subsection"><h3>lncRNA 可视化图</h3>'
            sec4 += figure_row(lncrna_plots, inline=inline_images)
            sec4 += "</div>"
        else:
            sec4 += '<p class="placeholder">lncRNA 分析图片未找到。</p>'

    sec4 += section_close()

    # ----- SECTION 5: ERCC Analysis (conditional) -----
    sec5_hidden = not ercc_enabled
    sec5 = section_open("ercc", "ERCC 定量分析", 5, hidden=sec5_hidden)
    if not ercc_enabled:
        sec5 += '<p class="placeholder">⚠ ERCC 分析未启用（ercc.use_ercc: false）。</p>'
    else:
        ercc_params = config.get("ercc_params", {})
        sec5 += '<div class="subsection"><h3>ERCC 参数</h3>'
        sec5 += param_grid([
            ("稀释倍数", ercc_params.get("dilution", "N/A")),
            ("ERCC 加入量 (µL)", ercc_params.get("ercc_ul", "N/A")),
            ("每细胞 RNA 量 (pg)", ercc_params.get("rna_pg_per_cell", "N/A")),
            ("总 RNA (µg)", ercc_params.get("total_rna_ug", "N/A")),
        ])
        sec5 += "</div>"

        for er in ercc_results:
            sec5 += f'<div class="subsection"><h3>样本：{er["sample"]}</h3>'
            # Standard curve (PDF — show download link)
            if er["curve_pdf"]:
                rel = Path(er["curve_pdf"]).name
                sec5 += (
                    f'<p>📄 ERCC 标准曲线（PDF）：<a href="{er["curve_pdf"]}" '
                    f'target="_blank">{rel}</a>'
                    f'<span style="color:#78909c;font-size:.8em"> （PDF 文件，点击下载）</span></p>'
                )
            else:
                sec5 += '<p class="placeholder">ERCC 标准曲线 PDF 未找到。</p>'

            mpc_h, mpc_r = read_mpc_table(er["mpc_table"])
            if mpc_h:
                sec5 += '<p style="margin-top:10px;font-weight:600">MPC 表格：</p>'
                sec5 += data_table_html(mpc_h, mpc_r)
            else:
                sec5 += '<p class="placeholder">MPC 表格未找到。</p>'
            sec5 += "</div>"

    sec5 += section_close()

    # ----- SECTION 6: DESeq2 Differential Expression (conditional) -----
    sec6_hidden = not deseq2_enabled
    sec6 = section_open("deseq2", "DESeq2 差异表达分析", 6, hidden=sec6_hidden)
    if not deseq2_enabled:
        sec6 += '<p class="placeholder">⚠ DESeq2 差异分析未启用（analysis.deseq2: false）。</p>'
    else:
        sec6 += '<div class="subsection"><h3>差异分析参数</h3>'
        sec6 += param_grid([
            ("比较组", f"{contrast[0]} vs {contrast[1]}"),
            ("FDR 阈值", config["analysis"].get("fdr_threshold", 0.05)),
            ("|log₂FC| 阈值", config["analysis"].get("lfc_threshold", 1)),
        ])
        sec6 += "</div>"

        # Volcano plot
        sec6 += '<div class="subsection"><h3>可视化图</h3>'
        img_pairs = []
        if volcano_png.is_file():
            img_pairs.append((str(volcano_png), "火山图 / Volcano Plot"))
        if go_png.is_file():
            img_pairs.append((str(go_png), "GO 富集点图 / GO Dotplot"))
        if img_pairs:
            sec6 += figure_row(img_pairs, inline=inline_images)
        else:
            sec6 += '<p class="placeholder">差异分析图片未找到。</p>'
        sec6 += "</div>"

        # DESeq2 result table preview
        deseq2_h, deseq2_r = parse_csv_file(deseq2_csv)
        if deseq2_h:
            sec6 += '<div class="subsection"><h3>DESeq2 结果（前 20 行）</h3>'
            sec6 += data_table_html(deseq2_h, deseq2_r, max_rows=20)
            sec6 += "</div>"
        else:
            sec6 += '<p class="placeholder">DESeq2 结果文件未找到。</p>'

    sec6 += section_close()

    # ----- SECTION 7: Output File Index -----
    def fstatus(p):
        return "✅" if Path(p).is_file() else "❌"

    file_items = [
        (f"{output_dir}/results/multiqc_report.html", "MultiQC 质控汇总报告"),
        (f"{output_dir}/results/featureCounts.gene_counts.symbol.tsv", "featureCounts 基因计数矩阵"),
    ]
    if use_rsem:
        file_items.append((f"{output_dir}/results/RSEM.gene_tpm.symbol.tsv", "RSEM TPM 表达矩阵"))
    if deseq2_enabled:
        file_items += [
            (f"{output_dir}/results/{contrast_str}.deseq2_results.csv", "DESeq2 差异分析结果"),
            (f"{output_dir}/results/{contrast_str}.volcano_plot.png", "火山图"),
            (f"{output_dir}/results/{contrast_str}.GO_enrichment.csv", "GO 富集分析结果"),
            (f"{output_dir}/results/{contrast_str}.GO_dotplot.png", "GO 富集点图"),
            (f"{output_dir}/lncRNA_analysis_output/tables/QC_summary.csv", "lncRNA 分析质控摘要"),
            (f"{output_dir}/lncRNA_analysis_output/tables/DEG_annotated_lncRNA_only.csv", "lncRNA 差异分析注释结果"),
            (f"{output_dir}/lncRNA_analysis_output/tables/SIG_lncRNA_padj0.05_lfc1.csv", "显著差异 lncRNA 列表"),
        ]
    if ercc_enabled:
        for s in samples:
            file_items += [
                (f"{output_dir}/results/{s}_ERCC_standard_curve.pdf", f"ERCC 标准曲线 ({s})"),
                (f"{output_dir}/results/{s}_mpc.txt", f"ERCC MPC 表 ({s})"),
            ]

    sec7 = section_open("file-index", "结果文件索引", 7)
    sec7 += '<ul class="file-index">'
    for fpath, desc in file_items:
        status = fstatus(fpath)
        rpath = fpath.replace(str(output_dir) + "/", "")
        sec7 += (
            f'<li>{status} '
            f'<span class="file-path">{rpath}</span>'
            f'<span class="file-desc">— {desc}</span>'
            f'</li>'
        )
    sec7 += "</ul>"
    sec7 += section_close()

    # ----- Assemble HTML -----
    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<style>{CSS}</style>
</head>
<body>

<div class="report-header">
  <h1>🧬 {title}</h1>
  <div class="meta">
    <span>📁 项目：{project}</span>
    <span>👤 分析人：{author}</span>
    <span>📅 生成时间：{now}</span>
    <span>🔬 样本数：{len(samples)}</span>
    <span>🆚 对照组：{contrast[0]} vs {contrast[1]}</span>
  </div>
</div>

<nav class="nav-bar">
{nav_links}
</nav>

<div class="content">
{sec1}
{sec2}
{sec3}
{sec4}
{sec5}
{sec6}
{sec7}
</div>

<div class="report-footer">
  lncRNA-seq Processing Pipeline Report &nbsp;|&nbsp; Generated by generate_report.py &nbsp;|&nbsp; {now}
</div>

</body>
</html>
"""
    return html


# ---------------------------------------------------------------------------
# Build Markdown report
# ---------------------------------------------------------------------------

def md_table(headers, rows):
    """Generate a Markdown table."""
    header_row = "| " + " | ".join(str(h) for h in headers) + " |"
    sep_row = "| " + " | ".join("---" for _ in headers) + " |"
    data_rows = "\n".join(
        "| " + " | ".join(str(c) for c in row) + " |"
        for row in rows
    )
    return f"{header_row}\n{sep_row}\n{data_rows}"


def build_markdown(config, output_dir):
    """Return the complete Markdown string."""
    report_cfg = config.get("report", {})
    title = report_cfg.get("title", "lncRNA-seq Analysis Report")
    author = report_cfg.get("author", "—")
    project = report_cfg.get("project", "—")
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    samples = config["samples"]
    contrast = config["analysis"].get("contrast", ["Treated", "Control"])
    ercc_enabled = bool(config.get("ercc", {}).get("use_ercc", False))
    deseq2_enabled = bool(config.get("analysis", {}).get("deseq2", True))
    use_rsem = bool(config.get("analysis", {}).get("use_rsem", False))
    contrast_str = f"{contrast[0]}_vs_{contrast[1]}"

    qc_rows = build_qc_table(config, output_dir)

    lines = [
        f"# {title}",
        "",
        f"- **项目 (Project):** {project}",
        f"- **分析人 (Author):** {author}",
        f"- **生成时间 (Generated):** {now}",
        f"- **样本数 (Samples):** {len(samples)}",
        f"- **比较对照 (Contrast):** {contrast[0]} vs {contrast[1]}",
        "",
        "---",
        "",
        "## 1. 项目与参数概览",
        "",
        "### 1.1 样本信息",
        "",
    ]

    # Sample table
    sample_rows = [[s, info.get("condition", "—"),
                    "Paired-end" if info.get("R2") else "Single-end"]
                   for s, info in samples.items()]
    lines.append(md_table(["样本名称", "实验组", "测序类型"], sample_rows))
    lines.append("")

    ref_cfg = config.get("ref", {})
    params_cfg = config.get("params", {})
    lines += [
        "### 1.2 关键参数",
        "",
        f"| 参数 | 值 |",
        f"| --- | --- |",
        f"| 参考基因组 (FASTA) | `{Path(ref_cfg.get('fasta', 'N/A')).name}` |",
        f"| 注释文件 (GTF) | `{Path(ref_cfg.get('gtf', 'N/A')).name}` |",
        f"| 比较对照 | {contrast[0]} vs {contrast[1]} |",
        f"| FDR 阈值 | {config['analysis'].get('fdr_threshold', 'N/A')} |",
        f"| |log\u2082FC| 阈值 | {config['analysis'].get('lfc_threshold', 'N/A')} |",
        f"| RSEM 定量 | {'✅ 启用' if use_rsem else '❌ 未启用（featureCounts）'} |",
        f"| ERCC 质控 | {'✅ 启用' if ercc_enabled else '❌ 未启用'} |",
        f"| DESeq2 分析 | {'✅ 启用' if deseq2_enabled else '❌ 未启用'} |",
        "",
        "---",
        "",
        "## 2. 数据质控",
        "",
        "### 2.1 质控指标汇总（三线表）",
        "",
    ]

    qc_headers = ["样本名称", "实验组", "原始 Reads (M)", "过滤后 Reads (M)",
                  "GC%", "Q20%", "Q30%", "比对率", "唯一比对率"]
    qc_data_rows = [
        [r["sample"], r["condition"], r["raw_reads"], r["clean_reads"],
         r["gc_pct"], r["q20_pct"], r["q30_pct"],
         r["mapping_rate"], r["unique_mapping_rate"]]
        for r in qc_rows
    ]
    lines.append(md_table(qc_headers, qc_data_rows))
    lines += ["", "> 💡 详细质控报告：`results/multiqc_report.html`", ""]

    # Section 3
    lines += [
        "---",
        "",
        "## 3. 比对与定量结果",
        "",
    ]
    multiqc_dir = Path(output_dir) / "results" / "multiqc_report_data"
    star_headers, star_rows = parse_tsv(multiqc_dir / "multiqc_star.txt")
    if star_headers and star_rows:
        display_cols = ["Sample", "total_reads", "uniquely_mapped",
                        "uniquely_mapped_percent", "multimapped_percent",
                        "unmapped_tooshort_percent"]
        idx = [i for i, h in enumerate(star_headers) if h in display_cols]
        fh = [star_headers[i] for i in idx]
        fr = [[row[i] if i < len(row) else "" for i in idx]
              for row in star_rows if row and row[0] in samples]
        if fr:
            lines.append("### 3.1 STAR 比对统计")
            lines.append("")
            lines.append(md_table(fh, fr))
            lines.append("")

    cnt_ok = (Path(output_dir) / "results" / "featureCounts.gene_counts.symbol.tsv").is_file()
    lines += [
        "### 3.2 定量结果文件",
        "",
        f"- featureCounts 计数：`results/featureCounts.gene_counts.symbol.tsv` {'✅' if cnt_ok else '❌ 未产出'}",
    ]
    if use_rsem:
        tpm_ok = (Path(output_dir) / "results" / "RSEM.gene_tpm.symbol.tsv").is_file()
        lines.append(f"- RSEM TPM 矩阵：`results/RSEM.gene_tpm.symbol.tsv` {'✅' if tpm_ok else '❌ 未产出'}")
    lines.append("")

    # Section 4: lncRNA
    lines += ["---", "", "## 4. lncRNA 相关结果汇总", ""]
    if not deseq2_enabled:
        lines += ["> ⚠ lncRNA 下游分析未启用（`analysis.deseq2: false`）。", ""]
    else:
        qs_h, qs_r = parse_csv_file(Path(output_dir) / "lncRNA_analysis_output" / "tables" / "QC_summary.csv")
        if qs_h and qs_r:
            lines.append("### 4.1 lncRNA 质控摘要")
            lines.append("")
            lines.append(md_table(qs_h, qs_r))
            lines.append("")

        lnc_h, lnc_r = parse_csv_file(Path(output_dir) / "lncRNA_analysis_output" / "tables" / "DEG_annotated_lncRNA_only.csv")
        if lnc_h and lnc_r:
            lines.append("### 4.2 lncRNA 差异分析结果（前 10 行）")
            lines.append("")
            lines.append(md_table(lnc_h, lnc_r[:10]))
            lines.append("")

        lncrna_plots = collect_lncrna_plots(output_dir)
        if lncrna_plots:
            lines.append("### 4.3 lncRNA 可视化图")
            lines.append("")
            for p, cap in lncrna_plots:
                rel = Path(p).relative_to(Path(output_dir))
                lines.append(f"![{cap}]({rel})")
            lines.append("")

    # Section 5: ERCC
    lines += ["---", "", "## 5. ERCC 定量分析", ""]
    if not ercc_enabled:
        lines += ["> ⚠ ERCC 分析未启用（`ercc.use_ercc: false`）。", ""]
    else:
        ercc_results = collect_ercc_results(config, output_dir)
        for er in ercc_results:
            lines.append(f"### 5. 样本 {er['sample']}")
            if er["curve_pdf"]:
                lines.append(f"- 标准曲线 PDF：`{er['curve_pdf']}`")
            mpc_h, mpc_r = read_mpc_table(er["mpc_table"])
            if mpc_h:
                lines.append(md_table(mpc_h, mpc_r))
            lines.append("")

    # Section 6: DESeq2
    lines += ["---", "", "## 6. DESeq2 差异表达分析", ""]
    if not deseq2_enabled:
        lines += ["> ⚠ DESeq2 差异分析未启用（`analysis.deseq2: false`）。", ""]
    else:
        deseq2_csv_path = Path(output_dir) / "results" / f"{contrast_str}.deseq2_results.csv"
        volcano_path = Path(output_dir) / "results" / f"{contrast_str}.volcano_plot.png"
        go_path = Path(output_dir) / "results" / f"{contrast_str}.GO_dotplot.png"
        lines += [
            f"比较：**{contrast[0]}** vs **{contrast[1]}**",
            "",
        ]
        if volcano_path.is_file():
            rel = volcano_path.relative_to(Path(output_dir))
            lines.append(f"![火山图]({rel})")
        if go_path.is_file():
            rel = go_path.relative_to(Path(output_dir))
            lines.append(f"![GO 富集点图]({rel})")
        lines.append("")
        deseq2_h, deseq2_r = parse_csv_file(deseq2_csv_path)
        if deseq2_h and deseq2_r:
            lines.append("### 6.1 差异分析结果（前 10 行）")
            lines.append("")
            lines.append(md_table(deseq2_h, deseq2_r[:10]))
            lines.append("")

    # Section 7: File index
    lines += [
        "---",
        "",
        "## 7. 结果文件索引",
        "",
        "| 状态 | 文件路径 | 说明 |",
        "| --- | --- | --- |",
    ]
    file_items = [
        (f"{output_dir}/results/multiqc_report.html", "MultiQC 质控汇总报告"),
        (f"{output_dir}/results/featureCounts.gene_counts.symbol.tsv", "featureCounts 基因计数矩阵"),
    ]
    if use_rsem:
        file_items.append((f"{output_dir}/results/RSEM.gene_tpm.symbol.tsv", "RSEM TPM 表达矩阵"))
    if deseq2_enabled:
        file_items += [
            (f"{output_dir}/results/{contrast_str}.deseq2_results.csv", "DESeq2 差异分析结果"),
            (f"{output_dir}/results/{contrast_str}.volcano_plot.png", "火山图"),
            (f"{output_dir}/results/{contrast_str}.GO_enrichment.csv", "GO 富集分析结果"),
            (f"{output_dir}/results/{contrast_str}.GO_dotplot.png", "GO 富集点图"),
            (f"{output_dir}/lncRNA_analysis_output/tables/QC_summary.csv", "lncRNA 质控摘要"),
            (f"{output_dir}/lncRNA_analysis_output/tables/DEG_annotated_lncRNA_only.csv", "lncRNA 差异分析注释结果"),
            (f"{output_dir}/lncRNA_analysis_output/tables/SIG_lncRNA_padj0.05_lfc1.csv", "显著差异 lncRNA 列表"),
        ]
    if ercc_enabled:
        for s in samples:
            file_items += [
                (f"{output_dir}/results/{s}_ERCC_standard_curve.pdf", f"ERCC 标准曲线 ({s})"),
                (f"{output_dir}/results/{s}_mpc.txt", f"ERCC MPC 表 ({s})"),
            ]

    for fpath, desc in file_items:
        status = "✅" if Path(fpath).is_file() else "❌"
        rpath = fpath.replace(str(output_dir) + "/", "")
        lines.append(f"| {status} | `{rpath}` | {desc} |")

    lines += ["", "---", "", f"*Report generated by lncRNA-seq Pipeline · {now}*", ""]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate lncRNA-seq pipeline HTML + Markdown report."
    )
    parser.add_argument("--config",  required=True, help="Path to config.yaml")
    parser.add_argument("--outdir",  required=True, help="Pipeline output directory")
    parser.add_argument("--report",  required=True, help="Report output directory")
    parser.add_argument("--no-inline-images", action="store_true",
                        help="Disable base64 image inlining (HTML will reference paths)")
    args = parser.parse_args()

    config_path = Path(args.config)
    output_dir = Path(args.outdir)
    report_dir = Path(args.report)
    inline_images = not args.no_inline_images

    # Load config
    if not config_path.is_file():
        print(f"ERROR: Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    config = load_yaml(str(config_path))

    # Allow config override of inline_images
    if config.get("report", {}).get("inline_images") is False:
        inline_images = False

    # Create report directory
    report_dir.mkdir(parents=True, exist_ok=True)

    # Generate HTML
    print("Generating HTML report...")
    try:
        html_content = build_html(config, str(output_dir), inline_images=inline_images)
        html_path = report_dir / "report.html"
        html_path.write_text(html_content, encoding="utf-8")
        print(f"  ✅ HTML report: {html_path}")
    except Exception as e:
        print(f"  ❌ HTML generation failed: {e}", file=sys.stderr)
        import traceback; traceback.print_exc()

    # Generate Markdown
    print("Generating Markdown report...")
    try:
        md_content = build_markdown(config, str(output_dir))
        md_path = report_dir / "report.md"
        md_path.write_text(md_content, encoding="utf-8")
        print(f"  ✅ Markdown report: {md_path}")
    except Exception as e:
        print(f"  ❌ Markdown generation failed: {e}", file=sys.stderr)
        import traceback; traceback.print_exc()

    print("Done.")


if __name__ == "__main__":
    main()
