"""Module 7/8: report — run summary, gene-set stats, completeness.

Quality-first single-page English report:
  - run summary: parameters + input evidence manifest
  - annotation overview: gene/transcript structure distributions (3 SVG charts)
  - completeness: BUSCO (proteins mode, longest isoform per gene) + OMArk
    (HOG consistency / contamination; needs an OMAmer DB via --omark-db,
    section omitted entirely when no DB is given)
  - evidence support: confidence tiers + goldset origin + intrinsic flags
  - software versions: probed live from the resolved executables
  - overall quality assessment: per-dimension scores + verdict
"""
from __future__ import annotations

import json
import logging
import re
import shutil
from collections import Counter
from pathlib import Path

from ..utils import ensure_dir, lazy, run, which
from .sanity import _parse_gff3, _read_genome

log = logging.getLogger("annomerge")

# ---------------------------------------------------------------- gene stats
def gene_set_stats(gff3: Path) -> dict:
    """Parse final annotation -> gene/transcript structure distributions."""
    genes, txs, order = _parse_gff3(Path(gff3))
    isoforms = Counter()              # transcripts per gene
    exons_per_tx = Counter()
    introns_per_gene = Counter()
    cds_len = []                      # CDS length per transcript (bp)
    tx_len = []                       # exon-span length per transcript (bp)
    n_coding = 0
    for tid in order:
        tx = txs[tid]
        isoforms[tx.gid] += 1
        exons_per_tx[len(tx.exons)] += 1
        if tx.cds:
            n_coding += 1
            cds_len.append(sum(c.end - c.start + 1 for c in tx.cds))
        if tx.exons:
            tx_len.append(sum(x.end - x.start + 1 for x in tx.exons))
    for gid in set(tx.gid for tx in txs.values()):
        members = [t for t in txs.values() if t.gid == gid]
        # introns of the gene = introns of its longest-CDS transcript
        cand = max(members, key=lambda t: sum(c.end - c.start + 1
                                              for c in t.cds) if t.cds else 0)
        intr = _tx_introns(cand)
        introns_per_gene[len(intr)] += 1
    def dist(c: Counter) -> dict:
        return {str(k): c[k] for k in sorted(c)}
    isoform_count_dist = Counter(isoforms.values())   # #txs per gene -> #genes
    return {
        "n_genes": len(set(tx.gid for tx in txs.values())) or len(genes),
        "n_transcripts": len(order),
        "n_coding_transcripts": n_coding,
        "isoforms_per_gene": dist(isoform_count_dist),
        "exons_per_transcript": dist(exons_per_tx),
        "introns_per_gene": dist(introns_per_gene),
        "cds_len_bp": _hist(cds_len),
        "transcript_len_bp": _hist(tx_len),
        "single_exon_transcripts": exons_per_tx.get(1, 0),
    }


def _tx_introns(tx) -> list:
    """Intron intervals (s,e) from exon chain of a transcript."""
    if len(tx.exons) < 2:
        return []
    ex = sorted(tx.exons, key=lambda f: f.start)
    return [(ex[i].end + 1, ex[i + 1].start - 1) for i in range(len(ex) - 1)]


def _hist(vals: list, bins: int = 40) -> dict:
    """Compact numeric summary: quartiles + histogram for charting."""
    if not vals:
        return {"n": 0}
    vs = sorted(vals)
    n = len(vs)
    def q(p):
        return vs[min(n - 1, int(p * n))]
    lo, hi = vs[0], vs[-1]
    edges = [lo + (hi - lo) * i / bins for i in range(bins + 1)]
    counts = [0] * bins
    for v in vs:
        i = min(bins - 1, int((v - lo) / (hi - lo + 1e-9) * bins))
        counts[i] += 1
    return {"n": n, "min": lo, "max": hi, "mean": round(sum(vs) / n, 1),
            "q25": q(0.25), "median": q(0.5), "q75": q(0.75),
            "hist_edges": [round(e, 1) for e in edges], "hist_counts": counts}


# ------------------------------------------------------- longest isoform faa
def longest_isoform_proteins(gff3: Path, genome_fa: Path, out_faa: Path) -> Path:
    """One protein per gene: longest CDS transcript, translated (table 1)."""
    if lazy([out_faa], desc="longest-isoform proteins"):
        return out_faa
    from .sanity import _translate_cds
    genome = _read_genome(Path(genome_fa))
    _, txs, order = _parse_gff3(Path(gff3))
    best: dict = {}
    for tid in order:
        tx = txs[tid]
        if not tx.cds:
            continue
        L = sum(c.end - c.start + 1 for c in tx.cds)
        if tx.gid not in best or L > best[tx.gid][0]:
            best[tx.gid] = (L, tx)
    with open(out_faa, "w") as fh:
        for gid, (_, tx) in sorted(best.items()):
            prot = _translate_cds(genome, tx).rstrip("*")
            fh.write(f">{gid}\n")
            for i in range(0, len(prot), 60):
                fh.write(prot[i:i + 60] + "\n")
    log.info("longest-isoform proteins: %d genes -> %s", len(best), out_faa)
    return out_faa


# ------------------------------------------------------------------ BUSCO
def busco_proteins(faa: Path, lineage: Path, out_dir: Path,
                   threads: int) -> dict:
    """BUSCO proteins-mode (offline, local lineage). Returns parsed summary."""
    bdir = ensure_dir(Path(out_dir))
    # BUSCO -o <name> writes <name>/short_summary.*.txt under the run dir
    summaries = list(bdir.glob("busco_run/short_summary*.txt"))
    if not summaries:
        run(f"cd {bdir} && PYTHONPATH= {which('busco')} -i {Path(faa).resolve()} -m proteins "
            f"-l {Path(lineage).resolve()} -o busco_run --offline -c {threads} -f",
            shell=True, log_file=bdir / "busco.log", desc="busco (report)")
        summaries = list(bdir.glob("busco_run/short_summary*.txt"))
    if not summaries:
        return {"ok": False, "error": "no BUSCO summary produced"}
    txt = summaries[0].read_text()
    out = {"ok": True, "summary_file": str(summaries[0])}
    m = re.search(r"# The lineage dataset is:\s*(\S+)", txt)
    if m:
        out["lineage"] = m.group(1)
    m = re.search(r"# BUSCO version is:\s*(\S+)", txt)
    if m:
        out["busco_version"] = m.group(1)
    # numbers come FIRST (tab-separated): "15\tComplete BUSCOs (C)"
    for key, pat in [("complete", r"Complete BUSCOs"),
                     ("single", r"Complete and single-copy BUSCOs"),
                     ("duplicated", r"Complete and duplicated BUSCOs"),
                     ("fragmented", r"Fragmented BUSCOs"),
                     ("missing", r"Missing BUSCOs"),
                     ("total", r"Total BUSCO groups searched")]:
        mm = re.search(r"\n\s*(\d+)\s*" + pat, txt)
        if mm:
            out[key] = int(mm.group(1))
    if "complete" in out and "total" in out and out["total"]:
        out["complete_pct"] = round(100.0 * out["complete"] / out["total"], 1)
        if "single" in out:
            out["single_pct"] = round(100.0 * out["single"] / out["total"], 1)
        if "duplicated" in out:
            out["duplicated_pct"] = round(100.0 * out["duplicated"] / out["total"], 1)
    return out


# ------------------------------------------------------------------ OMArk
def omark_proteins(faa: Path, db: Path, out_dir: Path, threads: int) -> dict:
    """OMArk quality assessment (needs an OMAmer database, --omark-db).

    Degrades gracefully: {"ok": False, "error": ...} when the DB or the
    omark/omamer installation is missing — the report then shows a skip note.
    """
    from .omark import run_omark
    return run_omark(Path(faa), Path(db), Path(out_dir), threads)


# ------------------------------------------------------------- evidence
def evidence_distribution(confidence_tsv: Path) -> dict:
    """confidence.tsv -> tier counts + goldset origin + intrinsic flags."""
    if not confidence_tsv or not Path(confidence_tsv).exists():
        return {}
    tiers, gold, flags = Counter(), 0, Counter()
    with open(confidence_tsv) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 18 or f[0] == "gene_id":
                continue
            tiers[f[6]] += 1
            if f[7] == "yes":
                gold += 1
            for fl in (f[17] or "").split():
                if fl and fl != ".":
                    flags[fl] += 1
    return {"tiers": dict(tiers), "goldset_origin": gold,
            "intrinsic_flags": dict(flags)}


# ---------------------------------------------------------------- versions
def _probe_version(cmd: list) -> str | None:
    """Run --version/-version/-v; parse a version token from stdout+stderr."""
    import subprocess
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        text = p.stdout + "\n" + p.stderr
    except Exception:
        return None
    for line in text.splitlines()[:8]:
        line = line.strip()
        if not line:
            continue
        m = re.search(r"(\d+\.\d+(?:\.\d+)?(?:[-_][0-9A-Za-z.]+)?)", line)
        if m:
            return m.group(1)
    return None


def _pkg_version(dist: str) -> str | None:
    try:
        from importlib.metadata import version
        return version(dist)
    except Exception:
        return None


def software_versions() -> dict:
    """Versions of the tools the pipeline invoked, probed from the resolved
    executables (fixed: gffcompare/hisat2 print to stderr, EVM/TransDecoder
    have no --version flag, BUSCO --version cannot run on read-only mounts)."""
    from ..deps import home
    rows: list[tuple[str, str]] = []
    from .. import __version__
    rows.append(("AnnoMerge", __version__))
    rows.append(("Python", f"{__import__('sys').version_info.major}."
                           f"{__import__('sys').version_info.minor}."
                           f"{__import__('sys').version_info.micro}"))
    for label, exe in [("HISAT2", "hisat2"), ("StringTie", "stringtie"),
                       ("SAMtools", "samtools"), ("miniprot", "miniprot"),
                       ("gffcompare", "gffcompare"), ("gffread", "gffread"),
                       ("DIAMOND", "diamond")]:
        try:
            v = _probe_version([which(exe), "--version"])
        except Exception:
            v = None
        if not v:
            try:
                v = _probe_version([which(exe), "-version"])
            except Exception:
                v = None
        if not v:
            try:
                v = _probe_version([which(exe), "-v"])
            except Exception:
                v = None
        rows.append((label, v or "not found"))
    # TransDecoder: no --version; latest "## vX.Y.Z" from the bundled Changelog
    td_v = None
    try:
        td_dir = home("transdecoder")
        changelog = td_dir / "Changelog.txt"
        if changelog.exists():
            m = re.search(r"^##\s*v?(\d+\.\d+(?:\.\d+)?)",
                          changelog.read_text(errors="ignore"), re.M)
            if m:
                td_v = m.group(1)
        if not td_v:
            td_v = f"{td_dir.name} (install dir)"
    except Exception:
        td_v = "not found"
    rows.append(("TransDecoder", td_v))
    # EVidenceModeler ships no --version. Version sources, in priority order:
    #   1. the resolved install directory name (conda: opt/evidencemodeler-2.1.0)
    #   2. the top-level dir of the distribution tarball (EVidenceModeler-X.Y.Z)
    evm_v = None
    try:
        evm_home = home("evm")
        m = re.search(r"evidencemodeler[-_](\d+\.\d+(?:\.\d+)?)",
                      str(evm_home), re.I)
        if m:
            evm_v = m.group(1)
    except Exception:
        evm_v = None
    if not evm_v:
        try:
            import tarfile
            evm_home = home("evm")
            for tg in sorted(evm_home.parent.glob("evm*.tar.gz")):
                try:
                    with tarfile.open(tg) as tf:
                        names = tf.getnames()
                except Exception:
                    continue
                m = re.search(r"EVidenceModeler-(\d+\.\d+(?:\.\d+)?)",
                              "\n".join(names))
                if m:
                    evm_v = m.group(1)
                    break
        except Exception:
            evm_v = None
    rows.append(("EVidenceModeler", evm_v or "unknown"))
    # BUSCO: `busco --version` writes logs and dies on read-only mounts;
    # derive the version from the resolved installation path instead.
    busco_v = None
    try:
        m = re.search(r"busco[/\\](\d+\.\d+(?:\.\d+)?)", which("busco"))
        if m:
            busco_v = m.group(1)
    except Exception:
        pass
    rows.append(("BUSCO", busco_v or _pkg_version("busco") or "unknown"))
    rows.append(("OMArk", _pkg_version("omark") or "not installed"))
    rows.append(("OMAmer", _pkg_version("omamer") or "not installed"))
    return {"rows": rows}


# ------------------------------------------------------- quality assessment
def quality_assessment(ctx: dict) -> dict:
    """Per-dimension scores + an overall verdict (quality-first summary)."""
    dims: list[dict] = []
    busco = ctx.get("busco") or {}
    if busco.get("ok"):
        dims.append({"dimension": "Completeness (BUSCO complete)",
                     "score": busco.get("complete_pct"),
                     "detail": f"{busco.get('complete', '?')}/{busco.get('total', '?')} "
                               f"BUSCOs complete (S {busco.get('single_pct')}%, "
                               f"D {busco.get('duplicated_pct')}%)"})
    omark = ctx.get("omark") or {}
    if omark.get("ok"):
        dims.append({"dimension": "Consistency (OMArk)",
                     "score": omark.get("consistent_pct"),
                     "detail": f"{omark.get('consistent', '?')}/{omark.get('total', '?')} "
                               f"genes consistent with HOGs; contamination "
                               f"{omark.get('contamination_pct', '?')}%"})
    # NOTE: the overall verdict averages exactly the three user-facing
    # dimensions (evidence support / BUSCO / OMArk)
    ev = ctx.get("evidence") or {}
    gs = ctx.get("gene_set") or {}
    n_tx = gs.get("n_transcripts") or 0
    strong = (ev.get("tiers") or {}).get("strong", 0)
    if n_tx and strong:
        dims.append({"dimension": "Evidence-backed (strong-tier transcripts)",
                     "score": round(100.0 * strong / n_tx, 1),
                     "detail": f"{strong}/{n_tx} transcripts strong"})
    # verdict averages exactly the three headline dims collected above
    scores = [d["score"] for d in dims if isinstance(d.get("score"), (int, float))]
    overall = round(sum(scores) / len(scores), 1) if scores else None
    if overall is None:
        label = "n/a"
    elif overall >= 95:
        label = "Excellent"
    elif overall >= 90:
        label = "Good"
    elif overall >= 75:
        label = "Fair"
    else:
        label = "Poor"
    return {"dims": dims, "overall": overall, "label": label}


# --------------------------------------------------------------------- html
_TEMPLATE = (Path(__file__).resolve().parent.parent / "assets"
             / "report_template.html").read_text()
_WORKFLOW_SVG = (Path(__file__).resolve().parent.parent / "assets"
                 / "workflow_diagram.svg")
def render_html(ctx: dict, out_html: Path,
                meta: dict | None = None) -> Path:
    """Render the single-page report."""
    from jinja2 import Environment, BaseLoader
    from markupsafe import Markup
    env = Environment(loader=BaseLoader(), autoescape=True)
    tpl = env.from_string(_TEMPLATE)
    diagram = Markup(_WORKFLOW_SVG.read_text()) if _WORKFLOW_SVG.exists() else None
    html = tpl.render(meta=meta or ctx.get("meta", {}),
                      workflow_diagram=diagram,
                      gene_set=ctx.get("gene_set", {}),
                      busco=ctx.get("busco", {}),
                      omark=ctx.get("omark", {}),
                      evidence=ctx.get("evidence", {}),
                      versions=ctx.get("versions", {}),
                      assessment=ctx.get("assessment", {}),
                      integ=ctx.get("integ"),
                      workflow=ctx.get("workflow", {}),
                      tags=ctx.get("tags", {}),
                      final_score=ctx.get("final_score"),
                      orthodb=ctx.get("orthodb"))
    Path(out_html).write_text(html)
    log.info("report html -> %s", out_html)
    return Path(out_html)


# --------------------------------------------------------- inline SVG charts
def _esc(s: str) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def svg_bars(rows: list, maxv: float = 100.0, width: int = 520,
             bar_h: int = 18, color: str = "#4C72B0") -> str:
    """Horizontal percentage bar chart as inline SVG (no JS, offline-safe)."""
    if not rows:
        return ""
    label_w = max(80, min(220, max(len(str(r[0])) for r in rows) * 7 + 12))
    h = (bar_h + 6) * len(rows) + 10
    plot_w = width - label_w - 70
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{h}" '
             f'font-family="Helvetica,Arial,sans-serif" font-size="12">']
    y = 5
    for label, v in rows:
        w = max(1.0, plot_w * min(float(v), maxv) / maxv)
        parts.append(f'<text x="{label_w - 6}" y="{y + bar_h - 5}" '
                     f'text-anchor="end">{_esc(label)}</text>')
        parts.append(f'<rect x="{label_w}" y="{y}" width="{w:.1f}" '
                     f'height="{bar_h}" fill="{color}" rx="2"/>')
        parts.append(f'<text x="{label_w + w + 6:.1f}" y="{y + bar_h - 5}">'
                     f'{float(v):.1f}%</text>')
        y += bar_h + 6
    parts.append("</svg>")
    return "".join(parts)


def svg_histogram(vals: list, bins: int = 20, width: int = 520,
                  height: int = 170, color: str = "#4C72B0",
                  title: str = "") -> str:
    """Vertical histogram of values in [0,1] as inline SVG."""
    vals = [float(v) for v in vals if v is not None]
    if not vals:
        return ""
    lo, hi = min(vals), max(vals)
    span = (hi - lo) or 1e-9
    counts = [0] * bins
    for v in vals:
        i = min(bins - 1, int((v - lo) / span * bins))
        counts[i] += 1
    top = max(counts) or 1
    pad_l, pad_b, pad_t = 44, 30, 22
    plot_w, plot_h = width - pad_l - 12, height - pad_b - pad_t
    bw = plot_w / bins
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
             f'font-family="Helvetica,Arial,sans-serif" font-size="11">']
    if title:
        parts.append(f'<text x="{pad_l}" y="14" font-size="12" font-weight="bold">'
                     f'{_esc(title)}</text>')
    for i, c in enumerate(counts):
        bh = plot_h * c / top
        x = pad_l + i * bw
        parts.append(f'<rect x="{x:.1f}" y="{pad_t + plot_h - bh:.1f}" '
                     f'width="{max(bw - 1.5, 1):.1f}" height="{bh:.1f}" '
                     f'fill="{color}"/>')
    # axes + x range labels
    parts.append(f'<line x1="{pad_l}" y1="{pad_t + plot_h}" x2="{width - 12}" '
                 f'y2="{pad_t + plot_h}" stroke="#888"/>')
    parts.append(f'<text x="{pad_l}" y="{height - 10}">{lo:.2f}</text>')
    parts.append(f'<text x="{width - 12}" y="{height - 10}" text-anchor="end">{hi:.2f}</text>')
    parts.append(f'<text x="{pad_l - 6}" y="{pad_t + 10}" text-anchor="end">{top}</text>')
    parts.append(f'<text x="{pad_l - 6}" y="{pad_t + plot_h}" text-anchor="end">0</text>')
    parts.append("</svg>")
    return "".join(parts)


# ------------------------------------------------- workflow evidence figures
def _gauss_smooth_2d(grid: list, sigma: float = 2.0) -> list:
    """Separable Gaussian smoothing of a 2D count grid (rows x cols)."""
    import math
    radius = max(2, int(3 * sigma))
    kern = [math.exp(-0.5 * (k / sigma) ** 2) for k in range(-radius, radius + 1)]
    s = sum(kern)
    kern = [k / s for k in kern]
    nr, nc = len(grid), len(grid[0])

    def blur(src, horizontal: bool) -> list:
        out = [[0.0] * nc for _ in range(nr)]
        for r in range(nr):
            for c in range(nc):
                acc = 0.0
                for k in range(-radius, radius + 1):
                    if horizontal:
                        cc = c + k
                        v = src[r][cc] if 0 <= cc < nc else 0
                    else:
                        rr = r + k
                        v = src[rr][c] if 0 <= rr < nr else 0
                    acc += kern[k + radius] * v
                out[r][c] = acc
        return out

    return blur(blur(grid, True), False)


def svg_density_marginal(xs: list, ys: list, x_label: str, y_label: str,
                         x_log: bool = False, xbins: int = 150, ybins: int = 90,
                         width: int = 620, height: int = 470,
                         title: str = "") -> str:
    """NanoPlot-style SMOOTH density scatter with marginal KDE curves.

    Central panel: Gaussian-smoothed 2D density field with NanoPlot's YlGnBu
    gradient (cream -> green -> teal -> dark blue). Marginals: filled smooth
    density curves in NanoPlot teal-green. Pure inline SVG, no JS."""
    import math
    pts = [(float(x), float(y)) for x, y in zip(xs, ys)
           if x is not None and y is not None and float(x) > 0]
    if not pts:
        return ""
    xs_v = [p[0] for p in pts]
    ys_v = [p[1] for p in pts]
    if x_log:
        xs_v = [math.log10(v) for v in xs_v]
    x_lo, x_hi = min(xs_v), max(xs_v)
    y_lo, y_hi = min(ys_v), max(ys_v)
    x_sp = (x_hi - x_lo) or 1e-9
    y_sp = (y_hi - y_lo) or 1e-9
    grid = [[0] * xbins for _ in range(ybins)]
    for x, y in pts:
        xx = math.log10(x) if x_log else x
        i = min(xbins - 1, int((xx - x_lo) / x_sp * xbins))
        j = min(ybins - 1, int((y - y_lo) / y_sp * ybins))
        grid[j][i] += 1
    sm = _gauss_smooth_2d(grid, sigma=2.2)
    vmax = max(max(r) for r in sm) or 1.0

    # geometry: top marginal sits ABOVE the main panel; both start below the
    # optional title. Fixes the previous overflow where curves were drawn
    # upward from the main panel's top edge and escaped the canvas.
    pad_l = 64
    title_h = 26 if title else 0
    top_h, right_w, gap = 95, 115, 8
    main_t = title_h + top_h + gap            # main panel top edge
    main_w = width - pad_l - right_w - gap - 14
    main_h = height - main_t - 30
    cw, ch = main_w / xbins, main_h / ybins
    NANO_GREEN = "#4cb391"                    # NanoPlot signature teal-green

    # NanoPlot YlGnBu-style gradient stops (low -> high density)
    _STOPS = [(0.00, (255, 255, 204)), (0.25, (161, 218, 180)),
              (0.50, (65, 182, 196)), (0.75, (44, 127, 184)),
              (1.00, (37, 52, 148))]

    def dens_color(v: float) -> str:
        if v <= 1e-4 * vmax:
            return "#ffffff"
        t = math.log10(1 + 9 * v / vmax)          # 0..1, log compressed
        for (t0, c0), (t1, c1) in zip(_STOPS, _STOPS[1:]):
            if t <= t1:
                f = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
                r = int(c0[0] + (c1[0] - c0[0]) * f)
                g = int(c0[1] + (c1[1] - c0[1]) * f)
                b = int(c0[2] + (c1[2] - c0[2]) * f)
                return f"#{r:02x}{g:02x}{b:02x}"
        r, g, b = _STOPS[-1][1]
        return f"#{r:02x}{g:02x}{b:02x}"

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
             f'height="{height}" font-family="Helvetica,Arial,sans-serif" font-size="11">']
    if title:
        parts.append(f'<text x="{pad_l}" y="17" font-size="13" font-weight="bold">'
                     f'{_esc(title)}</text>')
    # main smoothed density field (y axis: high values on top)
    for j in range(ybins):
        row = sm[j]
        y0 = main_t + (ybins - 1 - j) * ch
        for i in range(xbins):
            v = row[i]
            if v <= 1e-4 * vmax:
                continue
            parts.append(f'<rect x="{pad_l + i * cw:.2f}" y="{y0:.2f}" '
                         f'width="{cw + 0.05:.2f}" height="{ch + 0.05:.2f}" '
                         f'fill="{dens_color(v)}"/>')
    parts.append(f'<rect x="{pad_l}" y="{main_t}" width="{main_w:.1f}" '
                 f'height="{main_h:.1f}" fill="none" stroke="#94a3b8"/>')

    def smooth_1d(counts: list, sigma: float = 2.0) -> list:
        return _gauss_smooth_2d([counts], sigma=sigma)[0]

    def filled_curve(vals: list, orientation: str, base: float, span: float,
                     amp: float, color: str) -> str:
        n = len(vals)
        mx = max(vals) or 1.0
        pts_s = []
        for i, v in enumerate(vals):
            frac = i / (n - 1)
            if orientation == "top":
                px = pad_l + frac * span
                py = base - amp * v / mx
            else:  # right
                px = base + amp * v / mx
                py = main_t + main_h - frac * span
            pts_s.append(f"{px:.1f},{py:.1f}")
        if orientation == "top":
            d = (f"M {pad_l:.1f},{base:.1f} L " + " L ".join(pts_s) +
                 f" L {pad_l + span:.1f},{base:.1f} Z")
        else:
            d = (f"M {base:.1f},{main_t + main_h:.1f} L " + " L ".join(pts_s) +
                 f" L {base:.1f},{main_t:.1f} Z")
        return (f'<path d="{d}" fill="{color}" fill-opacity="0.55" '
                f'stroke="{color}" stroke-width="1.2"/>')

    # top marginal: smoothed x-density curve; baseline = main panel top edge,
    # curve grows upward into its own band -> never leaves the canvas
    xcnt = [sum(grid[j][i] for j in range(ybins)) for i in range(xbins)]
    xsm = smooth_1d(xcnt, sigma=2.5)
    parts.append(filled_curve(xsm, "top", main_t, main_w, top_h - 4, NANO_GREEN))
    # right marginal: smoothed y-density curve
    ycnt = [sum(row) for row in grid]
    ysm = smooth_1d(ycnt, sigma=2.5)
    rx = pad_l + main_w + gap
    parts.append(filled_curve(ysm, "right", rx, main_h, right_w - 4, NANO_GREEN))

    # axes labels
    def tick(v, log):
        if log:
            p = int(round(v))
            return f"10^{p}" if p < 4 else f"1e{p}"
        return f"{v:.2f}"
    y_base = main_t + main_h
    for frac in (0.0, 0.5, 1.0):
        vx = x_lo + frac * x_sp
        parts.append(f'<text x="{pad_l + frac * main_w:.1f}" y="{y_base + 14}" '
                     f'text-anchor="middle">{tick(vx, x_log)}</text>')
    parts.append(f'<text x="{pad_l + main_w / 2:.1f}" y="{height - 4}" '
                 f'text-anchor="middle" font-size="12">{_esc(x_label)}</text>')
    for frac in (0.0, 0.5, 1.0):
        vy = y_lo + frac * y_sp
        parts.append(f'<text x="{pad_l - 5}" y="{y_base - frac * main_h + 4:.1f}" '
                     f'text-anchor="end">{vy:.2f}</text>')
    parts.append(f'<text x="14" y="{main_t + main_h / 2:.1f}" text-anchor="middle" '
                 f'font-size="12" transform="rotate(-90 14 {main_t + main_h / 2:.1f})">'
                 f'{_esc(y_label)}</text>')
    parts.append(f'<text x="{pad_l + main_w:.1f}" y="{title_h + 12:.1f}" '
                 f'text-anchor="end" font-size="10" fill="#64748b">'
                 f'n = {len(pts)}</text>')
    parts.append("</svg>")
    return "".join(parts)


# ------------------------------------------------- workflow evidence figures
def rna_alignment_stats(out: Path) -> list:
    """Per-sample HISAT2 overall alignment rate from 02.rnaseq/bam/*.hisat2.log."""
    rows = []
    for logf in sorted(Path(out).glob("02.rnaseq/bam/*.hisat2.log")):
        m = re.search(r"([\d.]+)% overall alignment rate",
                      logf.read_text(errors="ignore"))
        if m:
            rows.append({"sample": logf.name.replace(".hisat2.log", ""),
                         "rate": float(m.group(1))})
    return rows


def protein_alignment_stats(out: Path) -> dict:
    """Homolog protein alignment: coverage/identity distributions + mnp filter."""
    res: dict = {"qcov": [], "ident": [], "qlen": [], "n_mapped": 0, "n_raw": 0,
                 "mnp_rows": [], "mnp_removed": 0}
    tsv = Path(out) / "03.homolog" / "alignments.tsv"
    if tsv.exists():
        seen: set = set()
        with open(tsv) as fh:
            header = next(fh).rstrip("\n").split("\t")
            ci = {n: i for i, n in enumerate(header)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < len(header):
                    continue
                q = f[ci["query"]]
                if q in seen:          # one record per protein for the distribution
                    continue
                seen.add(q)
                res["qcov"].append(float(f[ci["qcov"]]))
                res["ident"].append(float(f[ci["ident"]]))
                res["qlen"].append(int(f[ci["qlen"]]))
        res["n_mapped"] = len(seen)
    raw_paf = Path(out) / "03.homolog" / "miniprot.raw.paf"
    if raw_paf.exists():
        queries: set = set()
        with open(raw_paf) as fh:
            for line in fh:
                queries.add(line.split("\t", 1)[0])
        res["n_raw"] = len(queries)
    mnp = Path(out) / "03.homolog" / "mnp_filter.tsv"
    if mnp.exists():
        with open(mnp) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == 2 and f[0] != "reason" and f[1].isdigit():
                    res["mnp_rows"].append((f[0], int(f[1])))
                    res["mnp_removed"] += int(f[1])
    return res


def goldset_stats(out: Path) -> dict:
    """Gold set: StringTie transcripts in, TransDecoder ORFs out, protein length."""
    gdir = Path(out) / "04.goldset"
    res: dict = {}
    t_fa, g_fa = gdir / "transcripts.fa", gdir / "goldset.faa"
    if t_fa.exists():
        res["n_transcripts"] = sum(1 for l in open(t_fa) if l.startswith(">"))
    if g_fa.exists():
        lens, cur = [], 0
        for l in open(g_fa):
            if l.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
            else:
                cur += len(l.strip())
        if cur:
            lens.append(cur)
        res["n_proteins"] = len(lens)
        if lens:
            lens.sort()
            res["median_protein_aa"] = lens[len(lens) // 2]
            res["mean_protein_aa"] = round(sum(lens) / len(lens), 1)
    return res


def _read_kv(tsv: Path) -> dict:
    out: dict = {}
    if tsv.exists():
        for line in open(tsv):
            f = line.rstrip("\n").split("\t")
            if len(f) == 2 and f[0] != "key":
                out[f[0]] = f[1]
    return out


# ------------------------------------------------- final-annotation scoring
def final_annotation_score(out: Path, final_gff3: Path, genome_fa: Path,
                           threads: int, ctx: dict,
                           rna_junctions=None, stringtie_gtf=None,
                           homologs_faa=None, referee_ids=None) -> dict | None:
    """Score the FINAL annotation with the same P formula as the predictors
    (fixed default weights). B_c / O_s reuse the BUSCO / OMArk results already
    computed for this report; R_s / H_s / C_s are computed on the final GFF3
    and cached under 08.assessment/final_scoring."""
    from .scoring import (compute_rs, compute_hs, compute_cs, aggregate,
                          extract_referee_faa, load_referee_loci,
                          extract_proteins)
    work = ensure_dir(Path(out) / "08.assessment" / "final_scoring")
    sub: dict = {}
    try:
        pred_faa = extract_proteins(Path(genome_fa), Path(final_gff3),
                                    work / "proteins.faa")
        if rna_junctions and stringtie_gtf and Path(str(rna_junctions)).exists():
            r = compute_rs(Path(final_gff3), Path(str(rna_junctions)),
                           Path(str(stringtie_gtf)), work, 0.5, 0.35, 0.15)
            sub["R_s"] = r["R_s"]
        else:
            sub["R_s"] = None
        ref_faa = None
        if homologs_faa and referee_ids and Path(str(referee_ids)).exists():
            ref_faa = extract_referee_faa(Path(str(homologs_faa)),
                                          Path(str(referee_ids)),
                                          work / "referee.faa")
        if ref_faa:
            h = compute_hs(pred_faa, ref_faa, None, work / "hs", threads, 0.9)
            sub["H_s"] = h["Hs"]
            alignments = Path(out) / "03.homolog" / "alignments.tsv"
            if alignments.exists() and h.get("_hits_cache"):
                loci = load_referee_loci(alignments)
                c = compute_cs(pred_faa, h["_hits_cache"], loci)
                sub["C_s"] = c["C_s"]
            else:
                sub["C_s"] = None
        else:
            sub["H_s"] = sub["C_s"] = None
        busco = ctx.get("busco") or {}
        sub["B_c"] = (busco["complete_pct"] / 100.0) \
            if busco.get("ok") and busco.get("complete_pct") is not None else None
        omark = ctx.get("omark") or {}
        sub["O_s"] = (float(omark["consistent_pct"]) / 100.0) \
            if omark.get("ok") and omark.get("consistent_pct") is not None else None
        w, eff = aggregate(sub, {"R_s": 0.45, "H_s": 0.25, "B_c": 0.15,
                                 "O_s": 0.10, "C_s": 0.15})
        log.info("final annotation P=%s subscores=%s",
                 f"{w:.4f}" if w is not None else "NA",
                 {k: (f"{v:.3f}" if isinstance(v, float) else "NA")
                  for k, v in sub.items()})
        return {"P": w, "subscores": sub,
                "effective_weights": {k: round(v, 3) for k, v in eff.items() if v > 0}}
    except Exception as e:
        log.warning("final annotation scoring failed: %s", e)
        return None


# --------------------------------------------------------------- tag table
def tag_overview(final_gff3: Path) -> dict:
    """Count every support tag written onto the final GFF3 (transcript-level
    for Tier/RNA/protein/denovo/flags, gene-level for BUSCO/OMArk)."""
    tiers, flags, denovo = Counter(), Counter(), Counter()
    rna_yes = prot_yes = busco_sup = omark_sup = n_tx = 0
    with open(final_gff3) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = f[8]
            if f[2] == "mRNA":
                n_tx += 1
                m = re.search(r"Tier=([^;]+)", attrs)
                if m:
                    tiers[m.group(1)] += 1
                m = re.search(r"RNA_support=([\d.]+)", attrs)
                if m and float(m.group(1)) > 0:
                    rna_yes += 1
                m = re.search(r"Protein_support=([\d.]+)", attrs)
                if m and float(m.group(1)) > 0:
                    prot_yes += 1
                m = re.search(r"denovo_support=([^;]+)", attrs)
                if m:
                    for kv in m.group(1).split(","):
                        k, _, v = kv.partition(":")
                        try:
                            if float(v) > 0:
                                denovo[k] += 1
                        except ValueError:
                            pass
                m = re.search(r"Flags=([^;]+)", attrs)
                if m:
                    for fl in m.group(1).split(","):
                        fl = fl.strip()
                        if fl and fl != ".":
                            flags[fl] += 1
            elif f[2] == "gene":
                busco_sup += "BUSCO_support=" in attrs
                omark_sup += "OMArk_support=" in attrs
    return {"n_transcripts": n_tx, "tiers": dict(tiers),
            "rna_supported": rna_yes, "protein_supported": prot_yes,
            "denovo_supported": dict(denovo), "intrinsic_flags": dict(flags),
            "busco_supported_genes": busco_sup, "omark_supported_genes": omark_sup}


# ------------------------------------------------------------- orchestrate
def _summary_meta(run_meta: dict | None) -> dict:
    """Run-summary display fields: version, timing, full command (one
    parameter per line)."""
    from .. import __version__
    import shlex as _shlex
    m = dict(run_meta or {})
    m.setdefault("version", __version__)
    # full command, one parameter per line, for the quoted command block
    cmd = str(m.get("command") or "").strip()
    if cmd:
        try:
            args = _shlex.split(cmd)
        except ValueError:
            args = cmd.split()
        m["cmd_head"] = args[0] if args else ""
        m["cmd_args"] = args[1:] if len(args) > 1 else []
    else:
        m["cmd_head"], m["cmd_args"] = "", []
    return m


def _integration_stats(out: Path, final_gff3: Path) -> dict | None:
    """Integration-pipeline figures for report sections 2.1/2.2:
    predictor weights (scoring/evm_weights.tsv), post-EVM gene/transcript
    counts, and how many genes/transcripts augment recovered on top of EVM.
    None when this is not an integration run."""
    wtab = Path(out) / "05.score" / "evm_weights.tsv"
    evm_gff = Path(out) / "06.merge" / "evm.sanity.gff3"
    if not wtab.exists() or not evm_gff.exists():
        return None
    pred_score: dict = {}
    with open(wtab) as fh:
        for line in fh:
            if line.startswith("predictor"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 7 and f[0]:
                pred_score[f[0]] = f[1]
    # full EVM evidence manifest (evm/weights.txt): evidence type, source,
    # calibrated score (de novo predictors only) and the EVM weight used
    kind_label = {"ABINITIO_PREDICTION": "Ab initio prediction",
                  "TRANSCRIPT": "Transcript assembly",
                  "PROTEIN": "Homologous protein"}
    ev_rows = []
    wt = Path(out) / "06.merge" / "weights.txt"
    if wt.exists():
        with open(wt) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and f[0]:
                    ev_rows.append({"kind": kind_label.get(f[0], f[0]),
                                    "evidence": f[1],
                                    "score": pred_score.get(f[1], "n/a"),
                                    "evm_weight": f[2]})
    try:
        evm_gs = gene_set_stats(evm_gff)
        fin_gs = gene_set_stats(final_gff3)
    except Exception:
        return None
    # per-predictor scoring transparency: effective sub-block weights (after
    # NA renormalization / adaptive coverage scaling) and sub-blocks that
    # could not be scored — so users see exactly what each P was built on
    scoring_notes = []
    score_root = Path(out) / "05.score"
    for name in pred_score:
        dj = score_root / name / "score_detail.json"
        if not dj.exists():
            continue
        try:
            d = json.loads(dj.read_text())
        except Exception:
            continue
        eff = {k: round(v, 3) for k, v in
               (d.get("effective_weights") or {}).items() if v > 0}
        missing = {k: v for k, v in d.items()
                   if k.endswith("_note") and d.get("subscores", {}).get(k[:-5]) is None}
        cov = {k[:-12]: round(v, 3) for k, v in
               (d.get("informative_coverage") or {}).items() if v < 1.0}
        scoring_notes.append({"predictor": name, "effective_weights": eff,
                              "missing": missing, "coverage": cov})
    return {
        "evidence_rows": ev_rows,
        "scoring_notes": scoring_notes,
        "evm_genes": evm_gs["n_genes"],
        "evm_transcripts": evm_gs["n_transcripts"],
        "added_genes": max(0, fin_gs["n_genes"] - evm_gs["n_genes"]),
        "added_transcripts": max(0, fin_gs["n_transcripts"] - evm_gs["n_transcripts"]),
    }


def _elapsed_str(start_iso: str, end_dt) -> str:
    """Human-readable run duration from an ISO start timestamp."""
    try:
        from datetime import datetime
        start = datetime.fromisoformat(str(start_iso))
        secs = max(0, int((end_dt - start).total_seconds()))
        h, rem = divmod(secs, 3600)
        mnt, sec = divmod(rem, 60)
        if h:
            return f"{h}h {mnt}m {sec}s"
        if mnt:
            return f"{mnt}m {sec}s"
        return f"{sec}s"
    except Exception:
        return "n/a"


def final_sequences(final_gff3: Path, genome_fa: Path, out_dir: Path) -> dict:
    """09.result deliverables: annomerge.gff (full annotation),
    annomerge.pep / annomerge.cds (one representative isoform per gene:
    longest CDS)."""
    import shutil
    res = {}
    gff_out = Path(out_dir) / "annomerge.gff"
    shutil.copy2(final_gff3, gff_out)
    res["gff"] = gff_out
    # representative isoform per gene = longest CDS
    _, txs, order = _parse_gff3(Path(final_gff3))
    best: dict = {}
    for tid in order:
        tx = txs[tid]
        if not tx.cds:
            continue
        L = sum(c.end - c.start + 1 for c in tx.cds)
        if tx.gid not in best or L > best[tx.gid][0]:
            best[tx.gid] = (L, tid)
    keep = {tid for _, tid in best.values()}
    gffread = which("gffread")
    tmp_pep = Path(out_dir) / ".annomerge.all.pep.fa"
    tmp_cds = Path(out_dir) / ".annomerge.all.cds.fa"
    run([gffread, "-y", str(tmp_pep), "-x", str(tmp_cds),
         "-g", str(genome_fa), str(final_gff3)],
        desc="gffread pep/cds extraction")
    for tmp, name in ((tmp_pep, "annomerge.pep"), (tmp_cds, "annomerge.cds")):
        out_fa = Path(out_dir) / name
        n = 0
        with open(tmp) as fh, open(out_fa, "w") as out:
            take = False
            for ln in fh:
                if ln.startswith(">"):
                    rid = ln[1:].split()[0]
                    take = rid in keep
                    if take:
                        n += 1
                        out.write(f">{rid}\n")
                elif take:
                    out.write(ln)
        res[name.split(".", 1)[1]] = out_fa
        log.info("%s: %d representative sequences -> %s", name, n, out_fa)
    tmp_pep.unlink(missing_ok=True)
    tmp_cds.unlink(missing_ok=True)
    return res


def run_report(out: Path, final_gff3: Path, genome_fa: Path,
               final_d: Path, busco_lineage: Path | None,
               threads: int,
               run_meta: dict | None = None,
               omark_db: Path | None = None,
               rna_junctions=None, stringtie_gtf=None,
               homologs_faa=None, referee_ids=None) -> dict:
    """Collect every report input; writes 09.result/report.html + stats.json
    and 08.assessment/BUSCO (+OMArk). Returns context."""
    d = ensure_dir(Path(out) / "09.result")
    assess = ensure_dir(Path(out) / "08.assessment")
    ctx: dict = {"meta": _summary_meta(run_meta)}
    log.info("module 7/8 [report]: collecting stats")
    ctx["gene_set"] = gene_set_stats(final_gff3)
    ctx["evidence"] = evidence_distribution(assess / "confidence.tsv")
    faa = longest_isoform_proteins(final_gff3, genome_fa,
                                   d / "longest_isoforms.faa")
    if busco_lineage and Path(busco_lineage).exists():
        ctx["busco"] = busco_proteins(faa, busco_lineage, assess / "busco",
                                      threads)
    else:
        ctx["busco"] = {"ok": False, "error": "no BUSCO lineage given"}
    if omark_db and Path(omark_db).exists():
        ctx["omark"] = omark_proteins(faa, omark_db, assess / "omark", threads)
    else:
        # no usable OMAmer database (not given or file missing): leave the
        # OMArk section out entirely — no placeholder, no skip note
        if omark_db:
            log.info("omark db %s not found; OMArk section omitted", omark_db)
        ctx["omark"] = None
    ctx["versions"] = software_versions()
    ctx["assessment"] = quality_assessment(ctx)
    # integration-only sections 2.4-2.6 (None for evaluate runs)
    ctx["integ"] = _integration_stats(out, final_gff3)
    # workflow evidence figures (section 2)
    wf = {
        "mode": (run_meta or {}).get("mode"),
        "rna_align": rna_alignment_stats(out),
        "protein": protein_alignment_stats(out),
        "goldset": goldset_stats(out),
        "augment": _read_kv(Path(out) / "07.augment" / "augment_stats.tsv"),
    }
    from markupsafe import Markup
    wf["rna_chart"] = Markup(svg_bars(
        [(r["sample"], r["rate"]) for r in wf["rna_align"]]))
    wf["qcov_hist"] = Markup(svg_histogram(
        wf["protein"].get("qcov") or [], title="Query coverage (fraction of protein aligned)"))
    wf["ident_hist"] = Markup(svg_histogram(
        wf["protein"].get("ident") or [], title="Alignment identity"))
    wf["density_plot"] = Markup(svg_density_marginal(
        wf["protein"].get("qcov") or [], wf["protein"].get("ident") or [],
        x_label="Query coverage", y_label="Alignment identity"))
    ctx["workflow"] = wf
    # support-tag census (section 3.4)
    ctx["tags"] = tag_overview(final_gff3)
    # orthology-database support tags (gene level): runs only when a BUSCO
    # lineage and/or an OMAmer DB was provided AND produced results; writes
    # BUSCO_support/BUSCO_group/OMArk_support/OMArk_HOG onto gene rows of the
    # final GFF3 so downstream GFF3 consumers (incl. the evidence browser)
    # see the same verdict. Counts feed report section 3.4.
    from .confidence import annotate_orthodb_support
    od = annotate_orthodb_support(
        final_gff3, assess,
        busco_lineage=str(busco_lineage) if busco_lineage else None,
        omark_db=str(omark_db) if omark_db else None)
    ctx["orthodb"] = {
        "applied": od["applied"],
        "busco_supported": od["busco_supported"],
        "omark_supported": od["omark_supported"],
        "busco_lineage": Path(busco_lineage).name if busco_lineage else None,
        "omark_db": Path(omark_db).name if omark_db else None,
    }
    # section 3.0: score the final annotation with the predictor P formula
    ctx["final_score"] = final_annotation_score(
        out, final_gff3, genome_fa, threads, ctx,
        rna_junctions=rna_junctions, stringtie_gtf=stringtie_gtf,
        homologs_faa=homologs_faa, referee_ids=referee_ids)
    # finish time + elapsed, for section 1 (Run Summary)
    from datetime import datetime
    _end = datetime.now()
    ctx["meta"]["finished_at"] = _end.isoformat(timespec="seconds")
    ctx["meta"]["elapsed"] = _elapsed_str(ctx["meta"].get("timestamp", ""), _end)
    final_sequences(final_gff3, genome_fa, d)
    html = render_html(ctx, d / "report.html", None)
    with open(d / "stats.json", "w") as fh:
        json.dump(ctx, fh, ensure_ascii=False, indent=2)
    log.info("report stats -> %s", d / "stats.json")
    return {"report_dir": d, "ctx": ctx, "html": html}
