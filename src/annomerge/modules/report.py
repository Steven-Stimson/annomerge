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
        run(f"cd {bdir} && PYTHONPATH= {which('busco')} -i {faa} -m proteins "
            f"-l {lineage} -o busco_run --offline -c {threads} -f",
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
    # EVidenceModeler ships no --version; the release number lives in the
    # top-level directory name of the distribution tarball (EVidenceModeler-X.Y.Z)
    evm_v = None
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
def render_html(ctx: dict, out_html: Path,
                meta: dict | None = None) -> Path:
    """Render the single-page report."""
    from jinja2 import Environment, BaseLoader
    env = Environment(loader=BaseLoader(), autoescape=True)
    tpl = env.from_string(_TEMPLATE)
    html = tpl.render(meta=meta or ctx.get("meta", {}),
                      gene_set=ctx.get("gene_set", {}),
                      busco=ctx.get("busco", {}),
                      omark=ctx.get("omark", {}),
                      evidence=ctx.get("evidence", {}),
                      versions=ctx.get("versions", {}),
                      assessment=ctx.get("assessment", {}),
                      integ=ctx.get("integ"))
    Path(out_html).write_text(html)
    log.info("report html -> %s", out_html)
    return Path(out_html)


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
    return {
        "evidence_rows": ev_rows,
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
               omark_db: Path | None = None) -> dict:
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
    # integration-only sections 2.1/2.2 (None for evaluate runs)
    ctx["integ"] = _integration_stats(out, final_gff3)
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
