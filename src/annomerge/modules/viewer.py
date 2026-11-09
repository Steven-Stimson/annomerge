"""Evidence browser (igv.js) — auto-generated after every merge run.

Produces, under <out>/viewer/:
  data/                              genome.fa(.fai), final.gff3, goldset.gff3,
                                     protein.gff3, stringtie.gtf, junctions.bed,
                                     denovo_<name>.gff3, coverage bedgraph
                                     (optional), genes.json, tracks.json
  AnnoMerge_evidence_browser.html    single-file offline browser (all data
                                     inlined as dataURIs, double-click to open)

The igv.js bundle and the page template live in the package assets
(src/annomerge/assets/), so the browser works without network access.
"""
from __future__ import annotations

import base64
import json
import logging
import os
import re
import shutil
from pathlib import Path

from ..utils import ensure_dir, run, which

log = logging.getLogger("annomerge")

_ASSETS = Path(__file__).resolve().parent.parent / "assets"

_PALETTE = ["#d62728", "#e07b39", "#2ca02c", "#9467bd", "#8c564b",
            "#17becf", "#bcbd22", "#7f7f7f"]

_TIER_RANK = {"strong": 0, "moderate": 1, "weak": 2}

_FLAG_MAP = {
    "single_exon": "single-exon",
    "noncanonical_splice": "noncanonical-splice",
    "neighbour_overlap": "neighbour-overlap",
    "short_orf": "short-ORF",
}
_NOTE_MAP = {
    "no_stop": "no-stop",
    "non_atg_start": "non-ATG-start",
    "extended_stop": "extended-stop",
}


def _link_or_copy(src: Path, dst: Path):
    if dst.exists() and dst.stat().st_size == src.stat().st_size:
        return
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def _datauri(path: Path) -> str:
    with open(path, "rb") as fh:
        return "data:text/plain;base64," + base64.b64encode(fh.read()).decode()


# ------------------------------------------------------------- genes.json
def _parse_final_gff3(final_gff3: Path) -> tuple:
    """Per-gene structural counts + gene order of appearance (GFF order).

    Returns (info, order_index, tx_cds_exon, tx_cds_iv):
      info[gid] = {tx, cds, exon}; order_index[gid] = position among gene
      features; tx_cds_exon[tid] = [cds_n, exon_n]; tx_cds_iv[tid] = sorted
      [(start, end)] CDS intervals of the transcript.
    """
    info: dict = {}
    order: list = []
    tx_parent: dict = {}
    tx_ce: dict = {}
    tx_cds_iv: dict = {}
    with open(final_gff3) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] == "":
                continue
            ftype, attrs = p[2], p[8]
            av = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    av[k] = v
            if ftype == "gene":
                gid = av.get("ID")
                if gid and gid not in info:
                    info[gid] = {"tx": 0, "cds": 0, "exon": 0}
                    order.append(gid)
            elif ftype == "mRNA":
                tid = av.get("ID")
                gid = (av.get("Parent") or "").split(",")[0]
                if tid:
                    tx_parent[tid] = gid
                    tx_ce.setdefault(tid, [0, 0])
                    tx_cds_iv.setdefault(tid, [])
                if gid in info:
                    info[gid]["tx"] += 1
            elif ftype in ("CDS", "exon"):
                for parent in (av.get("Parent") or "").split(","):
                    if parent in tx_ce:
                        tx_ce[parent][0 if ftype == "CDS" else 1] += 1
                    gid = tx_parent.get(parent)
                    if gid in info:
                        info[gid]["cds" if ftype == "CDS" else "exon"] += 1
                    if ftype == "CDS" and parent in tx_cds_iv:
                        s, e = int(p[3]), int(p[4])
                        tx_cds_iv[parent].append((s, e))
    for iv in tx_cds_iv.values():
        iv.sort()
    return info, {g: i for i, g in enumerate(order)}, tx_ce, tx_cds_iv


def _denovo_gene_coverage(gff3: Path, genes: list,
                          tx_cds_iv: dict | None = None) -> tuple:
    """Per-transcript and per-gene fraction of CDS covered by any feature of
    this de novo prediction GFF3.  CDS-basis (NOT whole-span): UTR-less
    predictors (e.g. annevo, CDS-only) stay comparable with UTR-inclusive
    ones (e.g. helixer); whole-span coverage would systematically penalise
    the former.  Gene-level value = coverage of the UNION of all isoforms'
    CDS intervals (overlaps counted once).  Single pass over the file;
    records grouped by chrom.  Returns ({transcript_id: frac},
    {gene_id: frac})."""
    by_chrom: dict = {}
    for g in genes:
        entry = {"gene": g["gene"], "iso": []}
        for iso in g.get("isoforms", []):
            cds = (tx_cds_iv or {}).get(iso["tid"])
            if not cds:  # no CDS parsed: fall back to whole span
                cds = [(iso["start"], iso["end"])]
            entry["iso"].append({"tid": iso["tid"], "cds": cds})
        by_chrom.setdefault(g["seqid"], []).append(entry)
    intervals: dict = {}
    with open(gff3) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] == "":
                continue
            try:
                s, e = int(p[3]), int(p[4])
            except ValueError:
                continue
            intervals.setdefault(p[0], []).append((s, e))

    def _union(ivs: list) -> list:
        iv = []
        for s, e in sorted(ivs):
            if iv and s <= iv[-1][1] + 1:
                iv[-1] = (iv[-1][0], max(iv[-1][1], e))
            else:
                iv.append((s, e))
        return iv

    def _covered(cds: list, iv: list) -> float:
        """fraction of the (sorted, disjoint) CDS intervals covered by iv."""
        total = sum(e - s + 1 for s, e in cds)
        if not total:
            return 0.0
        cov = 0
        j = 0
        for s, e in cds:
            while j < len(iv) and iv[j][1] < s:
                j += 1
            for ps, pe in iv[j:]:
                if ps > e:
                    break
                cov += min(pe, e) - max(ps, s) + 1
        return round(cov / total, 3)

    tx_out: dict = {}
    gene_out: dict = {}
    for chrom, entries in by_chrom.items():
        iv = _union(intervals.get(chrom, []))
        for entry in entries:
            gene_cds: list = []
            for t in entry["iso"]:
                if not iv:
                    tx_out[t["tid"]] = 0.0
                    continue
                tx_out[t["tid"]] = _covered(t["cds"], iv)
                gene_cds.extend(t["cds"])
            if not iv:
                gene_out[entry["gene"]] = 0.0
                continue
            gene_out[entry["gene"]] = _covered(_union(gene_cds), iv)
    return tx_out, gene_out


def build_genes_json(confidence_tsv: Path, out_json: Path, *,
                     final_gff3: Path | None = None,
                     denovo_gff3: dict | None = None) -> int:
    """Aggregate confidence.tsv -> one viewer record per gene."""
    genes: dict = {}
    with open(confidence_tsv) as fh:
        for ln in fh:
            if ln.startswith("#") or ln.startswith("gene_id"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 18:
                continue
            (gid, _tid, seqid, strand, start, end, tier, gold,
             rna_frac, rna_basis, cov_frac, depth, rna_reached,
             prot_frac, prot_basis, prot_reached, sanity, flags) = p[:18]
            g = genes.setdefault(gid, {
                "gene": gid, "seqid": seqid, "strand": strand,
                "start": int(start), "end": int(end), "tier": tier,
                "iso": 0, "rna": 0.0, "prot": 0.0, "cov": 0.0, "depth": 0.0,
                "goldset": "no", "rna_reached": "no", "prot_reached": "no",
                "rna_basis": "", "prot_basis": "",
                "flags": set(), "notes": set(), "isoforms": [],
            })
            g["iso"] += 1
            g["start"] = min(g["start"], int(start))
            g["end"] = max(g["end"], int(end))
            g["rna"] = max(g["rna"], float(rna_frac))
            g["prot"] = max(g["prot"], float(prot_frac))
            g["cov"] = max(g["cov"], float(cov_frac))
            g["depth"] = max(g["depth"], float(depth))
            if gold == "yes":
                g["goldset"] = "yes"
            if rna_reached == "yes":
                g["rna_reached"] = "yes"
            if prot_reached == "yes":
                g["prot_reached"] = "yes"
            if float(rna_frac) >= g["rna"]:
                g["rna_basis"] = rna_basis
            if float(prot_frac) >= g["prot"]:
                g["prot_basis"] = prot_basis
            if flags != ".":
                for fl in flags.split(","):
                    g["flags"].add(fl)
            if sanity != ".":
                for nt in sanity.split(","):
                    g["notes"].add(nt)
            # per-transcript record (metrics are computed per transcript)
            tx_marks = []
            if flags != ".":
                tx_marks += [_FLAG_MAP[f] for f in sorted(set(flags.split(",")))
                             if f in _FLAG_MAP]
            if sanity != ".":
                tx_marks += [_NOTE_MAP[n] for n in sorted(set(sanity.split(",")))
                             if n in _NOTE_MAP]
            g["isoforms"].append({
                "tid": _tid, "start": int(start), "end": int(end),
                "tier": tier,
                "rna": f"{float(rna_frac):.3f}", "rna_basis": rna_basis,
                "rna_reached": rna_reached,
                "prot": f"{float(prot_frac):.3f}", "prot_basis": prot_basis,
                "prot_reached": prot_reached,
                "cov": f"{float(cov_frac):.3f}", "depth": float(depth),
                "gold": gold, "marks": ",".join(tx_marks),
                "cds": 0, "exon": 0,
            })

    def gid_num(g: str) -> int:
        m = re.search(r"(\d+)$", g)
        return int(m.group(1)) if m else 0

    out = []
    for gid, g in genes.items():
        marks = [_FLAG_MAP[f] for f in sorted(g["flags"]) if f in _FLAG_MAP]
        marks += [_NOTE_MAP[n] for n in sorted(g["notes"]) if n in _NOTE_MAP]
        out.append({
            "gene": gid,
            "seqid": g["seqid"],
            "start": g["start"],
            "end": g["end"],
            "locus": f'{g["seqid"]}:{g["start"]}-{g["end"]}',
            "strand": g["strand"],
            "tier": g["tier"],
            "tx": g["iso"],
            "rna": f'{g["rna"]:.3f}',
            "prot": f'{g["prot"]:.3f}',
            "cov": f'{g["cov"]:.3f}',
            "depth": round(g["depth"]),
            "goldset": g["goldset"],
            "rna_reached": g["rna_reached"],
            "prot_reached": g["prot_reached"],
            "rna_basis": g["rna_basis"],
            "prot_basis": g["prot_basis"],
            "isoforms": g["isoforms"],
            "marks": ",".join(marks),
            "flag_desc": "; ".join(_FLAG_MAP[f] for f in sorted(g["flags"])
                                   if f in _FLAG_MAP),
            "note_desc": "; ".join(_NOTE_MAP[n] for n in sorted(g["notes"])
                                   if n in _NOTE_MAP),
        })
    # structural counts + GFF order + CDS intervals from the final annotation
    tx_cds_iv: dict = {}
    if final_gff3 and Path(final_gff3).exists():
        info, order_idx, tx_ce, tx_cds_iv = _parse_final_gff3(Path(final_gff3))
        for r in out:
            st = info.get(r["gene"])
            if st:
                r["tx"] = st["tx"]
                r["cds"] = st["cds"]
                r["exon"] = st["exon"]
            r["gff_order"] = order_idx.get(r["gene"], 10 ** 9)
            for iso in r["isoforms"]:
                ce = tx_ce.get(iso["tid"])
                if ce:
                    iso["cds"], iso["exon"] = ce
    else:
        for r in out:
            r["cds"] = 0
            r["exon"] = 0
            r["gff_order"] = 10 ** 9
    # de novo evidence coverage per predictor (read from existing GFF3s);
    # per transcript AND per gene (gene = union of all isoforms' CDS)
    for name, path in (denovo_gff3 or {}).items():
        if path and Path(path).exists():
            tx_map, gene_map = _denovo_gene_coverage(Path(path), out, tx_cds_iv)
            for r in out:
                for iso in r["isoforms"]:
                    iso.setdefault("denovo", {})[name] = \
                        tx_map.get(iso["tid"], 0.0)
                r.setdefault("denovo", {})[name] = \
                    gene_map.get(r["gene"], 0.0)

    out.sort(key=lambda r: (_TIER_RANK.get(r["tier"], 9), gid_num(r["gene"])))
    with open(out_json, "w") as fh:
        json.dump(out, fh, ensure_ascii=False)
    log.info("viewer genes.json: %d genes", len(out))
    return len(out)


# ------------------------------------------------------------ tracks.json
def _split_protein_tracks(pgff: Path, data_dir: Path,
                          set_weights: dict) -> list:
    """Split protein evidence GFF3 into one file/track per homolog set.

    Multi-set runs: 06.merge/protein.gff3 lines carry miniprot_<label>
    sources -> one track per set, ordered by relatedness weight desc, the
    closest set keeping the legacy cyan. Single-set runs: one legacy
    "Protein (miniprot)" track. Returns [{set, weight, file, color}].
    """
    by_src: dict = {}
    with open(pgff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            by_src.setdefault(f[1], []).append(line)
    tracks = []
    multi = any(s.startswith("miniprot_") for s in by_src)
    for src in by_src:
        label = src[len("miniprot_"):] if src.startswith("miniprot_") else "all"
        fname = f"protein_{label}.gff3"
        with open(data_dir / fname, "w") as out:
            out.write("##gff-version 3\n")
            out.writelines(by_src[src])
        w = set_weights.get(label)
        tracks.append({"set": label, "weight": round(w, 4) if w else None,
                       "file": fname})
    tracks.sort(key=lambda t: -(t["weight"] or 1.0))
    for rank, t in enumerate(tracks):
        # closest set keeps the legacy cyan; others rotate the palette
        t["color"] = "#17becf" if rank == 0 else _PALETTE[rank % len(_PALETTE)]
    if not multi and tracks:
        tracks[0]["color"] = "#17becf"
    log.info("viewer protein tracks: %s",
             ", ".join(f"{t['set']}(w={t['weight']})" for t in tracks))
    return tracks


def build_tracks_json(preds: dict, details: dict, data_dir: Path,
                      out_json: Path, protein_tracks: list | None = None,
                      coverage: dict | None = None) -> list:
    """denovo predictor tracks, colored by W_raw (desc). Copies each GFF3
    into data_dir as denovo_<name>.gff3. Also records the per-set protein
    tracks and the per-sample/merged RNA coverage metadata."""
    denovo = []
    for name, path in preds.items():
        w = 0.0
        d = (details or {}).get(name) or {}
        if isinstance(d.get("W_raw"), (int, float)):
            w = round(float(d["W_raw"]), 3)
        dst = data_dir / f"denovo_{name}.gff3"
        if Path(path).exists() and not dst.exists():
            _link_or_copy(Path(path), dst)
        denovo.append({"pred": name, "w": w, "file": f"denovo_{name}.gff3"})
    denovo.sort(key=lambda d: -d["w"])
    for rank, d in enumerate(denovo):
        d["color"] = _PALETTE[rank % len(_PALETTE)]
    with open(out_json, "w") as fh:
        json.dump({"denovo": denovo, "protein": protein_tracks or [],
                   "coverage": coverage or {}},
                  fh, ensure_ascii=False, indent=1)
    return denovo


# ---------------------------------------------------------------- coverage
def coverage_bedgraph(bam: Path, out_bedgraph: Path, threads: int) -> Path | None:
    """Whole-genome 50-bp-bin coverage bedGraph.

    Preferred: deeptools bamCoverage.  Fallback: samtools depth -a (per-base,
    then averaged into 50-bp bins) so the coverage track appears even without
    deeptools.  Any failure -> None and the browser renders without the track.
    """
    if out_bedgraph.exists() and out_bedgraph.stat().st_size > 0:
        return out_bedgraph
    import shutil
    py = None
    env_dir = os.environ.get("DEEPTOOLS")
    if env_dir:
        cand = Path(env_dir) / "env/bin/python"
        if cand.exists():
            py = str(cand)
    bamcov = shutil.which("bamCoverage")
    if py is not None or bamcov is not None:
        try:
            if py:
                run([py, str(Path(py).parent / "bamCoverage"), "-b", str(bam),
                     "-o", str(out_bedgraph), "--outFileFormat", "bedgraph",
                     "--binSize", "50", "--numberOfProcessors", str(threads)],
                    log_file=out_bedgraph.parent / "coverage.log",
                    desc="viewer coverage bedgraph")
            else:
                run([bamcov, "-b", str(bam), "-o", str(out_bedgraph),
                     "--outFileFormat", "bedgraph", "--binSize", "50",
                     "--numberOfProcessors", str(threads)],
                    log_file=out_bedgraph.parent / "coverage.log",
                    desc="viewer coverage bedgraph")
            return out_bedgraph if out_bedgraph.exists() else None
        except Exception as e:
            log.warning("viewer: bamCoverage failed (%s); trying samtools", e)
    # fallback: samtools depth -a -> 50-bp bin means
    sam = shutil.which("samtools")
    if sam is None:
        log.warning("viewer: no bamCoverage and no samtools; coverage track skipped")
        return None
    try:
        raw = out_bedgraph.parent / "coverage.depth.tsv"
        run([sam, "depth", "-a", "-@", str(threads), str(bam)],
            log_file=raw, desc="viewer samtools depth")
        # bin per-base depth into 50-bp windows per reference
        bins: dict = {}
        with open(raw) as fh:
            for line in fh:
                ref, pos, dp = line.split("\t")
                key = (ref, int(pos) // 50)
                s, n = bins.get(key, (0, 0))
                bins[key] = (s + int(dp), n + 1)
        with open(out_bedgraph, "w") as out:
            for (ref, b), (s, n) in sorted(bins.items(),
                                           key=lambda x: (x[0][0], x[0][1])):
                out.write(f"{ref}\t{b * 50}\t{(b + 1) * 50}\t"
                          f"{s / n:.1f}\n")
        raw.unlink(missing_ok=True)
    except Exception as e:
        log.warning("viewer: coverage bedgraph failed (%s); track skipped", e)
        return None
    return out_bedgraph if out_bedgraph.exists() else None


def sample_coverage(bams: list, data_dir: Path, threads: int) -> dict:
    """Per-sample + merged 50-bp-bin coverage bedgraphs.

    One bedgraph per BAM (coverage_<sample>.bedgraph) plus, for >=2 BAMs, a
    per-bin summed coverage_all.bedgraph. Returns coverage metadata
    {"merged": {file, max, mean, n_samples, name?}, "samples": [{name, file,
    max, mean}]} or {} when no BAM produced a bedgraph. Track y-axes are left
    to igv per-track autoscale; max/mean are recorded here for the track
    labels so different sequencing depths stay comparable.
    """
    samples = []
    merged: dict = {}          # (ref, start) -> [end, depth_sum]
    used_names: set = set()
    for bam in bams:
        bam = Path(bam)
        if not bam.exists():
            continue
        base = bam.name.split(".")[0].replace(".sorted", "")
        name, i = base, 2
        while name in used_names:
            name = f"{base}_{i}"
            i += 1
        used_names.add(name)
        bg = data_dir / f"coverage_{name}.bedgraph"
        if not coverage_bedgraph(bam, bg, threads):
            log.warning("viewer: coverage failed for %s; sample track skipped", bam)
            continue
        mx = tot = 0.0
        n = 0
        vals = []
        with open(bg) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                v = float(f[3])
                vals.append((f[0], int(f[1]), int(f[2]), v))
                mx = max(mx, v)
                tot += v
                n += 1
        for ref, s, e, v in vals:
            rec = merged.setdefault((ref, s), [e, 0.0])
            rec[1] += v
        samples.append({"name": name, "file": bg.name,
                        "max": round(mx, 1),
                        "mean": round(tot / n, 1) if n else 0.0})
    if not samples:
        return {}
    if len(samples) > 1:
        mbg = data_dir / "coverage_all.bedgraph"
        mx = tot = 0.0
        n = 0
        with open(mbg, "w") as out:
            for (ref, s), (e, dsum) in sorted(merged.items()):
                out.write(f"{ref}\t{s}\t{e}\t{dsum:.1f}\n")
                mx = max(mx, dsum)
                tot += dsum
                n += 1
        merged_info = {"file": mbg.name, "max": round(mx, 1),
                       "mean": round(tot / n, 1) if n else 0.0,
                       "n_samples": len(samples)}
    else:
        merged_info = dict(samples[0])
        merged_info["n_samples"] = 1
    log.info("viewer coverage: %d sample(s), merged max %.0fx mean %.1fx",
             len(samples), merged_info["max"], merged_info["mean"])
    return {"merged": merged_info, "samples": samples}


_COV_SAMPLE_PALETTE = ["#31a354", "#756bb1", "#e6550d", "#636363",
                       "#dd1c77", "#2c7fb8"]


# ---------------------------------------------------------- offline bundle
def _inline_script_vars(files: dict, genes_json: str, tracks_cfg: dict,
                        fai_txt: str) -> str:
    return ("<script>\nconst DATA_URI = " + json.dumps(files) + ";\n"
            "const GENES_JSON = " + json.dumps(genes_json) + ";\n"
            "const TRACKS_JSON = " +
            json.dumps(json.dumps(tracks_cfg, ensure_ascii=False)) + ";\n"
            "const FAI_TXT = " + json.dumps(fai_txt) + ";\n"
            # lite reference: chromsizes served from a blob URL (pure
            # coordinates, no base sequence; fetch blob: works on file://)
            "const REF_URL = URL.createObjectURL(new Blob([FAI_TXT], "
            "{type: \"text/plain\"}));\n</script>\n")


def render_offline_html(data_dir: Path, out_html: Path, *,
                        ref_name: str, fallback_locus: str,
                        coverage: dict | None,
                        full: bool = False) -> Path:
    """Pack the browser into a single offline HTML (data inlined as dataURIs).

    full=False (default, lite): the genome sequence is NOT embedded; the
    reference is pure coordinates (chromsizes via blob URL) so the file
    stays small even for large genomes.  full=True additionally embeds the
    genome FASTA as a dataURI so base-level sequence is shown in igv.
    coverage: sample_coverage() metadata (merged + per-sample bedgraphs).
    """
    html = (_ASSETS / "browser_template.html").read_text(encoding="utf-8")
    igv_js = (_ASSETS / "igv.min.js").read_text(encoding="utf-8")

    # 1) inline igv.min.js
    html = html.replace('<script src="js/igv.min.js"></script>',
                        '<script>\n' + igv_js + '\n</script>')

    # 2) data files -> dataURI
    tracks_cfg = json.load(open(data_dir / "tracks.json"))
    files: dict = {}
    for name in ("final.gff3", "goldset.gff3", "protein.gff3", "stringtie.gtf",
                 "junctions.bed"):
        p = data_dir / name
        if p.exists():
            files[name] = _datauri(p)
    if full:
        # lite mode skips the genome sequence entirely
        files["genome.fa"] = _datauri(data_dir / "genome.fa")
    # merged + per-sample RNA coverage bedgraphs
    covcfg = tracks_cfg.get("coverage") or {}
    cov_files = []
    if covcfg.get("merged"):
        cov_files.append(covcfg["merged"]["file"])
    for s in covcfg.get("samples", []):
        cov_files.append(s["file"])
    for f in dict.fromkeys(cov_files):
        p = data_dir / f
        if p.exists():
            files[f] = _datauri(p)
    # cytoband from fai (full mode only; chromsizes reference cannot take it)
    if full:
        cytoband = []
        for ln in (data_dir / "genome.fa.fai").read_text().splitlines():
            p = ln.split("\t")
            if len(p) >= 2:
                cytoband.append(f"{p[0]}\t0\t{p[1]}\t{p[0]}\tgneg")
        files["cytoband.txt"] = ("data:text/plain;base64,"
                                 + base64.b64encode(("\n".join(cytoband) + "\n")
                                                    .encode()).decode())
    for d in tracks_cfg.get("denovo", []):
        p = data_dir / d["file"]
        if p.exists():
            files[d["file"]] = _datauri(p)
    # per-set protein evidence tracks (multi-set homolog runs)
    for t in tracks_cfg.get("protein", []):
        p = data_dir / t["file"]
        if p.exists():
            files[t["file"]] = _datauri(p)
    genes_json = (data_dir / "genes.json").read_text().replace("</", "<\\/")
    fai_txt = (data_dir / "genome.fa.fai").read_text()

    # 3) url references -> dataURI lookups
    for name in files:
        html = html.replace(f'url: "data/{name}"', f'url: DATA_URI["{name}"]')
    html = html.replace('url: "data/" + d.file', 'url: DATA_URI[d.file]')

    # 4) reference block: drop the id (an id makes igv consult its built-in
    #    genome registry and fail with "Unknown genome id" offline)
    html = html.replace('name: "$REF_NAME$"', f'name: "{ref_name}"')
    if full:
        html = html.replace(
            'fastaURL: "data/genome.fa",\n      indexURL: "data/genome.fa.fai",\n      cytobandURL: "data/cytoband.txt"',
            'fastaURL: DATA_URI["genome.fa"],\n      '
            'cytobandURL: DATA_URI["cytoband.txt"]')
    else:
        # lite: pure-coordinate reference (chromsizes via blob URL).
        # igv only picks its chromsizes loader when format is set explicitly;
        # without it the fai text is misparsed as FASTA and genome init fails.
        # No cytobandURL here — the chromsizes loader lacks getChromosomeName
        # and crashes when a cytoband is present.
        html = html.replace(
            'fastaURL: "data/genome.fa",\n      indexURL: "data/genome.fa.fai",\n      cytobandURL: "data/cytoband.txt"',
            'fastaURL: REF_URL,\n      format: "chromsizes"')
    html = html.replace("$FALLBACK_LOCUS$", fallback_locus)

    # 5) coverage tracks: merged depth over ALL RNA samples on top, then one
    #    wig per sample (>=2 samples). Y-axes are per-track autoscale — no
    #    fixed upper bound, each track adapts to its own depth range; the
    #    observed max depth is shown in the track label for comparability.
    def _wig(name: str, fname: str, color: str) -> str:
        return ('  { name: ' + json.dumps(name) + ', type: "wig",\n'
                '    format: "bedgraph",\n'
                '    url: DATA_URI[' + json.dumps(fname) + '], color: '
                + json.dumps(color) + ',\n'
                '    height: 90, min: 0, autoscale: true }')

    cov_tracks = []
    m = covcfg.get("merged") or {}
    multi = (m.get("n_samples") or 1) > 1
    if m and m.get("file") in files:
        if multi:
            label = f'RNA coverage (all {m["n_samples"]} samples, ~{m["max"]:,.0f}x max)'
        else:
            label = f'RNA coverage: {m.get("name", "sample")} (~{m["max"]:,.0f}x max)'
        cov_tracks.append(_wig(label, m["file"], "#3182bd"))
    if multi:
        for rank, s in enumerate(covcfg.get("samples", [])):
            if s["file"] not in files:
                continue
            color = _COV_SAMPLE_PALETTE[rank % len(_COV_SAMPLE_PALETTE)]
            label = f'RNA coverage: {s["name"]} (~{s["max"]:,.0f}x max)'
            cov_tracks.append(_wig(label, s["file"], color))
    html = html.replace("/*__COVERAGE_TRACK__*/",
                        ",\n".join(cov_tracks) + ("," if cov_tracks else ""))

    # 5b) protein evidence tracks: one per homolog set (relatedness-weighted
    #     multi-set runs) or the single legacy track
    protein_js = []
    for t in tracks_cfg.get("protein", []):
        if t["file"] not in files:
            continue
        if t["set"] == "all":
            name = "Protein (miniprot)"
        else:
            wtxt = f", w={float(t['weight']):.2f}" if t.get("weight") else ""
            name = f"Protein (miniprot: {t['set']}{wtxt})"
        protein_js.append(
            '  { name: ' + json.dumps(name) + ', type: "annotation", format: "gff3",\n'
            '    url: DATA_URI[' + json.dumps(t["file"]) + '], color: '
            + json.dumps(t["color"]) + ', height: 120,\n'
            '    visibilityWindow: -1 }')
    html = html.replace("/*__PROTEIN_TRACKS__*/", ",\n".join(protein_js))

    # 6) genes/tracks/fai fetches -> inlined constants
    html = html.replace('const resp = await fetch("data/genes.json");\n  genes = await resp.json();',
                        'genes = JSON.parse(GENES_JSON);')
    html = html.replace('const resp = await fetch("data/tracks.json");\n  const cfg = await resp.json();',
                        'const cfg = JSON.parse(TRACKS_JSON);')
    html = html.replace('const resp = await fetch("data/genome.fa.fai");\n  const text = await resp.text();',
                        'const text = FAI_TXT;')

    # 7) inject DATA_URI / GENES_JSON / TRACKS_JSON / FAI_TXT
    inject = _inline_script_vars(files, genes_json, tracks_cfg, fai_txt)
    html = html.replace('<div id="topbar">', inject + '<div id="topbar">')

    # 8) offline identity
    html = html.replace("<title>AnnoMerge Evidence Browser (igv.js)</title>",
                        "<title>AnnoMerge Evidence Browser (offline single-file)</title>")

    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    leftover = [l for l in html.splitlines()
                if 'url: "data/' in l or 'fetch("data/' in l]
    log.info("viewer offline html -> %s (%.1f MB, %d leftover data refs)",
             out_html, out_html.stat().st_size / 1048576, len(leftover))
    if leftover:
        for l in leftover[:5]:
            log.warning("viewer leftover ref: %s", l.strip())
    return out_html


# ------------------------------------------------------------------ driver
def run_viewer(out: Path, vi, results: dict, preds: dict, details: dict,
               final_gff3: Path, final_d: Path, integ_d: Path,
               threads: int, full: bool = False) -> dict:
    """Assemble data/ and the offline single-file evidence browser.

    full=False (default): lite browser, genome sequence not embedded.
    full=True: also embed the FASTA so igv shows base-level sequence.
    """
    vdir = ensure_dir(Path(out) / "09.result" / "viewer")
    data = ensure_dir(vdir / "data")

    # 1) genome (+ fai)
    genome = Path(vi.genome)
    _link_or_copy(genome, data / "genome.fa")
    fai = data / "genome.fa.fai"
    if not fai.exists():
        if (genome.parent / (genome.name + ".fai")).exists():
            _link_or_copy(genome.parent / (genome.name + ".fai"), fai)
        else:
            run([which("samtools"), "faidx", str(data / "genome.fa")],
                desc="viewer samtools faidx")

    # 2) annotation + evidence tracks
    _link_or_copy(Path(final_gff3), data / "final.gff3")
    for cand in (final_d / "goldset.genome.gff3", integ_d / "goldset.genome.gff3"):
        if cand.exists():
            _link_or_copy(cand, data / "goldset.gff3")
            break
    pgff = None
    for cand in (integ_d / "protein.gff3",
                 Path((results.get("homolog") or {}).get("gff") or "")):
        if cand and cand.exists():
            pgff = cand
            break
    # protein evidence split into one track per homolog set (multi-set runs:
    # 06.merge/protein.gff3 carries miniprot_<label> sources). Single-set
    # runs keep one legacy "Protein (miniprot)" track.
    set_weights = (results.get("homolog") or {}).get("set_weights") or {}
    protein_tracks = _split_protein_tracks(pgff, data, set_weights) if pgff else []
    rna = results.get("rna") or {}
    st_gtf = rna.get("stringtie_gtf")
    if st_gtf and Path(st_gtf).exists():
        _link_or_copy(Path(st_gtf), data / "stringtie.gtf")
    junc = rna.get("junctions_bed")
    if junc and Path(junc).exists():
        _link_or_copy(Path(junc), data / "junctions.bed")

    # coverage bedgraphs: one per BAM + a merged track over ALL samples;
    # must run BEFORE build_tracks_json so its metadata lands in tracks.json
    cov = sample_coverage(rna.get("bams") or [], data, threads)

    # 3) tracks.json first: it copies denovo_<name>.gff3 into data/, which
    #    build_genes_json needs on disk to compute de novo evidence coverage
    build_tracks_json(preds, details, data, data / "tracks.json",
                      protein_tracks=protein_tracks,
                      coverage=cov)
    build_genes_json(Path(out) / "08.assessment" / "confidence.tsv",
                     data / "genes.json",
                     final_gff3=data / "final.gff3",
                     denovo_gff3={n: data / f"denovo_{n}.gff3"
                                  for n in (preds or {})})

    # 4) coverage bedgraphs (already computed above, before tracks.json)
    # 5) pack the single-file offline browser
    first_chrom = (data / "genome.fa.fai").read_text().split("\t")[0]
    chrom_len = (data / "genome.fa.fai").read_text().split("\t")[1]
    fallback = f"{first_chrom}:1-{min(int(chrom_len), 200000)}"
    html = render_offline_html(
        data, vdir / "AnnoMerge_evidence_browser.html",
        ref_name=genome.name, fallback_locus=fallback,
        coverage=cov, full=full)
    return {"dir": vdir, "data": data, "html": html, "full": full}
