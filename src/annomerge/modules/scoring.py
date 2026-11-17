"""Module 3: predictor global scoring.

P = w_r * R_s + w_h * H_s'' + w_b * B_s + w_o * O_s
                          (defaults w_r=0.45, w_h=0.25, w_b=0.10, w_o=0.20)

Each sub-block maps to exactly ONE evidence source (no double counting):
  R_s: the RNA evidence bundle — (1-d)*(alpha*J_s* + beta*E_sp + gamma*M_e)
       + d*GoldMatch, with d = --rna-gold-weight (default 0.25)
      J_s*: junction support rate (multi-intron transcripts only)
      E_sp: exon precision, gffcompare vs StringTie transcripts
      M_e : single-exon bias correction vs RNA single-exon fraction
      GoldMatch: full-length protein match rate vs the RNA-assembled gold set
  H_s'': structure-aware homolog support — fraction of predicted proteins that
      are BOTH full-length bidirectional matches vs the homolog referee set
      AND structurally clean (no chimeric/split anomaly from the same
      miniprot referee hits). Pure homolog evidence; the RNA-derived gold set
      is mixed in ONLY as a rescue when the referee set has fewer than
      --hs-referee-min proteins (near-reference too distant).
  B_s : BUSCO complete recall (proteins mode, local lineage, offline)
  O_s : OMArk Consistent fraction (needs --omark-db; skipped otherwise)

Degradation: missing evidence -> NA + weight renormalization (w_i/(1-w_missing));
fewer than two usable sub-blocks -> refuse to score.
"""
from __future__ import annotations

import logging
import math
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from ..utils import ensure_dir, lazy, run, which

log = logging.getLogger("annomerge")


# ---------------------------------------------------------------- GFF parsing

@dataclass
class Transcript:
    tid: str
    seqid: str
    strand: str
    exons: list = field(default_factory=list)   # (start, end) 1-based
    cds: list = field(default_factory=list)     # (start, end) 1-based

    @property
    def introns(self) -> list[tuple[int, int]]:
        """Introns from CDS chain (1-based, inclusive)."""
        ex = sorted(self.cds)
        return [(a[1] + 1, b[0] - 1) for a, b in zip(ex, ex[1:])]

    @property
    def exon_introns(self) -> list[tuple[int, int]]:
        ex = sorted(self.exons)
        return [(a[1] + 1, b[0] - 1) for a, b in zip(ex, ex[1:])]


def parse_gff(gff: Path) -> dict[str, Transcript]:
    """Lightweight GFF3/GTF parser -> {transcript_id: Transcript}."""
    txs: dict[str, Transcript] = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feat = f[2]
            if feat not in ("exon", "CDS"):
                continue
            attrs = f[8]
            m = re.search(r'Parent=([^;]+)', attrs) or re.search(r'transcript_id "([^"]+)"', attrs)
            if not m:
                continue
            tid = m.group(1)
            tx = txs.get(tid)
            if tx is None:
                tx = txs[tid] = Transcript(tid=tid, seqid=f[0], strand=f[6])
            (tx.exons if feat == "exon" else tx.cds).append((int(f[3]), int(f[4])))
    return txs


def load_junctions(bed: Path) -> set[tuple[str, int, int]]:
    """junctions.bed -> {(seqid, start_1based, end_1based)}."""
    out = set()
    with open(bed) as fh:
        for line in fh:
            c, s, e = line.split("\t")[:3]
            out.add((c, int(s) + 1, int(e)))
    return out


# ---------------------------------------------------------------- R_s

def junction_support(txs: dict[str, Transcript], junc: set) -> tuple[float | None, dict]:
    """J_s*: fraction of predictor introns supported by RNA junctions.
    Only transcripts with >= 1 intron are considered."""
    per_tx = {}
    n_sup = n_tot = 0
    for tid, tx in txs.items():
        ivs = tx.introns
        if not ivs:
            continue
        hits = sum(1 for (s, e) in ivs if (tx.seqid, s, e) in junc)
        per_tx[tid] = (hits, len(ivs))
        n_sup += hits
        n_tot += len(ivs)
    js = n_sup / n_tot if n_tot else None
    return js, {"supported_introns": n_sup, "total_introns": n_tot,
                "multi_intron_tx": len(per_tx), "per_tx": per_tx}


def exon_precision(pred_gff: Path, ref_gtf: Path, out_dir: Path) -> tuple[float | None, Path]:
    """E_sp: exon-level precision vs RNA transcripts (pure Python).

    gffcompare is still run for auxiliary stats, but its default behaviour drops
    unstranded single-exon reference transcripts (StringTie cannot strand-assign
    single-exon transcripts without junction reads), which silently empties the
    referee set in low-splicing species. This implementation is strand-agnostic:
      - internal exons: require exact boundary match with a reference exon
      - terminal exons (first/last of transcript): >=50% reciprocal overlap
    """
    prefix = out_dir / "gffcmp"
    stats = Path(str(prefix) + ".stats")
    if not lazy([stats], desc="gffcompare"):
        run([which("gffcompare"), "-r", str(ref_gtf), "-o", str(prefix), str(pred_gff)],
            log_file=out_dir / "gffcompare.log", desc="gffcompare")

    pred = parse_gff(pred_gff)
    ref = parse_gff(ref_gtf)
    ref_exact: set[tuple[str, int, int]] = set()
    ref_ivals: dict[str, list[tuple[int, int]]] = {}
    for tx in ref.values():
        for s, e in tx.exons:
            ref_exact.add((tx.seqid, s, e))
            ref_ivals.setdefault(tx.seqid, []).append((s, e))
    for v in ref_ivals.values():
        v.sort()

    def overlapped(seqid, s, e) -> bool:
        for rs, re_ in ref_ivals.get(seqid, []):
            if rs > e:
                break
            ov = min(e, re_) - max(s, rs) + 1
            if ov > 0 and ov >= 0.5 * min(e - s + 1, re_ - rs + 1):
                return True
        return False

    n_match = n_tot = 0
    for tx in pred.values():
        exons = sorted(tx.exons)
        for i, (s, e) in enumerate(exons):
            n_tot += 1
            terminal = (i == 0 or i == len(exons) - 1)
            if terminal:
                n_match += overlapped(tx.seqid, s, e)
            else:
                n_match += (tx.seqid, s, e) in ref_exact
    prec = n_match / n_tot if n_tot else None
    log.info("Esp: %d/%d exons matched RNA reference", n_match, n_tot)
    return prec, stats


def single_exon_fraction(txs: dict[str, Transcript]) -> float:
    if not txs:
        return 0.0
    return sum(1 for tx in txs.values() if len(tx.exon_introns) == 0) / len(txs)


def compute_rs(pred_gff: Path, junctions_bed: Path, ref_gtf: Path,
               out_dir: Path, alpha: float, beta: float, gamma: float,
               gold_match: float | None = None,
               gold_weight: float = 0.0) -> dict:
    """R_s: the RNA evidence bundle. Js/Esp/Me are transcript-level RNA
    agreement; GoldMatch (full-length protein match vs the RNA-assembled
    gold set) is the protein-level confirmation of the same RNA evidence —
    it belongs here (not in H_s) so every scoring dimension maps to one
    independent evidence source. gold_weight=0 (or gold_match=None) keeps
    the legacy Js/Esp/Me-only behaviour."""
    txs = parse_gff(pred_gff)
    junc = load_junctions(junctions_bed)
    js, js_detail = junction_support(txs, junc)
    esp, stats_path = exon_precision(pred_gff, ref_gtf, out_dir)
    fse_pred = single_exon_fraction(txs)
    fse_rna = single_exon_fraction(parse_gff(ref_gtf))
    me = None
    if fse_rna > 0:
        me = max(0.0, 1.0 - abs(fse_pred - fse_rna) / fse_rna)
    parts = {"Js": (js, alpha), "Esp": (esp, beta), "Me": (me, gamma)}
    if gold_match is not None and gold_weight > 0:
        parts["GoldMatch"] = (gold_match, gold_weight)
    avail = [(v, w) for v, w in parts.values() if v is not None and w > 0]
    rs = None
    if avail:
        wsum = sum(w for _, w in avail)
        rs = sum(v * w for v, w in avail) / wsum
    return {
        "R_s": rs, "Js": js, "Esp": esp, "Me": me, "GoldMatch": gold_match,
        "fse_pred": round(fse_pred, 4), "fse_rna": round(fse_rna, 4),
        "Js_detail": {k: v for k, v in js_detail.items() if k != "per_tx"},
        "gffcompare_stats": str(stats_path),
    }


# ---------------------------------------------------------------- H_s''

def extract_referee_faa(homologs: Path, referee_ids: Path, out_faa: Path) -> Path:
    """Subset homolog FASTA to referee ids."""
    if lazy([out_faa], desc="extract referee faa"):
        return out_faa
    ensure_dir(out_faa.parent)
    keep = set(l.strip() for l in open(referee_ids) if l.strip())
    n = 0
    with open(homologs) as fh, open(out_faa, "w") as out:
        write = False
        for line in fh:
            if line.startswith(">"):
                write = line[1:].split()[0] in keep
                n += write
            if write:
                out.write(line)
    log.info("referee faa: %d proteins", n)
    return out_faa


def _diamond_blastp(query: Path, db_faa: Path, out_tsv: Path, threads: int,
                    db_dir: Path) -> Path:
    if lazy([out_tsv], desc=f"diamond vs {db_faa.stem}"):
        return out_tsv
    db = db_dir / (db_faa.stem + ".dmnd")
    if not lazy([Path(str(db) + ".dmnd")], desc="diamond makedb"):
        run([which("diamond"), "makedb", "--in", str(db_faa), "-d", str(db)],
            log_file=db_dir / f"{db_faa.stem}.makedb.log", desc="diamond makedb")
    run([which("diamond"), "blastp", "-q", str(query), "-d", str(db),
         "-p", str(threads), "--outfmt", "6", "qseqid", "sseqid", "pident",
         "length", "qlen", "slen", "qstart", "qend", "sstart", "send",
         "evalue", "qcovhsp", "scovhsp", "--max-target-seqs", "10",
         "--evalue", "1e-5", "-o", str(out_tsv)],
        log_file=db_dir / f"blastp_{db_faa.stem}.log", desc=f"diamond blastp vs {db_faa.stem}")
    return out_tsv


def full_length_match_rate(blast_tsv: Path, min_cov: float = 0.9,
                           subject_weights: dict | None = None) -> tuple[float, set, dict, dict]:
    """Fraction of query proteins with >=1 hit satisfying bidirectional coverage.

    subject_weights: {subject_id: weight} — relatedness weights of homolog
    sets (0.5-2.0). When given, each matched query contributes
    weight(subject)/2 (so support from the closest possible set = 1.0 and
    the score stays in [0,1]); unmatched contribute 0.

    Returns (rate, matched_qids, hits, qweight): qweight maps each matched
    query to its weight contribution (1.0 when unweighted), so callers can
    re-score any SUBSET of matched proteins with the same weighting (used by
    the structure-aware joint H_s).
    """
    qtotal = set()
    matched = set()
    qweight: dict[str, float] = {}
    hits: dict[str, list] = {}
    with open(blast_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, s = f[0], f[1]
            qcov, scov = float(f[11]) / 100.0, float(f[12]) / 100.0
            qtotal.add(q)
            hits.setdefault(q, []).append({
                "subject": s, "pident": float(f[2]),
                "qstart": int(f[6]), "qend": int(f[7]),
                "sstart": int(f[8]), "send": int(f[9]),
                "qlen": int(f[4]), "slen": int(f[5]),
                "qcov": qcov, "scov": scov,
            })
            if qcov >= min_cov and scov >= min_cov and q not in matched:
                matched.add(q)
                qweight[q] = (subject_weights.get(s, 1.0) / 2.0
                              if subject_weights else 1.0)
    if not qtotal:
        rate = None
    else:
        rate = sum(qweight.values()) / len(qtotal)
    return rate, matched, hits, qweight


def compute_hs(pred_faa: Path, referee_faa: Path, gold_faa: Path | None,
               out_dir: Path, threads: int, full_cov: float,
               subject_weights: dict | None = None) -> dict:
    """Full-length match rates vs the referee set AND vs the RNA-assembled
    gold set, returned separately (Hs_homolog / Hs_gold). The two rates are
    NOT blended here anymore: goldset is RNA-derived evidence, so the caller
    (score_one_predictor) routes Hs_gold to the R_s bundle, or back into H_s
    only as a rescue when the referee set is too small to be informative
    (--hs-referee-min) — never both, to avoid double-counting.

    _matched/_qweight expose the matched ids + per-protein weight so the
    caller can compute the structure-aware joint rate (matched AND clean)."""
    hs_dir = ensure_dir(out_dir)
    out = {}
    tsv = _diamond_blastp(pred_faa, referee_faa, hs_dir / "pred_vs_referee.tsv",
                          threads, hs_dir)
    hs_hom, matched, hits, qweight = full_length_match_rate(
        tsv, full_cov, subject_weights)
    out["Hs_homolog"] = hs_hom
    out["n_fullmatch_referee"] = len(matched)
    out["_matched"] = matched
    out["_qweight"] = qweight

    hs_gold = None
    if gold_faa and Path(gold_faa).exists():
        tsv_g = _diamond_blastp(pred_faa, gold_faa, hs_dir / "pred_vs_goldset.tsv",
                                threads, hs_dir)
        hs_gold, matched_g, _, _ = full_length_match_rate(tsv_g, full_cov)
        out["Hs_gold"] = hs_gold
        out["n_fullmatch_gold"] = len(matched_g)

    out["_hits_cache"] = tsv  # reused by the structural-anomaly scan
    return out


# ---------------------------------------------------------------- B_s

def compute_bc(pred_faa: Path, busco_path: Path, out_dir: Path, threads: int) -> dict:
    """B_s: BUSCO proteins-mode complete recall, local lineage, offline."""
    bdir = ensure_dir(out_dir)
    summary_candidates = list(bdir.glob("busco_out/short_summary*.txt"))
    if not summary_candidates:
        run(f"cd {bdir} && PYTHONPATH= {which('busco')} -i {Path(pred_faa).resolve()} -m proteins "
            f"-l {Path(busco_path).resolve()} -o busco_out --offline -c {threads} -f",
            shell=True, log_file=bdir / "busco.log", desc="busco (proteins)")
        summary_candidates = list(bdir.glob("busco_out/short_summary*.txt"))
    if not summary_candidates:
        raise RuntimeError(f"BUSCO summary not found under {bdir}")
    txt = summary_candidates[0].read_text()
    m = re.search(r"C:([\d.]+)%", summary_candidates[0].read_text())
    bc = float(m.group(1)) / 100.0 if m else None
    out = {"B_s": bc, "busco_summary": str(summary_candidates[0])}
    # informative coverage: fraction of predicted proteins assigned to any
    # BUSCO group (non-Missing rows of full_table.tsv)
    tables = sorted(bdir.glob("busco_out/run_*/full_table.tsv"))
    if tables:
        assessed: set = set()
        with open(tables[0]) as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                f = ln.rstrip("\n").split("\t")
                if len(f) >= 3 and f[1] != "Missing" and f[2]:
                    assessed.add(f[2])
        n_pred = sum(1 for l in open(pred_faa) if l.startswith(">"))
        out["n_busco_assessed"] = len(assessed)
        out["B_s_informative"] = len(assessed) / n_pred if n_pred else None
    return out
def build_isoform_file(gff: Path, out_path: Path) -> Path | None:
    """Build the omark -i isoform file from a predictor GFF:
    one line per gene, ';'-separated ids of ALL its isoform proteins
    (gffread -y protein ids == transcript ids). omark then picks the
    best-scoring isoform per gene instead of counting isoforms as
    duplicates. Returns None when no gene/mRNA structure is found."""
    if lazy([out_path], desc="isoform file"):
        return out_path
    genes: dict[str, list] = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("mRNA", "transcript"):
                continue
            attrs = f[8]
            m_id = re.search(r"ID=([^;]+)", attrs)
            if m_id:
                tid = m_id.group(1)
                m_par = re.search(r"Parent=([^;,]+)", attrs)
                gid = m_par.group(1) if m_par else tid
            else:
                m_tid = re.search(r'transcript_id "([^"]+)"', attrs)
                m_gid = re.search(r'gene_id "([^"]+)"', attrs)
                if not m_tid:
                    continue
                tid = m_tid.group(1)
                gid = m_gid.group(1) if m_gid else tid
            tids = genes.setdefault(gid, [])
            if tid not in tids:
                tids.append(tid)
    if not genes:
        return None
    ensure_dir(out_path.parent)
    with open(out_path, "w") as out:
        for gid, tids in genes.items():
            out.write(";".join(tids) + "\n")
    return out_path


# ---------------------------------------------------------------- O_s
def compute_os(pred_faa: Path, omark_db: Path, out_dir: Path, threads: int,
               isoform_file: Path | None = None,
               extra_opts: str | None = None) -> dict:
    """O_s: OMArk Consistent fraction (A: line of the .sum), 0..1.

    Runs omamer search + omark per predictor protein set; cached under
    out_dir (a *.sum file reuses the previous result; changing
    isoform_file/extra_opts invalidates the cache).
    """
    from .omark import run_omark
    ctx = run_omark(Path(pred_faa), Path(omark_db), ensure_dir(out_dir),
                    threads, isoform_file=isoform_file,
                    extra_opts=extra_opts)
    if not ctx.get("ok"):
        raise RuntimeError(f"OMArk failed for {pred_faa}: {ctx.get('error')}")
    pct = ctx.get("consistent_pct")
    os_ = float(pct) / 100.0 if pct is not None else None
    out = {"O_s": os_, "omark_consistent_pct": pct,
           "omark_contamination_pct": ctx.get("contamination_pct"),
           "omark_inconsistent": ctx.get("inconsistent"),
           "omark_summary": ctx.get("summary_file")}
    # informative coverage: Consistent + Inconsistent + Contamination all
    # mean the database gave a verdict (Unknown = no information)
    verdicts = [ctx.get("consistent"), ctx.get("inconsistent"),
                ctx.get("contamination")]
    unknown = ctx.get("unknown")
    if all(isinstance(x, int) for x in verdicts) and isinstance(unknown, int):
        total = sum(verdicts) + unknown
        out["O_s_informative"] = sum(verdicts) / total if total else None
    return out


# ---------------------------------------------------------------- C_s

def compute_cs(pred_faa: Path, blast_tsv: Path, referee_loci: dict[str, tuple],
               min_seg_cov: float = 0.25, locus_gap: int = 1_000_000,
               subject_weights: dict | None = None) -> dict:
    """C_s = 1 - r_anomaly.  Anomalies from homolog-hit patterns:
      chimeric: one predicted protein hits >=2 referee proteins whose genomic
                loci are non-adjacent, each covering a distinct query segment
      split   : >=2 predicted proteins cover non-overlapping segments of the
                same referee protein
    blast_tsv: pred_vs_referee.tsv (outfmt with qstart/qend/sstart/send/qlen/slen).
    referee_loci: {referee_id: (seqid, start, end)} from homolog alignments.
    subject_weights: relatedness weights of homolog sets; an anomaly supported
    only by distant sets counts less (weight/2 normalization, in [0,1]).
    """
    def sw(s: str) -> float:
        return (subject_weights.get(s, 1.0) / 2.0) if subject_weights else 1.0
    by_query: dict[str, dict[str, dict]] = {}
    by_subject: dict[str, dict[str, dict]] = {}
    n_pred = 0
    with open(pred_faa) as fh:
        n_pred = sum(1 for l in fh if l.startswith(">"))
    with open(blast_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, s = f[0], f[1]
            rec = {
                "qs": min(int(f[6]), int(f[7])), "qe": max(int(f[6]), int(f[7])),
                "ss": min(int(f[8]), int(f[9])), "se": max(int(f[8]), int(f[9])),
                "qlen": int(f[4]), "slen": int(f[5]), "pident": float(f[2]),
            }
            cur = by_query.setdefault(q, {}).get(s)
            if cur is None or (rec["qe"] - rec["qs"]) > (cur["qe"] - cur["qs"]):
                by_query[q][s] = rec
            cur2 = by_subject.setdefault(s, {}).get(q)
            if cur2 is None or (rec["se"] - rec["ss"]) > (cur2["se"] - cur2["ss"]):
                by_subject[s][q] = rec

    def overlap(a0, a1, b0, b1):
        return max(0, min(a1, b1) - max(a0, b0))

    chimeric: dict[str, float] = {}
    split_genes: dict[str, float] = {}
    for q, subs in by_query.items():
        # distinct subjects covering distinct query segments
        good = [(s, r) for s, r in subs.items()
                if (r["qe"] - r["qs"]) / r["qlen"] >= min_seg_cov and s in referee_loci]
        for i in range(len(good)):
            for j in range(i + 1, len(good)):
                (s1, r1), (s2, r2) = good[i], good[j]
                seg_ov = overlap(r1["qs"], r1["qe"], r2["qs"], r2["qe"])
                if seg_ov > 0.3 * min(r1["qe"] - r1["qs"], r2["qe"] - r2["qs"]):
                    continue  # same segment -> alternative homolog, not chimera
                l1, l2 = referee_loci[s1], referee_loci[s2]
                if l1[0] != l2[0] or abs(l1[1] - l2[1]) > locus_gap:
                    w = min(sw(s1), sw(s2))
                    if w > chimeric.get(q, 0.0):
                        chimeric[q] = w
    for s, qs in by_subject.items():
        good = [(q, r) for q, r in qs.items()
                if (r["se"] - r["ss"]) / r["slen"] >= min_seg_cov]
        for i in range(len(good)):
            for j in range(i + 1, len(good)):
                (q1, r1), (q2, r2) = good[i], good[j]
                seg_ov = overlap(r1["ss"], r1["se"], r2["ss"], r2["se"])
                if seg_ov < 0.1 * max(r1["se"] - r1["ss"], r2["se"] - r2["ss"]):
                    for q in (q1, q2):
                        if sw(s) > split_genes.get(q, 0.0):
                            split_genes[q] = sw(s)

    anomalous = dict(chimeric)
    for q, w in split_genes.items():
        if w > anomalous.get(q, 0.0):
            anomalous[q] = w
    if subject_weights:
        r_anom = min(1.0, sum(anomalous.values()) / n_pred) if n_pred else None
    else:
        r_anom = len(anomalous) / n_pred if n_pred else None
    return {
        "C_s": (1 - r_anom) if r_anom is not None else None,
        "r_anomaly": r_anom,
        "n_chimeric": len(chimeric), "n_split": len(split_genes),
        "chimeric_ids": sorted(chimeric), "split_ids": sorted(split_genes),
    }


def load_referee_loci(alignments_tsv: Path) -> dict[str, tuple]:
    loci = {}
    with open(alignments_tsv) as fh:
        header = next(fh).rstrip("\n").split("\t")
        # column lookup by name (robust to schema additions like 'set')
        ci = {name: i for i, name in enumerate(header)}
        i_t, i_s, i_e = ci.get("target", 7), ci.get("tstart", 8), ci.get("tend", 9)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_t, i_s, i_e):
                continue
            loci[f[0]] = (f[i_t], int(f[i_s]), int(f[i_e]))
    return loci


# ---------------------------------------------------------------- aggregation

def aggregate(subscores: dict[str, float | None],
              weights: dict[str, float],
              coverage: dict[str, float] | None = None) -> tuple[float | None, dict[str, float]]:
    """Weighted sum with NA renormalization. Refuse if <2 sub-blocks available.

    coverage (optional, adaptive mode): per-sub-block informative coverage in
    [0,1] — the fraction of predicted proteins the block actually gave a
    usable verdict on (e.g. OMArk non-Unknown fraction, BUSCO proteins hit).
    Effective weight = w * coverage, then renormalized over available blocks.
    """
    avail = {k: (v, weights[k]) for k, v in subscores.items()
             if v is not None and weights.get(k, 0) > 0}
    eff = {k: 0.0 for k in weights}
    if len(avail) < 2:
        return None, eff
    cov = coverage or {}
    wsum = sum(w * cov.get(k, 1.0) for k, (_, w) in avail.items())
    for k, (v, w) in avail.items():
        eff[k] = w * cov.get(k, 1.0) / wsum
    return sum(v * eff[k] for k, (v, _) in avail.items()), eff


def map_weights(scores: dict[str, float | None], mapping: str = "auto",
                min_score: float = 0.35) -> dict[str, int | None]:
    """Map P (prediction score) to EVM prediction weights 1..5.
    >=2 predictors: relative scaling (linear: min-max; rank: by rank order).
    auto (default): rank when >=3 qualified predictors, linear otherwise —
    rank spreads ANY number of inputs across the full 1..5 weight space,
    so adding predictors never compresses everyone into the same bucket.
    Single predictor: absolute mapping round(1 + 4*W). Below min_score -> None (excluded).
    """
    valid = {k: v for k, v in scores.items() if v is not None}
    out: dict[str, int | None] = {k: None for k in scores}
    if not valid:
        return out
    for k, v in valid.items():
        if v < min_score:
            out[k] = None  # hard drop
    good = {k: v for k, v in valid.items() if v >= min_score}
    if not good:
        return out
    if len(good) == 1:
        for k, v in good.items():
            out[k] = max(1, min(5, round(1 + 4 * v)))
        return out
    if mapping == "auto":
        mapping = "rank" if len(good) >= 3 else "linear"
    vals = sorted(good.values())
    if mapping == "rank":
        ranks = {k: sorted(vals).index(v) for k, v in good.items()}
        n = len(good) - 1
        for k, v in good.items():
            out[k] = max(1, min(5, round(1 + 4 * ranks_frac(ranks[k], n))))
    else:  # linear min-max
        lo, hi = vals[0], vals[-1]
        for k, v in good.items():
            frac = (v - lo) / (hi - lo) if hi > lo else 1.0
            out[k] = max(1, min(5, round(1 + 4 * frac)))
    return out


def ranks_frac(r: int, n: int) -> float:
    return r / n if n > 0 else 1.0


# ---------------------------------------------------------------- driver

def extract_proteins(genome: Path, gff: Path, out_faa: Path) -> Path:
    """gffread -y: translate CDS to proteins."""
    if not lazy([out_faa], desc=f"gffread -y {gff.name}"):
        run([which("gffread"), "-g", str(genome), "-y", str(out_faa), str(gff)],
            desc=f"gffread proteins {gff.name}")
    return out_faa


def score_one_predictor(name: str, gff: Path, vi, out_root: Path, threads: int,
                        evidence: dict, cfg: dict) -> dict:
    """Score a single predictor across all available sub-blocks."""
    pdir = ensure_dir(out_root / name)
    log.info("module 3 [%s]: scoring %s", name, gff)
    ablate = set(cfg.get("ablate", []))

    pred_faa = extract_proteins(vi.genome, gff, pdir / "proteins.faa")

    detail: dict = {"predictor": name, "gff": str(gff)}
    sub: dict[str, float | None] = {}

    # relatedness weights of homolog sets: protein id -> weight (multi-set
    # homolog input only; single set stays unweighted legacy behavior)
    hom_ev = evidence.get("homolog") or {}
    subject_weights = cfg.get("subject_weights")
    if subject_weights is None and len(hom_ev.get("set_weights") or {}) > 1:
        mem = hom_ev.get("membership") or {}
        subject_weights = {p: hom_ev["set_weights"].get(lab, 1.0)
                           for p, lab in mem.items()}
    if subject_weights:
        detail["H_s_weighted"] = True

    # ---- H_s'' (structure-aware: full-length AND anomaly-free; pure
    # homolog-protein evidence. The gold set is RNA-derived, so its match
    # rate is routed to the R_s bundle — EXCEPT when the referee set is too
    # small to be informative, in which case it rescues H_s instead.
    # Never both: goldset is counted in exactly one dimension.)
    goldmatch_for_rs = None
    if "homolog" in ablate:
        sub["H_s"] = None
        detail["H_s_note"] = "ablated"
    else:
        referee_faa = extract_referee_faa(vi.homologs, evidence["homolog"]["referee_ids"],
                                          out_root / "_shared" / "referee.faa")
        n_referee = sum(1 for ln in open(referee_faa) if ln.startswith(">"))
        ref_min = cfg.get("hs_referee_min", 1000)
        rescue = n_referee < ref_min
        detail["n_referee"] = n_referee
        gold_faa = evidence.get("goldset", {}).get("goldset_faa") if evidence.get("goldset") else None
        h = compute_hs(pred_faa, referee_faa, gold_faa, pdir / "hs", threads,
                       cfg["fulllength_cov"], subject_weights=subject_weights)
        detail.update({k: v for k, v in h.items() if not k.startswith("_")})
        detail["_hits_cache"] = h["_hits_cache"]
        # joint criterion: from the full-length matched proteins, drop the
        # structurally anomalous ones (chimeric / split) and recompute the
        # (relatedness-weighted) match rate on the remaining clean set
        if "structure" in ablate:
            detail["H_s_note"] = "structure screen ablated: plain full-length rate"
            hs_hom_joint = h["Hs_homolog"]
        else:
            loci = load_referee_loci(evidence["homolog"]["alignments"])
            c = compute_cs(pred_faa, h["_hits_cache"], loci,
                           subject_weights=subject_weights)
            detail.update({k: v for k, v in c.items()
                           if not isinstance(v, list) and k != "C_s"})
            detail["anomalous_ids"] = sorted(set(c["chimeric_ids"]) | set(c["split_ids"]))
            anomalous = set(c["chimeric_ids"]) | set(c["split_ids"])
            n_pred = sum(1 for l in open(pred_faa) if l.startswith(">"))
            clean = h["_matched"] - anomalous
            if not n_pred:
                hs_hom_joint = None
            else:
                hs_hom_joint = sum(h["_qweight"].get(q, 1.0) for q in clean) / n_pred
            detail["Hs_homolog_joint"] = hs_hom_joint
            detail["n_fullmatch_clean"] = len(clean)
        hs_gold = h.get("Hs_gold")
        if rescue and hs_gold is not None:
            gw = cfg["hs_gold_weight"]
            if hs_hom_joint is not None:
                sub["H_s"] = (1 - gw) * hs_hom_joint + gw * hs_gold
            else:
                sub["H_s"] = hs_gold
            detail["H_s_note"] = (f"goldset rescue: referee set too small "
                                  f"({n_referee} < {ref_min} proteins)")
            detail["goldset_routed_to"] = "H_s"
        else:
            sub["H_s"] = hs_hom_joint
            if hs_gold is not None:
                goldmatch_for_rs = hs_gold
                detail["goldset_routed_to"] = "R_s"

    # ---- R_s (RNA evidence bundle: junction support + exon precision +
    # single-exon bias + gold-set protein match, all RNA-derived)
    if "rna" in ablate or not evidence.get("rna"):
        sub["R_s"] = None
        detail["R_s_note"] = "ablated" if "rna" in ablate else "no RNA evidence"
    else:
        r = compute_rs(gff, evidence["rna"]["junctions_bed"],
                       evidence["rna"]["stringtie_gtf"], pdir,
                       cfg["rna_alpha"], cfg["rna_beta"], cfg["rna_gamma"],
                       gold_match=goldmatch_for_rs,
                       gold_weight=cfg.get("rna_gold_weight", 0.25))
        detail.update(r)
        sub["R_s"] = r["R_s"]

    # ---- B_s
    if "conserve" in ablate or not vi.busco:
        sub["B_s"] = None
        detail["B_s_note"] = "ablated" if "conserve" in ablate else "no conservation set"
    else:
        b = compute_bc(pred_faa, vi.busco, pdir / "busco", threads)
        detail.update(b)
        sub["B_s"] = b["B_s"]

    # ---- O_s
    if "omark" in ablate or not cfg.get("omark_db"):
        sub["O_s"] = None
        detail["O_s_note"] = "ablated" if "omark" in ablate else "no --omark-db"
    else:
        try:
            iso_file = build_isoform_file(gff, pdir / "isoforms.txt")
            o = compute_os(pred_faa, Path(cfg["omark_db"]), pdir / "omark",
                           threads, isoform_file=iso_file,
                           extra_opts=cfg.get("omark_options"))
            detail.update(o)
            sub["O_s"] = o["O_s"]
        except Exception as e:
            log.warning("module 3 [%s]: O_s unavailable: %s", name, e)
            sub["O_s"] = None
            detail["O_s_note"] = f"omark failed: {e}"

    weights = {"R_s": cfg["w_rna"], "H_s": cfg["w_homolog"],
               "B_s": cfg["w_conserve"], "O_s": cfg.get("w_omark", 0.0)}
    coverage = None
    if cfg.get("weight_mode") == "adaptive":
        # adaptive: scale each database-driven sub-block by its informative
        # coverage (fraction of proteins with a usable verdict); R_s/H_s
        # are global rate metrics over the whole prediction -> 1.0
        coverage = {"R_s": 1.0, "H_s": 1.0,
                    "B_s": detail.get("B_s_informative") or 0.0,
                    "O_s": detail.get("O_s_informative") or 0.0}
        detail["informative_coverage"] = coverage
    w, eff = aggregate(sub, weights, coverage)
    detail["P"] = w
    detail["subscores"] = sub
    detail["effective_weights"] = eff
    # structural anomalies (r_anomaly / n_chimeric / n_split) are kept as
    # informational metrics only: EVM entry is gated solely by P via
    # map_weights(min_score); the <2-usable-sub-blocks guard below is the
    # only hard reject (a score cannot be computed at all then)
    if w is None:
        detail["hard_reject"] = "fewer than 2 usable sub-blocks"

    import json
    with open(pdir / "score_detail.json", "w") as fh:
        import json as _json
        _json.dump({k: v for k, v in detail.items() if not k.startswith("_")},
                   fh, indent=2, ensure_ascii=False, default=str)
    log.info("module 3 [%s] P=%s subscores=%s", name,
             f"{w:.4f}" if w is not None else "NA",
             {k: (f"{v:.3f}" if v is not None else "NA") for k, v in sub.items()})
    return detail
