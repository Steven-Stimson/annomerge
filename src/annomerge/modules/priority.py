"""Module 4/5 (priority mode): goldset-verbatim + de novo gap-filling merge.

Design doc v2.2 §1.2 — priority merge (Paludamentum-style, single-de novo
predictor scenario). Instead of EVM weighted consensus:
  1. goldset (lifted TransDecoder ORFs) kept VERBATIM.
  2. de novo transcripts overlapping a goldset locus (same seqid+strand, any
     exon-vs-locus-span overlap) are DISCARDED — goldset wins its loci.
  3. only de novo at goldset-free loci are kept (gap filling).
  4. cross-predictor structural dedup: identical CDS or exon set collapses
     to the highest-scored predictor's copy.
  5. locus arbitration: a locus contested by several predictors (different
     structures, same seqid+strand, overlapping spans) is OWNED by the
     highest-scoring predictor (W_raw, tie-break R_s); its transcripts there
     are kept, lower scorers' are dropped. Exact (W_raw, R_s) ties keep every
     tied predictor's transcripts and are flagged in the stats.
Output is a gene/mRNA/exon/CDS GFF3 in the same shape as the EVM output, so
module 5b (sanity filter) and module 6 (confidence tiering) run unchanged.

Overlap semantics (documented): a goldset *locus* is its gene span [start,end]
on (seqid, strand). A de novo transcript is dropped iff any of its exons (or
CDS, if it has no exons) overlaps that span. This is the literal "任意重叠"
reading of the design; nested genes inside a goldset intron are therefore also
dropped (goldset owns the whole locus footprint).
"""
from __future__ import annotations

import bisect
import logging
from pathlib import Path

from ..utils import ensure_dir, lazy
from .evm import lift_goldset
from .sanity import _parse_gff3

log = logging.getLogger("annomerge")


# ---------------------------------------------------------------- interval idx
def _locus_index(genes: dict) -> dict:
    """{gid: Feat} -> {(seqid, strand): (sorted_starts, prefix_max_end)}."""
    by_strand: dict[tuple[str, str], list] = {}
    for g in genes.values():
        by_strand.setdefault((g.seqid, g.strand), []).append((g.start, g.end))
    index = {}
    for key, ivals in by_strand.items():
        ivals.sort()
        starts = [s for s, _ in ivals]
        pmax, m = [], 0
        for _, e in ivals:
            m = max(m, e)
            pmax.append(m)
        index[key] = (starts, pmax)
    return index


def _overlaps(index: dict, seqid: str, strand: str, s: int, e: int) -> bool:
    """True iff [s,e] overlaps any locus span on (seqid, strand)."""
    key = (seqid, strand)
    if key not in index:
        return False
    starts, pmax = index[key]
    k = bisect.bisect_right(starts, e) - 1   # last locus with start <= e
    return k >= 0 and pmax[k] >= s


# ---------------------------------------------------------------- merge
def _sc(scores: dict | None, name: str) -> float:
    """Predictor W_raw (0.0 if unknown)."""
    if not scores:
        return 0.0
    v = scores.get(name, (0.0, 0.0))[0]
    return v if isinstance(v, (int, float)) else 0.0


def _rs(scores: dict | None, name: str) -> float:
    """Predictor R_s, tie-breaker after W_raw (0.0 if unknown)."""
    if not scores:
        return 0.0
    v = scores.get(name, (0.0, 0.0))[1]
    return v if isinstance(v, (int, float)) else 0.0


def _cds_exon_keys(tx):
    ck = (tx.mrna.seqid, tx.mrna.strand,
          tuple(sorted((c.start, c.end) for c in tx.cds)))
    ek = (tx.mrna.seqid, tx.mrna.strand,
          tuple(sorted((x.start, x.end) for x in tx.exons)))
    return ck, ek


def _arbitrate_loci(records: list) -> tuple:
    """Group records into overlap-connected loci per (seqid, strand); the
    highest-scoring predictor OWNS each locus (all its tx there are kept).
    Different-structure-same-locus tx from lower-scoring predictors are
    dropped. Tie at (W_raw, R_s) -> keep every tied predictor's tx, flagged.

    Returns (kept_records, stats). Local evidence override (design doc d):
    intentionally not implemented yet — global score is the only arbiter.
    """
    by_strand: dict = {}
    for r in records:
        foot = r["tx"].exons if r["tx"].exons else r["tx"].cds
        s = min(f.start for f in foot)
        e = max(f.end for f in foot)
        by_strand.setdefault((r["tx"].mrna.seqid, r["tx"].mrna.strand), []
                             ).append((s, e, r))
    kept: list = []
    stats = {"loci_total": 0, "loci_single_predictor": 0,
             "loci_arbitrated": 0, "loci_arbitrated_drop": 0,
             "loci_tie_kept_both": 0}
    for key in by_strand:
        items = sorted(by_strand[key], key=lambda x: (x[0], x[1]))
        # connected components over span overlap (chained loci merge)
        groups, cur, cur_end = [], [], None
        for s, e, r in items:
            if cur and s > cur_end:
                groups.append(cur)
                cur = []
            cur.append((s, e, r))
            cur_end = e if cur_end is None else max(cur_end, e)
        if cur:
            groups.append(cur)
        for g in groups:
            stats["loci_total"] += 1
            names_in = {x[2]["name"] for x in g}
            if len(names_in) == 1:
                stats["loci_single_predictor"] += 1
                kept.extend(x[2] for x in g)
                continue
            top = max((x[2]["score"], x[2]["rs"]) for x in g)
            winners = {x[2]["name"] for x in g
                       if (x[2]["score"], x[2]["rs"]) == top}
            if len(winners) == 1:
                stats["loci_arbitrated"] += 1
                w = winners.pop()
                kept.extend(x[2] for x in g if x[2]["name"] == w)
                stats["loci_arbitrated_drop"] += sum(1 for x in g
                                                     if x[2]["name"] != w)
            else:
                stats["loci_tie_kept_both"] += 1
                kept.extend(x[2] for x in g if x[2]["name"] in winners)
    return kept, stats


def run_priority_merge(vi, out_dir: Path, preds: dict,
                       stringtie_gtf: Path, goldset_dir: Path,
                       scores: dict | None = None) -> dict:
    """preds: {name: gff3_path}; scores: {name: (W_raw, R_s)} or None.
    Returns {merged_gff3, goldset_gff, dir, stats}.

    goldset transcripts are kept verbatim. de novo transcripts at goldset-free
    loci go through: (1) cross-predictor structural dedup — identical CDS or
    exon set collapses to the highest-scored copy; (2) per-locus arbitration —
    at a locus contested by >1 predictor, the highest-scoring predictor owns
    the locus and its transcripts there are kept (different structures from
    lower scorers are dropped; exact (W_raw, R_s) ties keep all tied, flagged).
    """
    d = ensure_dir(Path(out_dir))
    log.info("module 4/5 [priority]: %d predictors + goldset", len(preds))

    # 1. lift goldset to genome coords (same artifact EVM mode builds)
    goldset_gff = d / "goldset.genome.gff3"
    if goldset_dir and (goldset_dir / "transcripts.fa.transdecoder.gff3").exists():
        lift_goldset(goldset_dir, stringtie_gtf, d, goldset_gff,
                     genome_fa=vi.genome)
    else:
        goldset_gff = None
        log.warning("priority merge: no goldset; de novo kept as-is")

    g_genes, g_txs, _ = _parse_goldset(goldset_gff) if goldset_gff else ({}, {}, [])
    locus_idx = _locus_index(g_genes)

    # 2. per predictor: drop goldset-overlapping tx; collect the rest
    collected: list = []   # records: {name, tx, gid, score, rs}
    stats = {"goldset_genes": len(g_genes), "goldset_txs": len(g_txs)}
    for name, gff in preds.items():
        _, txs, order = _parse_gff3(Path(gff))
        n_overlap = 0
        for tid in order:
            tx = txs[tid]
            foot = tx.exons if tx.exons else tx.cds
            if not foot:                       # nothing to place on the genome
                continue
            if any(_overlaps(locus_idx, f.seqid, tx.mrna.strand,
                             f.start, f.end) for f in foot):
                n_overlap += 1
                continue
            collected.append({"name": name, "tx": tx, "gid": tx.gid,
                              "score": _sc(scores, name),
                              "rs": _rs(scores, name)})
        stats[f"{name}_total"] = len(order)
        stats[f"{name}_goldset_overlap_drop"] = n_overlap

    # 3. cross-predictor structural dedup: identical CDS/exon set collapses
    #    to the highest-scored copy (same gene model — keep the better source)
    collected.sort(key=lambda r: (-r["score"], -r["rs"]))
    seen_cds, seen_exon, n_dedup = set(), set(), 0
    after_dedup: list = []
    for r in collected:
        ck, ek = _cds_exon_keys(r["tx"])
        if (r["tx"].cds and ck in seen_cds) or (r["tx"].exons and ek in seen_exon):
            n_dedup += 1
            continue
        if r["tx"].cds:
            seen_cds.add(ck)
        if r["tx"].exons:
            seen_exon.add(ek)
        after_dedup.append(r)
    stats["cross_pred_dedup_drop"] = n_dedup

    # 4. locus arbitration (different structures, same locus -> best predictor)
    kept_records, arb_stats = _arbitrate_loci(after_dedup)
    stats.update(arb_stats)
    kept_txs = [(r["name"], r["tx"], r["gid"]) for r in kept_records]
    for name in preds:
        stats[f"{name}_kept"] = sum(1 for r in kept_records if r["name"] == name)
        log.info("priority [%s]: total=%d goldset-overlap=%d kept=%d",
                 name, stats[f"{name}_total"],
                 stats[f"{name}_goldset_overlap_drop"], stats[f"{name}_kept"])

    # 3. write merged GFF3: goldset verbatim + de novo regrouped per gene
    merged = d / "priority.merged.gff3"
    if not lazy([merged], desc="priority merge"):
        _write_merged(merged, goldset_gff, g_genes, g_txs, kept_txs)

    # 4. stats table
    stsv = d / "priority_stats.tsv"
    with open(stsv, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in stats.items():
            fh.write(f"{k}\t{v}\n")
    log.info("priority merge stats -> %s", stsv)
    return {"merged_gff3": merged, "goldset_gff": goldset_gff,
            "dir": d, "stats": stats}


def _parse_goldset(gff: Path):
    """goldset genome gff3 -> (genes, txs, order) with goldset-only filter."""
    genes, txs, order = _parse_gff3(Path(gff))
    genes = {g: f for g, f in genes.items() if f.source == "goldset"}
    txs = {t: x for t, x in txs.items() if x.mrna.source == "goldset"}
    order = [t for t in order if t in txs]
    return genes, txs, order


def _write_merged(merged: Path, goldset_gff: Path | None,
                  g_genes: dict, g_txs: dict, kept_txs: list) -> None:
    with open(merged, "w") as out:
        out.write("##gff-version 3\n")
        # goldset verbatim (gene/mRNA/exon/CDS lines as lifted)
        if goldset_gff:
            for line in open(goldset_gff):
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 9 and f[1] == "goldset" and f[2] in \
                        ("gene", "mRNA", "transcript", "exon", "CDS"):
                    out.write(line)
        # de novo regrouped per original gene, re-IDed with predictor prefix
        gene_order: list = []
        gene_txs: dict = {}
        for name, tx, ogid in kept_txs:
            if ogid not in gene_txs:
                gene_order.append((name, ogid))
                gene_txs[ogid] = []
            gene_txs[ogid].append((name, tx))
        n = 0
        for name, ogid in gene_order:
            n += 1
            gid = f"{name}.g{n}"
            txs_g = gene_txs[ogid]
            feats = [f for _, tx in txs_g for f in (tx.exons or tx.cds)]
            gs = min(f.start for f in feats)
            ge = max(f.end for f in feats)
            gstrand = txs_g[0][1].mrna.strand
            gseqid = txs_g[0][1].mrna.seqid
            out.write(f"{gseqid}\t{name}\tgene\t{gs}\t{ge}\t.\t{gstrand}\t.\t"
                      f"ID={gid};Name={ogid}\n")
            for i, (_, tx) in enumerate(txs_g, 1):
                tid = f"{gid}.t{i}"
                out.write(f"{tx.mrna.seqid}\t{name}\tmRNA\t{tx.mrna.start}\t"
                          f"{tx.mrna.end}\t.\t{tx.mrna.strand}\t.\t"
                          f"ID={tid};Parent={gid}\n")
                for j, x in enumerate(sorted(tx.exons, key=lambda f: f.start), 1):
                    out.write(f"{x.seqid}\t{name}\texon\t{x.start}\t{x.end}\t.\t"
                              f"{x.strand}\t.\tID={tid}.exon{j};Parent={tid}\n")
                for c in tx.cds_tx_order():
                    out.write(f"{c.seqid}\t{name}\tCDS\t{c.start}\t{c.end}\t.\t"
                              f"{c.strand}\t{c.phase}\tID={tid}.cds;Parent={tid}\n")
