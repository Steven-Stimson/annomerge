"""Module 5b: sanity filter of the integrated annotation (post-EVM first pass).

Port of Paludamentum's sanity_filter_gff3.py rules (design doc v2.2 §2.2):

  remove (hard errors):
    no_cds            transcript without CDS
    cds_length_mod3   CDS length (less phase of first CDS) not a multiple of 3
    cds_structure     CDS segments overlap / on different strands or seqids
    internal_stop     stop codon before the last codon
  repair (instead of remove):
    cds_outside_exon  CDS segment outside its transcript's exons: exons are
                      extended (or a new exon created) to cover the CDS, then
                      overlapping exons are merged and renumbered. Removal is
                      only the last resort if coverage still cannot be fixed.
  keep with note:
    extended_stop     CDS has no terminal stop but the next 3 bases are one:
                      extend last CDS segment (and its exon) by 3 bp
    no_stop           CDS ends without stop codon (3' partial)
    non_atg_start     first codon is not ATG

Genes whose transcripts are all removed are removed. Gene/transcript spans are
recomputed from children after edits. Output: filtered GFF3 + sanity_filter.tsv
(gene_id, transcript_id, action, reason).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

from ..utils import lazy

log = logging.getLogger("annomerge")

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class Feat:
    seqid: str
    source: str
    ftype: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attrs: str


@dataclass
class Tx:
    tid: str
    gid: str
    mrna: Feat
    exons: list = field(default_factory=list)   # list[Feat]
    cds: list = field(default_factory=list)     # list[Feat]

    def cds_tx_order(self) -> list:
        """CDS segments in transcription order (5'->3' of the transcript)."""
        segs = sorted(self.cds, key=lambda f: f.start)
        if self.mrna.strand == "-":
            segs.reverse()
        return segs


def _parse_gff3(path: Path):
    """Return (genes: {gid: Feat}, txs: {tid: Tx}, order: [tid])."""
    genes: dict[str, Feat] = {}
    txs: dict[str, Tx] = {}
    order: list[str] = []
    tx2gene: dict[str, str] = {}

    def attrs_get(a: str, key: str) -> str | None:
        for kv in a.split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                if k.strip() == key:
                    return v.strip()
        return None

    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        feat = Feat(f[0], f[1], f[2], int(f[3]), int(f[4]), f[5], f[6], f[7], f[8])
        if feat.ftype == "gene":
            gid = attrs_get(feat.attrs, "ID")
            if gid:
                genes[gid] = feat
        elif feat.ftype in ("mRNA", "transcript"):
            tid = attrs_get(feat.attrs, "ID")
            gid = attrs_get(feat.attrs, "Parent") or ""
            if tid:
                txs[tid] = Tx(tid, gid, feat)
                order.append(tid)
                tx2gene[tid] = gid
        elif feat.ftype == "exon":
            tid = attrs_get(feat.attrs, "Parent")
            if tid in txs:
                txs[tid].exons.append(feat)
        elif feat.ftype == "CDS":
            tid = attrs_get(feat.attrs, "Parent")
            if tid in txs:
                txs[tid].cds.append(feat)
    return genes, txs, order


def _read_genome(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def _subseq(genome: dict[str, str], seqid: str, s: int, e: int, strand: str) -> str:
    """1-based closed interval -> sequence (strand-aware)."""
    seq = genome.get(seqid, "")
    s = max(s, 1)
    e = min(e, len(seq))
    if s > e:
        return ""
    sub = seq[s - 1:e]
    if strand == "-":
        comp = str.maketrans("ACGTN", "TGCAN")
        sub = sub.translate(comp)[::-1]
    return sub


def _translate_cds(genome, tx: Tx) -> str:
    """Translate the CDS in transcription order, honouring the first phase."""
    segs = tx.cds_tx_order()
    if not segs:
        return ""
    first_phase = int(segs[0].phase) if segs[0].phase in "012" else 0
    nuc = "".join(
        _subseq(genome, f.seqid, f.start, f.end, tx.mrna.strand) for f in segs)
    nuc = nuc[first_phase:]
    return "".join(CODON_TABLE.get(nuc[i:i + 3], "X")
                   for i in range(0, len(nuc) - 2, 3))


def _structure_problem(tx: Tx) -> str | None:
    if not tx.cds:
        return "no_cds"
    segs = sorted(tx.cds, key=lambda f: f.start)
    first = tx.cds_tx_order()[0]
    phase = int(first.phase) if first.phase in "012" else 0
    if (sum(f.end - f.start + 1 for f in segs) - phase) % 3:
        return "cds_length_mod3"
    if len({f.seqid for f in segs}) > 1 or len({f.strand for f in segs}) > 1 \
            or segs[0].strand not in "+-" \
            or segs[0].seqid != tx.mrna.seqid or segs[0].strand != tx.mrna.strand:
        return "cds_structure"
    if any(b.start <= a.end for a, b in zip(segs, segs[1:])):
        return "cds_structure"
    # CDS must lie inside exons
    exons = sorted((f.start, f.end) for f in tx.exons)
    for f in segs:
        if not any(s <= f.start and f.end <= e for s, e in exons):
            return "cds_outside_exon"
    return None


def _repair_cds_outside_exon(tx: Tx) -> bool:
    """Extend (or create) exons so every CDS segment is covered by an exon.

    CDS-outside-exon here is almost always an artefact of the augment step's
    UTR backfill shifting exon boundaries, not a bad gene model: dropping the
    transcript (and with it, often the whole gene) would silently lose valid
    consensus models. Repair instead: extend the nearest overlapping exon to
    cover the offending CDS segment, synthesise an exon for an uncovered
    segment, then merge overlapping exons and renumber them. Returns True if
    every CDS segment ends up inside an exon.
    """
    if not tx.cds:
        return False
    for seg in tx.cds:
        overlapping = [ex for ex in tx.exons
                       if ex.seqid == seg.seqid and ex.strand == seg.strand
                       and not (ex.end < seg.start or ex.start > seg.end)]
        if overlapping:
            target = min(overlapping, key=lambda ex: (abs(ex.start - seg.start)
                                                      + abs(ex.end - seg.end)))
            target.start = min(target.start, seg.start)
            target.end = max(target.end, seg.end)
        else:
            tx.exons.append(Feat(seg.seqid, seg.source, "exon",
                                 seg.start, seg.end, seg.score, seg.strand,
                                 ".", f"Parent={tx.tid}"))
    # merge exons that overlap or are directly adjacent
    tx.exons.sort(key=lambda f: (f.seqid, f.start))
    merged: list[Feat] = []
    for ex in tx.exons:
        if (merged and merged[-1].seqid == ex.seqid
                and ex.start <= merged[-1].end + 1):
            merged[-1].start = min(merged[-1].start, ex.start)
            merged[-1].end = max(merged[-1].end, ex.end)
        else:
            merged.append(ex)
    tx.exons = merged
    # renumber exon IDs in coordinate order
    for i, ex in enumerate(tx.exons, 1):
        rest = [kv for kv in ex.attrs.split(";")
                if not kv.startswith("ID=") and kv]
        ex.attrs = f"ID={tx.tid}.exon{i}" + (";" + ";".join(rest) if rest else "")
    # verify coverage
    for seg in tx.cds:
        if not any(ex.seqid == seg.seqid and ex.start <= seg.start
                   and seg.end <= ex.end for ex in tx.exons):
            return False
    return True


def sanity_filter(gff_in: Path, genome_fa: Path, gff_out: Path,
                  report_tsv: Path, keep_only: bool = False) -> dict:
    """Filter the annotation. Returns stats dict. Idempotent via lazy()."""
    if lazy([gff_out, report_tsv], desc="sanity filter"):
        return {}
    genome = _read_genome(genome_fa)
    genes, txs, order = _parse_gff3(gff_in)

    actions: list[tuple[str, str, str, str]] = []  # gid, tid, action, reason
    removed: set[str] = set()
    extended: list[tuple[Tx, str]] = []            # (tx, "extended_stop")

    for tid in order:
        tx = txs[tid]
        prob = _structure_problem(tx)
        if prob == "cds_outside_exon" and _repair_cds_outside_exon(tx):
            actions.append((tx.gid, tid, "repaired", "cds_outside_exon"))
            prob = _structure_problem(tx)  # re-check after repair
        if prob:
            removed.add(tx.tid)
            actions.append((tx.gid, tid, "removed", prob))
            continue
        prot = _translate_cds(genome, tx)
        # internal stop: any stop before the last codon
        if "*" in prot[:-1]:
            removed.add(tid)
            actions.append((tx.gid, tid, "removed", "internal_stop"))
            continue
        if prot and prot[-1] == "*":
            pass  # proper terminal stop: nothing to do
        else:
            # try extending the 3'-most CDS segment by 3 bp
            last = tx.cds_tx_order()[-1]
            seq = genome.get(tx.mrna.seqid, "")
            if tx.mrna.strand == "+":
                nxt = seq[last.end:last.end + 3].upper()   # bases end+1..end+3
            else:
                nxt = seq[max(last.start - 4, 0):last.start - 1].upper()
                nxt = nxt.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
            if len(nxt) == 3 and nxt in STOP_CODONS:
                extended.append((tx, "extended_stop"))
                actions.append((tx.gid, tid, "extended", "extended_stop"))
            else:
                actions.append((tx.gid, tid, "kept", "no_stop"))
        # note non-ATG start (first codon in transcription order)
        if prot and not prot.startswith("M"):
            actions.append((tx.gid, tid, "kept", "non_atg_start"))

    # apply extensions (find covering exon by ORIGINAL coords before mutating)
    for tx, _reason in extended:
        last = tx.cds_tx_order()[-1]
        orig_s, orig_e = last.start, last.end
        cover = next((ex for ex in tx.exons
                      if ex.start <= orig_s and orig_e <= ex.end), None)
        if tx.mrna.strand == "+":
            last.end += 3
            if cover is not None:
                cover.end = max(cover.end, last.end)
        else:
            last.start -= 3
            if cover is not None:
                cover.start = min(cover.start, last.start)

    # recompute spans; drop genes without surviving transcripts
    gene_has_tx: dict[str, list[str]] = {}
    for tid in order:
        if tid not in removed:
            gene_has_tx.setdefault(txs[tid].gid, []).append(tid)
    removed_genes = [g for g in genes if g not in gene_has_tx]
    for gid in removed_genes:
        actions.append((gid, ".", "removed", "all_transcripts_removed"))

    from collections import Counter
    stats = Counter(r for _, _, _, r in actions)
    log.info("sanity filter: %d transcripts removed, %d extended, genes removed %d; %s",
             len(removed), len(extended), len(removed_genes), dict(stats))

    if keep_only:
        # report only: write input unchanged is caller's business; here just report
        with open(report_tsv, "w") as fh:
            for k, v in sorted(stats.items()):
                fh.write(f"# {k}: {v}\n")
            fh.write("gene_id\ttranscript_id\taction\treason\n")
            for a in actions:
                fh.write("\t".join(a) + "\n")
        gff_out.write_text(open(gff_in).read())
        return dict(stats)

    # write filtered gff3
    with open(gff_out, "w") as out:
        out.write("##gff-version 3\n")
        written_genes: set[str] = set()
        for tid in order:
            if tid in removed:
                continue
            tx = txs[tid]
            g = genes.get(tx.gid)
            if g is None:
                continue
            # recompute transcript span from exons
            if tx.exons:
                tx.mrna.start = min(f.start for f in tx.exons)
                tx.mrna.end = max(f.end for f in tx.exons)
            if tx.gid not in written_genes:
                gtxs = gene_has_tx[tx.gid]
                g.start = min(txs[t].mrna.start for t in gtxs)
                g.end = max(txs[t].mrna.end for t in gtxs)
                out.write(_fmt(g))
                written_genes.add(tx.gid)
            out.write(_fmt(tx.mrna))
            for f in sorted(tx.exons, key=lambda x: x.start):
                out.write(_fmt(f))
            for f in sorted(tx.cds, key=lambda x: x.start):
                out.write(_fmt(f))

    with open(report_tsv, "w") as fh:
        for k, v in sorted(stats.items()):
            fh.write(f"# {k}: {v}\n")
        fh.write("gene_id\ttranscript_id\taction\treason\n")
        for a in actions:
            fh.write("\t".join(a) + "\n")
    return dict(stats)


def _fmt(f: Feat) -> str:
    return (f"{f.seqid}\t{f.source}\t{f.ftype}\t{f.start}\t{f.end}\t{f.score}"
            f"\t{f.strand}\t{f.phase}\t{f.attrs}\n")
