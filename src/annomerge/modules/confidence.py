"""Module 6: coverage-aware evidence tiering of the final annotation.

Design doc v2.2 §2.2. Tiering, NOT deletion: absence of evidence is not
evidence of error. Only where evidence actually reached a locus may missing
support count against a transcript.

Per transcript (coding):
  rna support      multi-exon: fraction of CDS introns exactly matched by an
                   RNA junction (seqid, start, end; strand-aware when the
                   junction has a strand). single-exon: fraction of CDS bases
                   with read depth >= min_depth (evidence can only be coverage).
  protein support  fraction of CDS segments overlapped >=1 bp by a miniprot
                   alignment block (match_part), same seqid, strand-aware.
  coverage context RNA "reached" the locus iff >= min_cov_frac of the
                   transcript span has depth >= min_depth (from the BAMs);
                   protein "reached" iff any miniprot alignment exists on the
                   same seqid (protein evidence is sparse by nature, so the
                   bar is deliberately low) — else its absence is neutral.
  goldset origin   CDS structure identical to a goldset model -> strong.

Tiers:
  strong        goldset-origin, or (rna_support AND protein_support)
  moderate      exactly one of rna/protein support
  weak          evidence reached but does not support (negative conclusion)
  uncalibrated  evidence did not reach (neutral; source predictor annotated
                by the caller/report, score NA)
  flagged       transcript carries a sanity-filter note (no_stop /
                non_atg_start / extended_stop) — independent column, not a
                tier override

RNA support definition: introns all matched (multi-exon) or CDS coverage
fraction >= cov_support_min (single-exon). Protein support: >= prot_support_min
of CDS segments covered.

Intrinsic warning flags (do not affect tier): short_orf (< min_orf_aa),
noncanonical_splice (intron not GT-AG/GC-AG), single_exon, neighbour_overlap
(transcript overlaps another transcript on the same strand).
"""
from __future__ import annotations

import logging
import re
import shutil
from collections import Counter
from pathlib import Path

from ..utils import lazy
from .sanity import _parse_gff3, _read_genome, _subseq

log = logging.getLogger("annomerge")


# ---------------------------------------------------------------- evidence io

def load_junctions(bed: Path) -> dict:
    """junctions.bed -> {(seqid, s, e): strand} (1-based, e = intron end)."""
    out: dict[tuple[str, int, int], str] = {}
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            strand = f[5] if len(f) >= 6 and f[5] in "+-" else "."
            out[(f[0], int(f[1]) + 1, int(f[2]))] = strand
    return out


def load_protein_blocks(gff: Path) -> dict[str, list]:
    """miniprot EVM-style protein gff -> {seqid: [(s, e, strand, prot)]} (match_part).

    The protein id (Name=/Target=) is kept so support can be weighted by the
    protein's homolog-set relatedness."""
    blocks: dict[str, list] = {}
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("match_part", "CDS"):
            continue
        m = re.search(r"Name=([^;]+)", f[8]) or re.search(r"Target=(\S+)", f[8])
        prot = m.group(1) if m else ""
        blocks.setdefault(f[0], []).append((int(f[3]), int(f[4]), f[6], prot))
    for v in blocks.values():
        v.sort()
    return blocks


def load_goldset_cds_keys(gff: Path | None) -> set:
    """goldset genome gff3 -> set of CDS structure keys for origin matching."""
    if not gff or not Path(gff).exists():
        return set()
    _, txs, _ = _parse_gff3(Path(gff))
    keys = set()
    for tx in txs.values():
        if tx.cds:
            keys.add(_cds_key(tx))
    return keys


def _cds_key(tx) -> tuple:
    return (tx.mrna.seqid, tx.mrna.strand,
            tuple(sorted((f.start, f.end) for f in tx.cds)))


# ---------------------------------------------------------------- coverage

class BedCovRunner:
    """Region coverage via one `samtools depth -aa -b <super_regions>` pass.

    bedcov scales badly with many small regions (per-region overhead), so
    instead: merge all query intervals (locus spans + single-exon CDS) into
    per-seqid super-regions, stream per-base depth once in C, and attribute
    each base to its transcript via bisect in Python.
    """

    def __init__(self, bams: list, min_depth: int = 5):
        self.bams = bams
        self.min_depth = min_depth
        self.regions: list[tuple[str, int, int, str]] = []  # 1-based, tag
        self.results: dict[str, tuple[float, float]] = {}
        self._ran = False

    def add_locus(self, tag: str, seqid: str, s: int, e: int) -> None:
        if s <= e:
            self.regions.append((seqid, s, e, tag))

    def add_cds(self, tag: str, seqid: str, segs: list) -> None:
        for i, (s, e) in enumerate(segs):
            self.regions.append((seqid, s, e, f"{tag}#cds{i}"))

    def run(self) -> None:
        if self._ran:
            return
        self._ran = True
        if not self.bams or shutil.which("samtools") is None:
            return
        import subprocess
        import tempfile
        # group intervals per seqid and merge into super-regions
        by_seq: dict[str, list] = {}
        for seqid, s, e, tag in self.regions:
            by_seq.setdefault(seqid, []).append((s, e))
        super_bed: dict[str, list] = {}  # seqid -> merged (s,e)
        for seqid, ivs in by_seq.items():
            ivs.sort()
            merged: list[list] = []
            for s, e in ivs:
                if merged and s <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            super_bed[seqid] = merged
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tf:
            for seqid, ivs in super_bed.items():
                for s, e in ivs:
                    tf.write(f"{seqid}\t{s - 1}\t{e}\n")
            bed_path = tf.name
        cmd = ["samtools", "depth", "-aa", "-b", bed_path,
               *[str(b) for b in self.bams]]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, bufsize=1 << 20)
        Path(bed_path).unlink(missing_ok=True)
        if proc.returncode != 0:
            log.warning("samtools depth failed (rc=%d)", proc.returncode)
            return
        # attribute each reported base to the interval(s) covering it.
        # NOTE: keep duplicate intervals (a single-exon transcript's CDS
        # segment equals its locus span) — one base may serve several tags.
        from bisect import bisect_right
        starts: dict[str, list] = {}
        iv_map: dict[str, list] = {}
        for seqid, ivs in by_seq.items():
            srt = sorted(ivs)  # duplicates kept on purpose
            starts[seqid] = [s for s, _ in srt]
            iv_map[seqid] = srt
        self.results = self._attribute(proc.stdout, starts, iv_map)

    def _attribute(self, stdout, starts, iv_map):
        from bisect import bisect_right
        tag_of: dict[tuple, list] = {}
        for seqid, s, e, t in self.regions:
            tag_of.setdefault((seqid, s, e), []).append(t)
        acc: dict[str, list] = {}
        for line in stdout.decode().splitlines():
            f = line.split("\t")
            if len(f) < 3:
                continue
            seqid, pos, d = f[0], int(f[1]), int(f[2])
            st = starts.get(seqid)
            if st is None:
                continue
            i = bisect_right(st, pos) - 1
            if i < 0:
                continue
            # add to every interval containing pos (sorted, may overlap)
            ivs = iv_map[seqid]
            j = i
            while j < len(ivs) and ivs[j][0] <= pos:
                s, e = ivs[j]
                if pos <= e:
                    for tag in tag_of.get((seqid, s, e), []):
                        a = acc.setdefault(tag, [0, 0, 0])
                        a[1] += 1
                        a[2] += d
                        if d >= self.min_depth:
                            a[0] += 1
                j += 1
        out: dict[str, tuple[float, float]] = {}
        for tag, (ge, tot, dsum) in acc.items():
            out[tag] = ((ge / tot) if tot else 0.0,
                        (dsum / tot) if tot else 0.0)
        # per-transcript CDS rollup: mean of segment cov fracs weighted by len
        cds_acc: dict[str, list] = {}
        for tag, (ge, tot, dsum) in acc.items():
            if "#cds" in tag:
                tid = tag.rsplit("#cds", 1)[0]
                a = cds_acc.setdefault(tid, [0, 0])
                a[0] += ge
                a[1] += tot
        for tid, (ge, tot) in cds_acc.items():
            out[f"{tid}#cds"] = ((ge / tot) if tot else 0.0, 0.0)
        return out

    def locus(self, tag: str) -> tuple[float, float]:
        return self.results.get(tag, (0.0, 0.0))

    def cds(self, tag: str) -> float:
        return self.results.get(f"{tag}#cds", (0.0,))[0]


# ---------------------------------------------------------------- tiering

def tier_transcript(rna_ok: bool | None, prot_ok: bool,
                    rna_reached: bool, prot_reached: bool,
                    gold: bool, notes: list[str]) -> str:
    """rna_ok None = no basis (single-exon without coverage data)."""
    if gold:
        return "strong"
    if rna_ok and prot_ok:
        return "strong"
    if rna_ok or prot_ok:
        return "moderate"
    # no support: negative only where evidence actually reached
    neg_rna = rna_reached and rna_ok is False
    neg_prot = prot_reached and not prot_ok
    if neg_rna or neg_prot:
        return "weak"
    return "uncalibrated"


def confidence_annotation(gff_in: Path, genome_fa: Path, junctions_bed: Path | None,
                          protein_gff: Path | None, bams: list[Path],
                          goldset_gff: Path | None, sanity_tsv: Path | None,
                          out_tsv: Path,
                          min_depth: int = 5, min_cov_frac: float = 0.5,
                          cov_support_min: float = 0.8,
                          prot_support_min: float = 0.5,
                          min_orf_aa: int = 50,
                          protein_weights: dict | None = None) -> dict:
    """Annotate every coding transcript with evidence support and a tier.

    protein_weights: {protein_id: relatedness weight} — when given, the
    protein support fraction is weighted: each CDS segment hit contributes
    weight(best covering protein)/2 (closest set = 1.0), so support from
    closer homolog relatives counts more in the final assessment.
    """
    if lazy([out_tsv], desc="confidence annotation"):
        return {}
    genome = _read_genome(genome_fa)
    _, txs, order = _parse_gff3(Path(gff_in))
    junc = load_junctions(junctions_bed) if junctions_bed and Path(junctions_bed).exists() else {}
    pblocks = load_protein_blocks(protein_gff) if protein_gff and Path(protein_gff).exists() else {}
    gold_keys = load_goldset_cds_keys(goldset_gff)

    # sanity notes per transcript
    notes_map: dict[str, list[str]] = {}
    if sanity_tsv and Path(sanity_tsv).exists():
        for line in open(sanity_tsv):
            if line.startswith("#") or line.startswith("gene_id"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[2] in ("kept", "extended"):
                notes_map.setdefault(f[1], []).append(f[3])

    # neighbour overlap (same strand, any other transcript)
    by_seq: dict[str, list] = {}
    for tid in order:
        tx = txs[tid]
        by_seq.setdefault((tx.mrna.seqid, tx.mrna.strand), []).append(tx)
    overlap_tids: set[str] = set()
    for group in by_seq.values():
        group.sort(key=lambda t: t.mrna.start)
        for a, b in zip(group, group[1:]):
            if b.mrna.start <= a.mrna.end:
                overlap_tids.add(a.tid)
                overlap_tids.add(b.tid)

    has_cov = bool(bams) and shutil.which("samtools") is not None
    cov = BedCovRunner(bams, min_depth) if has_cov else None

    rows: list[str] = []
    tier_count: Counter = Counter()
    n_rna_reached = 0

    # pass 1: register every coverage query (locus span + single-exon CDS)
    planned: list[tuple[str, list, bool]] = []  # tid, cds_segs, single_exon
    for tid in order:
        tx = txs[tid]
        if not tx.cds:
            continue
        cds_segs = sorted((f.start, f.end) for f in tx.cds)
        planned.append((tid, cds_segs, len(cds_segs) == 1))
        if cov is not None:
            cov.add_locus(tid, tx.mrna.seqid, tx.mrna.start, tx.mrna.end)
            if len(cds_segs) == 1:
                cov.add_cds(tid, tx.mrna.seqid, cds_segs)
    if cov is not None:
        log.info("confidence: bedcov over %d regions...", len(cov.regions))
        cov.run()

    # pass 2: tier every transcript
    for tid, cds_segs, single_exon in planned:
        tx = txs[tid]
        m = tx.mrna
        introns = [(cds_segs[i][1] + 1, cds_segs[i + 1][0] - 1)
                   for i in range(len(cds_segs) - 1)
                   if cds_segs[i][1] + 1 <= cds_segs[i + 1][0] - 1]
        cov_frac, mean_depth = cov.locus(tid) if cov else (0.0, 0.0)
        rna_reached = bool(junc) and cov_frac >= min_cov_frac
        if rna_reached:
            n_rna_reached += 1
        if introns:
            hits = 0
            for (s, e) in introns:
                st = junc.get((m.seqid, s, e))
                if st is not None and (st == "." or st == m.strand):
                    hits += 1
            rna_frac = hits / len(introns)
            rna_ok: bool | None = hits == len(introns)
            rna_frac_out = f"{rna_frac:.3f}"
            intron_note = f"{hits}/{len(introns)}"
        elif single_exon and cov is not None:
            cds_cov = cov.cds(tid)
            rna_ok = cds_cov >= cov_support_min
            rna_frac_out = f"{cds_cov:.3f}"
            rna_reached = rna_reached or cov_frac >= min_cov_frac
            intron_note = "coverage"
        else:
            rna_ok = None
            rna_frac_out = "NA"
            intron_note = "coverage"

        # ---- protein support (relatedness-weighted when protein_weights given)
        blocks = pblocks.get(m.seqid, [])
        wsum = 0.0
        n_hit = 0
        for (s, e) in cds_segs:
            best_w = 0.0
            for bs, be, bst, bprot in blocks:
                if bs <= e and s <= be and (bst == m.strand or bst == "."):
                    if protein_weights:
                        w = protein_weights.get(bprot, 1.0) / 2.0
                        best_w = max(best_w, w)
                    else:
                        best_w = 1.0
                        break
            if best_w > 0.0:
                n_hit += 1
                wsum += best_w
        if protein_weights:
            prot_frac = wsum / len(cds_segs)
        else:
            prot_frac = n_hit / len(cds_segs)
        prot_ok = prot_frac >= prot_support_min
        prot_reached = bool(blocks)

        # ---- goldset origin
        gold = _cds_key(tx) in gold_keys

        # ---- intrinsic flags
        flags: list[str] = []
        cds_len = sum(e - s + 1 for s, e in cds_segs)
        if cds_len // 3 < min_orf_aa:
            flags.append("short_orf")
        if len(cds_segs) == 1:
            flags.append("single_exon")
        for (s, e) in introns:
            donor = _subseq(genome, m.seqid, s, s + 1, m.strand)
            acceptor = _subseq(genome, m.seqid, e - 1, e, m.strand)
            if (donor, acceptor) not in (("GT", "AG"), ("GC", "AG")):
                flags.append("noncanonical_splice")
                break
        if tid in overlap_tids:
            flags.append("neighbour_overlap")

        notes = notes_map.get(tid, [])
        tier = tier_transcript(rna_ok, prot_ok, rna_reached, prot_reached,
                               gold, notes)
        tier_count[tier] += 1
        rows.append("\t".join([
            tx.gid, tid, m.seqid, m.strand, str(m.start), str(m.end),
            tier, "yes" if gold else "no",
            rna_frac_out, intron_note if introns else "coverage",
            f"{cov_frac:.3f}", f"{mean_depth:.1f}",
            "yes" if rna_reached else "no",
            f"{prot_frac:.3f}", f"{n_hit}/{len(cds_segs)}",
            "yes" if prot_reached else "no",
            ",".join(notes) if notes else ".",
            ",".join(flags) if flags else ".",
        ]))

    with open(out_tsv, "w") as fh:
        fh.write("# tiers: " + ", ".join(f"{k}={v}" for k, v in
                                         sorted(tier_count.items())) + "\n")
        fh.write("# rna_reached_loci: %d/%d\n" % (n_rna_reached,
                                                  tier_count.total()))
        fh.write("gene_id\ttranscript_id\tseqid\tstrand\tstart\tend\ttier\t"
                 "goldset_origin\trna_support_frac\trna_basis\tlocus_cov_frac\t"
                 "locus_mean_depth\trna_reached\tprotein_support_frac\t"
                 "protein_basis\tprotein_reached\tsanity_notes\tintrinsic_flags\n")
        fh.write("\n".join(rows) + "\n")
    log.info("confidence: %s", dict(tier_count))
    return dict(tier_count)
