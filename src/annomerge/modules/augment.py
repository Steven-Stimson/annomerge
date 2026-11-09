"""Module: augment — EVM-anchored isoform resolution + two-tier UTR backfill.

Design increment v2.3 (AnnoMerge_isoform模式设计增量_v2.3.md):
  L1 EVM loci     anchor = EVM consensus CDS; StringTie-compatible donors
                  adopted as isoforms (CDS=EVM, exon chain/UTR=donor);
                  differently-spliced donors contribute their own best goldset
                  ORF as AS isoform; fusion donors naturally fail and drop.
  L2 goldset-only transcripts not overlapping any EVM gene -> genes by
                  overlap clustering (best ORF per transcript).
  L3 de novo     at loci covered by neither -> W_raw arbitration (reuse
                  priority._arbitrate_loci).
  UTR backfill   two-tier for UTR-less transcripts: tier1 StringTie donors
                  (no ORF required), tier2 de novo donors; <=5000bp/side.

Matching discipline: per-transcript, NEVER per-gene longest-donor.
Internal-invariant: a donor's introns within the CDS span must equal the CDS
introns exactly — CDS internal structure is never rewritten by backfill.
"""
from __future__ import annotations
from ..utils import inject_utr_features

import logging
import re
from collections import Counter
from pathlib import Path

from ..utils import ensure_dir, lazy
from .priority import _arbitrate_loci, _locus_index, _overlaps
from .sanity import _parse_gff3, _read_genome, _translate_cds

log = logging.getLogger("annomerge")

TX_RE = re.compile(r"^(?P<tx>.+)\.p(?P<n>\d+)$")
MAX_EXTEND = 5000
JITTER_MERGE_BP = 50          # donors closer than this on BOTH ends collapse


# ------------------------------------------------------------------ helpers
def _cds_introns(tx) -> tuple:
    segs = sorted((c.start, c.end) for c in tx.cds)
    return tuple((segs[i][1] + 1, segs[i + 1][0] - 1) for i in range(len(segs) - 1))


def _cds_span(tx) -> tuple[int, int]:
    return min(c.start for c in tx.cds), max(c.end for c in tx.cds)


def _chain(exons) -> tuple:
    return tuple(sorted((x.start, x.end) for x in exons))


def _utr_parts(tx) -> tuple[list, list]:
    """(five_prime, three_prime) UTR intervals from exon-CDS difference."""
    if not tx.cds or not tx.exons:
        return [], []
    cds = sorted(tx.cds, key=lambda c: c.start)
    exons = sorted(tx.exons, key=lambda x: x.start)
    plus = tx.mrna.strand != "-"
    cs, ce = cds[0].start, cds[-1].end
    five: list = []
    three: list = []
    for x in exons:
        if x.end < cs:
            (five if plus else three).append((x.start, x.end))
        elif x.start > ce:
            (three if plus else five).append((x.start, x.end))
        else:
            if x.start < cs:
                (five if plus else three).append((x.start, cs - 1))
            if x.end > ce:
                (three if plus else five).append((ce + 1, x.end))
    return five, three


def _donor_covers(donor_introns: tuple, donor_span, cds_introns, cds_span) -> bool:
    """CDS-structure invariant: every donor intron that TOUCHES the CDS span
    must be exactly one of the CDS introns (fully-inside-only checking misses
    introns straddling cs/ce, which would cut the first/last CDS segment);
    donor span must also cover the full CDS."""
    ds, de = donor_span
    cs, ce = cds_span
    if ds > cs or de < ce:
        return False
    cds_set = tuple(cds_introns)
    for iv in donor_introns:
        if iv[1] < cs or iv[0] > ce:
            continue                  # entirely outside the CDS span
        if tuple(iv) not in cds_set:
            return False
    return True


def _neighbor_window(cores, seqid, cs, ce, max_extend):
    """UTR extension window [lo,hi] around CDS [cs,ce], clipped at the nearest
    other gene's CDS core (any strand) and at max_extend per side."""
    lo, hi = cs - max_extend, ce + max_extend
    for (s, e) in cores.get(seqid, []):
        if s <= cs and e >= ce:
            continue                      # self/container (shouldn't happen)
        if e < cs:
            lo = max(lo, e + 1)
        elif s > ce:
            hi = min(hi, s - 1)
        elif s >= cs and s <= ce:
            lo = max(lo, s - 1)
        elif e >= cs and e <= ce:
            hi = min(hi, e + 1)
    return lo, hi


class _Gene:
    __slots__ = ("seqid", "strand", "isoforms", "src", "name")

    def __init__(self, seqid, strand, src, name=""):
        self.seqid = seqid
        self.strand = strand
        self.isoforms: list = []      # list[Tx]
        self.src = src                # L1 | L2 | L3
        self.name = name

    def span(self):
        lo = min(x.start for tx in self.isoforms for x in (tx.exons or tx.cds))
        hi = max(x.end for tx in self.isoforms for x in (tx.exons or tx.cds))
        return lo, hi


# ------------------------------------------------------------ goldset parse
def _group_goldset(goldset_gff: Path) -> dict:
    """goldset genome gff3 -> {stringtie_txkey: {exons, strand, seqid, orfs}}.
    Every ORF mRNA carries the full StringTie exon chain; ORFs on the same
    transcript share it. orfs: {p_number: Tx}."""
    _, txs, order = _parse_gff3(Path(goldset_gff))
    groups: dict = {}
    for tid in order:
        tx = txs[tid]
        m = TX_RE.match(tid)
        key, n = (m.group("tx"), int(m.group("n"))) if m else (tid, 1)
        g = groups.setdefault(key, {
            "seqid": tx.mrna.seqid, "strand": tx.mrna.strand,
            "exons": tx.exons, "orfs": {}})
        g["orfs"][n] = tx
    return groups


def _best_orf(group) -> object | None:
    """Best ORF of a StringTie transcript: longest CDS, tie -> lower p number
    (TransDecoder p1 is the representative)."""
    best, best_len = None, -1
    for n in sorted(group["orfs"]):
        tx = group["orfs"][n]
        L = sum(c.end - c.start + 1 for c in tx.cds)
        if L > best_len:
            best, best_len = tx, L
    return best


def _chain_introns(exons) -> tuple:
    ex = sorted(exons, key=lambda x: x.start)
    return tuple((ex[i].end + 1, ex[i + 1].start - 1) for i in range(len(ex) - 1))


# ------------------------------------------------------------- main entry
def run_augment(vi, out_dir: Path, evm_gff3: Path, goldset_gff: Path,
                stringtie_gtf: Path, preds: dict, scores: dict | None,
                max_extend: int = MAX_EXTEND,
                jitter_merge_bp: int = JITTER_MERGE_BP,
                dedup_mode: str = "geometric",
                junctions_bed: Path | None = None,
                tpm_ratio: float = 0.1,
                min_tpm: float = 0.1,
                junc_min_reads: int = 3,
                pwm_gate: float | None = None) -> dict:
    """dedup_mode: geometric = legacy ±jitter_merge_bp span clustering;
    evidence = EviAnn-style voting — splice-PWM gate on non-anchor introns,
    junction-read corroboration of unique UTR introns, and TPM-ratio test
    for terminal-span duplicates (no bp-distance threshold anywhere)."""
    from .splicepwm import train_from_transcripts, DEFAULT_GATE
    if pwm_gate is None:
        pwm_gate = DEFAULT_GATE
    d = ensure_dir(Path(out_dir))
    stats: Counter = Counter()
    _, evm_txs, evm_order = _parse_gff3(Path(evm_gff3))
    groups = _group_goldset(goldset_gff)
    genome = _read_genome(vi.genome)
    st_tpms = _parse_stringtie_tpm(stringtie_gtf)
    # evidence-mode machinery (built lazily; cheap when geometric)
    pwm = train_from_transcripts(genome, [g["orfs"][n] for g in groups.values()
                                          for n in g["orfs"]]) \
        if dedup_mode == "evidence" else None
    junc = _load_junctions(junctions_bed, junc_min_reads) \
        if (dedup_mode == "evidence" and junctions_bed) else None
    stats["dedup_mode"] = dedup_mode

    used_donors: set = set()
    genes: list = []
    # CDS cores of all final genes (anchor CDS per L1 gene, ORF CDS for L2/L3)
    # used for the neighbor-stop rule; built incrementally (genes don't overlap)
    cores: dict = {}

    def _add_core(gene):
        cs = min(c.start for tx in gene.isoforms for c in tx.cds)
        ce = max(c.end for tx in gene.isoforms for c in tx.cds)
        cores.setdefault(gene.seqid, []).append((cs, ce))

    # ============ Phase 1: EVM anchor stubs + register ALL anchor cores.
    # The neighbor-stop index must be complete BEFORE any UTR extension, or
    # early genes extend into later genes' CDS (order-dependent overlaps).
    l1_genes: dict = {}
    for gid in evm_order:
        anc = evm_txs[gid]
        if not anc.cds:
            continue
        gene = _Gene(anc.mrna.seqid, anc.mrna.strand, "L1", name=gid)
        gene.isoforms.append(anc)
        l1_genes[gid] = gene
        stats["l1_genes"] += 1
        _add_core(gene)

    # ============ Phase 2: L2 goldset genes at EVM-absent loci ONLY
    l2_records = []
    for key, g in groups.items():
        gspan = (min(x.start for x in g["exons"]), max(x.end for x in g["exons"]))
        if any(s <= gspan[1] and e >= gspan[0] for (s, e) in cores.get(g["seqid"], [])):
            # donor on an EVM-anchored locus: routed to L1 adoption in Phase 4,
            # never an L2 gene (the EVM anchor owns its locus)
            stats["l2_routed_to_l1"] += 1
            continue
        orf = _best_orf(g)
        if orf is None or not orf.cds:
            stats["l2_no_orf_skip"] += 1
            continue
        l2_records.append({"name": "goldset", "tx": orf, "gid": key,
                           "score": 1.0, "rs": 0.0})
    kept_l2, l2_arb = _arbitrate_loci(l2_records)   # single name -> pure clustering
    stats["l2_genes"] = l2_arb["loci_total"]
    l2_clusters: dict = {}
    for r in kept_l2:
        l2_clusters.setdefault(r["gid"], []).append(r["tx"])
    l2_gene_list: list = []
    for key, members in l2_clusters.items():
        tx0 = members[0]
        gene = _Gene(tx0.mrna.seqid, tx0.mrna.strand, "L2", name=key)
        gene.isoforms.extend(members)
        stats["l2_isoforms"] += len(members)
        l2_gene_list.append(gene)
        _add_core(gene)

    # ============ Phase 3: L3 de novo gap-fill at uncovered loci
    final_idx = _locus_index({f"g{i}": type("F", (), {"seqid": g.seqid, "strand": g.strand,
                                                      "start": g.span()[0], "end": g.span()[1]})()
                             for i, g in enumerate(list(l1_genes.values()) + l2_gene_list)})
    l3_records = []
    for name, gff in preds.items():
        _, txs, order = _parse_gff3(Path(gff))
        for tid in order:
            tx = txs[tid]
            foot = tx.exons if tx.exons else tx.cds
            if not foot:
                continue
            hit = False
            for f in foot:
                if _overlaps(final_idx, f.seqid, tx.mrna.strand, f.start, f.end):
                    hit = True
                    break
            if hit:
                stats["l3_overlap_drop"] += 1
                continue
            l3_records.append({"name": name, "tx": tx, "gid": tx.gid,
                               "score": (scores or {}).get(name, (0.0, 0.0))[0],
                               "rs": (scores or {}).get(name, (0.0, 0.0))[1]})
    kept_l3, l3_arb = _arbitrate_loci(l3_records)
    stats.update({f"l3_{k}": v for k, v in l3_arb.items()})
    l3_groups: dict = {}
    for r in kept_l3:
        l3_groups.setdefault((r["name"], r["gid"]), []).append(r["tx"])
    l3_gene_list: list = []
    for (name, ogid), members in l3_groups.items():
        tx0 = members[0]
        gene = _Gene(tx0.mrna.seqid, tx0.mrna.strand, "L3", name=f"{name}:{ogid}")
        gene.isoforms.extend(members)
        stats["l3_genes"] += 1
        stats["l3_isoforms"] += len(members)
        l3_gene_list.append(gene)
        _add_core(gene)

    # ============ Phase 4: L1 isoform adoption (complete core index)
    gene_donor_tpm_max: dict = {}     # gid -> max TPM over all donors evaluated
    for gid, gene in l1_genes.items():
        anc = gene.isoforms[0]
        cs, ce = _cds_span(anc)
        ai = _cds_introns(anc)
        utr_cands = []                     # (key, group, dspan) awaiting jitter merge
        gene_donor_keys: list = []
        for key, g in groups.items():
            if g["strand"] != anc.mrna.strand or g["seqid"] != anc.mrna.seqid:
                continue
            dspan = (min(x.start for x in g["exons"]), max(x.end for x in g["exons"]))
            if dspan[0] > ce or dspan[1] < cs:      # span disjoint
                continue
            gene_donor_keys.append(key)
            din = _chain_introns(g["exons"])
            if _donor_covers(din, dspan, ai, (cs, ce)):
                utr_cands.append((key, g, dspan, din))
                continue
            # AS branch: donor's own best ORF
            orf = _best_orf(g)
            if orf is None:
                stats["donor_no_orf_skip"] += 1
                continue
            cds_seg = tuple(sorted((c.start, c.end) for c in orf.cds))
            if cds_seg == tuple(sorted((c.start, c.end) for c in anc.cds)):
                continue                      # identical to anchor, no value
            if "*" in _translate_cds(genome, orf)[:-1]:
                stats["as_orf_invalid_skip"] += 1
                continue
            if not any(c.start <= ce and c.end >= cs for c in orf.cds):
                stats["donor_no_cds_overlap"] += 1
                continue
            # reject if the AS exon chain overlaps another gene's CDS core
            clash = False
            for x in orf.exons:
                for (s, e) in cores.get(x.seqid, []):
                    if x.start <= e and x.end >= s and not (s <= cs and e >= ce):
                        clash = True
                        break
                if clash:
                    break
            if clash:
                stats["as_core_clash_skip"] += 1
                continue
            if dedup_mode == "evidence" and not _donor_gates_pass(
                    g, din, ai, cs, ce, genome, pwm, pwm_gate, junc, stats, "as"):
                continue
            orf.mrna.attrs = f"donor_tpm={st_tpms.get(key, 0.0)}"
            gene.isoforms.append(orf)
            used_donors.add(key)
            stats["isoform_as"] += 1
        # near-duplicate (assembly jitter) collapse BEFORE isoform building
        if dedup_mode == "evidence":
            kept, n_merged = _dedup_evidence(utr_cands, st_tpms, genome, pwm,
                                             pwm_gate, junc, tpm_ratio, min_tpm,
                                             (cs, ce), stats)
        else:
            kept, n_merged = _merge_jitter(utr_cands, st_tpms, jitter_merge_bp)
        stats["donor_jitter_merged"] += n_merged
        gene_donor_tpm_max[gid] = max((st_tpms.get(k, 0.0) for k in gene_donor_keys),
                                      default=0.0)
        lo, hi = _neighbor_window(cores, anc.mrna.seqid, cs, ce, max_extend)
        for key, g, _, _ in kept:
            iso = _build_utr_isoform(anc, g, lo, hi, tpm=st_tpms.get(key, 0.0))
            if iso is not None:
                gene.isoforms.append(iso)
                used_donors.add(key)
                stats["isoform_utr"] += 1
        # collapse exact duplicate isoforms (identical exon chain + CDS)
        seen, uniq = set(), []
        for tx in gene.isoforms:
            k = (_chain(tx.exons), tuple(sorted((c.start, c.end) for c in tx.cds)))
            if k not in seen:
                seen.add(k)
                uniq.append(tx)
        if len(uniq) != len(gene.isoforms):
            stats["isoform_exact_dup_collapsed"] += len(gene.isoforms) - len(uniq)
            gene.isoforms = uniq
        stats["isoforms_total"] += len(uniq)

    genes.extend(l1_genes.values())
    genes.extend(l2_gene_list)
    genes.extend(l3_gene_list)

    # ---------------- UTR backfill (tier1 StringTie, tier2 de novo)
    st_donors = _parse_stringtie_chains(stringtie_gtf)
    _backfill_utrs(genes, st_donors, preds, scores, genome, stats, max_extend, cores)

    # -------- collapse same-exon-chain ORF shadows (AFTER backfill) ---------
    # Backfill rebuilds the anchor's exon chain from the same StringTie pool
    # that produced the AS isoforms, so chains only become final here.  An
    # identical exon chain with a different CDS is ONE transcript with
    # ORF-calling ambiguity, not a biological isoform (EviAnn one-model-per-
    # locus semantics at the ORF level).  Keep the longest-CDS representative;
    # ties keep the earliest (anchor has priority).
    for gene in genes:
        by_chain: dict = {}
        for tx in gene.isoforms:
            by_chain.setdefault((_chain(tx.exons), tx.mrna.strand), []).append(tx)
        if not any(len(m) > 1 for m in by_chain.values()):
            continue
        orig = {id(t): i for i, t in enumerate(gene.isoforms)}
        merged = []
        for members in by_chain.values():
            if len(members) == 1:
                merged.extend(members)
            else:
                members = sorted(members,
                                 key=lambda t: (-sum(c.end - c.start + 1
                                                     for c in t.cds),
                                                orig[id(t)]))
                merged.append(members[0])
                stats["isoform_orf_shadow_collapsed"] += len(members) - 1
        merged.sort(key=lambda t: orig[id(t)])       # anchor stays first
        gene.isoforms = merged

    # ---------------- fold near-anchor UTR duplicates (AFTER backfill) ------
    # Backfill extends the anchor from the same donor pool that produced the
    # adopted UTR isoforms, so the real jitter test is against the FINAL
    # anchor chain: an isoform with the same intron signature as the anchor
    # carries no splicing diversity.  geometric mode folds those within
    # +-jitter_merge_bp on both ends; evidence mode folds those whose donor
    # TPM fails the TPM-ratio test against the gene donor pool (no bp distance).
    TPM_RE = re.compile(r"donor_tpm=([\d.eE+-]+)")
    for gene in l1_genes.values():
        if len(gene.isoforms) < 2:
            continue
        anc = gene.isoforms[0]
        ai_anc = _chain_introns(anc.exons)
        if dedup_mode == "evidence":
            thr_tpm = max(min_tpm,
                          tpm_ratio * gene_donor_tpm_max.get(gene.name, 0.0))
        kept = [anc]
        for tx in gene.isoforms[1:]:
            if _chain_introns(tx.exons) != ai_anc:
                kept.append(tx)
                continue
            if dedup_mode == "evidence":
                m = TPM_RE.search(getattr(tx.mrna, "attrs", "") or "")
                t = float(m.group(1)) if m else 0.0
                fold = t < thr_tpm
            else:
                fold = (abs(tx.mrna.start - anc.mrna.start) < jitter_merge_bp
                        and abs(tx.mrna.end - anc.mrna.end) < jitter_merge_bp)
            if fold:
                if (tx.mrna.end - tx.mrna.start) > (anc.mrna.end - anc.mrna.start):
                    anc.exons = tx.exons
                    anc.mrna.start = tx.mrna.start
                    anc.mrna.end = tx.mrna.end
                stats["isoform_folded_into_anchor"] += 1
            else:
                kept.append(tx)
        gene.isoforms = kept
    # final isoform counts AFTER the fold pass
    stats["isoforms_total"] = sum(len(g.isoforms) for g in genes)
    iso_dist = Counter(len(g.isoforms) for g in genes)
    stats["isoforms_per_gene"] = dict(sorted(iso_dist.items()))

    # ---------------- write
    out_gff3 = d / "augment.gff3"
    if not lazy([out_gff3], desc="augment merge"):
        _collapse_isoforms(genes, stats)
        log.info("augment: %d isoform(s) collapsed (identical CDS, longest-UTR kept)", stats.get("isoforms_collapsed", 0))
        _write(out_gff3, genes, stats)
        n5, n3 = inject_utr_features(out_gff3, out_gff3)
        log.info("augment: UTR rows added: %d 5' + %d 3'", n5, n3)
    stsv = d / "augment_stats.tsv"
    with open(stsv, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in sorted(stats.items()):
            fh.write(f"{k}\t{v}\n")
    log.info("augment stats -> %s", stsv)
    return {"augment_gff3": out_gff3, "dir": d, "stats": dict(stats)}


def _build_utr_isoform(anc, group, lo, hi, tpm: float = 0.0):
    """New Tx: CDS = anchor's, exon chain = donor's clipped to [lo,hi].
    donor_tpm is carried in mRNA attrs for the evidence-mode fold test."""
    from .sanity import Tx, Feat
    cs, ce = _cds_span(anc)
    exons = []
    for x in sorted(group["exons"], key=lambda f: f.start):
        s, e = max(x.start, lo), min(x.end, hi)
        if s <= e:
            exons.append(Feat(x.seqid, "goldset", "exon", s, e, ".", x.strand, ".", ""))
    if not exons or exons[0].start > cs or exons[-1].end < ce:
        return None
    mrna = Feat(anc.mrna.seqid, "aug", "mRNA", exons[0].start, exons[-1].end,
                ".", anc.mrna.strand, ".",
                f"utr_source=stringtie;donor_tpm={tpm}")
    return Tx("", "", mrna, exons, list(anc.cds))


def _parse_stringtie_chains(gtf: Path) -> dict:
    """stringtie gtf -> {tid: (seqid, strand, [(s,e),...])} (donors need no ORF)."""
    import re as _re
    chains: dict = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = _re.search(r'transcript_id "([^"]+)"', f[8])
        if m:
            chains.setdefault(m.group(1), (f[0], f[6], []))[2].append((int(f[3]), int(f[4])))
    return chains


def _parse_stringtie_tpm(gtf: Path) -> dict:
    """stringtie gtf -> {transcript_id: TPM} (missing -> 0.0)."""
    tpms: dict = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        t = re.search(r'TPM "([^"]+)"', f[8])
        if m:
            try:
                tpms[m.group(1)] = float(t.group(1)) if t else 0.0
            except ValueError:
                tpms[m.group(1)] = 0.0
    return tpms


def _load_junctions(bed: Path, min_reads: int) -> dict:
    """junctions.bed (chrom,start,end,n_samples,n_reads,strand; 0-based
    half-open from pysam) -> {(seqid, strand, s, e)} in GFF 1-based inclusive
    coords with n_reads >= min_reads; strand '.' is stored both ways so either
    transcription direction corroborates it."""
    junc: set = set()
    with open(bed) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6 or int(c[4]) < min_reads:
                continue
            s, e = int(c[1]) + 1, int(c[2])      # BED 0-based half-open -> GFF 1-based inclusive
            if c[5] == ".":
                junc.add((c[0], "+", s, e))
                junc.add((c[0], "-", s, e))
            else:
                junc.add((c[0], c[5], s, e))
    return junc


def _utr_introns(din: tuple, cs: int, ce: int) -> tuple:
    """Donor introns lying entirely outside the CDS span (UTR-region introns);
    introns straddling the CDS boundary cannot exist here (donor_cover check
    guarantees any intron touching the CDS span IS a CDS intron)."""
    return tuple(iv for iv in din if iv[1] < cs or iv[0] > ce)


def _donor_gates_pass(g, din, ai, cs, ce, genome, pwm, pwm_gate, junc, stats,
                      prefix: str) -> bool:
    """Evidence gates shared by UTR-candidate and AS-branch donors (EviAnn
    fix_suspect_introns + RNA-intron-DB vote, reimplemented):
      PWM gate   — every non-anchor intron must have a splice-site score
                   strictly above pwm_gate (weak site = assembly artifact);
      junction vote — if ALL of a donor's non-anchor introns are absent from
                   the junction-read set, the donor is a spurious splicing ->
                   drop.  Donors with no non-anchor introns always pass."""
    novel = tuple(iv for iv in din if tuple(iv) not in tuple(ai))
    if pwm is not None:
        for s, e in novel:
            sc = pwm.score_intron(genome, g["seqid"], g["strand"], s, e)
            if sc is not None and sc <= pwm_gate:
                stats[f"{prefix}_pwm_gate_drop"] += 1
                return False
    if novel and junc is not None and \
            not any((g["seqid"], g["strand"], s, e) in junc for s, e in novel):
        stats[f"{prefix}_junction_vote_drop"] += 1
        return False
    return True


def _dedup_evidence(cands, tpms, genome, pwm, pwm_gate, junc,
                    tpm_ratio, min_tpm, cds_span, stats) -> tuple[list, int]:
    """EviAnn-style evidence voting over UTR donor candidates, replacing the
    geometric ±bp clustering.  Three mechanisms, no distance threshold:
      1. shared evidence gates (PWM + junction vote) per donor;
      2. cluster by UTR-intron signature (terminal-span-only variants land in
         one cluster);
      3. TPM-ratio test inside each cluster: keep the representative plus any
         donor with TPM >= max(min_tpm, tpm_ratio * rep_TPM) — a donor that
         is both structurally identical to and much less abundant than the
         representative is an assembly shadow, not an isoform.
    cands: [(key, group, (ds,de), dintrons), ...] with the anchor CDS span."""
    cs, ce = cds_span
    n0 = len(cands)
    if not cands:
        return [], 0
    # gate 1: shared evidence gates per donor
    gated = []
    for key, g, dspan, din in cands:
        if _donor_gates_pass(g, din, (), cs, ce, genome, pwm, pwm_gate, junc,
                             stats, "donor"):
            gated.append((key, g, dspan, din))
    # gates 2+3: cluster by UTR-intron signature, TPM-ratio test inside
    by_sig: dict = {}
    for c in gated:
        sig = _utr_introns(c[3], cs, ce)
        by_sig.setdefault(sig, []).append(c)
    kept = []
    for sig, members in by_sig.items():
        rep = max(members, key=lambda c: (tpms.get(c[0], 0.0),
                                          c[2][1] - c[2][0]))
        rep_tpm = tpms.get(rep[0], 0.0)
        for c in members:
            t = tpms.get(c[0], 0.0)
            if c is rep or t >= max(min_tpm, tpm_ratio * rep_tpm):
                kept.append(c)
            else:
                stats["donor_tpm_ratio_merged"] += 1
    return kept, n0 - len(kept)


def _merge_jitter(cands: list, tpms: dict, thr: int) -> tuple[list, int]:
    """Collapse near-duplicate UTR donors BEFORE isoform building: donors whose
    span starts AND ends both lie within thr bp of a cluster representative
    AND whose internal intron chains are identical are assembly jitter of one
    biological transcript — keep the highest-TPM representative. Donors with
    a distinct UTR intron structure always survive (real UTR-alt-splicing).
    cands: [(key, group, (ds,de), dintrons), ...] -> (kept, n_merged)."""
    if thr <= 0 or len(cands) < 2:
        return list(cands), 0
    cands = sorted(cands, key=lambda c: (c[2][0], c[2][1]))
    clusters: list[list] = []
    for c in cands:
        for cl in clusters:
            rs, re_ = cl[0][2]
            if (abs(c[2][0] - rs) < thr and abs(c[2][1] - re_) < thr
                    and c[3] == cl[0][3]):
                cl.append(c)
                break
        else:
            clusters.append([c])
    kept = []
    for cl in clusters:
        best = max(cl, key=lambda c: (tpms.get(c[0], 0.0),
                                      c[2][1] - c[2][0]))
        kept.append(best)
    return kept, len(cands) - len(kept)


def _backfill_utrs(genes, st_donors, preds, scores, genome, stats, max_extend, cores):
    """For UTR-less transcripts extend exons from tier1 StringTie then tier2
    de novo donors (internal introns must match CDS exactly)."""
    for gene in genes:
        for tx in gene.isoforms:
            five, three = _utr_parts(tx)
            if five or three or not tx.cds:
                continue
            cs, ce = _cds_span(tx)
            ai = _cds_introns(tx)
            best = None
            # tier1: StringTie chains
            for tid, (sq, st, chain) in st_donors.items():
                if sq != tx.mrna.seqid or (st not in (".", tx.mrna.strand)):
                    continue
                span = (min(s for s, _ in chain), max(e for _, e in chain))
                din = tuple(sorted((chain[i][1] + 1, chain[i + 1][0] - 1)
                                   for i in range(len(chain) - 1)))
                if not _donor_covers(din, span, ai, (cs, ce)):
                    continue
                ext = (cs - span[0]) + (span[1] - ce)
                if best is None or ext > best[0]:
                    best = (ext, "stringtie", chain)
            # tier2: de novo chains
            if best is None:
                cand = []
                for name, gff in preds.items():
                    _, txs, order = _parse_gff3(Path(gff))
                    for tid in order:
                        dtx = txs[tid]
                        if dtx.mrna.seqid != tx.mrna.seqid or dtx.mrna.strand != tx.mrna.strand:
                            continue
                        foot = dtx.exons if dtx.exons else dtx.cds
                        if not foot:
                            continue
                        span = (min(f.start for f in foot), max(f.end for f in foot))
                        din = _chain_introns(foot)
                        if not _donor_covers(din, span, ai, (cs, ce)):
                            continue
                        ext = (cs - span[0]) + (span[1] - ce)
                        sc = (scores or {}).get(name, (0.0, 0.0))
                        cand.append((ext, sc, name,
                                     [(f.start, f.end) for f in sorted(foot, key=lambda f: f.start)]))
                if cand:
                    cand.sort(key=lambda c: (-c[1][0], -c[1][1], -c[0]))
                    best = (cand[0][0], f"denovo:{cand[0][2]}", cand[0][3])
            if best is None:
                stats["utr_none"] += 1
                continue
            ext, src, chain = best
            lo, hi = _neighbor_window(cores, tx.mrna.seqid, cs, ce, max_extend)
            from .sanity import Feat
            new_exons = [Feat(tx.mrna.seqid, "aug", "exon", max(s, lo), min(e, hi),
                              ".", tx.mrna.strand, ".", "")
                         for s, e in sorted(chain) if max(s, lo) <= min(e, hi)]
            # guard: rebuilt chain must be non-empty and still cover the CDS
            if not new_exons or new_exons[0].start > cs or new_exons[-1].end < ce:
                stats["utr_backfill_failed"] += 1
                continue
            tx.exons = new_exons
            tx.mrna.start = new_exons[0].start
            tx.mrna.end = new_exons[-1].end
            tx.mrna.attrs = f"utr_source={src.split(':')[0]}"
            stats[f"utr_from_{src.split(':')[0]}"] += 1


def _collapse_isoforms(genes: list, stats: Counter) -> None:
    """Collapse isoforms with identical CDS structure within a gene (same
    strand + same CDS intervals => same intron chain).  The representative
    is the isoform with the longest total UTR (UTR evidence reaches
    furthest); ties break to the higher donor TPM.  Recorded in
    stats['isoforms_collapsed']."""
    for gene in genes:
        if len(gene.isoforms) <= 1:
            continue
        groups: dict = {}
        for tx in gene.isoforms:
            key = (tx.mrna.strand,
                   tuple(sorted((c.start, c.end) for c in tx.cds)))
            groups.setdefault(key, []).append(tx)
        if len(groups) == len(gene.isoforms):
            continue
        keep = []
        for txs in groups.values():
            if len(txs) == 1:
                keep.append(txs[0])
                continue

            def _utr_total(t):
                return sum(e - s + 1 for part in _utr_parts(t)
                           for (s, e) in part)

            def _tpm(t):
                m = re.search(r"donor_tpm=([\d.eE+-]+)",
                              getattr(t.mrna, "attrs", "") or "")
                return float(m.group(1)) if m else 0.0

            txs.sort(key=lambda t: (_utr_total(t), _tpm(t)), reverse=True)
            keep.append(txs[0])
            stats["isoforms_collapsed"] += len(txs) - 1
        order = {id(t): i for i, t in enumerate(gene.isoforms)}
        keep.sort(key=lambda t: order[id(t)])
        gene.isoforms = keep


def _write(out: Path, genes: list, stats: Counter) -> None:
    iso_dist = Counter()
    with open(out, "w") as fh:
        fh.write("##gff-version 3\n")
        for i, gene in enumerate(genes, 1):
            gid = f"AUG.g{i}"
            gs, ge = gene.span()
            fh.write(f"{gene.seqid}\tAnnoMerge\tgene\t{gs}\t{ge}\t.\t{gene.strand}\t.\t"
                     f"ID={gid};Name={gene.name};layer={gene.src}\n")
            iso_dist[len(gene.isoforms)] += 1
            for j, tx in enumerate(gene.isoforms, 1):
                tid = f"{gid}.t{j}"
                src = (getattr(tx.mrna, "attrs", "") or "")
                extra = f";{src}" if "utr_source" in src else ""
                fh.write(f"{tx.mrna.seqid}\tAnnoMerge\tmRNA\t{tx.mrna.start}\t{tx.mrna.end}\t.\t"
                         f"{tx.mrna.strand}\t.\tID={tid};Parent={gid}{extra}\n")
                for k, x in enumerate(sorted(tx.exons, key=lambda f: f.start), 1):
                    fh.write(f"{x.seqid}\tAnnoMerge\texon\t{x.start}\t{x.end}\t.\t"
                             f"{x.strand}\t.\tID={tid}.exon{k};Parent={tid}\n")
                for c in tx.cds_tx_order():
                    fh.write(f"{c.seqid}\tAnnoMerge\tCDS\t{c.start}\t{c.end}\t.\t"
                             f"{c.strand}\t{c.phase}\tID={tid}.cds;Parent={tid}\n")
                five, three = _utr_parts(tx)
                for k, (s, e) in enumerate(five, 1):
                    fh.write(f"{tx.mrna.seqid}\tAnnoMerge\tfive_prime_UTR\t{s}\t{e}\t.\t"
                             f"{tx.mrna.strand}\t.\tID={tid}.utr5.{k};Parent={tid}\n")
                for k, (s, e) in enumerate(three, 1):
                    fh.write(f"{tx.mrna.seqid}\tAnnoMerge\tthree_prime_UTR\t{s}\t{e}\t.\t"
                             f"{tx.mrna.strand}\t.\tID={tid}.utr3.{k};Parent={tid}\n")
    stats["isoforms_per_gene"] = dict(sorted(iso_dist.items()))
