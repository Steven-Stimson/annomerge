"""Splice-site PWM scoring (EviAnn-style evidence voting, reimplemented).

EviAnn (compute_junction_scores_bed.pl + score_transcripts_with_hmms.pl)
trains 0/1/2-order Markov models on all observed splice junctions and gates
transcripts by a junction score threshold.  We use an equivalent
position-weight matrix with log2-lift-over-random scoring:

  donor window    : 3 exonic + 6 intronic nt around the 5' splice site
  acceptor window : 20 intronic + 3 exonic nt around the 3' splice site

An intron score is min(donor_score, acceptor_score) — the weakest site
gates.  Canonical GT-AG sites score strongly positive; shifted or
mis-assembled junctions score near zero or negative, which makes the gate
threshold species-agnostic (see _calibrate for the default).

Windows are always scored in transcription orientation; minus-strand
introns are reverse-complemented before scoring.
"""
from __future__ import annotations
import logging
log = logging.getLogger("annomerge")

DON_EX, DON_IN = 3, 6          # donor: 3 exonic + 6 intronic
ACC_IN, ACC_EX = 20, 3         # acceptor: 20 intronic + 3 exonic
PSEUDO = 0.5                   # pseudocount per base per position
_COMP = str.maketrans("ACGT", "TGCA")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _window(genome: dict, seqid: str, start0: int, end0: int) -> str | None:
    """1/2-open genomic window, clamped at scaffold bounds; None if empty."""
    seq = genome.get(seqid)
    if seq is None:
        return None
    s, e = max(0, start0), min(len(seq), end0)
    if s >= e:
        return None
    return seq[s:e].upper()


def donor_acceptor_seqs(genome: dict, seqid: str, strand: str,
                        i_s: int, i_e: int) -> tuple[str | None, str | None]:
    """Transcription-oriented (donor, acceptor) window sequences for a
    1-based inclusive intron [i_s, i_e]."""
    g0, g1 = i_s - 1, i_e            # 0-based half-open
    if strand != "-":
        don = _window(genome, seqid, g0 - DON_EX, g0 + DON_IN)
        acc = _window(genome, seqid, g1 - ACC_IN, g1 + ACC_EX)
    else:
        # transcription runs high->low: donor at the genomic high end
        don = _window(genome, seqid, g1 - DON_IN, g1 + DON_EX)
        acc = _window(genome, seqid, g0 - ACC_EX, g0 + ACC_IN)
        if don is not None:
            don = _revcomp(don)
        if acc is not None:
            acc = _revcomp(acc)
    return don, acc


class SplicePWM:
    """Position-weight matrices trained on observed introns."""

    def __init__(self):
        self.don_counts: list[dict[str, float]] = [dict.fromkeys("ACGT", 0.0)
                                                   for _ in range(DON_EX + DON_IN)]
        self.acc_counts: list[dict[str, float]] = [dict.fromkeys("ACGT", 0.0)
                                                   for _ in range(ACC_IN + ACC_EX)]
        self.n = 0

    # ---------------------------------------------------------------- train
    def add_intron(self, genome: dict, seqid: str, strand: str,
                   i_s: int, i_e: int) -> float | None:
        """Train on one intron; returns its score (None if windows clipped)."""
        don, acc = donor_acceptor_seqs(genome, seqid, strand, i_s, i_e)
        if don is None or acc is None or len(don) < DON_EX + DON_IN \
                or len(acc) < ACC_IN + ACC_EX:
            return None
        self.n += 1
        for i, b in enumerate(don):
            self.don_counts[i][b] += 1.0
        for i, b in enumerate(acc):
            self.acc_counts[i][b] += 1.0
        return self.score_windows(don, acc)

    # ---------------------------------------------------------------- score
    def _site_score(self, counts: list[dict[str, float]], seq: str) -> float:
        n = self.n
        if n == 0:
            return 0.0
        total = {}
        for pos in counts:
            for b in "ACGT":
                total[b] = total.get(b, 0.0) + pos[b]
        # position-specific log2 lift over random (0.25)
        sc = 0.0
        for i, b in enumerate(seq):
            if i >= len(counts) or b not in "ACGT":
                continue
            p = (counts[i][b] + PSEUDO) / (sum(counts[i].values()) + 4 * PSEUDO)
            sc += __import__("math").log2(p / 0.25)
        return sc

    def score_windows(self, don: str, acc: str) -> float:
        return min(self._site_score(self.don_counts, don),
                   self._site_score(self.acc_counts, acc))

    def score_intron(self, genome: dict, seqid: str, strand: str,
                     i_s: int, i_e: int) -> float | None:
        don, acc = donor_acceptor_seqs(genome, seqid, strand, i_s, i_e)
        if don is None or acc is None or len(don) < DON_EX + DON_IN \
                or len(acc) < ACC_IN + ACC_EX:
            return None
        return self.score_windows(don, acc)


def train_from_transcripts(genome: dict, transcripts: list) -> SplicePWM:
    """Train on a list of Tx-like objects (need .exons with start/end and
    .mrna.seqid/.mrna.strand).  Introns are derived from exon chains."""
    pwm = SplicePWM()
    n_in, n_scored = 0, 0
    for tx in transcripts:
        ex = sorted(tx.exons, key=lambda x: x.start)
        for i in range(len(ex) - 1):
            n_in += 1
            if pwm.add_intron(genome, tx.mrna.seqid, tx.mrna.strand,
                              ex[i].end + 1, ex[i + 1].start - 1) is not None:
                n_scored += 1
    log.info("splice PWM trained on %d introns (%d scored)", n_in, n_scored)
    return pwm


def _calibrate(pwm: SplicePWM, genome: dict,
               introns: list[tuple[str, str, int, int]]) -> dict:
    """Score observed introns and their +1bp-shifted versions; report the
    separation (used to justify the default gate threshold)."""
    import math
    real, shifted = [], []
    for seqid, strand, s, e in introns:
        r = pwm.score_intron(genome, seqid, strand, s, e)
        if r is not None:
            real.append(r)
        # shift donor+acceptor by +1bp: same intron length, offset sites
        sh = pwm.score_intron(genome, strand and seqid, strand, s + 1, e + 1) \
            if strand != "-" else \
            pwm.score_intron(genome, seqid, strand, s - 1, e - 1)
        if sh is not None:
            shifted.append(sh)
    def q(v, p):
        v = sorted(v)
        return v[min(len(v) - 1, max(0, int(p * len(v))))] if v else float("nan")
    return {"n": len(real), "real_min": q(real, 0.0), "real_p5": q(real, 0.05),
            "real_med": q(real, 0.5), "shift_p50": q(shifted, 0.5),
            "shift_p95": q(shifted, 0.95)}


# Default gate: a non-anchor intron must score strictly above this.
# Calibrated on S. cerevisiae deep data (scripts/calibrate_splice_pwm.py):
# 413 real introns vs their 1bp-shifted versions — shifted sites never
# scored above -5; the real-intron tail below -3 (2.2%) is exactly the
# suspicious/non-canonical set the gate exists to drop.
DEFAULT_GATE = -3.0
