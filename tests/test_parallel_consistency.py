#!/usr/bin/env python3
"""Parallel-consistency + Iso-Seq platform-layer smoke test for augment.py.

Builds a fully synthetic two-contig dataset (EVM anchors, goldset donors,
StringTie GTF, Iso-Seq GTF, one de novo predictor) and asserts:

  1. threads=1 vs threads=8 produce BYTE-IDENTICAL augment.gff3 and
     augment_stats.tsv (fork-parallel Phase 4 / backfill / platform flagging
     must not change results), for geometric x evidence mode, each with and
     without an Iso-Seq pool.
  2. With Iso-Seq: screening demotes the internal-priming artifact, tier0
     backfill fires (utr_from_isoseq), AS FSM support is reported.
  3. Iso-Seq on/off outputs differ (the platform layer actually does work).

Run:  python3 tests/test_parallel_consistency.py
"""
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from annomerge.modules.augment import run_augment  # noqa: E402


# ------------------------------------------------------------ synthetic data
def _build_genome(fa: Path):
    import random
    rng = random.Random(42)
    chr1 = list("".join(rng.choice("ACGT") for _ in range(40000)))
    chr2 = list("".join(rng.choice("ACGT") for _ in range(35000)))

    def force(seq, pos1, s):                     # 1-based start
        seq[pos1 - 1:pos1 - 1 + len(s)] = list(s)

    def write_cds(seq, segs):                    # valid ORF across the segments
        total = sum(e - s + 1 for s, e in segs)
        assert total % 3 == 0, total
        codons = "".join(["ATG"] + ["GCT"] * (total // 3 - 2) + ["TAA"])
        i = 0
        for s, e in segs:
            force(seq, s, codons[i:i + (e - s + 1)])
            i += e - s + 1

    # CDS first (splice-site forcing wins on overlaps; checked: no stops)
    write_cds(chr1, [(10001, 10200), (10501, 10702)])   # g1 anchor / d1-d3
    write_cds(chr1, [(10001, 10300), (10601, 10702)])   # d4 AS ORF
    write_cds(chr1, [(30051, 30200), (30501, 30752)])   # g3 (UTR-less)
    write_cds(chr2, [(20001, 20200), (20501, 20702)])   # g2 / d5
    write_cds(chr2, [(25001, 25402)])                   # g4 (UTR-less, 1-exon)

    def splice(seq, a, b):                       # intron a..b (1-based)
        force(seq, a, "GT")
        force(seq, b - 1, "AG")

    for a, b in [(10201, 10500), (9701, 9800), (10301, 10600),
                 (30201, 30500)]:
        splice(chr1, a, b)
    splice(chr2, 20201, 20500)
    # internal-priming trap: A-rich window right after iso.2's claimed TES
    force(chr2, 30501, "AAAAAAAAAAAAAAAAAAAA")

    with open(fa, "w") as fh:
        fh.write(">chr1\n" + "".join(chr1) + "\n")
        fh.write(">chr2\n" + "".join(chr2) + "\n")


def _g(ftype, s, e, attrs, seqid, strand="+", src="gold", phase="."):
    return (f"{seqid}\t{src}\t{ftype}\t{s}\t{e}\t.\t{strand}\t{phase}\t"
            f"{attrs}\n")


def _mrna_block(tid, exons, cds, seqid, strand="+", src="gold"):
    out = [_g("mRNA", exons[0][0], exons[-1][1], f"ID={tid}", seqid,
              strand, src)]
    for s, e in exons:
        out.append(_g("exon", s, e, f"Parent={tid}", seqid, strand, src))
    for s, e in cds:
        out.append(_g("CDS", s, e, f"Parent={tid}", seqid, strand, src))
    return out


def _gtf(chains: dict, meta: dict) -> str:
    out = []
    for tid, (sq, exons) in chains.items():
        tpm, cov = meta.get(tid, (0.0, 1.0))
        attr = (f'gene_id "{tid.split(".")[0]}"; transcript_id "{tid}"; '
                f'cov "{cov}"; TPM "{tpm}";')
        s0, e0 = exons[0][0], exons[-1][1]
        out.append(f"{sq}\tStringTie\ttranscript\t{s0}\t{e0}\t.\t+\t.\t{attr}\n")
        for s, e in exons:
            out.append(f"{sq}\tStringTie\texon\t{s}\t{e}\t.\t+\t.\t{attr}\n")
    return "".join(out)


def _write_inputs(wd: Path):
    """Returns dict of all input paths."""
    fa = wd / "genome.fa"
    _build_genome(fa)

    # ---- EVM anchors
    evm = wd / "evm.gff3"
    lines = []
    anchors = {
        "evm.TU.1": ("chr1", [(9901, 10200), (10501, 10800)],
                     [(10001, 10200), (10501, 10702)]),
        "evm.TU.2": ("chr2", [(19801, 20200), (20501, 20800)],
                     [(20001, 20200), (20501, 20702)]),
        # UTR-less anchors (exons == CDS) to exercise backfill
        "evm.TU.3": ("chr1", [(30051, 30200), (30501, 30752)],
                     [(30051, 30200), (30501, 30752)]),
        "evm.TU.4": ("chr2", [(25001, 25402)], [(25001, 25402)]),
    }
    for tid, (sq, ex, cd) in anchors.items():
        lines += _mrna_block(tid, ex, cd, sq, src="evm")
    evm.write_text("##gff-version 3\n" + "".join(lines))

    # ---- goldset (donors WITH ORFs; key = MSTRG.x.y)
    gold = wd / "goldset.genome.gff3"
    g1cds = [(10001, 10200), (10501, 10702)]
    donors = {
        "MSTRG.1.1.p1": ("chr1", [(9601, 10200), (10501, 11000)], g1cds),
        "MSTRG.1.2.p1": ("chr1", [(9611, 10200), (10501, 10990)], g1cds),
        "MSTRG.1.3.p1": ("chr1", [(9401, 9700), (9801, 10200),
                                  (10501, 11000)], g1cds),
        "MSTRG.1.4.p1": ("chr1", [(9901, 10300), (10601, 10800)],
                         [(10001, 10300), (10601, 10702)]),      # AS donor
        "MSTRG.2.1.p1": ("chr2", [(19601, 20200), (20501, 21050)],
                         [(20001, 20200), (20501, 20702)]),
    }
    lines = []
    for tid, (sq, ex, cd) in donors.items():
        lines += _mrna_block(tid, ex, cd, sq)
    gold.write_text("##gff-version 3\n" + "".join(lines))

    # ---- stringtie gtf (TPM per donor + one ORF-less backfill donor)
    st = wd / "stringtie.gtf"
    tpms = {"MSTRG.1.1": (10.0, 5.0), "MSTRG.1.2": (0.5, 1.0),
            "MSTRG.1.3": (8.0, 4.0), "MSTRG.1.4": (5.0, 3.0),
            "MSTRG.2.1": (6.0, 3.0), "ST.100": (2.0, 1.5)}
    chains = {
        "MSTRG.1.1": ("chr1", [(9601, 10200), (10501, 11000)]),
        "MSTRG.1.2": ("chr1", [(9611, 10200), (10501, 10990)]),
        "MSTRG.1.3": ("chr1", [(9401, 9700), (9801, 10200),
                               (10501, 11000)]),
        "MSTRG.1.4": ("chr1", [(9901, 10300), (10601, 10800)]),
        "MSTRG.2.1": ("chr2", [(19601, 20200), (20501, 21050)]),
        "ST.100":    ("chr1", [(29801, 30200), (30501, 31000)]),
    }
    st.write_text(_gtf(chains, tpms))

    # ---- isoseq gtf
    iso = wd / "isoseq.stringtie.gtf"
    iso_chains = {
        "ISO.1": ("chr1", [(9501, 10200), (10501, 11100)]),   # FSM of d1/d2
        "ISO.2": ("chr2", [(30001, 30500)]),                  # priming trap
        "ISO.3": ("chr1", [(9750, 10300), (10601, 11050)]),   # FSM of d4
        "ISO.4": ("chr1", [(29501, 30200), (30501, 31200)]),  # tier0 for g3
    }
    iso_meta = {"ISO.1": (3.0, 5.0), "ISO.2": (1.0, 2.0),
                "ISO.3": (2.0, 3.0), "ISO.4": (2.0, 4.0)}
    iso.write_text(_gtf(iso_chains, iso_meta))

    # ---- one de novo predictor (tier2 donor for g4)
    dn = wd / "denovo1.gff3"
    dn.write_text("##gff-version 3\n" + "".join(
        _mrna_block("dn.1", [(24501, 25600)], [], "chr2", src="denovo1")))

    return {"genome": fa, "evm": evm, "gold": gold, "st": st,
            "iso": iso, "preds": {"denovo1": dn}}


# ------------------------------------------------------------------ harness
def _run(inputs, mode, isoseq, threads, outdir):
    vi = SimpleNamespace(genome=str(inputs["genome"]))
    return run_augment(
        vi, outdir, inputs["evm"], inputs["gold"], inputs["st"],
        preds=inputs["preds"], scores=None,
        dedup_mode=mode, junctions_bed=None,
        isoseq_gtf=inputs["iso"] if isoseq else None,
        isoseq_tol_bp=25, isoseq_priority=True,
        coverage_bams=None, isoseq_fl_min=2.0, isoseq_cov_cliff=0.5,
        threads=threads)


def _read(outdir):
    d = Path(outdir)
    return ((d / "augment.gff3").read_bytes(),
            (d / "augment_stats.tsv").read_text())


def main():
    failures = []

    def check(name, cond, extra=""):
        print(("PASS" if cond else "FAIL"), name, extra)
        if not cond:
            failures.append(name)

    with tempfile.TemporaryDirectory() as tmp:
        inputs = _write_inputs(Path(tmp))
        for mode in ("geometric", "evidence"):
            for isoseq in (False, True):
                tag = f"{mode}/isoseq={'on' if isoseq else 'off'}"
                o1 = Path(tmp) / f"o_{mode}_{int(isoseq)}_t1"
                o8 = Path(tmp) / f"o_{mode}_{int(isoseq)}_t8"
                _run(inputs, mode, isoseq, 1, o1)
                _run(inputs, mode, isoseq, 8, o8)
                g1, s1 = _read(o1)
                g8, s8 = _read(o8)
                check(f"{tag}: gff3 byte-identical t1 vs t8", g1 == g8)
                check(f"{tag}: stats identical t1 vs t8", s1 == s8)
                if isoseq:
                    st = dict(l.split("\t") for l in s1.strip().splitlines()[1:]
                              if "\t" in l)
                    check(f"{tag}: priming artifact demoted",
                          int(st.get("isoseq_demoted_priming", 0)) >= 1,
                          st.get("isoseq_demoted_priming", "0"))
                    check(f"{tag}: tier0 backfill fired",
                          int(st.get("utr_from_isoseq", 0)) >= 1)
                    check(f"{tag}: AS FSM support reported",
                          int(st.get("as_fsm_supported", 0)) >= 1)
                    check(f"{tag}: platform-flagged goldset groups",
                          int(st.get("isoseq_goldset_groups", 0)) >= 2,
                          st.get("isoseq_goldset_groups", "0"))
        # platform layer must change the outcome
        g_off, _ = _read(Path(tmp) / "o_evidence_0_t8")
        g_on, _ = _read(Path(tmp) / "o_evidence_1_t8")
        check("isoseq on != off (platform layer effective)", g_off != g_on)

    print()
    if failures:
        print(f"{len(failures)} FAILURE(S)")
        sys.exit(1)
    print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
