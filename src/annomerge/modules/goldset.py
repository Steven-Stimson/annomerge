"""Module 2b - Gold set construction.

StringTie transcripts -> gffread sequences -> TransDecoder ORFs ->
DIAMOND validation against the homolog referee set -> goldset.faa.

A gold transcript carries BOTH RNA expression support AND homology-verifiable
protein-coding capacity; used as a high-confidence, species-matched referee
inside the scoring system only (never as a held-out validation metric).

Outputs (under <outdir>/goldset/):
  transcripts.fa              transcript sequences (gffread)
  transdecoder.longest_orfs.pep
  blastp.outfmt6              DIAMOND blastp vs homolog referee
  transdecoder.pep            TransDecoder final predictions (homology-retained)
  goldset.faa                 final gold protein set
  goldset_details.tsv         per-transcript ORF/homology evidence table
"""
from __future__ import annotations

import logging
from pathlib import Path

from ..utils import ensure_dir, lazy, run, which
from .validate import ValidatedInputs

log = logging.getLogger("annomerge")


def _read_fasta_ids(path: Path) -> dict[str, str]:
    ids = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                ids[name] = line[1:].rstrip("\n")
    return ids


# TransDecoder.Predict-only flags: passed to Predict but not LongOrfs
TD_PREDICT_ONLY = ("--single_best_only", "--single_best_ORF",
                   "--retain_long_orfs_length", "--retain_long_orfs_mode",
                   "--no_refine_starts")


def split_td_opts(opts: list[str] | None) -> tuple[list[str], list[str]]:
    """Split user TransDecoder options into (longorfs_extra, predict_extra).

    Flags in TD_PREDICT_ONLY go to Predict only; everything else goes to
    both stages (e.g. -S strand-specific, -G genetic code).
    """
    both, predict_only = [], []
    for t in opts or []:
        base = t.split("=", 1)[0]
        (predict_only if base in TD_PREDICT_ONLY else both).append(t)
    return both, predict_only


def run_goldset(vi: ValidatedInputs, out_dir: Path, threads: int,
                stringtie_gtf: Path | None, referee_faa: Path | None,
                min_orf_aa: int = 50,
                td_opts: list[str] | None = None) -> dict | None:
    lo_opts, pr_opts = split_td_opts(td_opts)
    """Build the gold set. Requires RNA transcripts AND homolog evidence."""
    if stringtie_gtf is None:
        log.warning("goldset: no RNA transcripts available -> skipped (auto-degraded)")
        return None
    gdir = ensure_dir(out_dir / "04.goldset")
    log.info("module 2b [goldset]: transcripts=%s", stringtie_gtf)

    # 1. transcript sequences
    tfa = gdir / "transcripts.fa"
    if not lazy([tfa], desc="gffread -w"):
        run([which("gffread"), "-w", str(tfa), "-g", str(vi.genome), str(stringtie_gtf)],
            desc="gffread extract transcripts")

    # 2. TransDecoder.LongOrfs
    td_dir = gdir / "transdecoder"
    orfs = gdir / "transcripts.fa.transdecoder_dir" / "longest_orfs.pep"
    if not lazy([orfs], desc="TransDecoder.LongOrfs"):
        run(f"cd {gdir} && {which('TransDecoder.LongOrfs')} -t {tfa} -m {min_orf_aa} "
            + " ".join(lo_opts),
            shell=True, log_file=gdir / "longorfs.log", desc="TransDecoder.LongOrfs")

    # 3. DIAMOND blastp: candidate ORFs vs homolog referee (fall back to full homologs)
    ref = referee_faa if referee_faa else vi.homologs
    dmnd_db = gdir / "referee.dmnd"
    blast_out = gdir / "blastp.outfmt6"
    if not lazy([blast_out], desc="diamond blastp"):
        if not lazy([Path(str(dmb) + ".dmnd") if (dmb := dmnd_db) else None], desc="diamond makedb"):
            run([which("diamond"), "makedb", "--in", str(ref), "-d", str(dmnd_db)],
                log_file=gdir / "diamond_makedb.log", desc="diamond makedb")
        run([which("diamond"), "blastp", "-q", str(orfs), "-d", str(dmnd_db),
             "-p", str(threads), "--outfmt", "6", "--max-target-seqs", "5",
             "--evalue", "1e-5", "-o", str(blast_out)],
            log_file=gdir / "diamond_blastp.log", desc="diamond blastp")

    # 4. TransDecoder.Predict with homology retention
    final_pep = gdir / "transcripts.fa.transdecoder.pep"
    if not lazy([final_pep], desc="TransDecoder.Predict"):
        run(f"cd {gdir} && {which('TransDecoder.Predict')} -t {tfa} "
            f"--retain_blastp_hits {blast_out} "
            + " ".join(lo_opts + pr_opts),
            shell=True, log_file=gdir / "predict.log", desc="TransDecoder.Predict")

    # 5. goldset.faa = TransDecoder final pep; details table
    gold = gdir / "goldset.faa"
    details = gdir / "goldset_details.tsv"
    if not lazy([gold, details], desc="goldset finalize"):
        pep_ids = _read_fasta_ids(final_pep)
        # homology hits per ORF
        hits: dict[str, int] = {}
        with open(blast_out) as fh:
            for line in fh:
                q = line.split("\t", 1)[0]
                hits[q] = hits.get(q, 0) + 1
        with open(gold, "w") as out, open(details, "w") as det:
            det.write("orf_id\thas_homology_hit\tn_blastp_hits\n")
            seqs = {}
            name, buf = None, []
            with open(final_pep) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if name:
                            seqs[name] = "".join(buf)
                        name = line[1:].split()[0]
                        buf = []
                    else:
                        buf.append(line.strip())
                if name:
                    seqs[name] = "".join(buf)
            for name, seq in seqs.items():
                has_hit = name in hits
                det.write(f"{name}\t{has_hit}\t{hits.get(name, 0)}\n")
                out.write(f">{name} homology_validated={has_hit}\n{seq}\n")
    n_gold = sum(1 for line in open(gold) if line.startswith(">"))
    log.info("module 2b [goldset] done: %d gold proteins", n_gold)
    return {"goldset_faa": gold, "details": details, "n_gold": n_gold}
