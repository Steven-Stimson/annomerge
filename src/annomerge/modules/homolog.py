"""Module 2 - homology branch: miniprot -> PAF parsing -> referee set (wi = cov x ident).

Outputs (under <outdir>/homolog/):
  miniprot.paf / miniprot.gff   raw miniprot alignments
  alignments.tsv                per-protein quality table (query, species, cov, ident, wi, ...)
  referee.tsv / referee.ids     proteins passing wi >= --referee-min-wi (default 0.6)
  species_stats.tsv             per-species alignment quality summary
"""
from __future__ import annotations

import logging
import re
from collections import defaultdict
from pathlib import Path

from ..utils import ensure_dir, lazy, run, which
from .validate import ValidatedInputs

log = logging.getLogger("annomerge")

SPECIES_RE = re.compile(r"\|species=([^\s|]+)")


def run_miniprot(genome: Path, proteins: Path, out_dir: Path, threads: int) -> tuple[Path, Path]:
    """Run miniprot once, producing PAF and GFF."""
    ensure_dir(out_dir)
    paf = out_dir / "miniprot.paf"
    gff = out_dir / "miniprot.gff"
    if not lazy([paf, gff], desc="miniprot"):
        run(f"{which('miniprot')} -t{threads} --paf {genome} {proteins} > {paf}",
            shell=True, log_file=out_dir / "miniprot.paf.log", desc="miniprot (paf)")
        run(f"{which('miniprot')} -t{threads} --gff {genome} {proteins} > {gff}",
            shell=True, log_file=out_dir / "miniprot.gff.log", desc="miniprot (gff)")
    return paf, gff


def parse_paf(paf: Path) -> list[dict]:
    """Parse miniprot PAF.

    cov   = aligned query length / query length   (PAF cols 2,3,4,7,8)
    ident = residue matches / alignment block length (PAF cols 10,11)
    wi    = cov * ident
    """
    rows = []
    with open(paf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
            tname, ts, te = f[5], int(f[7]), int(f[8])
            nmatch, alen = int(f[9]), int(f[10])
            mapq = int(f[11])
            cov = (qe - qs) / qlen if qlen else 0.0
            ident = nmatch / alen if alen else 0.0
            m = SPECIES_RE.search(qname)
            rows.append({
                "query": qname, "species": m.group(1) if m else "NA",
                "qlen": qlen, "qcov": round(cov, 4), "ident": round(ident, 4),
                "wi": round(cov * ident, 4), "target": tname,
                "tstart": ts, "tend": te, "mapq": mapq,
            })
    return rows


def best_per_query(rows: list[dict]) -> dict[str, dict]:
    """Keep the highest-wi alignment per query protein."""
    best: dict[str, dict] = {}
    for r in rows:
        if r["query"] not in best or r["wi"] > best[r["query"]]["wi"]:
            best[r["query"]] = r
    return best


def run_homolog_branch(vi: ValidatedInputs, out_dir: Path, threads: int,
                       referee_min_wi: float = 0.6) -> dict:
    hdir = ensure_dir(out_dir / "03.homolog")
    log.info("module 2 [homolog]: %s", vi.homologs)

    paf, gff = run_miniprot(vi.genome, vi.homologs, hdir, threads)

    aln_tsv = hdir / "alignments.tsv"
    ref_tsv = hdir / "referee.tsv"
    ref_ids = hdir / "referee.ids"
    stats_tsv = hdir / "species_stats.tsv"

    if lazy([aln_tsv, ref_tsv, ref_ids, stats_tsv], desc="homolog PAF parsing"):
        n_ref = sum(1 for _ in open(ref_ids))
        return {"paf": paf, "gff": gff, "alignments": aln_tsv,
                "referee_tsv": ref_tsv, "referee_ids": ref_ids,
                "species_stats": stats_tsv, "n_referee": n_ref}

    rows = parse_paf(paf)
    best = best_per_query(rows)

    with open(aln_tsv, "w") as fh:
        fh.write("query\tspecies\tqlen\tqcov\tident\twi\ttarget\ttstart\ttend\tmapq\n")
        for r in sorted(best.values(), key=lambda r: -r["wi"]):
            fh.write("\t".join(map(str, [r["query"], r["species"], r["qlen"], r["qcov"],
                                         r["ident"], r["wi"], r["target"],
                                         r["tstart"], r["tend"], r["mapq"]])) + "\n")

    referee = [r for r in best.values() if r["wi"] >= referee_min_wi]
    with open(ref_tsv, "w") as fh:
        fh.write("query\tspecies\tqcov\tident\twi\n")
        for r in sorted(referee, key=lambda r: -r["wi"]):
            fh.write(f"{r['query']}\t{r['species']}\t{r['qcov']}\t{r['ident']}\t{r['wi']}\n")
    ref_ids.write_text("\n".join(r["query"] for r in referee) + "\n")

    by_species: dict[str, list[dict]] = defaultdict(list)
    for r in best.values():
        by_species[r["species"]].append(r)
    with open(stats_tsv, "w") as fh:
        fh.write("species\tn_aligned\tn_referee\tmean_wi\tmedian_wi\n")
        for sp, rs in sorted(by_species.items(), key=lambda kv: -len(kv[1])):
            wis = sorted(r["wi"] for r in rs)
            n = len(wis)
            med = wis[n // 2] if n % 2 else (wis[n // 2 - 1] + wis[n // 2]) / 2
            n_ref = sum(1 for r in rs if r["wi"] >= referee_min_wi)
            fh.write(f"{sp}\t{n}\t{n_ref}\t{sum(wis)/n:.4f}\t{med:.4f}\n")

    log.info("module 2 [homolog] done: %d aligned proteins, %d referee (wi>=%.2f), %d species",
             len(best), len(referee), referee_min_wi, len(by_species))
    return {"paf": paf, "gff": gff, "alignments": aln_tsv, "referee_tsv": ref_tsv,
            "referee_ids": ref_ids, "species_stats": stats_tsv, "n_referee": len(referee)}
