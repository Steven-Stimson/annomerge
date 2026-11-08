"""Module 2 - homology branch: miniprot -> PAF parsing -> referee set (wi = cov x ident).
Supports multiple labeled homolog sets (-p label=path,...): each set is run
through miniprot separately, then alignments are merged transparently.
Per-set relatedness (median identity x median coverage over mapped proteins)
drives downstream weights (EVM protein evidence, H_s/C_s, confidence).
Outputs (under <outdir>/03.homolog/):
  miniprot.paf / miniprot.gff   raw miniprot alignments (all sets merged)
  alignments.tsv                per-protein quality table (query, set, species, cov, ident, wi, ...)
  referee.tsv / referee.ids     proteins passing wi >= --referee-min-wi (default 0.6)
  species_stats.tsv             per-species alignment quality summary
  sets.tsv                      per-set relatedness + weight table
  sets/<label>/                 per-set miniprot outputs
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


def _load_membership(vi: ValidatedInputs) -> dict[str, str]:
    """protein id -> set label. Single-set runs: all proteins map to the only label."""
    if vi.homolog_membership and Path(vi.homolog_membership).exists():
        m: dict[str, str] = {}
        with open(vi.homolog_membership) as fh:
            next(fh, None)  # header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    m[f[0]] = f[1]
        return m
    only = next(iter(vi.homolog_sets), None) if vi.homolog_sets else None
    return defaultdict(lambda: only or "homologs") if only else {}


def _median(xs: list[float]) -> float:
    xs = sorted(xs)
    n = len(xs)
    if not n:
        return 0.0
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def _run_sets_miniprot(vi: ValidatedInputs, hdir: Path, threads: int) -> tuple[Path, Path]:
    """Per-set miniprot runs, merged transparently into hdir/miniprot.{paf,gff}."""
    sroot = ensure_dir(hdir / "sets")
    n_sets = len(vi.homolog_sets)
    t_each = max(1, threads // n_sets)
    set_outputs = []
    for label, path in vi.homolog_sets.items():
        set_outputs.append(run_miniprot(vi.genome, path, ensure_dir(sroot / label), t_each))
    paf = hdir / "miniprot.paf"
    gff = hdir / "miniprot.gff"
    if not lazy([paf, gff], desc="merge per-set miniprot"):
        with open(paf, "w") as out:
            for spaf, _ in set_outputs:
                with open(spaf) as fh:
                    for line in fh:
                        if not line.startswith("#"):
                            out.write(line)
        # per-set miniprot feature IDs share one namespace (MP0_1, ...); prefix
        # them with the set label so merged IDs stay unique downstream (EVM
        # protein evidence derives match IDs from these)
        def _refix(line: str, label: str) -> str:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9:
                a = re.sub(r"ID=([^;\s]+)", rf"ID={label}__\1", f[8])
                a = re.sub(r"Parent=([^;\s]+)", rf"Parent={label}__\1", a)
                f[8] = a
                line = "\t".join(f) + "\n"
            return line
        with open(gff, "w") as out:
            out.write("##gff-version 3\n")
            for (label, _), (_, sgff) in zip(vi.homolog_sets.items(), set_outputs):
                with open(sgff) as fh:
                    for line in fh:
                        if not line.startswith("#"):
                            out.write(_refix(line, label))
    return paf, gff


def _relatedness(vi: ValidatedInputs, hdir: Path) -> tuple[Path, dict[str, float], dict[str, float]]:
    """Per-set relatedness table. Returns (sets_tsv, rel, weight).

    relatedness_k = median identity x median coverage over the set's mapped
    proteins (best alignment per protein). With >1 sets:
        weight_k = clamp(2 * rel_k / max_rel, 0.5, 2.0)
    Single set -> weight 1.0 (neutral, legacy behavior).
    """
    sets_tsv = hdir / "sets.tsv"
    if lazy([sets_tsv], desc="homolog set relatedness"):
        rel, weight = {}, {}
        with open(sets_tsv) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 7:
                    rel[f[0]] = float(f[5])
                    weight[f[0]] = float(f[6])
        return sets_tsv, rel, weight
    membership = _load_membership(vi)
    per_set_rows: dict[str, list[dict]] = {k: [] for k in vi.homolog_sets}
    only = next(iter(vi.homolog_sets)) if len(vi.homolog_sets) == 1 else None
    for label in vi.homolog_sets:
        spaf = hdir / "sets" / label / "miniprot.paf"
        if spaf.exists():
            per_set_rows[label] = list(best_per_query(parse_paf(spaf)).values())
        elif only:
            per_set_rows[label] = []
    if only and not per_set_rows[only]:  # single-set: parse merged paf
        per_set_rows[only] = list(best_per_query(parse_paf(hdir / "miniprot.paf")).values())
    rel: dict[str, float] = {}
    for label, rows in per_set_rows.items():
        rel[label] = round(_median([r["ident"] for r in rows]) *
                           _median([r["qcov"] for r in rows]), 6) if rows else 0.0
    if len(vi.homolog_sets) > 1:
        max_rel = max(rel.values()) or 1.0
        weight = {k: min(2.0, max(0.5, 2.0 * r / max_rel)) for k, r in rel.items()}
    else:
        weight = {k: 1.0 for k in rel}
    with open(sets_tsv, "w") as fh:
        fh.write("set\tn_proteins\tn_mapped\tmedian_ident\tmedian_qcov\trelatedness\tweight\n")
        for label, rows in per_set_rows.items():
            n_prot = 0
            with open(vi.homolog_sets[label]) as f2:
                n_prot = sum(1 for l in f2 if l.startswith(">"))
            fh.write(f"{label}\t{n_prot}\t{len(rows)}\t"
                     f"{_median([r['ident'] for r in rows]):.4f}\t"
                     f"{_median([r['qcov'] for r in rows]):.4f}\t"
                     f"{rel[label]:.6f}\t{weight[label]:.4f}\n")
    log.info("homolog set relatedness: %s",
             ", ".join(f"{k}: rel={rel[k]:.3f} w={weight[k]:.2f}" for k in rel))
    return sets_tsv, rel, weight


def run_homolog_branch(vi: ValidatedInputs, out_dir: Path, threads: int,
                       referee_min_wi: float = 0.6) -> dict:
    hdir = ensure_dir(out_dir / "03.homolog")
    log.info("module 2 [homolog]: %s (%d set(s))", vi.homologs, len(vi.homolog_sets))
    membership = _load_membership(vi)

    if len(vi.homolog_sets) > 1:
        paf, gff = _run_sets_miniprot(vi, hdir, threads)
    else:
        paf, gff = run_miniprot(vi.genome, vi.homologs, hdir, threads)

    sets_tsv, rel, weight = _relatedness(vi, hdir)

    aln_tsv = hdir / "alignments.tsv"
    ref_tsv = hdir / "referee.tsv"
    ref_ids = hdir / "referee.ids"
    stats_tsv = hdir / "species_stats.tsv"

    if lazy([aln_tsv, ref_tsv, ref_ids, stats_tsv], desc="homolog PAF parsing"):
        n_ref = sum(1 for _ in open(ref_ids))
        return {"paf": paf, "gff": gff, "alignments": aln_tsv,
                "referee_tsv": ref_tsv, "referee_ids": ref_ids,
                "species_stats": stats_tsv, "n_referee": n_ref,
                "sets_tsv": sets_tsv, "set_rel": rel, "set_weights": weight,
                "membership": membership}

    rows = parse_paf(paf)
    best = best_per_query(rows)

    with open(aln_tsv, "w") as fh:
        fh.write("query\tset\tspecies\tqlen\tqcov\tident\twi\ttarget\ttstart\ttend\tmapq\n")
        for r in sorted(best.values(), key=lambda r: -r["wi"]):
            fh.write("\t".join(map(str, [r["query"], membership.get(r["query"], "NA"),
                                         r["species"], r["qlen"], r["qcov"],
                                         r["ident"], r["wi"], r["target"],
                                         r["tstart"], r["tend"], r["mapq"]])) + "\n")

    referee = [r for r in best.values() if r["wi"] >= referee_min_wi]
    with open(ref_tsv, "w") as fh:
        fh.write("query\tset\tspecies\tqcov\tident\twi\n")
        for r in sorted(referee, key=lambda r: -r["wi"]):
            fh.write(f"{r['query']}\t{membership.get(r['query'], 'NA')}\t{r['species']}"
                     f"\t{r['qcov']}\t{r['ident']}\t{r['wi']}\n")
    ref_ids.write_text("\n".join(r["query"] for r in referee) + "\n")

    by_species: dict[str, list[dict]] = defaultdict(list)
    for r in best.values():
        by_species[r["species"]].append(r)
    with open(stats_tsv, "w") as fh:
        fh.write("species\tn_aligned\tn_referee\tmean_wi\tmedian_wi\n")
        for sp, rs in sorted(by_species.items(), key=lambda kv: -len(kv[1])):
            wis = sorted(r["wi"] for r in rs)
            n = len(wis)
            med = _median(wis)
            n_ref = sum(1 for r in rs if r["wi"] >= referee_min_wi)
            fh.write(f"{sp}\t{n}\t{n_ref}\t{sum(wis)/n:.4f}\t{med:.4f}\n")

    log.info("module 2 [homolog] done: %d aligned proteins, %d referee (wi>=%.2f), %d species",
             len(best), len(referee), referee_min_wi, len(by_species))
    return {"paf": paf, "gff": gff, "alignments": aln_tsv, "referee_tsv": ref_tsv,
            "referee_ids": ref_ids, "species_stats": stats_tsv, "n_referee": len(referee),
            "sets_tsv": sets_tsv, "set_rel": rel, "set_weights": weight,
            "membership": membership}
