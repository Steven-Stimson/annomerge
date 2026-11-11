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

from ..utils import ToolError, ensure_dir, lazy, run, which
from .validate import ValidatedInputs

log = logging.getLogger("annomerge")

SPECIES_RE = re.compile(r"\|species=([^\s|]+)")


REFEREE_MIN_WI = 0.6  # referee set threshold: wi = coverage x identity

# ---- miniprot alignment quality filters (vocabulary from mnpFormatter.pl) --
MNP_TOKENS = ("UnStop", "UnShift", "minIden", "minPos", "minSize",
              "maxSize", "maxRank", "R2Q", "Q2R", "topRank")


def parse_mnp_filter(spec: str) -> dict:
    """'--mnp-filter "UnStop,UnShift,Q2R=0.5,topRank=1"' -> filter dict.

    Token names mirror mnpFormatter.pl: UnStop/UnShift are switches (present
    = drop in-frame-stop / frameshift alignments); numeric thresholds default
    to 0 = disabled.  maxRank keeps only Rank<=N mappings (0 = off);
    topRank keeps the best N mappings per protein (0 = off); R2Q/Q2R bound
    the ratio reference-CDS-length / aligned-CDS-length from each side
    (reference lengths come from the homolog FASTA itself).
    """
    f = {"UnStop": False, "UnShift": False, "minIden": 0.0, "minPos": 0.0,
         "minSize": 0, "maxSize": 10**9, "maxRank": 0, "R2Q": 0.0,
         "Q2R": 0.0, "topRank": 0}
    for tok in (spec or "").split(","):
        tok = tok.strip()
        if not tok:
            continue
        key, _, val = tok.partition("=")
        key = key.strip()
        if key in ("UnStop", "UnShift"):
            f[key] = True
        elif key in f:
            f[key] = (float(val) if key in ("minIden", "minPos", "R2Q", "Q2R")
                      else int(float(val)))
        else:
            raise ValueError(
                f"unknown --mnp-filter token '{key}' "
                f"(allowed: {', '.join(MNP_TOKENS)})")
    return f


def _fasta_lens(path: Path | None) -> dict:
    lens: dict = {}
    if not path or not Path(path).exists():
        return lens
    name, n = None, 0
    for line in open(path):
        if line.startswith(">"):
            if name:
                lens[name] = n
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
    if name:
        lens[name] = n
    return lens


def filter_miniprot_alignments(raw_gff: Path, raw_paf: Path,
                               out_gff: Path, out_paf: Path,
                               ref_fa: Path | None, filt: dict) -> dict:
    """Apply mnp quality filters to one miniprot output pair (GFF+PAF).

    GFF records (mRNA + CDS/stop_codon children) are filtered on Rank /
    Identity / Positive / span / StopCodon / Frameshift / CDS-coverage
    ratios; topRank keeps the best N mappings per protein afterwards.
    PAF lines are dropped in sync (keyed by query+target+tstart).  Returns
    per-reason drop counts plus 'kept'."""
    lens = _fasta_lens(ref_fa)
    counts: dict = {r: 0 for r in
                    ("stop", "shift", "minIden", "minPos", "minSize",
                     "maxSize", "maxRank", "R2Q", "Q2R", "topRank", "kept")}
    records: list = []
    cur: dict | None = None
    with open(raw_gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                if cur:
                    records.append(cur)
                attrs = {}
                for kv in f[8].split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        attrs[k.strip()] = v.strip()
                cur = {"lines": [line], "attrs": attrs, "seqid": f[0],
                       "start": int(f[3]), "end": int(f[4]), "cds": 0}
            elif cur is not None and f[2] in ("CDS", "stop_codon"):
                cur["lines"].append(line)
                if f[2] == "CDS":
                    cur["cds"] += int(f[4]) - int(f[3]) + 1
    if cur:
        records.append(cur)

    def _metrics(rec):
        a = rec["attrs"]
        rank = int(a.get("Rank", 0) or 0)
        ident = float(a.get("Identity", 0) or 0)
        target = (a.get("Target", "") or "").split()[0]
        return rank, ident, target

    kept: list = []
    for rec in records:
        a = rec["attrs"]
        rank, ident, _t = _metrics(rec)
        pos = float(a.get("Positive", 0) or 0)
        stop = a.get("StopCodon", "0") not in ("", "0")
        shift = a.get("Frameshift", "0") not in ("", "0")
        span = rec["end"] - rec["start"] + 1
        target = (a.get("Target", "") or "").split()[0]
        reason = None
        if filt["UnStop"] and stop:
            reason = "stop"
        elif filt["UnShift"] and shift:
            reason = "shift"
        elif ident < filt["minIden"]:
            reason = "minIden"
        elif pos < filt["minPos"]:
            reason = "minPos"
        elif span < filt["minSize"]:
            reason = "minSize"
        elif span > filt["maxSize"]:
            reason = "maxSize"
        elif filt["maxRank"] and rank > filt["maxRank"]:
            reason = "maxRank"
        else:
            rl = lens.get(target)
            if rl and rec["cds"]:
                q2r = rec["cds"] / (3 * rl)
                r2q = (3 * rl) / rec["cds"]
                if filt["Q2R"] and q2r < filt["Q2R"]:
                    reason = "Q2R"
                elif filt["R2Q"] and r2q < filt["R2Q"]:
                    reason = "R2Q"
        if reason:
            counts[reason] += 1
        else:
            kept.append(rec)
    if filt["topRank"]:
        by_tgt: dict = {}
        for rec in kept:
            by_tgt.setdefault(_metrics(rec)[2], []).append(rec)
        kept2: list = []
        for tgt, rs in by_tgt.items():
            rs.sort(key=lambda r: (_metrics(r)[0], -_metrics(r)[1]))
            kept2.extend(rs[:filt["topRank"]])
            counts["topRank"] += max(0, len(rs) - filt["topRank"])
        kept = kept2
    counts["kept"] = len(kept)
    # PAF sync: miniprot PAF excludes the stop codon from its coords, so the
    # GFF mRNA start can differ from the PAF tstart by up to 3 bp (stop sits
    # at the transcription-start side).  Index kept starts per (query, chr)
    # and accept a match within 3 bp.
    from bisect import bisect_left
    groups: dict = {}
    for r in kept:
        groups.setdefault((_metrics(r)[2], r["seqid"]), []).append(r["start"] - 1)
    for g in groups.values():
        g.sort()

    def _paf_kept(c):
        g = groups.get((c[0], c[5]))
        if not g:
            return False
        t = int(c[7])
        i = bisect_left(g, t)
        for j in (i - 1, i):
            if 0 <= j < len(g) and abs(g[j] - t) <= 3:
                return True
        return False

    with open(out_gff, "w") as out:
        out.write("##gff-version 3\n")
        for rec in kept:
            out.writelines(rec["lines"])
    with open(raw_paf) as fh, open(out_paf, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) >= 9 and _paf_kept(c):
                out.write(line)
    return counts


def run_miniprot(genome: Path, proteins: Path, out_dir: Path, threads: int,
                 filt: dict | None = None,
                 extra_opts: list[str] | None = None) -> tuple[Path, Path]:
    """Run miniprot once, producing PAF and GFF.

    extra_opts are inserted into both invocations (protected: -t, --paf, --gff).
    """
    ensure_dir(out_dir)
    paf = out_dir / "miniprot.paf"
    gff = out_dir / "miniprot.gff"
    raw_paf = out_dir / "miniprot.raw.paf"
    raw_gff = out_dir / "miniprot.raw.gff"
    if not (raw_paf.exists() and raw_gff.exists()):
        for mode, raw in (("--paf", raw_paf), ("--gff", raw_gff)):
            cmd = [which("miniprot"), f"-t{threads}", mode,
                   *(extra_opts or []), str(genome), str(proteins)]
            with open(raw, "w") as fh:
                import subprocess
                log.info("RUN  miniprot (%s) %s", mode, Path(proteins).name)
                rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
            if rc != 0:
                raise ToolError(f"miniprot failed (rc={rc}): see {raw}")
    if not lazy([paf, gff], desc="miniprot mnp-filter"):
        counts = filter_miniprot_alignments(raw_gff, raw_paf, gff, paf,
                                            proteins, filt or {})
        drops = {k: v for k, v in counts.items() if k != "kept" and v}
        log.info("miniprot mnp-filter (%s): kept %d, dropped %s",
                 Path(proteins).name, counts["kept"], drops or 0)
        (out_dir / "mnp_filter.tsv").write_text(
            "reason\tcount\n" + "".join(f"{k}\t{v}\n"
                                          for k, v in sorted(counts.items())))
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


def _run_sets_miniprot(vi: ValidatedInputs, hdir: Path, threads: int,
                       filt: dict | None = None,
                       extra_opts: list[str] | None = None) -> tuple[Path, Path]:
    """Per-set miniprot runs, merged transparently into hdir/miniprot.{paf,gff}."""
    sroot = ensure_dir(hdir / "sets")
    n_sets = len(vi.homolog_sets)
    t_each = max(1, threads // n_sets)
    set_outputs = []
    for label, path in vi.homolog_sets.items():
        set_outputs.append(run_miniprot(vi.genome, path, ensure_dir(sroot / label),
                                        t_each, filt, extra_opts))
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
                       mnp_filter: str = "",
                       miniprot_opts: list[str] | None = None) -> dict:
    referee_min_wi = REFEREE_MIN_WI
    hdir = ensure_dir(out_dir / "03.homolog")
    log.info("module 2 [homolog]: %s (%d set(s))", vi.homologs, len(vi.homolog_sets))
    membership = _load_membership(vi)
    filt = parse_mnp_filter(mnp_filter)
    active = [k for k in ("UnStop", "UnShift", "minIden", "minPos", "minSize",
                          "maxRank", "R2Q", "Q2R", "topRank")
              if filt.get(k)] or filt.get("maxSize", 10**9) != 10**9
    if active:
        log.info("homolog mnp-filter active (%s): %s", mnp_filter, active)

    if len(vi.homolog_sets) > 1:
        paf, gff = _run_sets_miniprot(vi, hdir, threads, filt, miniprot_opts)
    else:
        paf, gff = run_miniprot(vi.genome, vi.homologs, hdir, threads, filt,
                                miniprot_opts)

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
