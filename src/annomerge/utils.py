"""Shared utilities: logging, subprocess wrapper, lazy (resumable) execution."""
from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from pathlib import Path

log = logging.getLogger("annomerge")


def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stderr,
        force=True,
    )


class ToolError(RuntimeError):
    """Raised when an external tool fails or is missing."""


def which(tool: str) -> str:
    """Resolve an external tool: dependency registry (env vars) first,
    PATH as fallback.  Import is local to avoid a module-level cycle."""
    from .deps import which as _deps_which
    return _deps_which(tool)


def run(cmd: list[str] | str, *, shell: bool = False, log_file: Path | None = None,
        desc: str = "", cwd: Path | str | None = None) -> None:
    """Run a command, raising ToolError on non-zero exit."""
    if desc:
        log.info("RUN  %s", desc)
    shown = cmd if isinstance(cmd, str) else " ".join(map(str, cmd))
    log.debug("cmd: %s", shown)
    stdout = subprocess.DEVNULL if log_file else None
    stderr = subprocess.STDOUT if log_file else None
    if log_file:
        with open(log_file, "ab") as fh:
            rc = subprocess.run(cmd, shell=shell, stdout=fh, stderr=subprocess.STDOUT,
                                cwd=cwd).returncode
    else:
        rc = subprocess.run(cmd, shell=shell, stdout=stdout, stderr=stderr,
                            cwd=cwd).returncode
    if rc != 0:
        raise ToolError(f"command failed (rc={rc}): {shown}" + (f" (log: {log_file})" if log_file else ""))


def lazy(outputs: list[Path] | Path, *, desc: str = "") -> bool:
    """Return True if all outputs already exist (step can be skipped)."""
    if not isinstance(outputs, list):
        outputs = [outputs]
    ok = all(Path(o).exists() and Path(o).stat().st_size > 0 for o in outputs)
    if ok and desc:
        log.info("SKIP %s (outputs exist)", desc)
    return ok


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def inject_utr_features(gff3_in: Path, gff3_out: Path) -> tuple[int, int]:
    """Add five_prime_UTR/three_prime_UTR rows for every mRNA that has both
    exon and CDS children (UTR = exon minus CDS coverage, strand-aware).

    Already-annotated files (any UTR rows present) are copied unchanged, so
    predictors like helixer keep their native UTR models. Extra UTR rows are
    appended (order-independent; consumers filter by feature type) and are
    inert for parsers that only look at gene/mRNA/CDS/exon.
    Returns (n_five_prime, n_three_prime) added.
    """
    lines = open(gff3_in).read().splitlines()
    if any(("\tfive_prime_UTR\t" in l) or ("\tthree_prime_UTR\t" in l)
           for l in lines):
        with open(gff3_out, "w") as out:
            out.write("\n".join(lines) + "\n")
        return (0, 0)
    mrna_meta: dict = {}   # tid -> (seqid, source, strand)
    exons: dict = {}       # tid -> [(s, e)]
    cds: dict = {}         # tid -> [(s, e)]
    for ln in lines:
        if ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 9 or f[0] == "":
            continue
        attrs = {}
        for kv in f[8].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                attrs[k] = v
        if f[2] == "mRNA":
            tid = attrs.get("ID")
            if tid:
                mrna_meta[tid] = (f[0], f[1], f[6])
        elif f[2] == "exon":
            for parent in (attrs.get("Parent") or "").split(","):
                if parent:
                    exons.setdefault(parent, []).append((int(f[3]), int(f[4])))
        elif f[2] == "CDS":
            for parent in (attrs.get("Parent") or "").split(","):
                if parent:
                    cds.setdefault(parent, []).append((int(f[3]), int(f[4])))
    added5 = added3 = 0
    utr_rows = []
    for tid, (seqid, source, strand) in mrna_meta.items():
        ex, cd = exons.get(tid), cds.get(tid)
        if not ex or not cd:
            continue
        cd.sort()
        cov = []
        for s_, e_ in cd:
            if cov and s_ <= cov[-1][1] + 1:
                cov[-1][1] = max(cov[-1][1], e_)
            else:
                cov.append([s_, e_])
        cds_min, cds_max = cov[0][0], cov[-1][1]
        for (s_, e_) in sorted(ex):
            segs = [(s_, e_)]
            for (cs, ce) in cov:
                nxt = []
                for (a, b) in segs:
                    if ce < a or cs > b:
                        nxt.append((a, b))
                    else:
                        if a < cs:
                            nxt.append((a, cs - 1))
                        if b > ce:
                            nxt.append((ce + 1, b))
                segs = nxt
            for (a, b) in segs:
                if strand == "-":
                    utr = ("five_prime_UTR" if a > cds_max
                           else "three_prime_UTR" if b < cds_min else None)
                else:
                    utr = ("five_prime_UTR" if b < cds_min
                           else "three_prime_UTR" if a > cds_min and a > cds_max
                           else None)
                if utr is None:
                    continue
                utr_rows.append((seqid, source, utr, a, b, ".", strand, ".",
                                 f"Parent={tid}"))
                if utr == "five_prime_UTR":
                    added5 += 1
                else:
                    added3 += 1
    with open(gff3_out, "w") as out:
        for ln in lines:
            out.write(ln + "\n")
        for row in utr_rows:
            out.write("\t".join(map(str, row)) + "\n")
    return added5, added3
