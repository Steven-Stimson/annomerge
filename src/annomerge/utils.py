"""Shared utilities: logging, subprocess wrapper, lazy (resumable) execution."""
from __future__ import annotations

import itertools
import logging
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import NamedTuple

log = logging.getLogger("annomerge")


def _fmt_dt(sec: float) -> str:
    """Format an elapsed time, e.g. 8s / 5m07s / 1h23m45s."""
    sec = int(round(sec))
    if sec < 60:
        return f"{sec}s"
    m, s = divmod(sec, 60)
    if m < 60:
        return f"{m}m{s:02d}s"
    h, m = divmod(m, 60)
    return f"{h}h{m:02d}m{s:02d}s"


class StepTracker:
    """English step-by-step progress logger.

    The pipeline is narrated as numbered big steps (``[STEP 02] rna - ...``)
    and, within each, compact one-line sub-step entries
    (``  [02|2/6] hisat2 DRR070508 (5m04s)``). Parallel fan-out is aggregated
    (``  [06|EVM] 3/7 done``) so worker output never interleaves into the main
    log. fail() pinpoints the failing step and its per-step log file.
    """

    def __init__(self) -> None:
        self.step_id = "00"
        self.step_name = "pipeline"
        self.sub_total: int | None = None
        self.sub_done = 0
        self._t0: float | None = None
        self._lock = threading.Lock()

    def big(self, sid: str, name: str, desc: str, total: int | None = None) -> None:
        """Start a big step (01-09); resets the sub-step counter."""
        with self._lock:
            self.step_id, self.step_name = sid, name
            self.sub_total, self.sub_done = total, 0
            self._t0 = time.monotonic()
        log.info("")
        log.info("[STEP %s] %s - %s", sid, name, desc)

    def sub(self, label: str, note: str = "") -> None:
        """One completed sub-step (main thread only, see run())."""
        with self._lock:
            self.sub_done += 1
            n, tot = self.sub_done, self.sub_total
        prog = f"{n}/{tot}" if tot else str(n)
        suffix = f" ({note})" if note else ""
        log.info("  [%s|%s] %s%s", self.step_id, prog, label, suffix)

    def skip(self, label: str) -> None:
        """A sub-step skipped because its outputs already exist (cached)."""
        with self._lock:
            self.sub_done += 1
            n, tot = self.sub_done, self.sub_total
        prog = f"{n}/{tot}" if tot else str(n)
        log.info("  [%s|%s] %s (cached)", self.step_id, prog, label)

    def parallel(self, label: str, done: int, total: int) -> None:
        """Aggregate progress for a parallel fan-out (e.g. EVM partitions)."""
        log.info("  [%s|%s] %d/%d done", self.step_id, label, done, total)

    def note(self, msg: str, *args) -> None:
        """A progress/info line inside the current big step (indented, no counter).
        Main thread only, so parallel workers do not interleave info lines."""
        if threading.current_thread() is threading.main_thread():
            log.info("  [%s] " + msg, self.step_id, *args)

    def big_done(self, note: str = "") -> None:
        dt = time.monotonic() - self._t0 if self._t0 else 0.0
        with self._lock:
            n, tot = self.sub_done, self.sub_total
        prog = f" - {n}/{tot} steps" if tot else ""
        extra = f" - {note}" if note else ""
        log.info("[STEP %s] %s done%s%s (%s)",
                 self.step_id, self.step_name, prog, extra, _fmt_dt(dt))

    def fail(self, label: str, detail: str = "") -> None:
        extra = f" - {detail}" if detail else ""
        log.error("[STEP %s] %s FAILED at: %s%s",
                  self.step_id, self.step_name, label, extra)


tracker = StepTracker()

# Per-step log directory: when set, run() redirects external-tool output here
# (one file per command) instead of letting it flood the main log.
_RUN_LOG_DIR: Path | None = None
_run_counter = itertools.count(1)


def set_run_log_dir(d: Path | str | None) -> None:
    """Set the directory that run() uses for auto per-command log files."""
    global _RUN_LOG_DIR
    _RUN_LOG_DIR = None if d is None else ensure_dir(Path(d))


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


def parse_extra_opts(spec: str | None, protected: tuple[str, ...], tool: str) -> list[str]:
    """shlex-split a user-supplied extra-options string and reject flags that
    the pipeline manages itself (e.g. 'hisat2 -x/-p/-S/-1/-2/-U/--dta').

    protected flags may not appear at all: not as a bare token, with '=',
    or glued to a value ('-p8' / '--threads=8' forms are rejected too).
    Returns [] for empty/None specs.
    """
    import shlex
    if not spec or not spec.strip():
        return []
    toks = shlex.split(spec)
    for tok in toks:
        for f in protected:
            base = tok.split("=", 1)[0]
            if base == f or (len(f) == 2 and tok.startswith(f) and len(tok) > 2):
                raise ToolError(
                    f"{tool}: '{tok}' overlaps a pipeline-managed option ({f}); "
                    f"remove it from the extra options. Protected: {', '.join(protected)}")
    return toks


class ToolOpts(NamedTuple):
    """Resolved extra options for one external tool: global flags plus
    per-target overrides (RNA sample name or homolog-set label)."""
    global_flags: list[str]
    per_target: dict[str, list[str]]

    def for_target(self, name: str | None) -> list[str]:
        """Global flags followed by the overrides for `name` (last wins)."""
        if name is None:
            return list(self.global_flags)
        return list(self.global_flags) + list(self.per_target.get(name, []))


def parse_targeted_opts(specs: list[str] | None, protected: tuple[str, ...],
                        tool: str) -> ToolOpts:
    """Parse repeatable --xxx-options values into a ToolOpts.

    Each value is either bare flags (global, e.g. "-k 3 --splice-flank=no")
    or a targeted override 'NAME=flags' where NAME is an RNA sample name or a
    homolog-set label (e.g. "isoseq_sub=-x splice:hq"). A value counts as
    targeted iff it contains '=' and the part before the first '=' is a bare
    name (no leading '-', no spaces), so flag forms like '--min-length=200'
    stay global. A bare '=' token in a global spec (e.g. the mistyped
    "iso = -x hq") is rejected up front instead of leaking into the tool
    command line. Protected-flag checks apply to every spec (global or
    targeted).
    """
    global_flags: list[str] = []
    per_target: dict[str, list[str]] = {}
    for spec in specs or []:
        spec = spec.strip()
        if not spec:
            continue
        head, sep, rest = spec.partition("=")
        if sep and head and not head.startswith("-") and " " not in head:
            name = head
            if not rest.strip():
                raise ToolError(f"{tool}: targeted option '{spec}' has no flags")
            if name in per_target:
                raise ToolError(f"{tool}: duplicate targeted options for '{name}'")
            per_target[name] = parse_extra_opts(rest, protected, tool)
        else:
            toks = parse_extra_opts(spec, protected, tool)
            if "=" in toks:
                raise ToolError(
                    f"{tool}: stray '=' in extra options '{spec}' — a targeted "
                    f"override must be written as NAME=flags with no space "
                    f"around '='")
            global_flags += toks
    return ToolOpts(global_flags, per_target)


def which(tool: str) -> str:
    """Resolve an external tool: dependency registry (env vars) first,
    PATH as fallback.  Import is local to avoid a module-level cycle."""
    from .deps import which as _deps_which
    return _deps_which(tool)


def run(cmd: list[str] | str, *, shell: bool = False, log_file: Path | None = None,
        desc: str = "", cwd: Path | str | None = None, track: bool = True) -> None:
    """Run a command, raising ToolError on non-zero exit.

    The external tool's stdout/stderr is redirected to a per-command log file
    and never interleaved into the main log: the explicit ``log_file`` if given,
    otherwise an auto-named file under the run-log dir (see set_run_log_dir),
    otherwise discarded. Main-thread calls log a one-line sub-step entry via the
    StepTracker; worker threads stay silent so parallel fan-out (EVM partitions,
    OMArk prefetch, multi-sample) does not interleave — the caller aggregates it.
    """
    shown = cmd if isinstance(cmd, str) else " ".join(map(str, cmd))
    log.debug("cmd: %s", shown)
    if log_file is None and _RUN_LOG_DIR is not None:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", desc or "cmd").strip("_")[:60] or "cmd"
        log_file = Path(_RUN_LOG_DIR) / f"{next(_run_counter):03d}_{safe}.log"
    t0 = time.monotonic()
    if log_file:
        with open(log_file, "ab") as fh:
            rc = subprocess.run(cmd, shell=shell, stdout=fh,
                                stderr=subprocess.STDOUT, cwd=cwd).returncode
    else:
        rc = subprocess.run(cmd, shell=shell, stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT, cwd=cwd).returncode
    dt = time.monotonic() - t0
    if rc != 0:
        tracker.fail(desc or shown,
                     f"rc={rc}" + (f", log: {log_file}" if log_file else ""))
        raise ToolError(f"command failed (rc={rc}): {shown}" +
                        (f" (log: {log_file})" if log_file else ""))
    if track and desc and threading.current_thread() is threading.main_thread():
        tracker.sub(desc, _fmt_dt(dt))


def lazy(outputs: list[Path] | Path, *, desc: str = "") -> bool:
    """Return True if all outputs already exist (step can be skipped)."""
    if not isinstance(outputs, list):
        outputs = [outputs]
    ok = all(Path(o).exists() and Path(o).stat().st_size > 0 for o in outputs)
    if ok and desc and threading.current_thread() is threading.main_thread():
        tracker.skip(desc)
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
