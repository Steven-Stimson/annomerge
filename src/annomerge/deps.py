"""External dependency registry — short-name environment variable resolution.

Every external tool is resolved in this order:
  1. the dependency's short env vars (directory or executable file both OK;
     for a directory, <dir>/<exe> and <dir>/bin/<exe> are probed)
  2. PATH lookup (fallback, keeps env.sh working as before)

Example (user side):
    export EVM=/opt/EVidenceModeler
    export MINIPROT=/opt/miniprot/miniprot      # direct executable is fine too
    export TRANSDECODER=/opt/TransDecoder
    export HISAT2=/opt/hisat2-2.2.1
    annomerge merge -c          # verify all dependencies, then exit

Dir-type dependencies (EVM, TRANSDECODER) have no PATH-resolvable binary of
their own, so their env var must point at the installation directory;
TRANSDECODER falls back to the parent directory of TransDecoder.LongOrfs in
PATH if it is in PATH at all.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .utils import ToolError

log = __import__("logging").getLogger("annomerge")


# ---------------------------------------------------------------- registry
@dataclass(frozen=True)
class ExeDep:
    """A dependency reached through executable(s)."""
    key: str                          # canonical id, e.g. "hisat2"
    envs: tuple[str, ...]             # env vars, in priority order
    exes: tuple[str, ...]             # executables belonging to this dep
    hint: str = ""                    # install hint shown on failure
    alternatives: tuple[str, ...] = ()  # any one of these also satisfies the dep
                                         # (e.g. the combined `TransDecoder` wrapper
                                         #  when the split .LongOrfs/.Predict are absent)
    optional: bool = False              # missing -> WARN (yellow), never counted
                                         # as a hard failure (busco/omark add-ons)


@dataclass(frozen=True)
class HomeDep:
    """A directory-based dependency (perl-script suites without a binary)."""
    key: str
    envs: tuple[str, ...]
    markers: tuple[str, ...]          # relative paths that must exist inside
    path_fallback: tuple[str, ...] = ()  # executables in PATH to derive home from
    hint: str = ""


EXE_DEPS: tuple[ExeDep, ...] = (
    ExeDep("hisat2", ("HISAT2", "HISAT2_HOME"), ("hisat2", "hisat2-build"),
           "hisat2: conda install -c bioconda hisat2"),
    ExeDep("samtools", ("SAMTOOLS", "SAMTOOLS_HOME"), ("samtools",),
           "samtools: conda install -c bioconda samtools"),
    ExeDep("stringtie", ("STRINGTIE", "STRINGTIE_HOME"), ("stringtie",),
           "stringtie: conda install -c bioconda stringtie"),
    ExeDep("gffread", ("GFFREAD", "GFFREAD_HOME"), ("gffread",),
           "gffread/gffcompare (gffclib): conda install -c bioconda gffclib"),
    ExeDep("gffcompare", ("GFFCOMPARE", "GFFCOMPARE_HOME"), ("gffcompare",),
           "gffread/gffcompare (gffclib): conda install -c bioconda gffclib"),
    ExeDep("diamond", ("DIAMOND", "DIAMOND_HOME"), ("diamond",),
           "diamond: conda install -c bioconda diamond"),
    ExeDep("busco", ("BUSCO", "BUSCO_HOME"), ("busco",),
           "busco: conda install -c bioconda busco (5.7.x for odb10 lineage) "
           "(only needed with --busco-lineage)", optional=True),
    ExeDep("miniprot", ("MINIPROT", "MINIPROT_HOME"), ("miniprot",),
           "miniprot: conda install -c bioconda miniprot"),
    ExeDep("transdecoder", ("TRANSDECODER", "TRANSDECODER_HOME"),
           ("TransDecoder.LongOrfs", "TransDecoder.Predict"),
           "TransDecoder: conda install -c bioconda transdecoder",
           alternatives=("TransDecoder",)),
    ExeDep("perl", ("PERL",), ("perl",), "perl 5: system package"),
    ExeDep("omark", ("OMARK", "OMARK_HOME"), ("omark",),
           "omark: conda install -c bioconda omark (only needed with --omark-db)",
           optional=True),
    ExeDep("omamer", ("OMAMER", "OMAMER_HOME"), ("omamer",),
           "omamer: conda install -c bioconda omamer (only needed with --omark-db)",
           optional=True),
)

HOME_DEPS: tuple[HomeDep, ...] = (
    HomeDep("evm", ("EVM", "EVM_HOME"),
            ("EvmUtils/partition_EVM_inputs.pl",
             "EvmUtils/convert_EVM_outputs_to_GFF3.pl"),
            path_fallback=None,
            hint="EVidenceModeler: https://github.com/EVidenceModeler/EVidenceModeler "
                 "(export EVM=/path/to/EVidenceModeler)"),
    HomeDep("transdecoder", ("TRANSDECODER", "TRANSDECODER_HOME"),
            ("util/gtf_to_alignment_gff3.pl",),
            path_fallback=("TransDecoder.LongOrfs", "TransDecoder"),
            hint="TransDecoder: conda install -c bioconda transdecoder "
                 "(export TRANSDECODER=/path/to/TransDecoder)"),
)

_EXE_LOOKUP: dict[str, ExeDep] = {exe: dep for dep in EXE_DEPS for exe in dep.exes}
_HOME_LOOKUP: dict[str, HomeDep] = {dep.key: dep for dep in HOME_DEPS}

# ---------------------------------------------------------------- resolution
def _probe_dir(dirp: Path, exe: str) -> str | None:
    for cand in (dirp / exe, dirp / "bin" / exe, dirp / "util" / exe):
        if cand.is_file() and os.access(cand, os.X_OK):
            return str(cand)
    return None


def which(exe: str) -> str:
    """Resolve an external executable: env vars first, then PATH."""
    dep = _EXE_LOOKUP.get(exe)
    if dep:
        for var in dep.envs:
            val = os.environ.get(var)
            if not val:
                continue
            p = Path(val)
            if p.is_dir():
                hit = _probe_dir(p, exe)
                if hit:
                    return hit
                log.warning("%s=%s is a directory but %s not found inside "
                            "(probed %s, bin/ and util/)", var, val, exe, p)
            elif p.is_file() and os.access(p, os.X_OK):
                if p.name == exe or len(dep.exes) == 1:
                    return str(p)
                log.warning("%s=%s is a file but does not match executable %r",
                            var, val, exe)
            else:
                log.warning("%s=%s does not exist or is not executable", var, val)
    hit = shutil.which(exe)
    if hit:
        return hit
    # conda layout fallback: sometimes only a wrapper is on PATH (e.g.
    # bin/TransDecoder -> ../opt/transdecoder/TransDecoder) while the real
    # scripts (TransDecoder.LongOrfs/.Predict, util/) live next to the
    # wrapper target.  Follow a resolvable sibling/alternative of the same
    # dependency and probe its real directory.
    for anchor in (exe, *(dep.alternatives if dep else ())):
        wrapper = shutil.which(anchor)
        if not wrapper:
            continue
        real_dir = Path(wrapper).resolve().parent
        for sub in (".", "util", "bin"):
            cand = real_dir / sub / exe
            if cand.is_file() and os.access(cand, os.X_OK):
                return str(cand)
    tried = ", ".join(dep.envs) if dep else ""
    raise ToolError(
        f"external tool not found: {exe!r}. "
        + (f"Tried env vars: {tried} (unset); not in PATH either. {dep.hint}"
           if dep else f"Not in PATH."))


def home(key: str) -> Path:
    """Resolve a directory-based dependency (e.g. evm, transdecoder)."""
    dep = _HOME_LOOKUP[key]
    for var in dep.envs:
        val = os.environ.get(var)
        if not val:
            continue
        p = Path(val)
        if p.is_dir() and all((p / m).exists() for m in dep.markers):
            return p
        raise ToolError(
            f"{var}={val} is set but does not look like a valid {dep.key} "
            f"installation (missing: {', '.join(dep.markers)}). {dep.hint}")
    if dep.path_fallback:
        hits = [shutil.which(fb) for fb in dep.path_fallback]
        parents: list[Path] = []
        for hit in hits:
            if hit:
                parents.append(Path(hit).resolve().parent)          # exe dir
                parents.append(Path(hit).resolve().parent.parent)   # e.g. conda prefix root
        # conda layout: binaries in <env>/bin/, suite under <env>/opt/<name>-<ver>/
        for hit in hits:
            if hit:
                env_root = Path(hit).resolve().parent.parent
                parents.extend(sorted(env_root.glob("opt/*")))
        for p in parents:
            if all((p / m).exists() for m in dep.markers):
                return p
    raise ToolError(
        f"directory dependency {dep.key!r} not configured: set one of "
        f"{', '.join(dep.envs)} to its install directory. {dep.hint}")


def utility(key: str, relpath: str) -> Path:
    """Path to a script inside a home dependency (e.g. EvmUtils script)."""
    p = home(key) / relpath
    if not p.exists():
        raise ToolError(f"{key} utility missing: {p}")
    return p


# ---------------------------------------------------------------- checking
def _version_of(exe_path: str) -> str:
    for args in (["--version"], ["-version"], ["-v"]):
        try:
            env = dict(__import__("os").environ)
            if "busco" in exe_path:
                env.pop("PYTHONPATH", None)
            out = subprocess.run([exe_path] + args, capture_output=True,
                                 text=True, timeout=15, env=env)
            text = (out.stdout + out.stderr).strip()
            if text:
                lines = [l.strip() for l in text.splitlines() if l.strip()][:2]
                for line in lines:                       # prefer a dotted version token
                    m = re.search(r"(\d+\.\d+[\w.\-]*)", line)
                    if m:
                        return m.group(1)
                m = re.search(r"(\d+[\w.\-]*)", lines[0])
                return m.group(1) if m else lines[0][:40]
        except (OSError, subprocess.TimeoutExpired):
            continue
    return "?"


def _home_version(key: str, p: Path) -> str:
    """Best-effort version for a directory dependency (no --version binary)."""
    text = str(p)
    m = re.search(r"evidencemodeler[-_](\d+(?:\.\d+)+)", text, re.I)
    if m:
        return m.group(1)
    m = re.search(r"rans\w*ecoder[-_](\d+(?:\.\d+)+)", text, re.I)
    if m:
        return m.group(1)
    if key == "transdecoder":
        # version embedded in the scripts, e.g.  my $VERSION = "5.7.1";
        for script in ("TransDecoder.LongOrfs", "TransDecoder.Predict", "TransDecoder"):
            f = p / script
            if not f.is_file():
                continue
            try:
                head = "\n".join(f.read_text(errors="replace").splitlines()[:80])
            except OSError:
                continue
            m = re.search(r"\$VERSION\s*=\s*[\"']([\d.]+)", head)
            if m:
                return m.group(1)
    m = re.search(r"[-_](\d+\.\d+(?:\.\d+)?)\s*$", p.name)
    return m.group(1) if m else "?"


def _dist_version(key: str, exe_path: str | None) -> str:
    """Python-distribution tools (omark/omamer) report via importlib; CLI tools
    (and wrappers whose --version argparse would only print usage) via --version."""
    try:
        from importlib.metadata import version as _v
        return _v(key)
    except Exception:
        pass
    if exe_path and "TransDecoder" in Path(exe_path).name:
        # perl script with the version embedded: my $VERSION = "5.7.1";
        # (not every script carries it — the wrapper does, so probe siblings)
        cands = [Path(exe_path),
                 Path(exe_path).parent / "TransDecoder.Predict",
                 Path(exe_path).parent.parent / "TransDecoder"]
        for f in cands:
            try:
                head = f.read_text(errors="replace")[:4000]
            except OSError:
                continue
            m = re.search(r"\$VERSION\s*=\s*[\"']([\d.]+)", head)
            if m:
                return m.group(1)
    return _version_of(exe_path) if exe_path else "-"


def check_all() -> list[dict]:
    """Resolve every dependency; return rows for display (never raises)."""
    rows: list[dict] = []
    for dep in EXE_DEPS:
        resolved: dict[str, tuple[str, str]] = {}
        missing: list[str] = []
        for exe in dep.exes:
            try:
                path = which(exe)
                src = "PATH"
                for var in dep.envs:
                    v = os.environ.get(var)
                    if v and path.startswith(os.path.abspath(v)):
                        src = f"env:{var}"
                resolved[exe] = (path, src)
            except ToolError:
                missing.append(exe)
        if missing and dep.alternatives:
            # e.g. only the combined `TransDecoder` wrapper is installed
            for alt in dep.alternatives:
                try:
                    resolved[alt] = (which(alt), "PATH")
                    missing = []
                    break
                except ToolError:
                    continue
        first = next(iter(resolved.values()), (None, None))[0]
        rows.append({
            "dep": dep.key, "env": ",".join(dep.envs), "exes": dep.exes,
            "ok": not missing, "missing": missing, "optional": dep.optional,
            "path": first or "-", "src": next(iter(resolved.values()), (None, "-"))[1],
            "version": _dist_version(dep.key, first) if first else "-",
            "hint": dep.hint,
        })
    for dep in HOME_DEPS:
        try:
            p = home(dep.key)
            rows.append({"dep": dep.key, "env": ",".join(dep.envs), "exes": (),
                         "ok": True, "missing": [], "path": str(p),
                         "src": "env:" + next((v for v in dep.envs
                                               if os.environ.get(v)), "derived/PATH"),
                         "version": _home_version(dep.key, p), "hint": dep.hint})
        except ToolError as e:
            rows.append({"dep": dep.key, "env": ",".join(dep.envs), "exes": (),
                         "ok": False, "missing": [dep.key], "path": "-",
                         "src": "-", "version": "-", "hint": dep.hint})
    # optional OMArk support: ete3's NCBI taxonomy db.  omark needs it on the
    # first run, where it would otherwise try to download taxdump.tar.gz from
    # the NCBI FTP site — that fails behind proxies, so surface it here.
    taxa_db = Path.home() / ".etetoolkit" / "taxa.sqlite"
    if taxa_db.is_file() and taxa_db.stat().st_size > 0:
        rows.append({"dep": "ete3-taxonomy", "env": "HOME", "exes": (),
                     "ok": True, "missing": [], "optional": True,
                     "path": str(taxa_db), "src": "home",
                     "version": f"{taxa_db.stat().st_size / 1e6:.0f} MB", "hint": ""})
    else:
        rows.append({"dep": "ete3-taxonomy", "env": "HOME", "exes": (),
                     "ok": False, "missing": ["taxa.sqlite"], "optional": True,
                     "path": str(taxa_db), "src": "-", "version": "-",
                     "hint": "ete3 NCBI taxonomy db not found; on the first omark "
                             "run ete3 would download taxdump.tar.gz from the NCBI "
                             "FTP site (usually blocked by proxies). Build it "
                             "manually: download taxdump.tar.gz from "
                             "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ then run "
                             "python -c \"from ete3 import NCBITaxa; "
                             "NCBITaxa().update_taxonomy_database('taxdump.tar.gz')\""})
    # python package sanity — the actual runtime deps from pyproject.toml
    py = []
    for mod, dist in (("pysam", "pysam"), ("typer", "typer"),
                      ("jinja2", "jinja2"), ("yaml", "pyyaml")):
        try:
            m = __import__(mod)
            ver = getattr(m, "__version__", None)
            if not ver:
                try:
                    from importlib.metadata import version as _dist_ver
                    ver = _dist_ver(dist)
                except Exception:
                    ver = "?"
            py.append(f"{mod}={ver}")
        except ImportError:
            py.append(f"{mod}=MISSING")
    rows.append({"dep": "python-packages", "env": "python", "exes": (),
                 "ok": not any("MISSING" in s for s in py), "missing": [],
                 "path": "-", "src": "-", "version": "; ".join(py),
                 "hint": "pip install typer pyyaml pysam jinja2"})
    # deps registered both as ExeDep and HomeDep (transdecoder: entry points +
    # util/ script suite) display as a single merged row
    merged: dict[str, dict] = {}
    order: list[str] = []
    for r in rows:
        if r["dep"] in merged:
            m = merged[r["dep"]]
            m["ok"] = m["ok"] and r["ok"]
            m["missing"] = m["missing"] + r["missing"]
            if m["version"] in ("-", "?") and r["version"] not in ("-", "?"):
                m["version"] = r["version"]
            if m["path"] == "-" and r["path"] != "-":
                m["path"], m["src"] = r["path"], r["src"]
        else:
            merged[r["dep"]] = dict(r)
            order.append(r["dep"])
    return [merged[k] for k in order]


_GREEN = "\033[32m"
_RED = "\033[31m"
_YELLOW = "\033[33m"
_RESET = "\033[0m"


def format_table(rows: list[dict]) -> str:
    out = ["dependency       env var(s)              status  version   resolved from"]
    out.append("-" * 78)
    for r in rows:
        if r["ok"]:
            status = f"{_GREEN}OK{_RESET}"
        elif r.get("optional"):
            miss = f" ({','.join(r['missing'])})" if r["missing"] else ""
            status = f"{_YELLOW}WARN{miss}{_RESET}"
        else:
            miss = f" ({','.join(r['missing'])})" if r["missing"] else ""
            status = f"{_RED}MISSING{miss}{_RESET}"
        out.append(f"{r['dep']:<16} {r['env']:<23} {status:<28} "
                   f"{r['version'][:42]:<42} {r['src']}: {r['path']}")
        if not r["ok"]:
            out.append(f"    {_YELLOW}->{_RESET} {r['hint']}")
    n_ok = sum(1 for r in rows if r["ok"])
    n_bad = sum(1 for r in rows if not r["ok"] and not r.get("optional"))
    n_warn = sum(1 for r in rows if not r["ok"] and r.get("optional"))
    out.append("-" * 78)
    if n_bad:
        out.append(f"{n_ok}/{len(rows)} dependencies OK; "
                   f"{_RED}{n_bad} MISSING{_RESET}")
    else:
        out.append(f"{_GREEN}{n_ok}/{len(rows)} dependencies OK"
                   f" — ready to run.{_RESET}")
    if n_warn:
        out.append(f"{_YELLOW}{n_warn} optional dependency warning(s) — "
                   f"busco/omark features unavailable until installed.{_RESET}")
    return "\n".join(out)
