"""OMArk orthology-based annotation quality assessment (BUSCO complement).

Runs OMArk against an OMAmer database on the longest-isoform proteome and
parses the ``<prefix>.omark`` summary for the report.  Everything is cached:
if the summary already exists it is parsed again instead of re-running.

The OMAmer database (.h5) is NOT shipped with the software -- download from
https://omabrowser.org/oma/archives/ (e.g. LUCA.h5 or a kingdom database)
and pass it via ``--omark-db``.  When absent the report shows a skip note
instead of failing the run.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

from ..utils import ensure_dir, run, which

log = logging.getLogger("annomerge")


def _num(s: str) -> int | None:
    m = re.search(r"[\d,]+", s)
    return int(m.group(0).replace(",", "")) if m else None


def _pct(n: int | None, total: int | None) -> str:
    if not n or not total:
        return "0"
    return f"{100.0 * n / total:.1f}"


def parse_omark_summary(txt: str, db: Path, version: str) -> dict:
    """Parse OMArk's ``*.omark`` text summary (tolerant to format tweaks)."""
    out: dict = {"ok": True, "db": Path(db).name,
                 "omark_version": version}
    total = None
    for line in txt.splitlines():
        m = re.search(r"(?:total number of proteins|number of proteins assessed"
                      r"|total proteins)\D*([\d,]+)", line, re.I)
        if m and total is None:
            total = _num(line)
        m = re.search(r"Consistent(?: proteins)?\D*([\d,]+)", line, re.I)
        if m and "consistent" not in out:
            out["consistent"] = _num(line)
        m = re.search(r"Inconsistent(?: proteins)?\D*([\d,]+)", line, re.I)
        if m and "inconsistent" not in out:
            out["inconsistent"] = _num(line)
        m = re.search(r"Unknown(?: proteins)?\D*([\d,]+)", line, re.I)
        if m and "unknown" not in out:
            out["unknown"] = _num(line)
        m = re.search(r"duplicated.*?(?:\D*?)(\d+(?:\.\d+)?)\s*%", line, re.I)
        if m and "duplicated_pct" not in out:
            # OMArk: "Consistent: N (X%) [ ... duplicated: M (Y%) ... ]"
            dm = re.search(r"duplicated\D*([\d,]+)\s*\((\d+(?:\.\d+)?)\s*%\)",
                           line, re.I)
            if dm:
                out["duplicated"] = _num(dm.group(1))
                out["duplicated_pct"] = dm.group(2)
            else:
                out["duplicated_pct"] = m.group(1)
    if total is None:
        parts = [out.get(k) for k in ("consistent", "inconsistent", "unknown")]
        parts = [p for p in parts if p is not None]
        total = sum(parts) if parts else None
    if not out.get("consistent") and not out.get("inconsistent"):
        return {"ok": False, "error": "could not parse OMArk summary"}
    out["total"] = total
    out["consistent_pct"] = _pct(out.get("consistent"), total)
    out["inconsistent_pct"] = _pct(out.get("inconsistent"), total)
    out["contamination_pct"] = out.get("duplicated_pct") or \
        _pct(out.get("duplicated"), total)
    return out


def _resolve_exe(name: str) -> str | None:
    """Find an omark/omamer executable: PATH first, then the bin/ directory
    sitting next to an importable omamer package (pip --target style installs
    put the console scripts there without touching PATH)."""
    import shutil
    exe = shutil.which(name)
    if exe:
        return exe
    try:
        import omamer  # noqa: F401
        cand = Path(omamer.__file__).resolve().parent.parent / "bin" / name
        if cand.exists():
            return str(cand)
    except Exception:
        pass
    return None


def run_omark(faa: Path, db: Path, out_dir: Path, threads: int) -> dict:
    """Run (or reuse cached) OMArk; returns report-ready context dict."""
    db, out_dir = Path(db), Path(out_dir)
    if not db.exists():
        return {"ok": False,
                "error": f"OMAmer database not found: {db} "
                         "(download from https://omabrowser.org/oma/archives/)"}
    omark_exe, omamer_exe = _resolve_exe("omark"), _resolve_exe("omamer")
    if not omark_exe or not omamer_exe:
        return {"ok": False,
                "error": "omark/omamer not found (install: conda install -c "
                         "bioconda omark, or pip install omark omamer)"}
    bdir = ensure_dir(out_dir)
    summaries = list(bdir.glob("*.omark"))
    if not summaries:
        tsv = bdir / "omamer_hits.tsv"
        if not tsv.exists():
            run([omamer_exe, "search", "--db", str(db),
                 "--query", str(faa), "--out", str(tsv),
                 "--nthreads", str(threads)],
                desc="omamer search (OMArk)")
        run([omark_exe, "-f", str(tsv), "-d", str(db),
             "-o", str(bdir / "omark"), "--threads", str(threads)],
            desc="omark assessment")
        summaries = list(bdir.glob("*.omark"))
    if not summaries:
        return {"ok": False, "error": "no OMArk summary produced"}
    version = ""
    try:
        from importlib.metadata import version as _v
        version = _v("omark")
    except Exception:
        pass
    ctx = parse_omark_summary(summaries[0].read_text(), db, version)
    ctx["summary_file"] = str(summaries[0])
    return ctx
