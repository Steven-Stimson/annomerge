"""OMArk orthology-based annotation quality assessment (BUSCO complement).

Runs OMArk against an OMAmer database on the longest-isoform proteome and
parses the ``<prefix>.sum`` summary for the report.  Everything is cached:
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
    """Parse OMArk's ``*.sum`` text summary (tolerant to format tweaks).

    OMArk 0.3.x ``<prefix>.sum`` layout::

        #On the whole proteome, there are 5425 proteins
        A:5362[P:40,F:3],I:4[P:4,F:0],C:0[P:0,F:0],U:59
        A:98.84%[P:0.74%,F:0.06%],I:0.07%[...],C:0.00%[...],U:1.09%
        S:3265,D:361[U:55,E:306],M:17
        S:89.62%,D:9.91%[U:1.51%,E:8.40%],M:0.47%
    """
    out: dict = {"ok": True, "db": Path(db).name,
                 "omark_version": version}
    total = None
    # --- OMArk 0.3.x explicit format -------------------------------------
    m = re.search(r"there are\s+([\d,]+)\s+proteins", txt, re.I)
    if m:
        total = _num(m.group(1))
    m = re.search(r"^A:([\d,]+)\[P:[\d,]+,F:[\d,]+\],I:([\d,]+)\["
                  r"P:[\d,]+,F:[\d,]+\],C:([\d,]+)\["
                  r"P:[\d,]+,F:[\d,]+\],U:([\d,]+)", txt, re.M)
    if m:
        out["consistent"], out["inconsistent"] = _num(m.group(1)), _num(m.group(2))
        out["contamination"], out["unknown"] = _num(m.group(3)), _num(m.group(4))
    m = re.search(r"^A:([\d.]+)%\[P:[\d.]+%,F:[\d.]+%\],I:[\d.]+%\["
                  r"P:[\d.]+%,F:[\d.]+%\],C:([\d.]+)%\["
                  r"P:[\d.]+%,F:[\d.]+%\],U:([\d.]+)%", txt, re.M)
    if m:
        out["consistent_pct"], out["contamination_pct"] = m.group(1), m.group(2)
        out["unknown_pct"] = m.group(3)
    m = re.search(r"^S:([\d,]+),D:([\d,]+)\[U:[\d,]+,E:[\d,]+\],M:([\d,]+)",
                  txt, re.M)
    if m:
        out["single"], out["duplicated"], out["missing"] = \
            _num(m.group(1)), _num(m.group(2)), _num(m.group(3))
    m = re.search(r"^S:[\d.]+%?,D:([\d.]+)%", txt, re.M)
    if m:
        out["duplicated_pct"] = m.group(1)
    # --- generic fallbacks (older/newer phrasings) ------------------------
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
        m = re.search(r"duplicated\D*([\d,]+)\s*\((\d+(?:\.\d+)?)\s*%\)",
                      line, re.I)
        if m and "duplicated_pct" not in out:
            out["duplicated"] = _num(m.group(1))
            out["duplicated_pct"] = m.group(2)
    if total is None:
        parts = [out.get(k) for k in ("consistent", "inconsistent", "unknown")]
        parts = [p for p in parts if p is not None]
        total = sum(parts) if parts else None
    if not out.get("consistent") and not out.get("inconsistent"):
        return {"ok": False, "error": "could not parse OMArk summary"}
    out["total"] = total
    if "consistent_pct" not in out:
        out["consistent_pct"] = _pct(out.get("consistent"), total)
    if "inconsistent_pct" not in out:
        out["inconsistent_pct"] = _pct(out.get("inconsistent"), total)
    if "unknown_pct" not in out:
        out["unknown_pct"] = _pct(out.get("unknown"), total)
    # contamination = OMArk "Likely Contamination" (C:), NOT duplicated HOGs
    if "contamination_pct" not in out:
        out["contamination_pct"] = _pct(out.get("contamination"), total) \
            if out.get("contamination") is not None else "0"
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
    # omark 0.3.x 的汇总文件是 <prefix>.sum（不是 *.omark）
    summaries = list(bdir.glob("*.sum"))
    if not summaries:
        tsv = bdir / "omamer_hits.tsv"
        if not tsv.exists():
            # omamer 2.x CLI: -o/--out 无长形式, 线程参数是 -t (0..48)
            run([omamer_exe, "search", "-d", str(db),
                 "-q", str(faa), "-o", str(tsv),
                 "-t", str(threads)],
                desc="omamer search (OMArk)")
        # omark 0.3.x CLI: -o 是输出目录(非文件前缀), 无线程参数(-t 是 TAXID)
        run([omark_exe, "-f", str(tsv), "-d", str(db),
             "-o", str(bdir)],
            desc="omark assessment")
        summaries = list(bdir.glob("*.sum"))
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
