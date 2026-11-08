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
