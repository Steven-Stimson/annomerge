"""Parameter resolution: CLI > --config YAML > defaults; dump resolved_params.yaml."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import yaml

log = logging.getLogger("annomerge")


def resolve_params(cli_params: dict[str, Any], config_file: str | None) -> dict[str, Any]:
    """Merge a previous resolved_params.yaml (--config) under CLI params.

    Convention: cli_params carries None for options the user did not pass,
    so config values only fill the gaps; CLI always wins.
    """
    base: dict[str, Any] = {}
    if config_file:
        with open(config_file) as fh:
            base = yaml.safe_load(fh) or {}
        log.info("loaded config snapshot from %s", config_file)
    merged = {**base, **{k: v for k, v in cli_params.items() if v is not None}}
    return merged


def dump_params(params: dict[str, Any], out_dir: Path) -> Path:
    path = out_dir / "resolved_params.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(params, fh, allow_unicode=True, sort_keys=True)
    log.info("resolved params -> %s", path)
    return path
