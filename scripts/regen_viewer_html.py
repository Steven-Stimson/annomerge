#!/usr/bin/env python3
"""Regenerate the offline viewer HTML from an existing viewer data/ dir
(template-only change: no pipeline stage re-run).

Usage:
  python3 regen_viewer_html.py                # classic -> AnnoMerge_evidence_browser.html
  python3 regen_viewer_html.py npg            # theme -> AnnoMerge_evidence_browser_theme_<theme>.html
  python3 regen_viewer_html.py npg out.html   # theme -> explicit output name
"""
import sys
from pathlib import Path
sys.path.insert(0, "/work/shidishen/AnnoMerge/annomerge/src")
from annomerge.modules.viewer import render_offline_html

vdir = Path("/work/shidishen/AnnoMerge/annomerge/result/fullrun_yeast_v3/09.result/viewer")
data = vdir / "data"
theme = sys.argv[1] if len(sys.argv) > 1 else "npg"
if len(sys.argv) > 2:
    out = vdir / sys.argv[2]
elif theme == "npg":
    out = vdir / "AnnoMerge_evidence_browser.html"
else:
    out = vdir / f"AnnoMerge_evidence_browser_theme_{theme}.html"

fai = (data / "genome.fa.fai").read_text().split("\t")
fallback = f"{fai[0]}:1-{min(int(fai[1]), 200000)}"
html = render_offline_html(
    data, out,
    ref_name="GCF_000146045.2_R64_genomic.fna",
    fallback_locus=fallback, coverage=None, full=False, theme=theme)
print("REGEN_OK", theme, html, html.stat().st_size)
