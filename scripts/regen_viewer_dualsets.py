#!/usr/bin/env python3
"""Regenerate the evidence browser for an existing 01-09 layout run
(fullrun_dualsets) with per-homolog-set protein tracks."""
import json
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, "/work/shidishen/AnnoMerge/annomerge/src")
from annomerge.modules.viewer import run_viewer

OUT = Path("/work/shidishen/AnnoMerge/annomerge/result/fullrun_dualsets")
TD = Path("/work/shidishen/AnnoMerge/test_data")
SC = TD / "Saccharomyces_cerevisiae"

vi = SimpleNamespace(genome=SC / "genome/GCF_000146045.2_R64_genomic.fna")

# homolog results: gff + per-set relatedness weights
sets_tsv = OUT / "03.homolog/sets.tsv"
set_weights = {}
for line in sets_tsv.read_text().splitlines()[1:]:
    c = line.split("\t")
    set_weights[c[0]] = float(c[6])
results = {
    "homolog": {"gff": OUT / "03.homolog/miniprot.gff",
                "set_weights": set_weights},
    "rna": {"stringtie_gtf": OUT / "02.rnaseq/stringtie/stringtie.gtf",
            "junctions_bed": OUT / "02.rnaseq/junctions.bed",
            "bams": [OUT / "02.rnaseq/SRR31963760.sorted.bam"]},
}
preds = {"helixer": str(SC / "denovo/helixer.gff"),
         "annevo": str(SC / "denovo/annevo.norm.gff3")}
details = {}
for name in preds:
    p = OUT / f"05.score/{name}/score_detail.json"
    if p.exists():
        details[name] = json.load(open(p))

final_gff3 = OUT / "07.augment/augment.sanity.gff3"
final_d = OUT / "07.augment"
integ_d = OUT / "06.merge"

# drop the stale single-file protein track data (replaced by per-set files)
stale = OUT / "09.result/viewer/data/protein.gff3"
if stale.exists():
    stale.unlink()

vw = run_viewer(OUT, vi, results, preds, details,
                final_gff3, final_d, integ_d, threads=4)
print("VIEWER_REGEN_OK", vw["html"])
