#!/usr/bin/env python3
"""Regenerate the new-style report.html + evidence browser on the existing
fullrun_evi3 output (BUSCO/accuracy cached where possible)."""
import datetime
import os
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, "/work/shidishen/AnnoMerge/annomerge/src")

from annomerge.modules.report import run_report
from annomerge.modules.viewer import run_viewer

OUT = Path("/work/shidishen/AnnoMerge/annomerge/result/fullrun_evi3")
TD = Path("/work/shidishen/AnnoMerge/test_data")
SC = TD / "Saccharomyces_cerevisiae"
genome = SC / "genome/GCF_000146045.2_R64_genomic.fna"
final_d = OUT / "augment"
final_gff3 = final_d / "augment.sanity.gff3"
integ_d = OUT / "evm"
assert final_gff3.exists(), final_gff3

preds = {"helixer": str(SC / "denovo/helixer.gff"),
         "annevo": str(SC / "denovo/annevo.norm.gff3")}
details = {}
wtab = OUT / "scoring/evm_weights.tsv"
for line in wtab.read_text().splitlines()[1:]:
    c = line.split("\t")
    try:
        details[c[0]] = {"W_raw": float(c[1])}
    except Exception:
        details[c[0]] = {}
results = {
    "rna": {"bams": [str(OUT / "rna/bam/SRR31963760.sorted.bam")],
            "junctions_bed": str(OUT / "rna/junctions.bed"),
            "stringtie_gtf": str(OUT / "rna/stringtie/stringtie.gtf")},
    "homolog": {"gff": str(OUT / "homolog/miniprot.gff")},
}
vi = SimpleNamespace(genome=str(genome))

rep = run_report(
    OUT, final_gff3, genome, final_d,
    busco_lineage=TD / "db/busco/saccharomycetes_odb10",
    threads=4,
    omark_db=TD / "db/omark/Fungi.h5",  # not present yet -> informative skip note
    run_meta={"mode": "isoform",
              "command": (
                  "python3 annomerge merge "
                  f"-g {SC / 'genome/GCF_000146045.2_R64_genomic.fna'} "
                  f"-P helixer={SC / 'denovo/helixer.gff'},"
                  f"annevo={SC / 'denovo/annevo.norm.gff3'} "
                  f"-p {SC / 'homologs/homologs.faa'} "
                  f"-r {SC / 'rnaseq_deep/rnaseq_deep.list'} "
                  f"--busco {TD / 'db/busco/saccharomycetes_odb10'} "
                  f"-t 4 -o {OUT}"
              ),
              "threads": 4,
              "evm_weights": str(wtab),
              "homologs": results["homolog"]["gff"],
              "rna_bam": ",".join(results["rna"]["bams"]),
              "rna_stringtie": results["rna"]["stringtie_gtf"],
              "rna_junctions": results["rna"]["junctions_bed"],
              "busco_lineage": str(TD / "db/busco/saccharomycetes_odb10"),
              "omark_db": str(TD / "db/omark/Fungi.h5") + " (not present; skipped)",
              "predictors": ",".join(preds),
              "genome": genome.name, "outdir": str(OUT),
              "timestamp": datetime.datetime.now().isoformat(timespec="seconds")})
print("REPORT_HTML:", rep["html"])

vw = run_viewer(OUT, vi, results, preds, details,
                final_gff3, final_d, integ_d, threads=4)
print("VIEWER_HTML:", vw["html"])

# refresh the published 07.report copies
pub = OUT / "07.report"
pub.mkdir(exist_ok=True)
for src, dst in [(rep["html"], pub / "report.html"),
                 (OUT / "report/stats.json", pub / "stats.json")]:
    src, dst = Path(src), Path(dst)
    if dst.exists() and os.path.samefile(src, dst):
        print("SKIP copy (same file):", dst)
        continue
    shutil.copy2(src, dst)
print("PUBLISHED:", pub / "report.html")
print("REGEN_DONE")
