#!/usr/bin/env python3
"""Re-run 07.augment -> sanity -> renamer -> 08.confidence on an existing
01-06 layout (fullrun_dualsets), then regenerate the viewer. Used to apply
the isoform-collapse + Name-strip changes without re-running 01-06."""
import json
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

REPO = Path("/work/shidishen/AnnoMerge/annomerge")
sys.path.insert(0, str(REPO / "src"))
from annomerge.modules.augment import run_augment          # noqa: E402
from annomerge.modules.confidence import confidence_annotation  # noqa: E402
from annomerge.modules.renamer import normalize_ids        # noqa: E402
from annomerge.modules.sanity import sanity_filter         # noqa: E402
from annomerge.modules.viewer import run_viewer            # noqa: E402

OUT = REPO / "result/fullrun_dualsets"
TD = Path("/work/shidishen/AnnoMerge/test_data")
SC = TD / "Saccharomyces_cerevisiae"
vi = SimpleNamespace(genome=SC / "genome/GCF_000146045.2_R64_genomic.fna")

details, wraw = {}, {}
for name in ("helixer", "annevo"):
    d = json.load(open(OUT / f"05.score/{name}/score_detail.json"))
    details[name] = d
    wraw[name] = None if d.get("hard_reject") else d["W_raw"]
preds = {"helixer": str(SC / "denovo/helixer.gff"),
         "annevo": str(SC / "denovo/annevo.norm.gff3")}
scores = {n: (wraw.get(n) or 0.0,
              (details.get(n, {}).get("subscores") or {}).get("R_s") or 0.0)
          for n in preds}

for d in (OUT / "07.augment", OUT / "08.assessment"):
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)

aug = run_augment(
    vi, OUT / "07.augment",
    evm_gff3=OUT / "06.merge/evm.sanity.gff3",
    goldset_gff=OUT / "06.merge/goldset.genome.gff3",
    stringtie_gtf=OUT / "02.rnaseq/stringtie/stringtie.gtf",
    preds=preds, scores=scores,
    jitter_merge_bp=50, dedup_mode="geometric",
    junctions_bed=OUT / "02.rnaseq/junctions.bed",
    tpm_ratio=0.1, min_tpm=0.1, junc_min_reads=3, pwm_gate=-3.0)
print("AUG_OK", aug["augment_gff3"])

sanity_filter(aug["augment_gff3"], vi.genome,
              aug["dir"] / "augment.sanity.gff3",
              aug["dir"] / "sanity_filter.tsv")
print("SANITY_OK")
print("RENAME", normalize_ids(aug["dir"] / "augment.sanity.gff3"))

# homolog relatedness weights (multi-set) for confidence protein weighting
set_w, membership = {}, {}
for line in (OUT / "03.homolog/sets.tsv").read_text().splitlines()[1:]:
    c = line.split("\t")
    set_w[c[0]] = float(c[6])
mm = OUT / "03.homolog/homologs.merged.faa.sets.tsv"
if mm.exists():
    for line in mm.read_text().splitlines()[1:]:
        c = line.split("\t")
        if len(c) >= 2:
            membership[c[0]] = c[1]
protein_weights = ({p: set_w.get(lab, 1.0) for p, lab in membership.items()}
                   if len(set_w) > 1 else None)
conf_stats = confidence_annotation(
    aug["dir"] / "augment.sanity.gff3", vi.genome,
    OUT / "02.rnaseq/junctions.bed",
    OUT / "06.merge/protein.gff3",
    [OUT / "02.rnaseq/bam/SRR31963760.sorted.bam"],
    OUT / "06.merge/goldset.genome.gff3",
    aug["dir"] / "sanity_filter.tsv",
    OUT / "08.assessment/confidence.tsv",
    protein_weights=protein_weights)
print("CONF_OK", conf_stats)
