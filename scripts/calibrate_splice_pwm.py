#!/usr/bin/env python3
"""Calibrate the splice-PWM gate: train on deep-test StringTie introns and
compare real vs 1bp-shifted site scores."""
import sys
sys.path.insert(0, "/work/shidishen/AnnoMerge/annomerge/src")
from annomerge.modules.splicepwm import SplicePWM, _calibrate
from annomerge.modules.sanity import _read_genome

GENOME = "/work/shidishen/AnnoMerge/test_data/Saccharomyces_cerevisiae/genome/GCF_000146045.2_R64_genomic.fna"
GTF = "/work/shidishen/AnnoMerge/annomerge/result/deep_test/rna/stringtie/stringtie.gtf"

genome = _read_genome(GENOME)
# parse GTF exon chains -> introns
tx_exons = {}
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "exon":
            continue
        import re
        m = re.search(r'transcript_id "([^"]+)"', c[8])
        if not m:
            continue
        tx_exons.setdefault((c[0], c[6], m.group(1)), []).append((int(c[3]), int(c[4])))

pwm = SplicePWM()
introns = []
for (seqid, strand, tid), exons in tx_exons.items():
    exons.sort()
    for i in range(len(exons) - 1):
        iv = (seqid, strand, exons[i][1] + 1, exons[i + 1][0] - 1)
        introns.append(iv)
        pwm.add_intron(genome, *iv)

print(f"training introns: {len(introns)} (pwm.n={pwm.n})")
import json
print(json.dumps(_calibrate(pwm, genome, introns), indent=1))
