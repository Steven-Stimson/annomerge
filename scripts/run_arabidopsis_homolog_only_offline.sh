#!/bin/bash
# Arabidopsis homolog-only OFFLINE run (no RNA): homolog proteins + de novo
# predictors only. isoform mode auto-degrades to EVM-only integration
# (yellow notice); confidence tiers are protein-only (rna=NA); the viewer
# drops goldset/stringtie tracks. No BUSCO lineage (B_c ablated, as v2).
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
AT=/work/shidishen/AnnoMerge/test_data/Arabidopsis_thaliana
python3 -u annomerge merge \
  -g $AT/repeat/genome.fa \
  -P helixer=$AT/denovo/Arabidopsis_helixer.gff,annevo=$AT/denovo/arabidopsis.fa.ANNEVO.gff \
  -p $AT/homologs/homologs.fasta \
  --max-anomaly-rate 0.25 \
  -t 32 \
  -o /work/shidishen/AnnoMerge/annomerge/result/fullrun_arabidopsis_homolog_only
echo "FULLRUN_ARA_HOMOLOG_ONLY_EXIT=$?"
