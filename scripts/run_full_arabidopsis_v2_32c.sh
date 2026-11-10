#!/bin/bash
# Arabidopsis v2 (new pipeline): single homolog pool + mnp filters default,
# isoform collapse tol 5bp, Name-free IDs, UTR injection, tabix lazy viewer.
# No BUSCO lineage (B_c ablated).  32 threads.
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
AT=/work/shidishen/AnnoMerge/test_data/Arabidopsis_thaliana
python3 -u annomerge merge \
  -g $AT/repeat/genome.fa \
  -P helixer=$AT/denovo/Arabidopsis_helixer.gff,annevo=$AT/denovo/arabidopsis.fa.ANNEVO.gff \
  -p $AT/homologs/homologs.fasta \
  -r $AT/rnaseq/arabidopsis.list \
  -t 32 \
  -o /work/shidishen/AnnoMerge/annomerge/result/fullrun_arabidopsis_v2
echo "FULLRUN_ARA_V2_EXIT=$?"
