#!/bin/bash
# AnnoMerge full pipeline on Arabidopsis thaliana (TAIR10.1) at 32 threads.
# No BUSCO lineage for Brassicales locally -> B_c ablated (weights renormalize).
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
  -o /work/shidishen/AnnoMerge/annomerge/result/fullrun_arabidopsis
echo "FULLRUN_ARABIDOPSIS_EXIT=$?"
