#!/bin/bash
# Arabidopsis BUSCO v2 (new pipeline): same as v2 plus --busco
# embryophyta_odb10 (B_c subscore active in W_raw).  32 threads.
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
AT=/work/shidishen/AnnoMerge/test_data/Arabidopsis_thaliana
python3 -u annomerge merge \
  -g $AT/repeat/genome.fa \
  -P helixer=$AT/denovo/Arabidopsis_helixer.gff,annevo=$AT/denovo/arabidopsis.fa.ANNEVO.gff \
  -p $AT/homologs/homologs.fasta \
  -r $AT/rnaseq/arabidopsis.list \
  --busco /work/shidishen/AnnoMerge/test_data/db/busco/embryophyta_odb10 \
  -t 32 \
  -o /work/shidishen/AnnoMerge/annomerge/result/fullrun_arabidopsis_busco_v2
echo "FULLRUN_ARA_BUSCO_V2_EXIT=$?"
