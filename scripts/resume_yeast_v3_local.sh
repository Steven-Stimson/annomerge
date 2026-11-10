#!/bin/bash
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
TD=/work/shidishen/AnnoMerge/test_data
SC=$TD/Saccharomyces_cerevisiae
python3 -u annomerge merge \
  -g $SC/genome/GCF_000146045.2_R64_genomic.fna \
  -P helixer=$SC/denovo/helixer.gff,annevo=$SC/denovo/annevo.norm.gff3 \
  -p $SC/homologs_sets/pool.faa \
  -r $SC/rnaseq_deep/rnaseq_deep.list \
  --busco $TD/db/busco/saccharomycetes_odb10 \
  -t 4 \
  -o /work/shidishen/AnnoMerge/annomerge/result/fullrun_yeast_v3 2>&1
echo "YEAST_V3_RESUME_EXIT=$?"
