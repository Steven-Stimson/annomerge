#!/bin/bash
# AnnoMerge yeast validation run (v3): single protein pool (near+far merged,
# unlabeled -> pool mode), mnp filters at default (UnStop,UnShift,Q2R=0.5,
# topRank=1), isoform collapse tol 5bp, Name-free IDs, UTR injection,
# tabix lazy viewer.  Fresh full run at 32 threads.
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
  --mnp-filter "UnStop,UnShift,Q2R=0.5,topRank=1" \
  -t 32 \
  -o /work/shidishen/AnnoMerge/annomerge/result/fullrun_yeast_v3
echo "FULLRUN_YEAST_V3_EXIT=$?"
