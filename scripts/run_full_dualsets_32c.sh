#!/bin/bash
# AnnoMerge full pipeline (dual homolog sets) at 32 threads, fresh output
# dir — robustness/validity run for the isoform-collapse + Name-strip changes.
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
TD=/work/shidishen/AnnoMerge/test_data
SC=$TD/Saccharomyces_cerevisiae
python3 -u annomerge merge \
  -g $SC/genome/GCF_000146045.2_R64_genomic.fna \
  -P helixer=$SC/denovo/helixer.gff,annevo=$SC/denovo/annevo.norm.gff3 \
  -p near=$SC/homologs/homologs.1.faa,far=$SC/homologs_sets/ncrassa.faa \
  -r $SC/rnaseq_deep/rnaseq_deep.list \
  --busco $TD/db/busco/saccharomycetes_odb10 \
  -t 32 \
  -o /work/shidishen/AnnoMerge/annomerge/result/fullrun_dualsets_v2
echo "FULLRUN_V2_EXIT=$?"
