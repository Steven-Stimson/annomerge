#!/bin/bash
# AnnoMerge homolog-only validation run: NO RNA input (-r/-bam/--transcripts).
# Expected behavior: RNA branch + gold set auto-skipped, isoform mode degrades
# to EVM-only integration, and the evidence viewer shows the EVM result
# (final = 06.merge/evm.sanity.gff3).
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
TD=/work/shidishen/AnnoMerge/test_data
SC=$TD/Saccharomyces_cerevisiae
python3 -u annomerge merge \
  -g $SC/genome/GCF_000146045.2_R64_genomic.fna \
  -P helixer=$SC/denovo/helixer.gff,annevo=$SC/denovo/annevo.norm.gff3 \
  -p $SC/homologs_sets/pool.faa \
  --mnp-filter "UnStop,UnShift,Q2R=0.5,topRank=1" \
  --busco $TD/db/busco/saccharomycetes_odb10 \
  -t 32 \
  -o /work/shidishen/AnnoMerge/annomerge/result/homolog_only_yeast
echo "HOMOLOG_ONLY_EXIT=$?"
