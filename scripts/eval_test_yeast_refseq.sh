#!/bin/bash
# E2E smoke test: annomerge evaluate on the yeast RefSeq annotation.
# Evidence layers (02/03/04) are pre-seeded via symlinks from the finished
# fullrun_yeast_v3 merge run, so only the NEW evaluate-specific stages run:
# gffread normalization, sanity, goldset lift, confidence, tags, report, viewer.
set -x
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
SC=/work/shidishen/AnnoMerge/test_data/Saccharomyces_cerevisiae
E=result/eval_test_yeast_refseq
mkdir -p "$E"
ln -sfn ../fullrun_yeast_v3/02.rnaseq "$E/02.rnaseq"
ln -sfn ../fullrun_yeast_v3/03.homolog "$E/03.homolog"
ln -sfn ../fullrun_yeast_v3/04.goldset "$E/04.goldset"
python3 -m annomerge.cli evaluate \
  -g "$SC/genome/GCF_000146045.2_R64_genomic.fna" \
  -a "$SC/annotation/GCF_000146045.2_R64_genomic.gff" \
  -p "$SC/homologs_sets/pool.faa" \
  -P "helixer=$SC/denovo/helixer.gff,annevo=$SC/denovo/annevo.norm.gff3" \
  --bam result/fullrun_yeast_v3/02.rnaseq/bam/SRR31963760.sorted.bam \
  --transcripts result/fullrun_yeast_v3/02.rnaseq/stringtie/stringtie.gtf \
  --busco /work/shidishen/AnnoMerge/test_data/db/busco/saccharomycetes_odb10 \
  -t 4 -o "$E"
EXIT=$?
echo "EVAL_TEST_EXIT=$EXIT"
ls -la "$E" "$E/06.evaluate" "$E/08.assessment" "$E/09.result" 2>/dev/null
