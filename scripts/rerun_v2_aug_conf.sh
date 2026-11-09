#!/bin/bash
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
export ANNO_OUT=/work/shidishen/AnnoMerge/annomerge/result/fullrun_dualsets_v2
rm -f $ANNO_OUT/09.result/viewer/data/final.gff3*
python3 -u scripts/rerun_aug_conf.py
echo "RERUN_V2_EXIT=$?"
python3 -u scripts/regen_viewer_dualsets.py
echo "REGEN_V2_VIEWER_EXIT=$?"
