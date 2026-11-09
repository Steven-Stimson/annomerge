#!/bin/bash
# Re-run 07.augment -> 08.confidence with proper tool PATH, then regen viewer.
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
python3 -u scripts/rerun_aug_conf.py
echo "RERUN_AUG_CONF_EXIT=$?"
python3 -u scripts/regen_viewer_dualsets.py
echo "REGEN_VIEWER_EXIT=$?"
