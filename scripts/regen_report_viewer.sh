#!/bin/bash
# regenerate report.html + evidence browser on fullrun_evi3
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
python3 -u scripts/regen_report_viewer.py
echo "EXIT=$?"
