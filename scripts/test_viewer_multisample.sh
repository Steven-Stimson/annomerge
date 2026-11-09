#!/bin/bash
# Multi-BAM evidence-browser coverage-track test (offline DCS job).
set -uo pipefail
cd /work/shidishen/AnnoMerge/annomerge
source env.sh
python3 -u scripts/test_viewer_multisample.py
echo "MULTISAMPLE_VIEWER_TEST_EXIT=$?"
