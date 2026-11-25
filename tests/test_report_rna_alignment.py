#!/usr/bin/env python3
"""Report RNA alignment statistics for pre-aligned BAM inputs."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from annomerge.modules.report import rna_alignment_stats  # noqa: E402


class RnaAlignmentReportTests(unittest.TestCase):
    def test_prealigned_bam_stats_are_merged_into_report(self):
        with tempfile.TemporaryDirectory(prefix="annomerge_rna_report_") as tmp:
            bam_dir = Path(tmp) / "02.rnaseq" / "bam"
            bam_dir.mkdir(parents=True)
            (bam_dir / "sample.sorted.alignment.json").write_text(json.dumps({
                "total": 1000, "mapped": 957, "paired": 1000,
                "proper_pairs": 940, "rate": 95.7,
            }))
            rows = rna_alignment_stats(Path(tmp))
            self.assertEqual(rows, [{
                "sample": "sample", "rate": 95.7,
                "aligner": "pre-aligned BAM",
                "total": 1000, "mapped": 957,
                "paired": 1000, "proper_pairs": 940,
            }])

    def test_bam_counts_extend_hisat2_log_rows(self):
        with tempfile.TemporaryDirectory(prefix="annomerge_rna_report_") as tmp:
            bam_dir = Path(tmp) / "02.rnaseq" / "bam"
            bam_dir.mkdir(parents=True)
            (bam_dir / "sample.hisat2.log").write_text(
                "95.50% overall alignment rate\n")
            (bam_dir / "sample.sorted.alignment.json").write_text(json.dumps({
                "total": 200, "mapped": 191, "paired": 200,
                "proper_pairs": 188, "rate": 95.5,
            }))
            rows = rna_alignment_stats(Path(tmp))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["aligner"], "hisat2")
            self.assertEqual(rows[0]["rate"], 95.5)
            self.assertEqual(rows[0]["total"], 200)
            self.assertEqual(rows[0]["mapped"], 191)
            self.assertEqual(rows[0]["proper_pairs"], 188)


if __name__ == "__main__":
    unittest.main()
