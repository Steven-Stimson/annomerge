#!/usr/bin/env python3
"""Report RNA alignment statistics exclusively from BAM-derived JSON."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from markupsafe import Markup  # noqa: E402
from annomerge.modules.report import (render_html, rna_alignment_stats,  # noqa: E402
                                      svg_bars)


class RnaAlignmentReportTests(unittest.TestCase):
    def test_current_schema_keeps_metrics_but_renders_only_rate_chart(self):
        with tempfile.TemporaryDirectory(prefix="annomerge_rna_report_") as tmp:
            bam_dir = Path(tmp) / "02.rnaseq" / "bam"
            bam_dir.mkdir(parents=True)
            (bam_dir / "sample.sorted.alignment.json").write_text(json.dumps({
                "sample_name": "sample", "alignment_source": "pre-aligned BAM",
                "primary_reads_total": 1000, "primary_reads_mapped": 957,
                "primary_reads_unmapped": 43, "primary_reads_paired": 1000,
                "primary_reads_properly_paired": 940,
                "primary_alignment_rate_pct": 95.7,
            }))
            rows = rna_alignment_stats(Path(tmp))
            self.assertEqual(rows, [{
                "sample_name": "sample", "alignment_source": "pre-aligned BAM",
                "primary_reads_total": 1000, "primary_reads_mapped": 957,
                "primary_reads_unmapped": 43, "primary_reads_paired": 1000,
                "primary_reads_properly_paired": 940,
                "primary_alignment_rate_pct": 95.7,
            }])
            report = Path(tmp) / "report.html"
            chart = Markup(svg_bars([
                ("sample (pre-aligned BAM)", rows[0]["primary_alignment_rate_pct"])
            ]))
            render_html({"workflow": {"rna_align": rows, "rna_chart": chart,
                                      "protein": {"n_mapped": 0},
                                      "goldset": {}},
                         "versions": {"rows": []}}, report)
            html = report.read_text()
            section = html.split('<h3>2.1 RNA alignment</h3>', 1)[1].split(
                '<h3>2.2 Homolog protein alignment</h3>', 1)[0]
            self.assertIn("<svg", section)
            self.assertIn("95.7%", section)
            self.assertNotIn("<table", section)
            self.assertNotIn("Properly paired primary reads", section)
            self.assertNotIn("read-pair fragments", section)

    def test_previous_schema_is_mapped_and_hisat2_log_rate_ignored(self):
        with tempfile.TemporaryDirectory(prefix="annomerge_rna_report_") as tmp:
            bam_dir = Path(tmp) / "02.rnaseq" / "bam"
            bam_dir.mkdir(parents=True)
            (bam_dir / "sample.hisat2.log").write_text(
                "1.00% overall alignment rate\n")
            (bam_dir / "sample.sorted.alignment.json").write_text(json.dumps({
                "sample": "sample", "source": "hisat2",
                "primary_total": 200, "mapped": 191, "unmapped": 9,
                "paired": 200, "proper_pairs": 188,
                "alignment_rate": 95.5,
            }))
            row = rna_alignment_stats(Path(tmp))[0]
            self.assertEqual(row["alignment_source"], "hisat2")
            self.assertEqual(row["primary_alignment_rate_pct"], 95.5)
            self.assertEqual(row["primary_reads_total"], 200)
            self.assertEqual(row["primary_reads_mapped"], 191)
            self.assertEqual(row["primary_reads_unmapped"], 9)
            self.assertEqual(row["primary_reads_properly_paired"], 188)

    def test_oldest_json_needs_no_minimap2_log_or_bam_rescan(self):
        with tempfile.TemporaryDirectory(prefix="annomerge_rna_report_") as tmp:
            bam_dir = Path(tmp) / "02.rnaseq" / "bam"
            bam_dir.mkdir(parents=True)
            (bam_dir / "sample.minimap2.log").write_text("minimap2 log\n")
            (bam_dir / "sample.sorted.alignment.json").write_text(json.dumps({
                "total": 100, "mapped": 90, "paired": 0,
                "proper_pairs": 0, "rate": 90.0,
            }))
            rows = rna_alignment_stats(Path(tmp))
            self.assertEqual(rows, [{
                "sample_name": "sample", "alignment_source": "BAM",
                "primary_reads_total": 100, "primary_reads_mapped": 90,
                "primary_reads_unmapped": 10, "primary_reads_paired": 0,
                "primary_reads_properly_paired": 0,
                "primary_alignment_rate_pct": 90.0,
            }])


if __name__ == "__main__":
    unittest.main()
