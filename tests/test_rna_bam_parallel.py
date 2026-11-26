#!/usr/bin/env python3
"""Parallel BAM scanning must be identical to the serial implementation."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import pysam  # noqa: E402
from annomerge.modules.rna import scan_bams  # noqa: E402


def make_bam(path: Path) -> None:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        primary = pysam.AlignedSegment()
        primary.query_name = "primary"
        primary.query_sequence = "A" * 100
        primary.flag = 99  # paired + proper pair
        primary.reference_id = 0
        primary.reference_start = 99
        primary.mapping_quality = 60
        primary.cigartuples = [(0, 50), (3, 100), (0, 50)]
        primary.query_qualities = pysam.qualitystring_to_array("I" * 100)
        out.write(primary)

        supp = pysam.AlignedSegment()
        supp.query_name = "supplementary"
        supp.query_sequence = "A" * 100
        supp.flag = 2048
        supp.reference_id = 0
        supp.reference_start = 99
        supp.mapping_quality = 60
        supp.cigartuples = [(0, 50), (3, 100), (0, 50)]
        supp.query_qualities = pysam.qualitystring_to_array("I" * 100)
        out.write(supp)

        unmapped = pysam.AlignedSegment()
        unmapped.query_name = "unmapped"
        unmapped.query_sequence = "A" * 100
        unmapped.flag = 4
        unmapped.reference_id = -1
        unmapped.reference_start = -1
        unmapped.mapping_quality = 0
        unmapped.cigartuples = []
        unmapped.query_qualities = pysam.qualitystring_to_array("I" * 100)
        out.write(unmapped)


class ParallelBamScanTests(unittest.TestCase):
    def test_serial_parallel_and_legacy_cache_upgrade(self):
        with tempfile.TemporaryDirectory(prefix="annomerge_bam_parallel_") as tmp:
            root = Path(tmp)
            bams = []
            for i in range(3):
                bam = root / f"sample{i}.sorted.bam"
                make_bam(bam)
                bams.append(bam)
            platform = {b: "short" for b in bams}
            source = {b: "pre-aligned BAM" for b in bams}
            outputs = {}
            for threads in (1, 4):
                out_dir = root / f"out{threads}"
                beds, _ = scan_bams(bams, out_dir, platform, source, threads)
                outputs[threads] = beds
                self.assertEqual([p.name for p in beds],
                                 [f"sample{i}.sorted.junctions.bed" for i in range(3)])
                for bed in beds:
                    self.assertEqual(bed.read_text(),
                                     "chr1\t149\t249\t.\t2\t+\n")
                    stats_path = out_dir / bed.name.replace(
                        ".junctions.bed", ".alignment.json")
                    stats = json.loads(stats_path.read_text())
                    self.assertEqual(stats["alignment_source"], "pre-aligned BAM")
                    self.assertEqual(stats["primary_reads_total"], 2)
                    self.assertEqual(stats["primary_reads_mapped"], 1)
                    self.assertEqual(stats["primary_reads_unmapped"], 1)
                    self.assertEqual(stats["primary_reads_paired"], 1)
                    self.assertEqual(stats["primary_reads_properly_paired"], 1)
                    self.assertEqual(stats["primary_alignment_rate_pct"], 50.0)
                self.assertFalse(list(out_dir.glob("*.tmp.*")))
            for serial, parallel in zip(outputs[1], outputs[4]):
                self.assertEqual(serial.read_bytes(), parallel.read_bytes())
                serial_stats = serial.with_name(
                    serial.name.replace(".junctions.bed", ".alignment.json"))
                parallel_stats = parallel.with_name(
                    parallel.name.replace(".junctions.bed", ".alignment.json"))
                self.assertEqual(serial_stats.read_bytes(),
                                 parallel_stats.read_bytes())

            # Old two-field-era JSON is readable by the report, but the RNA
            # scanner must refresh it once to the complete new schema.
            legacy = outputs[4][0].with_name(
                outputs[4][0].name.replace(".junctions.bed", ".alignment.json"))
            legacy.write_text(json.dumps({"total": 2, "mapped": 1, "rate": 50.0}))
            scan_bams(bams, root / "out4", platform, source, 4)
            upgraded = json.loads(legacy.read_text())
            self.assertIn("sample_name", upgraded)
            self.assertIn("primary_reads_unmapped", upgraded)
            self.assertIn("primary_alignment_rate_pct", upgraded)


if __name__ == "__main__":
    unittest.main()
