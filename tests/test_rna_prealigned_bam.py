#!/usr/bin/env python3
"""Regression test for junction generation from pre-aligned BAM input."""
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import pysam  # noqa: E402
from annomerge.modules.rna import run_rna_branch  # noqa: E402
from annomerge.modules.validate import RnaSample  # noqa: E402


class PrealignedBamTests(unittest.TestCase):
    def test_bam_only_input_creates_junction_directory_and_beds(self):
        with tempfile.TemporaryDirectory(prefix="annomerge_bam_junction_") as tmp:
            root = Path(tmp)
            genome = root / "genome.fa"
            genome.write_text(">chr1\n" + "A" * 1000 + "\n")
            bam = root / "sample.sorted.bam"
            header = {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [{"SN": "chr1", "LN": 1000}],
            }
            with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
                read = pysam.AlignedSegment()
                read.query_name = "spliced_read"
                read.query_sequence = "A" * 100
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 99
                read.mapping_quality = 60
                read.cigartuples = [(0, 50), (3, 100), (0, 50)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                out.write(read)
            pysam.index(str(bam))

            sample = RnaSample(files=[bam], tag="bam", name_override="sample")
            vi = SimpleNamespace(
                genome=genome, rna_fastqs=[], rna_bams=[sample], rna_gtfs=[],
            )
            outdir = root / "output"

            def fake_assemble(units, gtf_out, threads, extra_opts=None):
                gtf_out.parent.mkdir(parents=True, exist_ok=True)
                gtf_out.write_text("##gff-version 3\n")
                return gtf_out

            with patch("annomerge.modules.rna.assemble", side_effect=fake_assemble):
                result = run_rna_branch(vi, outdir, threads=2)

            bam_dir = outdir / "02.rnaseq" / "bam"
            self.assertTrue(bam_dir.is_dir())
            per_sample = bam_dir / "sample.sorted.junctions.bed"
            self.assertEqual(per_sample.read_text(), "chr1\t149\t249\t.\t1\t+\n")
            self.assertEqual(Path(result["junctions_bed"]).read_text(),
                             "chr1\t149\t249\t1\t1\t+\n")
            self.assertEqual(Path(result["junctions_short_bed"]).read_text(),
                             "chr1\t149\t249\t1\t1\t+\n")
            self.assertIsNone(result["junctions_isoseq_bed"])


if __name__ == "__main__":
    unittest.main()
