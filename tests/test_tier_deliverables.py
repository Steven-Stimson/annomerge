#!/usr/bin/env python3
"""Tier-specific GFF/CDS/PEP deliverables."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from unittest.mock import patch

from annomerge.modules.report import _extract_sequences, tiered_sequences  # noqa: E402
from annomerge.utils import ToolError  # noqa: E402


class TierDeliverableTests(unittest.TestCase):
    def test_three_bundles_preserve_structure_and_sequences(self):
        with tempfile.TemporaryDirectory(prefix="annomerge_tiers_") as tmp:
            root = Path(tmp)
            genome = root / "genome.fa"
            genome.write_text(
                ">chr1\n" + "ATGAAATAA" + "A" * 20 + "ATGCCCTAA" +
                "A" * 20 + "ATGGGGTAA" + "\n")
            rows = ["##gff-version 3"]
            for i, (tid, tier, start) in enumerate((
                    ("t1", "strong", 1),
                    ("t2", "moderate", 30),
                    ("t3", "weak", 59)), 1):
                gid = f"g{i}"
                end = start + 8
                rows.append(f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gid}")
                rows.append(f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={tid};Parent={gid};Tier={tier}")
                rows.append(f"chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\tParent={tid}")
                rows.append(f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tParent={tid}")
            gff = root / "final.gff3"
            gff.write_text("\n".join(rows) + "\n")
            out = root / "09.result"
            manifest = tiered_sequences(gff, genome, out)
            self.assertEqual([r["tier"] for r in manifest],
                             ["strong", "strong_moderate", "weak"])
            self.assertEqual([r["transcripts"] for r in manifest], [1, 2, 1])
            self.assertEqual([r["genes"] for r in manifest], [1, 2, 1])
            for row in manifest:
                self.assertEqual(row["cds_sequences"], row["transcripts"])
                self.assertEqual(row["pep_sequences"], row["transcripts"])
                for key in ("gff", "cds", "pep"):
                    self.assertTrue(Path(row[key]).is_file())

            strong = (out / "by_tier/strong/annomerge.strong.gff").read_text()
            self.assertIn("ID=t1", strong)
            self.assertNotIn("ID=t2", strong)
            self.assertNotIn("ID=g2", strong)
            self.assertIn("Parent=t1", strong)

            sm = (out / "by_tier/strong_moderate/annomerge.strong_moderate.gff").read_text()
            self.assertIn("ID=t1", sm)
            self.assertIn("ID=t2", sm)
            self.assertNotIn("ID=t3", sm)

            weak_dir = out / "by_tier/weak"
            weak = (weak_dir / "annomerge.weak.gff").read_text()
            self.assertIn("ID=t3", weak)
            self.assertNotIn("ID=t1", weak)
            cds = (weak_dir / "annomerge.weak.cds").read_text()
            pep = (weak_dir / "annomerge.weak.pep").read_text()
            self.assertTrue(cds.startswith(">t3\nATGGGGTAA"))
            self.assertTrue(pep.startswith(">t3\nM"))
            self.assertEqual(len((out / "by_tier/tier_deliverables.tsv").read_text().splitlines()), 4)

    def test_gffread_retry_uses_writable_output_cwd(self):
        """Read-only genome/caller cwd must not control FASTA index creation."""
        with tempfile.TemporaryDirectory(prefix="annomerge_gffread_cwd_") as tmp:
            root = Path(tmp)
            out = root / "writable-output"
            out.mkdir()
            gff = root / "input.gff3"
            genome = root / "genome.fa"
            gff.write_text("##gff-version 3\n")
            genome.write_text(">chr1\nATG\n")
            calls = []

            def fake_run(cmd, **kwargs):
                calls.append(kwargs)
                if len(calls) == 1:
                    raise ToolError("simulated index failure")
                Path(cmd[2]).write_text(">t1\nM\n")
                Path(cmd[4]).write_text(">t1\nATG\n")

            with patch("annomerge.modules.report.which", return_value="gffread"), \
                    patch("annomerge.modules.report.run", side_effect=fake_run):
                counts = _extract_sequences(gff, genome, out, "test")
            self.assertEqual(counts, {"pep": 1, "cds": 1})
            self.assertEqual(len(calls), 2)
            self.assertTrue(all(c["cwd"] == out for c in calls))


if __name__ == "__main__":
    unittest.main()
