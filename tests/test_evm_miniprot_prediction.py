#!/usr/bin/env python3
"""Tests for miniprot as an EVM structural prediction source."""
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from annomerge.modules.evm import miniprot_to_prediction_gff, run_evm  # noqa: E402


def _rows(path):
    return [
        line.rstrip("\n").split("\t")
        for line in Path(path).read_text().splitlines()
        if line and not line.startswith("#")
    ]


class MiniprotPredictionTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="annomerge_mnp_pred_")
        self.root = Path(self._tmp.name)
        self.miniprot = self.root / "miniprot.gff3"
        self.miniprot.write_text(
            "##gff-version 3\n"
            "chr1\tminiprot\tmRNA\t101\t400\t.\t+\t.\tID=p1;Target=protA 1 99\n"
            "chr1\tminiprot\tCDS\t101\t400\t.\t+\t0\tParent=p1\n"
            "chr1\tminiprot\tstop_codon\t398\t400\t.\t+\t0\tParent=p1\n"
            "chr1\tminiprot\tmRNA\t1001\t1200\t.\t-\t.\tID=n1;Target=protB 1 50\n"
            "chr1\tminiprot\tCDS\t1101\t1200\t.\t-\t0\tParent=n1\n"
            "chr1\tminiprot\tCDS\t1001\t1050\t.\t-\t0\tParent=n1\n"
            "chr1\tminiprot\tstop_codon\t1001\t1003\t.\t-\t0\tParent=n1\n"
        )

    def tearDown(self):
        self._tmp.cleanup()

    def test_converter_builds_models_and_preserves_evm_cds_spans(self):
        out = self.root / "prediction.gff3"
        miniprot_to_prediction_gff(self.miniprot, out)
        rows = _rows(out)
        self.assertEqual({r[1] for r in rows}, {"miniprot"})
        self.assertEqual(sum(r[2] == "gene" for r in rows), 2)
        self.assertEqual(sum(r[2] == "mRNA" for r in rows), 2)
        self.assertFalse(any(r[2] == "stop_codon" for r in rows))

        plus_cds = [r for r in rows if r[2] == "CDS" and "Parent=miniprot.p1" in r[8]]
        plus_exon = [r for r in rows if r[2] == "exon" and "Parent=miniprot.p1" in r[8]]
        self.assertEqual([(int(r[3]), int(r[4])) for r in plus_cds], [(101, 400)])
        self.assertEqual([(int(r[3]), int(r[4])) for r in plus_exon], [(101, 400)])

        minus_cds = [r for r in rows if r[2] == "CDS" and "Parent=miniprot.n1" in r[8]]
        minus_exon = [r for r in rows if r[2] == "exon" and "Parent=miniprot.n1" in r[8]]
        self.assertEqual(
            sorted((int(r[3]), int(r[4])) for r in minus_cds),
            [(1001, 1050), (1101, 1200)],
        )
        self.assertEqual(
            sorted((int(r[3]), int(r[4])) for r in minus_exon),
            [(1001, 1050), (1101, 1200)],
        )

    def test_run_evm_combines_denovo_miniprot_and_goldset_with_hierarchy(self):
        genome = self.root / "genome.fa"
        genome.write_text(">chr1\n" + "A" * 2000 + "\n")
        denovo = self.root / "denovo.gff3"
        denovo.write_text(
            "##gff-version 3\n"
            "chr1\tdenovo\tgene\t101\t400\t.\t+\t.\tID=d1.gene\n"
            "chr1\tdenovo\tmRNA\t101\t400\t.\t+\t.\tID=d1;Parent=d1.gene\n"
            "chr1\tdenovo\texon\t101\t400\t.\t+\t.\tID=d1.exon;Parent=d1\n"
            "chr1\tdenovo\tCDS\t101\t400\t.\t+\t0\tParent=d1\n"
        )
        gold_dir = self.root / "goldset_input"
        gold_dir.mkdir()
        (gold_dir / "transcripts.fa.transdecoder.gff3").write_text("##gff-version 3\n")
        lifted_gold = self.root / "goldset.lifted.gff3"
        lifted_gold.write_text(
            "##gff-version 3\n"
            "chr1\tgoldset\tgene\t501\t800\t.\t+\t.\tID=g1.gene\n"
            "chr1\tgoldset\tmRNA\t501\t800\t.\t+\t.\tID=g1;Parent=g1.gene\n"
            "chr1\tgoldset\texon\t501\t800\t.\t+\t.\tID=g1.exon;Parent=g1\n"
            "chr1\tgoldset\tCDS\t501\t800\t.\t+\t0\tParent=g1\n"
        )
        fake_evm = self.root / "fake.evm.gff3"
        fake_evm.write_text(
            "##gff-version 3\n"
            "chr1\tEVM\tgene\t101\t400\t.\t+\t.\tID=e1\n"
        )
        outdir = self.root / "06.merge"
        vi = SimpleNamespace(genome=genome)
        with (
            patch("annomerge.modules.evm.lift_goldset", return_value=lifted_gold),
            patch("annomerge.modules.evm._evm_run_one",
                  return_value={"chr1": fake_evm}) as run_one,
        ):
            result = run_evm(
                vi, outdir, {"denovo": (denovo, 5)},
                stringtie_gtf=None, miniprot_gff=self.miniprot,
                goldset_dir=gold_dir, threads=1,
            )

        weights = {}
        for line in (outdir / "weights.txt").read_text().splitlines():
            kind, source, weight = line.split("\t")
            weights[(kind, source)] = int(weight)
        self.assertEqual(weights[("ABINITIO_PREDICTION", "denovo")], 5)
        self.assertEqual(weights[("ABINITIO_PREDICTION", "miniprot")], 6)
        self.assertEqual(weights[("ABINITIO_PREDICTION", "goldset")], 10)
        self.assertEqual(weights[("PROTEIN", "miniprot")], 6)
        self.assertGreater(weights[("ABINITIO_PREDICTION", "miniprot")],
                           weights[("ABINITIO_PREDICTION", "denovo")])
        self.assertLess(weights[("ABINITIO_PREDICTION", "miniprot")],
                        weights[("ABINITIO_PREDICTION", "goldset")])

        combined = _rows(outdir / "pred.all.gff3")
        self.assertEqual({r[1] for r in combined}, {"denovo", "miniprot", "goldset"})
        self.assertTrue((outdir / "miniprot.prediction.gff3").exists())
        self.assertTrue((outdir / "protein.gff3").exists())
        self.assertIn("miniprot", run_one.call_args.args[2])
        self.assertTrue(Path(result["evm_gff3"]).exists())


if __name__ == "__main__":
    unittest.main()
