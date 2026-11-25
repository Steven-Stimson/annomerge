#!/usr/bin/env python3
"""Regression tests for --protein list sheets and CLI-supplied species tags."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from annomerge.modules.validate import (  # noqa: E402
    InputError,
    normalize_homolog_species,
    parse_homolog_sets,
    validate,
)


class HomologInputTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="annomerge_protein_input_")
        self.root = Path(self._tmp.name)
        (self.root / "data").mkdir()
        (self.root / "lists").mkdir()
        self.a = self.root / "data" / "a.fa"
        self.b = self.root / "lists" / "b.fa"
        self.a.write_text(">p1 first protein\nMAAA\n")
        self.b.write_text(">p2\nMBBB\n")

    def tearDown(self):
        self._tmp.cleanup()

    def _sheet(self, text="sp1=../data/a.fa\nsp2=b.fa\n"):
        sheet = self.root / "lists" / "proteins.list"
        sheet.write_text(text)
        return sheet

    def _genome_and_prediction(self):
        genome = self.root / "genome.fa"
        pred = self.root / "pred.gff3"
        genome.write_text(">chr1\n" + "A" * 100 + "\n")
        pred.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=g1\n"
        )
        return genome, pred

    def test_list_expands_relative_paths_and_comments(self):
        sheet = self._sheet(
            "# species=protein FASTA\n\n"
            "sp1 = ../data/a.fa\n"
            "sp2=b.fa\n"
        )
        sets = parse_homolog_sets(str(sheet))
        self.assertEqual(list(sets), ["sp1", "sp2"])
        self.assertEqual(sets["sp1"], [self.a.resolve()])
        self.assertEqual(sets["sp2"], [self.b.resolve()])

    def test_list_rejects_missing_id_or_path(self):
        bad = self._sheet("a.fa\n")
        with self.assertRaisesRegex(InputError, "expected ID=path"):
            parse_homolog_sets(str(bad))
        bad.write_text("sp1=\n")
        with self.assertRaisesRegex(InputError, "non-empty ID=path"):
            parse_homolog_sets(str(bad))

    def test_list_rejects_duplicate_ids(self):
        sheet = self._sheet("sp1=../data/a.fa\nsp1=b.fa\n")
        with self.assertRaisesRegex(InputError, "duplicate homolog set label 'sp1'"):
            parse_homolog_sets(str(sheet))

    def test_species_normalization_preserves_source_and_existing_tags(self):
        tagged = self.root / "tagged.fa"
        original = ">p1 description\nMAAA\n>p2|species=old kept\nMBBB\n"
        tagged.write_text(original)
        out = self.root / "normalized.fa"
        normalize_homolog_species(tagged, "new", out)
        self.assertEqual(tagged.read_text(), original)
        self.assertEqual(
            out.read_text(),
            ">p1|species=new description\nMAAA\n"
            ">p2|species=old kept\nMBBB\n",
        )

    def test_validate_list_builds_weighted_sets_from_plain_headers(self):
        sheet = self._sheet()
        genome, pred = self._genome_and_prediction()
        vi = validate(
            genome=str(genome), homologs=str(sheet),
            predictors=f"denovo={pred}", workdir=str(self.root / "run"),
        )
        self.assertEqual(set(vi.homolog_sets), {"sp1", "sp2"})
        self.assertIsNotNone(vi.homolog_membership)
        self.assertEqual(
            vi.homologs.read_text(),
            ">p1|species=sp1 first protein\nMAAA\n"
            ">p2|species=sp2\nMBBB\n",
        )
        self.assertEqual(
            Path(vi.homolog_membership).read_text(),
            "protein\tset\np1|species=sp1\tsp1\np2|species=sp2\tsp2\n",
        )
        # Source files are read-only inputs and retain their plain headers.
        self.assertEqual(self.a.read_text(), ">p1 first protein\nMAAA\n")
        self.assertEqual(self.b.read_text(), ">p2\nMBBB\n")

    def test_single_inline_label_also_accepts_plain_headers(self):
        genome, pred = self._genome_and_prediction()
        vi = validate(
            genome=str(genome), homologs=f"species_a={self.a}",
            predictors=f"denovo={pred}", workdir=str(self.root / "single"),
        )
        self.assertEqual(list(vi.homolog_sets), ["species_a"])
        self.assertIn(">p1|species=species_a first protein\n", vi.homologs.read_text())

    def test_unlabeled_legacy_input_still_requires_species_headers(self):
        genome, pred = self._genome_and_prediction()
        with self.assertRaisesRegex(InputError, r"lack '\|species=Genus_species' tag"):
            validate(
                genome=str(genome), homologs=str(self.a),
                predictors=f"denovo={pred}", workdir=str(self.root / "legacy_bad"),
            )
        self.a.write_text(">p1|species=old\nMAAA\n")
        vi = validate(
            genome=str(genome), homologs=str(self.a),
            predictors=f"denovo={pred}", workdir=str(self.root / "legacy_ok"),
        )
        self.assertIn(">p1|species=old\n", vi.homologs.read_text())


if __name__ == "__main__":
    unittest.main()
