"""Regression: miniprot progress must never contaminate PAF/GFF data."""
import sys
from pathlib import Path

from annomerge.modules import homolog


def test_run_miniprot_separates_stderr_from_data(tmp_path, monkeypatch):
    fake = tmp_path / "fake_miniprot.py"
    fake.write_text(f"""#!{sys.executable}
import sys
mode = sys.argv[2]
if mode == '--paf':
    print('q1\\t100\\t0\\t100\\t+\\tchr1\\t1000\\t0\\t300\\t300\\t300\\t60')
else:
    print('##gff-version 3')
    print('chr1\\tminiprot\\tmRNA\\t1\\t300\\t300\\t+\\t.\\tID=MP1;Rank=1;Identity=1.0;Positive=1.0;Target=q1 1 100')
    print('chr1\\tminiprot\\tCDS\\t1\\t300\\t300\\t+\\t0\\tParent=MP1;Target=q1 1 100')
print('[M::worker_pipeline] mapped 1 sequences', file=sys.stderr)
""")
    fake.chmod(0o755)
    monkeypatch.setattr(homolog, "which", lambda _name: str(fake))
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "A" * 1000 + "\n")
    proteins = tmp_path / "proteins.fa"
    proteins.write_text(">q1\n" + "M" * 100 + "\n")
    out = tmp_path / "out"

    paf, gff = homolog.run_miniprot(
        genome, proteins, out, 1, homolog.parse_mnp_filter(""))

    assert paf.read_text().splitlines() == [
        "q1\t100\t0\t100\t+\tchr1\t1000\t0\t300\t300\t300\t60"]
    raw_gff = (out / "miniprot.raw.gff").read_text()
    assert "Target=q1 1 100" in raw_gff
    assert "worker_pipeline" not in raw_gff
    assert "worker_pipeline" not in (out / "miniprot.raw.paf").read_text()
    assert "mapped 1 sequences" in (out / "miniprot.raw.paf.log").read_text()
    assert "mapped 1 sequences" in (out / "miniprot.raw.gff.log").read_text()
