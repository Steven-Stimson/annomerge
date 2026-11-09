#!/usr/bin/env python3
"""Multi-BAM viewer test: merged + per-sample coverage tracks.

Stages a lightweight copy of the fullrun_dualsets output, regenerates the
evidence browser with THREE BAMs, and verifies:
  1. one non-empty coverage_<sample>.bedgraph per BAM, no name collisions
  2. coverage_all.bedgraph (per-bin summed depth across all samples)
  3. tracks.json carries the coverage metadata (n_samples=3)
  4. browser HTML shows 1 merged + 3 per-sample wig tracks, all autoscale
  5. merged depth at its peak bin == sum of per-sample depths at that bin
Exit code 0 = TEST_PASS.
"""
import json
import re
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

REPO = Path("/work/shidishen/AnnoMerge/annomerge")
sys.path.insert(0, str(REPO / "src"))
from annomerge.modules.viewer import run_viewer  # noqa: E402

OUT = REPO / "result/multisample_test"
SRC = REPO / "result/fullrun_dualsets"
TD = Path("/work/shidishen/AnnoMerge/test_data")
SC = TD / "Saccharomyces_cerevisiae"
BAMS = [
    SRC / "02.rnaseq/bam/SRR31963760.sorted.bam",
    SC / "rnaseq/SRR/ERR1157562.bam",
    SC / "rnaseq/SRR/ERR1157579.bam",
]

fails = []


def check(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg)
    if not ok:
        fails.append(msg)


def load_bg(path):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4:
            d[(f[0], int(f[1]))] = float(f[3])
    return d


def main():
    # ---- stage: copy only the subdirs run_viewer reads ----
    if OUT.exists():
        shutil.rmtree(OUT)
    for sub in ("06.merge", "07.augment", "08.assessment"):
        shutil.copytree(SRC / sub, OUT / sub)

    sets_tsv = SRC / "03.homolog/sets.tsv"
    set_weights = {}
    for line in sets_tsv.read_text().splitlines()[1:]:
        c = line.split("\t")
        set_weights[c[0]] = float(c[6])
    results = {
        "homolog": {"gff": SRC / "03.homolog/miniprot.gff",
                    "set_weights": set_weights},
        "rna": {"stringtie_gtf": SRC / "02.rnaseq/stringtie/stringtie.gtf",
                "junctions_bed": SRC / "02.rnaseq/junctions.bed",
                "bams": BAMS},
    }
    SC = TD / "Saccharomyces_cerevisiae"
    preds = {"helixer": str(SC / "denovo/helixer.gff"),
             "annevo": str(SC / "denovo/annevo.norm.gff3")}
    details = {}
    for name in preds:
        p = SRC / f"05.score/{name}/score_detail.json"
        if p.exists():
            details[name] = json.load(open(p))
    vi = SimpleNamespace(genome=SC / "genome/GCF_000146045.2_R64_genomic.fna")
    vw = run_viewer(OUT, vi, results, preds, details,
                    SRC / "07.augment/augment.sanity.gff3",
                    OUT / "07.augment", OUT / "06.merge", threads=4)
    data = vw["data"]
    html = open(vw["html"]).read()
    tcfg = json.load(open(data / "tracks.json"))

    # 1) per-sample bedgraphs
    cov = tcfg.get("coverage") or {}
    names = [s["name"] for s in cov.get("samples", [])]
    check(len(names) == 3 and len(set(names)) == 3,
          f"3 per-sample tracks with unique names (got {names})")
    for s in cov.get("samples", []):
        p = data / s["file"]
        check(p.exists() and p.stat().st_size > 0,
              f"per-sample bedgraph non-empty: {s['file']}")

    # 2) merged bedgraph
    m = cov.get("merged") or {}
    check(m.get("file") == "coverage_all.bedgraph"
          and (data / "coverage_all.bedgraph").exists(),
          "coverage_all.bedgraph exists for the 3-sample run")

    # 3) metadata
    check(m.get("n_samples") == 3, "tracks.json coverage.n_samples == 3")

    # 4) HTML tracks
    merged_labels = re.findall(r'name: "(RNA coverage \(all[^"]+)', html)
    sample_labels = re.findall(r'name: "RNA coverage: ([^"]+)', html)
    check(len(merged_labels) == 1, f"1 merged track label (got {merged_labels})")
    check(len(sample_labels) == 3, f"3 per-sample labels (got {sample_labels})")
    check(html.count("autoscale: true") >= 4, "all coverage tracks autoscale")
    check('min: 0' in html, "coverage y-axis starts at 0")

    # 5) merged == sum of samples at the merged peak bin
    if (data / "coverage_all.bedgraph").exists() and len(names) == 3:
        mall = load_bg(data / "coverage_all.bedgraph")
        peak_ref, peak_pos = max(mall, key=mall.get)
        per = [load_bg(data / f"coverage_{n}.bedgraph") for n in names]
        ssum = sum(d.get((peak_ref, peak_pos), 0.0) for d in per)
        check(abs(mall[(peak_ref, peak_pos)] - ssum) < 1.0,
              f"merged depth {mall[(peak_ref, peak_pos)]:.1f} == sum of samples "
              f"{ssum:.1f} at {peak_ref}:{peak_pos}")

    print()
    if fails:
        print(f"TEST_FAIL ({len(fails)} check(s))")
        sys.exit(1)
    print("TEST_PASS")


if __name__ == "__main__":
    main()
