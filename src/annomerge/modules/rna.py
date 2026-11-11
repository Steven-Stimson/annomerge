"""Module 2 - RNA branch: HISAT2 alignment -> samtools -> StringTie -> junctions.bed.

Lazy/resumable: existing non-empty outputs are skipped.
Outputs (under <outdir>/rna/):
  hisat2_index/            HISAT2 index
  bam/<sample>.sorted.bam  sorted+indexed BAM per sample
  stringtie/stringtie.gtf  merged assembly (all samples)
  junctions.bed            per-junction read support, merged over samples
"""
from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

import pysam

from ..utils import ensure_dir, lazy, run, which
from .validate import ValidatedInputs

log = logging.getLogger("annomerge")


def build_index(genome: Path, index_dir: Path, threads: int) -> Path:
    prefix = index_dir / "genome"
    marker = Path(str(prefix) + ".1.ht2")
    if not lazy([marker], desc="hisat2-build"):
        ensure_dir(index_dir)
        run([which("hisat2-build"), "-p", str(threads), str(genome), str(prefix)],
            desc="hisat2-build")
    return prefix


def align_sample(sample, index_prefix: Path, bam_dir: Path, threads: int,
                 extra_opts: list[str] | None = None) -> Path:
    """HISAT2 -> samtools sort -> index. Returns sorted BAM path."""
    ensure_dir(bam_dir)
    bam = bam_dir / f"{sample.name}.sorted.bam"
    if lazy([bam, Path(str(bam) + ".bai")], desc=f"align {sample.name}"):
        return bam
    sam = bam_dir / f"{sample.name}.sam"
    cmd = [which("hisat2"), "-x", str(index_prefix), "-p", str(threads),
           "--dta", "-S", str(sam)]  # --dta: downstream transcriptome assembly
    cmd += list(extra_opts or [])
    if len(sample.files) == 2:
        cmd += ["-1", str(sample.files[0]), "-2", str(sample.files[1])]
    else:
        cmd += ["-U", str(sample.files[0])]
    run(cmd, log_file=bam_dir / f"{sample.name}.hisat2.log", desc=f"hisat2 {sample.name}")
    run([which("samtools"), "sort", "-@", str(threads), "-o", str(bam), str(sam)],
        desc=f"samtools sort {sample.name}")
    sam.unlink()
    run([which("samtools"), "index", str(bam)], desc=f"samtools index {sample.name}")
    return bam


def extract_junctions(bam: Path, out_bed: Path, min_reads: int = 1) -> Path:
    """Extract splice junctions (CIGAR N blocks) with per-junction read counts."""
    if lazy([out_bed], desc=f"junctions {bam.stem}"):
        return out_bed
    counts: dict[tuple[str, int, int, str], int] = defaultdict(int)
    with pysam.AlignmentFile(str(bam), "rb") as af:
        for read in af.fetch(until_eof=True):
            if read.is_unmapped or not read.cigartuples:
                continue
            pos = read.reference_start
            strand = "-" if read.is_reverse else "+"
            for op, length in read.cigartuples:
                if op == 3:  # N: skipped region = intron
                    counts[(read.reference_name, pos, pos + length, strand)] += 1
                    pos += length
                elif op in (0, 7, 8):  # M, =, X
                    pos += length
                elif op == 2:  # D
                    pos += length
    with open(out_bed, "w") as fh:
        for (chrom, s, e, strand), n in sorted(counts.items()):
            if n >= min_reads:
                fh.write(f"{chrom}\t{s}\t{e}\t.\t{n}\t{strand}\n")
    log.info("junctions: %d (>= %d reads) from %s", sum(1 for v in counts.values() if v >= min_reads),
             min_reads, bam.name)
    return out_bed


def merge_junctions(beds: list[Path], out_bed: Path) -> Path:
    """Merge per-sample junction BEDs, summing read support; track sample count."""
    if lazy([out_bed], desc="merge junctions"):
        return out_bed
    agg: dict[tuple[str, int, int], list] = defaultdict(lambda: [0, set(), set()])
    for bed in beds:
        sample = bed.stem.replace(".junctions", "")
        with open(bed) as fh:
            for line in fh:
                chrom, s, e, _, n, strand = line.rstrip("\n").split("\t")[:6]
                rec = agg[(chrom, int(s), int(e))]
                rec[0] += int(n)
                rec[1].add(sample)
                rec[2].add(strand)
    with open(out_bed, "w") as fh:
        for (chrom, s, e), (n, samples, strands) in sorted(agg.items()):
            strand = next(iter(strands)) if len(strands) == 1 else "."
            fh.write(f"{chrom}\t{s}\t{e}\t{len(samples)}\t{n}\t{strand}\n")
    log.info("merged junctions: %d sites from %d samples", len(agg), len(beds))
    return out_bed


def assemble(bams: list[Path], gtf_out: Path, threads: int,
             extra_opts: list[str] | None = None,
             long_read: set[Path] | None = None) -> Path:
    """StringTie per sample, then stringtie --merge.

    extra_opts are inserted into every per-sample assembly call (protected:
    input BAM, -p, -o, -L — isoseq samples get -L automatically).
    The --merge step is pipeline-managed and not user-tunable.
    """
    if lazy([gtf_out], desc="stringtie"):
        return gtf_out
    sd = gtf_out.parent
    ensure_dir(sd)
    per_sample = []
    for bam in bams:
        g = sd / f"{bam.stem.replace('.sorted', '')}.gtf"
        if not lazy([g], desc=f"stringtie {bam.stem}"):
            cmd = [which("stringtie"), str(bam), "-p", str(threads), "-o", str(g)]
            if long_read and bam in long_read:
                cmd.append("-L")
            cmd += list(extra_opts or [])
            run(cmd, log_file=sd / f"{bam.stem}.stringtie.log", desc=f"stringtie {bam.stem}")
        per_sample.append(g)
    if len(per_sample) == 1:
        per_sample[0].rename(gtf_out)
    else:
        lst = sd / "gtf_list.txt"
        lst.write_text("\n".join(map(str, per_sample)) + "\n")
        run([which("stringtie"), "--merge", "-p", str(threads), "-o", str(gtf_out), str(lst)],
            desc="stringtie --merge")
    return gtf_out


def run_rna_branch(vi: ValidatedInputs, out_dir: Path, threads: int,
                   hisat2_opts: list[str] | None = None,
                   stringtie_opts: list[str] | None = None) -> dict:
    """Execute the RNA branch. Returns dict of key outputs."""
    rdir = ensure_dir(out_dir / "02.rnaseq")
    log.info("module 2 [rna]: %d fastq + %d bam + %d gtf inputs",
             len(vi.rna_fastqs), len(vi.rna_bams), len(vi.rna_gtfs))

    bams: list[Path] = []
    if vi.rna_fastqs:
        idx = build_index(vi.genome, rdir / "hisat2_index", threads)
        for s in vi.rna_fastqs:
            bams.append(align_sample(s, idx, rdir / "bam", threads,
                                     extra_opts=hisat2_opts))
    bams += [s.files[0] for s in vi.rna_bams]
    long_read = {s.files[0] for s in vi.rna_fastqs + vi.rna_bams if s.isoseq}

    beds = []
    for bam in bams:
        beds.append(extract_junctions(bam, rdir / "bam" / (Path(bam).stem + ".junctions.bed")))
    merged_bed = merge_junctions(beds, rdir / "junctions.bed") if beds else None

    gtf = None
    if bams:
        gtf = assemble(bams, rdir / "stringtie" / "stringtie.gtf", threads,
                       extra_opts=stringtie_opts, long_read=long_read)
    elif vi.rna_gtfs:
        gtf = vi.rna_gtfs[0]  # precomputed transcripts passed through

    result = {"bams": bams, "junctions_bed": merged_bed, "stringtie_gtf": gtf}
    log.info("module 2 [rna] done: bams=%d junctions=%s stringtie=%s",
             len(bams), merged_bed, gtf)
    return result
