"""Module 2 - RNA branch: HISAT2/minimap2 alignment -> StringTie -> junctions.bed.

Aligner is chosen per sample tag:
  fastq/fasta            HISAT2 (Illumina SE/PE)
  isoseq*/ont            minimap2 -x splice[:hq] (long reads)
  transcripts            minimap2 -x splice:hq -u f (assembled transcripts;
                         excluded from coverage tracks: their depth is 1x)
  mix                    HISAT2 (Illumina part) + minimap2 (long-read part),
                         assembled with stringtie --mix
  bam/bam_isoseq/bam_mix pre-aligned, passed through

Lazy/resumable: existing non-empty outputs are skipped.
Outputs (under <outdir>/02.rnaseq/):
  hisat2_index/            HISAT2 index
  bam/<sample>.sorted.bam  sorted+indexed BAM per sample (mix: .ilmn/.lr)
  stringtie/stringtie.gtf  merged assembly (all samples)
  stringtie/isoseq.stringtie.gtf  Iso-Seq-only merged assembly (when any
                           isoseq/ont/bam_isoseq sample exists; augment
                           tier0 donor pool)
  stringtie/short.stringtie.gtf   short-read-only merged assembly (same
                           condition; augment tier1 donor pool)
  junctions.bed            per-junction read support, merged over samples
"""
from __future__ import annotations

import logging
import shlex
import shutil
from collections import defaultdict
from pathlib import Path

import pysam

from ..utils import ToolError, ToolOpts, ensure_dir, lazy, run, tracker, which
from .validate import ValidatedInputs

log = logging.getLogger("annomerge")

_FASTA_SUF = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def _is_fasta(p: Path) -> bool:
    return p.name.lower().endswith(_FASTA_SUF)


def build_index(genome: Path, index_dir: Path, threads: int) -> Path:
    prefix = index_dir / "genome"
    marker = Path(str(prefix) + ".1.ht2")
    if not lazy([marker], desc="hisat2-build"):
        ensure_dir(index_dir)
        run([which("hisat2-build"), "-p", str(threads), str(genome), str(prefix)],
            desc="hisat2-build")
    return prefix


def align_sample(sample, index_prefix: Path, bam_dir: Path, threads: int,
                 extra_opts: list[str] | None = None,
                 max_intron: int | None = None,
                 min_intron: int | None = None,
                 files: list | None = None,
                 name: str | None = None) -> Path:
    """HISAT2 -> samtools sort -> index. Returns sorted BAM path.

    files/name override: used by mix samples (Illumina part only).  -f is
    derived from the input file extension, not the sample tag (the Illumina
    part of a mix sample may be fasta).
    """
    files = list(files or sample.files)
    name = name or sample.name
    ensure_dir(bam_dir)
    bam = bam_dir / f"{name}.sorted.bam"
    if lazy([bam, Path(str(bam) + ".bai")], desc=f"align {name}"):
        return bam
    sam = bam_dir / f"{name}.sam"
    cmd = [which("hisat2"), "-x", str(index_prefix), "-p", str(threads),
           "--dta", "-S", str(sam)]  # --dta: downstream transcriptome assembly
    if _is_fasta(files[0]):
        cmd.append("-f")
    if min_intron:
        cmd += ["--min-intronlen", str(min_intron)]
    if max_intron:
        cmd += ["--max-intronlen", str(max_intron)]
    cmd += list(extra_opts or [])
    if len(files) == 2:
        cmd += ["-1", str(files[0]), "-2", str(files[1])]
    else:
        cmd += ["-U", str(files[0])]
    run(cmd, log_file=bam_dir / f"{name}.hisat2.log", desc=f"hisat2 {name}")
    run([which("samtools"), "sort", "-m", "256M", "-@", str(threads),
         "-o", str(bam), str(sam)],
        desc=f"samtools sort {name}")
    sam.unlink()
    run([which("samtools"), "index", str(bam)], desc=f"samtools index {name}")
    return bam


def align_minimap2(sample, genome: Path, bam_dir: Path, threads: int,
                   max_intron: int = 500000,
                   extra_opts: list[str] | None = None,
                   files: list | None = None,
                   name: str | None = None,
                   preset: str | None = None) -> Path:
    """minimap2 (splice-aware) -> samtools sort -> index. Returns sorted BAM.

    Used for long-read samples (isoseq_*/ont), assembled transcripts and the
    long-read part of mix samples.  preset comes from the sample tag
    (splice:hq for HiFi/transcripts, splice for CLR/ONT); the long-read part
    of a mix defaults to plain 'splice' and can be overridden with
    --minimap2-options '-x splice:hq'.  -u f (transcript strand is known)
    is only applied to transcripts samples.
    """
    files = list(files or sample.files)
    name = name or sample.name
    ensure_dir(bam_dir)
    bam = bam_dir / f"{name}.sorted.bam"
    if lazy([bam, Path(str(bam) + ".bai")], desc=f"minimap2 {name}"):
        return bam
    preset = preset or sample.mm2_preset or "splice:hq"
    cmd = [which("minimap2"), "-a", "-x", preset,
           "-t", str(threads), "-G", str(max_intron)]
    if sample.tag == "transcripts":
        cmd += ["-u", "f"]       # assembled transcripts: strand is known
    cmd += [str(genome), str(files[0])]
    cmd += list(extra_opts or [])
    pipe = " ".join(shlex.quote(str(c)) for c in cmd) \
        + f" | {shlex.quote(which('samtools'))} sort -m 256M -@ {threads}" \
        + f" -o {shlex.quote(str(bam))} -"
    run(pipe, shell=True, log_file=bam_dir / f"{name}.minimap2.log",
        desc=f"minimap2 {name}")
    run([which("samtools"), "index", str(bam)], desc=f"samtools index {name}")
    return bam


def extract_junctions(bam: Path, out_bed: Path, min_reads: int = 1) -> Path:
    """Extract splice junctions (CIGAR N blocks) with per-junction read counts.

    Strand: minimap2 splice alignments carry the ts tag (transcript strand);
    everything else falls back to the alignment direction.
    """
    if lazy([out_bed], desc=f"junctions {bam.stem}"):
        return out_bed
    counts: dict[tuple[str, int, int, str], int] = defaultdict(int)
    with pysam.AlignmentFile(str(bam), "rb") as af:
        for read in af.fetch(until_eof=True):
            if read.is_unmapped or not read.cigartuples:
                continue
            pos = read.reference_start
            strand = read.get_tag("ts") if read.has_tag("ts") \
                else ("-" if read.is_reverse else "+")
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
    tracker.note("splice junctions: %d (>= %d reads) from %s", sum(1 for v in counts.values() if v >= min_reads),
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
    tracker.note("merged splice junctions: %d sites from %d samples", len(agg), len(beds))
    return out_bed


def _gtf_has_transcripts(gtf: Path) -> bool:
    """True when the GTF contains at least one transcript/mRNA feature."""
    if not gtf.exists():
        return False
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) > 2 and f[2] in ("transcript", "mRNA"):
                return True
    return False


def assemble(units: list[dict], gtf_out: Path, threads: int,
             extra_opts: ToolOpts | list | None = None) -> Path:
    """StringTie per sample, then stringtie --merge.

    units: [{"name": str, "bams": [Path], "mode": str, "opts": [str]?,
             "platform": str?, "lr_preset": str?}]
    with mode
      default   one BAM, plain stringtie
      longread  one BAM, stringtie -L
      mix       two BAMs (Illumina, long-read): stringtie --mix b1 b2
    "platform" is informational only here ("short" / "isoseq"); run_rna_branch
    uses it to build per-platform merged GTFs for the augment donor pools.
    "lr_preset" records the minimap2 splice preset of a long-read unit
    (splice:hq = HiFi, splice = CLR/ONT) for downstream boundary tolerance.
    optional "opts" are per-unit, pipeline-managed StringTie flags (e.g.
    transcripts units get relaxed thresholds: single-copy transcripts are
    1x coverage by nature and stringtie -L assembles nothing from them).
    NOTE override direction: StringTie honours the FIRST occurrence of a
    repeated flag (verified on 3.0.3: '-c 99 -c 1' applies 99) — opposite
    of minimap2/miniprot (last wins). "opts" is inserted BEFORE extra_opts,
    so the pipeline's per-unit thresholds win and a user --stringtie-options
    repeat of the same flag is silently ignored; this is intentional
    (protects the transcripts thresholds).
    extra_opts are appended to every per-sample call (protected: input BAMs,
    -p, -o, -L, --mix — those are pipeline-managed).
    Samples that assemble nothing are dropped before --merge (stringtie
    --merge dies on a transcript-less input GTF).
    """
    if lazy([gtf_out], desc="stringtie"):
        return gtf_out
    sd = gtf_out.parent
    ensure_dir(sd)
    per_sample = []
    for u in units:
        g = sd / f"{u['name']}.gtf"
        if not lazy([g], desc=f"stringtie {u['name']}"):
            if u["mode"] == "mix":
                cmd = [which("stringtie"), str(u["bams"][0]), str(u["bams"][1]),
                       "--mix", "-p", str(threads), "-o", str(g)]
            else:
                cmd = [which("stringtie"), str(u["bams"][0]),
                       "-p", str(threads), "-o", str(g)]
                if u["mode"] == "longread":
                    cmd.append("-L")
                cmd += list(u.get("opts") or [])
            if extra_opts is not None:
                opts = (extra_opts if isinstance(extra_opts, ToolOpts)
                        else ToolOpts(extra_opts, {}))
                cmd += opts.for_target(u["name"])
            run(cmd, log_file=sd / f"{u['name']}.stringtie.log",
                desc=f"stringtie {u['name']}")
        if _gtf_has_transcripts(g):
            per_sample.append(g)
        else:
            log.warning("stringtie %s: no transcripts assembled — excluded "
                        "from merge", u["name"])
    if not per_sample:
        log.warning("stringtie: no transcripts assembled in any sample; "
                    "writing an empty merged GTF")
        gtf_out.write_text("##gff-version 3\n")
        return gtf_out
    if len(per_sample) == 1:
        # copy, not rename: per-sample GTFs are shared cache for the
        # per-platform merges run_rna_branch may trigger afterwards
        shutil.copy(per_sample[0], gtf_out)
    else:
        lst = sd / "gtf_list.txt"
        lst.write_text("\n".join(map(str, per_sample)) + "\n")
        run([which("stringtie"), "--merge", "-p", str(threads), "-o", str(gtf_out), str(lst)],
            desc="stringtie --merge")
    return gtf_out


def _coerce_opts(v) -> ToolOpts:
    """Accept a ToolOpts (CLI path), a plain flag list, or None."""
    if v is None:
        return ToolOpts([], {})
    return v if isinstance(v, ToolOpts) else ToolOpts(list(v), {})


def _check_targets(opts: ToolOpts, valid: set[str], tool: str) -> None:
    """Reject targeted names that no invocation of this tool will see."""
    bad = sorted(set(opts.per_target) - valid)
    if bad:
        raise ToolError(
            f"{tool}: targeted option(s) for unknown sample(s) "
            f"{bad}; samples this tool runs on: {sorted(valid)}")


def run_rna_branch(vi: ValidatedInputs, out_dir: Path, threads: int,
                   hisat2_opts: ToolOpts | list | None = None,
                   stringtie_opts: ToolOpts | list | None = None,
                   minimap2_opts: ToolOpts | list | None = None,
                   max_intron: int = 500000,
                   min_intron: int = 20) -> dict:
    """Execute the RNA branch. Returns dict of key outputs.

    bams:            every aligned BAM (junction extraction + assembly)
    coverage_bams:   BAMs whose depth is real sequencing depth (excludes
                     'transcripts' samples — their depth is a meaningless
                     1x plate, and the long-read part of mix samples)
    """
    rdir = ensure_dir(out_dir / "02.rnaseq")
    tracker.note("rna inputs: %d to align + %d pre-aligned + %d gtf",
             len(vi.rna_fastqs), len(vi.rna_bams), len(vi.rna_gtfs))

    hisat2_opts = _coerce_opts(hisat2_opts)
    stringtie_opts = _coerce_opts(stringtie_opts)
    minimap2_opts = _coerce_opts(minimap2_opts)
    # targeted-option sanity: a NAME=flags entry must point at a sample this
    # tool actually runs on (mix counts for both hisat2 and minimap2; stringtie
    # units exist for every fastq/bam sample)
    fq = {s.name: s for s in vi.rna_fastqs}
    _check_targets(hisat2_opts,
                   {n for n, s in fq.items() if s.tag in ("fastq", "fasta", "mix")},
                   "hisat2")
    _check_targets(minimap2_opts,
                   {n for n, s in fq.items() if s.use_minimap2 or s.tag == "mix"},
                   "minimap2")
    _check_targets(stringtie_opts,
                   set(fq) | {s.name for s in vi.rna_bams},
                   "stringtie")

    bams: list[Path] = []
    coverage_bams: list[Path] = []
    units: list[dict] = []
    idx = None
    bam_dir = rdir / "bam"

    if any(s.tag in ("fastq", "fasta", "mix") for s in vi.rna_fastqs):
        idx = build_index(vi.genome, rdir / "hisat2_index", threads)

    for s in vi.rna_fastqs:
        if s.tag in ("fastq", "fasta"):
            b = align_sample(s, idx, bam_dir, threads,
                             extra_opts=hisat2_opts.for_target(s.name),
                             max_intron=max_intron, min_intron=min_intron)
            bams.append(b)
            coverage_bams.append(b)
            units.append({"name": s.name, "bams": [b], "mode": "default",
                          "platform": "short"})
        elif s.use_minimap2:
            b = align_minimap2(s, vi.genome, bam_dir, threads,
                               max_intron=max_intron,
                               extra_opts=minimap2_opts.for_target(s.name))
            bams.append(b)
            if s.tag == "transcripts":
                # assembled transcripts: 1x coverage per isoform by nature —
                # stringtie -L assembles nothing from them; use the default
                # (short-read) mode with relaxed thresholds instead
                units.append({"name": s.name, "bams": [b], "mode": "default",
                              "opts": ["-c", "1", "-j", "1"],
                              "platform": "short"})
            else:
                units.append({"name": s.name, "bams": [b], "mode": "longread",
                              "platform": "isoseq",
                              "lr_preset": s.mm2_preset or "splice:hq"})
                coverage_bams.append(b)   # real reads: depth is meaningful
        elif s.tag == "mix":
            b1 = align_sample(s, idx, bam_dir, threads,
                              extra_opts=hisat2_opts.for_target(s.name),
                              max_intron=max_intron, min_intron=min_intron,
                              files=s.illumina_files, name=f"{s.name}.ilmn")
            b2 = align_minimap2(s, vi.genome, bam_dir, threads,
                                max_intron=max_intron,
                                extra_opts=minimap2_opts.for_target(s.name),
                                files=s.longread_files, name=f"{s.name}.lr",
                                preset="splice")
            bams += [b1, b2]
            coverage_bams.append(b1)         # long-read part stays out of coverage
            units.append({"name": s.name, "bams": [b1, b2], "mode": "mix",
                          "platform": "short"})

    for s in vi.rna_bams:
        if s.tag == "bam_mix":
            bams += list(s.files)
            coverage_bams.append(s.files[0])
            units.append({"name": s.name, "bams": list(s.files), "mode": "mix",
                          "platform": "short"})
        else:
            bams.append(s.files[0])
            coverage_bams.append(s.files[0])
            is_lr = s.stringtie_longread
            units.append({"name": s.name, "bams": [s.files[0]],
                          "mode": "longread" if is_lr else "default",
                          # bam_isoseq carries no preset info; 'splice'
                          # (CLR/ONT, wider tolerance) is the safe default
                          "platform": "isoseq" if is_lr else "short",
                          "lr_preset": "splice" if is_lr else None})

    beds = []
    for bam in bams:
        beds.append(extract_junctions(bam, bam_dir / (Path(bam).stem + ".junctions.bed")))
    merged_bed = merge_junctions(beds, rdir / "junctions.bed") if beds else None

    gtf = None
    isoseq_gtf = None
    short_gtf = None
    if units:
        st_dir = rdir / "stringtie"
        gtf = assemble(units, st_dir / "stringtie.gtf", threads,
                       extra_opts=stringtie_opts)
        # per-platform merged GTFs for the augment donor pools: Iso-Seq
        # assemblies keep their identity (full-length boundary evidence)
        # instead of being dissolved into the pooled merge.  The combined
        # stringtie.gtf above is untouched — EVM/goldset/report/viewer see
        # exactly what they saw before.  When no Iso-Seq platform sample
        # exists both stay None and augment falls back to stringtie.gtf.
        iso_units = [u for u in units if u.get("platform") == "isoseq"]
        short_units = [u for u in units if u.get("platform") != "isoseq"]
        if iso_units:
            isoseq_gtf = assemble(iso_units, st_dir / "isoseq.stringtie.gtf",
                                  threads, extra_opts=stringtie_opts)
            if short_units:
                short_gtf = assemble(short_units,
                                     st_dir / "short.stringtie.gtf",
                                     threads, extra_opts=stringtie_opts)
            lr_presets = sorted({u.get("lr_preset") or "splice"
                                 for u in iso_units})
            tracker.note("isoseq donor pool: %d sample(s) -> %s "
                         "(presets: %s); short pool: %s",
                         len(iso_units), isoseq_gtf, ",".join(lr_presets),
                         short_gtf or "none")
    elif vi.rna_gtfs:
        gtf = vi.rna_gtfs[0]  # precomputed transcripts passed through

    result = {"bams": bams, "coverage_bams": coverage_bams,
              "junctions_bed": merged_bed, "stringtie_gtf": gtf,
              "isoseq_gtf": isoseq_gtf, "short_gtf": short_gtf}
    tracker.note("rna done: bams=%d (coverage=%d), splice_junctions=%s, stringtie=%s",
             len(bams), len(coverage_bams), merged_bed, gtf)
    return result
