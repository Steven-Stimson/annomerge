"""Module 1: input validation (fail-fast).

Checks, per CLI design v2.1:
  - genome FASTA exists; seqids extracted
  - -p homolog FASTA is MANDATORY; headers must carry species tags (>id|species=Genus_species)
  - -P predictor list parsing (name=path, comma-separated; unnamed -> filename stem)
  - -r RNA list file parsing (fastq/bam/bam_isoseq/gtf tags), mutually exclusive with --bam/--transcripts
  - --busco-lineage local BUSCO odb lineage (offline completeness QC)
  - GFF seqids must be a subset of genome seqids (fail-fast)
"""
from __future__ import annotations

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path
from ..utils import tracker
log = logging.getLogger("annomerge")


class InputError(ValueError):
    """Fatal input validation error (fail-fast)."""


# ------------------------------------------------------- input materialization
# Every user-supplied file/dir is linked under <workdir>/01.validate/input/
# (genome/, denovo/, proteins/, rna/<sample>/, busco/, omark/, annotation/)
# and every downstream module uses the LINKS. User input dirs may live on a
# read-only mount (e.g. /Files); tools that drop sidecar files next to their
# input (gffcompare .tmap, samtools .fai/.bai) then write into our writable
# tree instead of dying with 'Error creating file'. input/manifest.tsv records
# original -> link for provenance.
_SIDECAR_SUF = (".fai", ".gzi", ".csi")


def _cache_stale(cache: Path, srcs) -> bool:
    """Rebuild a concatenation cache when any source is newer than it, or when
    the sizes stop adding up (catches mtime-preserving replacements such as
    rsync -a)."""
    if not cache.exists():
        return True
    cm = cache.stat()
    if any(s.stat().st_mtime > cm.st_mtime for s in srcs):
        return True
    return sum(s.stat().st_size for s in srcs) != cm.st_size


def _fmt_suffix(p: Path) -> str:
    """Suffix used when naming a materialized link: keep compound extensions
    for gzipped inputs (x.gff3.gz -> .gff3.gz, not bare .gz)."""
    return "".join(p.suffixes[-2:]) if p.suffix == ".gz" else p.suffix


def _materialize(src: Path, dst: Path, manifest: list, kind: str) -> Path:
    """Symlink src -> dst (absolute target); copy when symlinks are unsupported.
    Reuse rules: an existing symlink resolving to src is kept; a regular copy
    from an earlier run is kept when its <dst>.src record (target path + mtime
    + size) still matches; anything else is replaced. Index sidecars
    (.fai/.gzi/.csi, and both .bam.bai/.bai conventions) ride along so tools
    needing an index can read (or rebuild) it in our tree."""
    import shutil
    src, dst = Path(src), Path(dst)
    target = src.resolve()
    dst.parent.mkdir(parents=True, exist_ok=True)
    src_rec = dst.with_name(dst.name + ".src")
    reuse = False
    if dst.is_symlink():
        reuse = dst.exists() and dst.resolve() == target
    elif dst.exists():
        # copy fallback from an earlier run: reuse when the recorded source
        # (path + mtime + size) is unchanged
        try:
            t, mt, sz = src_rec.read_text().rstrip("\n").split("\t")
            st = target.stat()
            reuse = (t == str(target)
                     and float(mt) == st.st_mtime and int(sz) == st.st_size)
        except (OSError, ValueError):
            reuse = False
    if not reuse:
        if dst.is_symlink() or dst.is_file():
            dst.unlink()
        elif dst.is_dir():
            shutil.rmtree(dst)
        try:
            dst.symlink_to(target, target_is_directory=target.is_dir())
            src_rec.unlink(missing_ok=True)
        except OSError:
            if target.is_dir():
                shutil.copytree(target, dst)
            else:
                shutil.copy2(target, dst)
            st = target.stat()
            src_rec.write_text(f"{target}\t{st.st_mtime}\t{st.st_size}\n")
    manifest.append((str(target), str(dst), kind))
    if target.is_file():
        for suf in _SIDECAR_SUF:
            side = target.with_name(target.name + suf)
            if side.exists():
                _materialize(side, dst.with_name(dst.name + suf),
                             manifest, "sidecar")
        if target.suffix == ".bam":
            for side, dsuf in ((target.with_name(target.name + ".bai"), ".bai"),
                               (target.with_suffix(".bai"), ".bai")):
                if side.exists():
                    d = (dst.with_name(dst.name + dsuf)
                         if side.name.endswith(".bam.bai")
                         else dst.with_suffix(dsuf))
                    _materialize(side, d, manifest, "sidecar")
    return dst


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def _require(path: str | Path, what: str) -> Path:
    p = Path(path)
    if not p.exists():
        raise InputError(f"{what} not found: {p}")
    return p


def fasta_seqids(path: Path, max_scan: int | None = None) -> list[str]:
    ids = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
                if max_scan and len(ids) >= max_scan:
                    break
    if not ids:
        raise InputError(f"no sequences found in FASTA: {path}")
    return ids


def parse_predictors(spec: str | None) -> dict[str, Path]:
    """Parse '-P aug=aug.gff3,braker=braker.gff3' -> {name: path}."""
    out: dict[str, Path] = {}
    if not spec:
        return out
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            name, p = item.split("=", 1)
        else:
            p = item
            name = Path(p).stem.split(".")[0]
        name = name.strip()
        if name in out:
            raise InputError(f"duplicate predictor name in -P: {name}")
        out[name] = _require(p.strip(), f"predictor GFF3 ({name})")
    return out


_LR_FASTQ_TAGS = ("isoseq", "isoseq_hifi", "isoseq_clr", "ont", "transcripts")
_KNOWN_TAGS = ("fastq", "fasta") + _LR_FASTQ_TAGS + \
              ("bam", "bam_isoseq", "mix", "bam_mix", "gtf")
_MM2_PRESET = {"isoseq": "splice:hq", "isoseq_hifi": "splice:hq",
               "isoseq_clr": "splice", "ont": "splice",
               "transcripts": "splice:hq"}


@dataclass
class RnaSample:
    files: list[Path]
    tag: str               # one of _KNOWN_TAGS
    name_override: str | None = None

    @property
    def name(self) -> str:
        if self.name_override:
            return self.name_override
        n = self.files[0].name
        for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam",
                    ".gtf", ".gff3", ".gff",
                    ".fasta.gz", ".fa.gz", ".fna.gz",
                    ".fasta", ".fa", ".fna"):
            if n.endswith(suf):
                n = n[: -len(suf)]
        for tail in ("_1", "_2", "_R1", "_R2"):
            if n.endswith(tail):
                n = n[: -len(tail)]
        return n

    # ---- derived behaviour ----------------------------------------------
    @property
    def is_longread(self) -> bool:
        return self.tag in _LR_FASTQ_TAGS + ("bam_isoseq",)

    @property
    def needs_alignment(self) -> bool:
        return self.tag in ("fastq", "fasta", "mix") + _LR_FASTQ_TAGS

    @property
    def use_minimap2(self) -> bool:
        return self.tag in _LR_FASTQ_TAGS

    @property
    def mm2_preset(self) -> str | None:
        return _MM2_PRESET.get(self.tag)

    @property
    def stringtie_longread(self) -> bool:        # stringtie -L downstream
        return self.tag in _LR_FASTQ_TAGS + ("bam_isoseq",)

    @property
    def is_mix(self) -> bool:
        return self.tag in ("mix", "bam_mix")

    @property
    def illumina_files(self) -> list[Path]:      # mix: Illumina part
        if self.tag == "mix":
            return self.files[:-1]
        if self.tag == "bam_mix":
            return self.files[:1]
        return self.files

    @property
    def longread_files(self) -> list[Path]:      # mix: long-read part
        if self.tag in ("mix", "bam_mix"):
            return self.files[-1:]
        return self.files if self.is_longread else []


_FQ_SUF = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
_FA_SUF = (".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna")
_GTF_SUF = (".gtf", ".gff3", ".gff")
_LIST_SUF = (".txt", ".list")


def _sniff_rna(path: Path) -> str:
    n = path.name.lower()
    if n.endswith(_FQ_SUF):
        return "fastq"
    if n.endswith(_FA_SUF):
        return "fasta"
    if n.endswith(".bam"):
        return "bam"
    if n.endswith(_GTF_SUF):
        return "gtf"
    raise InputError(
        f"cannot tell RNA input type from extension: {path} "
        "(expected .fastq/.fq[.gz], .fasta/.fa/.fna[.gz], .bam, "
        ".gtf/.gff[3], or a .txt/.list file)")


def _sample_from_files(raw_files: list[str], name: str | None,
                       ln_tag: str, tag: str | None = None) -> RnaSample:
    files = [_require(f.strip(), f"RNA file ({ln_tag})") for f in raw_files if f.strip()]
    if not files:
        raise InputError(f"{ln_tag}: no input files given")
    if len({f.name for f in files}) != len(files):
        raise InputError(
            f"{ln_tag}: duplicate file NAME within one sample "
            f"({sorted(f.name for f in files)}): files are materialized under "
            "input/rna/<sample>/ by name, so same-named files from different "
            "directories would silently overwrite each other — rename one of them")
    if tag is None:
        # inference: a lone fastq/fasta is Illumina SE (long-read types must
        # be requested with an explicit tag)
        kinds = {_sniff_rna(f) for f in files}
        if len(kinds) != 1:
            raise InputError(
                f"{ln_tag}: mixed input types in one sample "
                f"({', '.join(sorted(kinds))}); give each RNA sample its own "
                "-r (or use the 'mix' tag for Illumina + long-read)")
        k = kinds.pop()
        tag = k
        if k in ("fastq", "fasta") and len(files) > 2:
            raise InputError(
                f"{ln_tag}: {len(files)} files and no tag — is this a 'mix' "
                "sample (Illumina + long-read)? Add an explicit tag")
    else:
        tag = tag.lower()
        if tag not in _KNOWN_TAGS:
            raise InputError(f"{ln_tag}: unknown tag {tag!r} "
                             f"(known tags: {', '.join(_KNOWN_TAGS)})")
    if tag in ("fastq", "fasta"):
        if len(files) not in (1, 2):
            raise InputError(f"{ln_tag}: {tag} sample needs 1 (SE) or 2 (PE) files")
        bad = [f for f in files if _sniff_rna(f) != tag]
        if bad:
            raise InputError(f"{ln_tag}: tag is {tag} but file is not: {bad[0]}")
    elif tag in _LR_FASTQ_TAGS:
        if len(files) != 1:
            raise InputError(f"{ln_tag}: {tag} sample needs exactly 1 file")
        if _sniff_rna(files[0]) not in ("fastq", "fasta"):
            raise InputError(f"{ln_tag}: {tag} expects a fastq/fasta file, got {files[0]}")
    elif tag in ("bam", "bam_isoseq"):
        if len(files) != 1 or _sniff_rna(files[0]) != "bam":
            raise InputError(f"{ln_tag}: {tag} sample needs exactly 1 .bam file")
    elif tag == "mix":
        if len(files) not in (2, 3):
            raise InputError(f"{ln_tag}: mix sample needs Illumina (1-2 files) "
                             "+ long-read (1 file)")
        if any(_sniff_rna(f) not in ("fastq", "fasta") for f in files):
            raise InputError(f"{ln_tag}: mix expects fastq/fasta files only")
        if len({_sniff_rna(f) for f in files[:-1]}) != 1:
            raise InputError(f"{ln_tag}: mix Illumina part must be a single type")
    elif tag == "bam_mix":
        if len(files) != 2 or any(_sniff_rna(f) != "bam" for f in files):
            raise InputError(f"{ln_tag}: bam_mix needs exactly 2 .bam files")
    elif tag == "gtf":
        if any(_sniff_rna(f) != "gtf" for f in files):
            raise InputError(f"{ln_tag}: gtf sample expects .gtf/.gff3/.gff files")
    return RnaSample(files=files, tag=tag, name_override=name)


_INLINE_TAG_RE = __import__("re").compile(r":([A-Za-z_]+)$")


def _split_inline(entry: str):
    """'ID=f1,f2:tag' -> (id | None, [files], tag | None).
    The 'ID=' prefix is honoured only when '=' sits inside the first file
    token; ':tag' is honoured only as the LAST token and must be a known tag
    (so paths containing ':' or '=' are never misparsed)."""
    sid, filespec = None, entry
    head, sep, rest = entry.partition("=")
    if sep and head.strip() and "/" not in head and "," not in head:
        sid, filespec = head.strip(), rest
    tag = None
    m = _INLINE_TAG_RE.search(filespec)
    if m:
        tok = m.group(1).lower()
        if tok in _KNOWN_TAGS:
            tag = tok
            filespec = filespec[: m.start()]
        elif not Path(filespec).exists():
            raise InputError(
                f"unknown RNA sample tag :{m.group(1)} in {entry!r} "
                f"(known tags: {', '.join(_KNOWN_TAGS)})")
        # else: a real file whose name contains ':' -- keep the name as-is
    return sid, filespec.split(","), tag


def _parse_list_line(line: str, ln_tag: str):
    """EviAnn-style list line: '[ID=name] file [file ...] [tag]'.
    Files may also be comma-joined (legacy 'ID=1.fq,2.fq' form).
    NOTE: the sample ID must use the 'ID=' prefix (with '='); a bare
    leading token is always treated as a file name, never as an ID."""
    tokens = line.split()
    tag = None
    if tokens and tokens[-1].lower() in _KNOWN_TAGS:
        tag = tokens.pop(-1).lower()
    sid = None
    if tokens and "=" in tokens[0] and "/" not in tokens[0].split("=", 1)[0]:
        head = tokens.pop(0)
        sid, _, first = head.partition("=")
        sid = sid.strip() or None
        if first:
            tokens.insert(0, first)
    files: list[str] = []
    for tok in tokens:
        files.extend(t for t in tok.split(",") if t)
    return sid, files, tag


def parse_rna_entries(entries: list[str] | None) -> list[RnaSample]:
    """Parse -r/--rna entries.
    Each entry is either an inline sample or a list file:
      -r sample_1.fq,sample_2.fq      PE Illumina sample (inferred)
      -r S1=sample_1.fq,sample_2.fq   same, with an explicit sample ID
      -r se.fq                        SE Illumina sample (inferred)
      -r ccs.fq:isoseq_hifi           PacBio HiFi long-read sample
      -r ont.fq:ont                   ONT long-read sample
      -r tx.fa:transcripts            assembled transcripts (fasta)
      -r "M=il_1.fq,il_2.fq,lr.fq:mix"  one sample sequenced both ways
      -r aligned.bam                  pre-aligned BAM (Illumina)
      -r iso.bam:bam_isoseq           pre-aligned long-read BAM
      -r a.gtf,b.gtf                  assembled transcripts (GTF)
      -r samples.txt                  list file, one sample per line
                                      (EviAnn style: [ID=name] files... [tag])
    """
    samples: list[RnaSample] = []
    for entry in entries or []:
        entry = entry.strip()
        if not entry:
            continue
        if entry.lower().endswith(_LIST_SUF) and Path(entry).is_file():
            with open(entry) as fh:
                for ln, line in enumerate(fh, 1):
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    sid, files, tag = _parse_list_line(line, f"{entry}:{ln}")
                    samples.append(
                        _sample_from_files(files, sid,
                                           f"{entry}:{ln}" + (f" ({sid})" if sid else ""),
                                           tag=tag))
        else:
            sid, files, tag = _split_inline(entry)
            samples.append(_sample_from_files(files, sid, f"-r {entry}", tag=tag))
    seen: set[str] = set()
    for s in samples:
        if s.name in seen:
            raise InputError(f"duplicate RNA sample name: {s.name!r} "
                             "(name samples explicitly with -r ID=files "
                             "or a list file)")
        seen.add(s.name)
    return samples


_LABEL_RE = __import__("re").compile(r"^[A-Za-z0-9_.\-]+$")

def parse_homolog_sets(spec: str) -> "dict[str, list[Path]]":
    """Parse '-p id1=a.fa,b.fa,id2=c.pep' -> {label: [paths]}.

    A 'label=path' item starts a labeled set; every following UNLABELED
    item joins the most recent labeled set, so one label may span several
    files (-p ID1=1.1.fa,1.2.fa,ID2=2.1.fa,2.2.fa = two 2-file sets).
    Items seen before any label (or when no label is given at all) each
    form their own set named by filename stem (legacy behaviour).
    Labels must be unique and match [A-Za-z0-9_.-] (they become EVM
    source names / report columns).
    """
    out: "dict[str, list[Path]]" = {}
    cur: str | None = None   # most recent explicit label
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            label, p = item.split("=", 1)
            label = label.strip()
            if not _LABEL_RE.match(label):
                raise InputError(
                    f"invalid homolog set label '{label}': only [A-Za-z0-9_.-] allowed"
                )
            if label in out:
                raise InputError(
                    f"duplicate homolog set label '{label}': to give one set several "
                    f"files, list them after the label instead — e.g. {label}=a.fa,b.fa"
                )
            out[label] = []
            cur = label
        else:
            p = item
            if cur is None:
                # no label seen yet: one stem-named set per file (legacy)
                label = Path(p.strip()).stem.split(".")[0]
                if label in out:
                    raise InputError(
                        f"duplicate homolog set label '{label}': give each set an "
                        "explicit label, e.g. -p scer=scer.pep,sklu=sklu.pep"
                    )
                out[label] = []
            else:
                label = cur   # joins the most recent labeled set
        out[label].append(_require(p.strip(), f"homolog protein FASTA (set {label})"))
    if not out:
        raise InputError("homolog protein FASTA (-p) is empty")
    return out

def check_homolog_headers(path: Path, scan: int = 1000) -> None:
    """Require '>id|species=Genus_species' headers in homolog FASTA."""
    bad = 0
    with _open(path) as fh:
        for i, line in enumerate(fh):
            if i >= scan:
                break
            if line.startswith(">") and "|species=" not in line:
                bad += 1
    if bad:
        raise InputError(
            f"{path}: {bad} header(s) in first {scan} lack '|species=Genus_species' tag. "
            "Homolog headers must carry a species tag, e.g. '>protID|species=Genus_species'."
        )

def merge_homolog_fastas(sets: "dict[str, Path]", dest: Path) -> tuple[int, Path]:
    """Concatenate multiple labeled homolog protein FASTAs into one (cached at dest).

    Every source must pass check_homolog_headers. Also writes <dest>.sets.tsv
    mapping each protein id -> set label, so downstream consumers (miniprot
    per-set runs, weighted scoring) can recover set membership.
    Returns (total protein count, membership tsv path).
    Callers downstream (miniprot, referee) see a single transparent FASTA.
    """
    n = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    map_path = Path(str(dest) + ".sets.tsv")
    with open(dest, "w") as out, open(map_path, "w") as mp:
        mp.write("protein\tset\n")
        for label, p in sets.items():
            check_homolog_headers(p)
            with _open(p) as fh:
                for line in fh:
                    out.write(line)
                    if line.startswith(">"):
                        n += 1
                        mp.write(f"{line[1:].split()[0]}\t{label}\n")
    log.info("merged %d homolog FASTAs -> %s (%d proteins)", len(sets), dest, n)
    return n, map_path


def gff_seqids(path: Path, max_lines: int = 200000) -> set[str]:
    ids: set[str] = set()
    with _open(path) as fh:
        for i, line in enumerate(fh):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            ids.add(line.split("\t", 1)[0])
            if i >= max_lines:
                break
    return ids


def check_gff_seqids(gff: Path, genome_ids: set[str], label: str) -> None:
    ids = gff_seqids(gff)
    missing = ids - genome_ids
    if missing:
        raise InputError(
            f"{label}: {len(missing)} seqid(s) in GFF not found in genome FASTA "
            f"(e.g. {sorted(missing)[:3]}). Provide a seqid mapping or fix the annotation."
        )


@dataclass
class ValidatedInputs:
    genome: Path
    genome_ids: set[str]
    homologs: Path
    homolog_sets: dict = field(default_factory=dict)      # {label: Path}
    homolog_membership: Path | None = None                # protein -> set tsv
    predictors: dict[str, Path] = field(default_factory=dict)
    annotation: Path | None = None          # evaluate mode: -a
    rna_samples: list[RnaSample] = field(default_factory=list)
    transcripts: Path | None = None         # precomputed StringTie GTF
    busco: Path | None = None
    omark: Path | None = None               # OMAmer .h5 (linked into input/)
    input_dir: Path | None = None           # <outdir>/01.validate/input

    @property
    def has_rna(self) -> bool:
        return bool(self.rna_samples) or self.transcripts is not None

    @property
    def rna_fastqs(self) -> list[RnaSample]:
        return [s for s in self.rna_samples if s.needs_alignment]
    @property
    def rna_bams(self) -> list[RnaSample]:
        return [s for s in self.rna_samples
                if s.tag in ("bam", "bam_isoseq", "bam_mix")]

    @property
    def rna_gtfs(self) -> list[Path]:
        out = [s.files[0] for s in self.rna_samples if s.tag == "gtf"]
        if self.transcripts:
            out.append(self.transcripts)
        return out


def validate(
    genome: str,
    homologs: str,
    predictors: str | None = None,
    annotation: str | None = None,
    rna: list[str] | None = None,
    busco: str | None = None,
    omark_db: str | None = None,
    workdir: str | None = None,
) -> ValidatedInputs:
    """Fail-fast validation shared by merge and evaluate.

    -p accepts 'label=path' items or plain paths, comma-separated
    ('-p scer=scer.pep,sklu=sklu.pep'); plain items take their filename stem.
    Multiple files are header-checked individually and concatenated (cached
    under the input dir) so every downstream consumer keeps seeing one
    transparent FASTA; a protein->set membership table is kept alongside.

    Input materialization: with a workdir, every user-supplied file/dir is
    linked into <workdir>/01.validate/input/ (genome/, denovo/, proteins/,
    rna/<sample>/, busco/, omark/, annotation/) and ALL ValidatedInputs paths
    point at those links — user input dirs may be read-only mounts; tools
    that drop sidecars next to their input (gffcompare .tmap, samtools
    .fai/.bai) then write into our writable tree. input/manifest.tsv records
    original -> link.
    """
    tracker.note("input validation")
    inroot = Path(workdir) / "01.validate" / "input" if workdir else None
    if inroot is None:
        log.warning("no workdir: inputs not materialized — read-only input "
                    "dirs may break tools that write sidecars next to inputs")
    manifest: list = []

    g = _require(genome, "genome FASTA (-g)")
    gids = set(fasta_seqids(g))
    tracker.note("genome: %d sequences", len(gids))

    h_lists = parse_homolog_sets(homologs)
    # flatten multi-file labels into one cached FASTA per set (downstream
    # expects {label: Path}); single-file sets pass through unchanged
    h_sets: dict[str, Path] = {}
    for _lab, _paths in h_lists.items():
        if len(_paths) == 1:
            h_sets[_lab] = _paths[0]
        else:
            if not inroot:
                raise InputError(
                    "multi-file homolog sets require workdir to cache the merged file"
                )
            _mp = inroot / "proteins" / f"homolog_set.{_lab}.faa"
            if _cache_stale(_mp, _paths):
                _mp.parent.mkdir(parents=True, exist_ok=True)
                with open(_mp, "w") as _out:
                    for _p in _paths:
                        _out.write(_p.read_text())
            manifest.append(("+".join(str(_p.resolve()) for _p in _paths),
                             str(_mp), "generated(proteins-merge)"))
            log.info("homolog set %s: %d files merged -> %s",
                     _lab, len(_paths), _mp.name)
            h_sets[_lab] = _mp
    # option B: multiple UNLABELED files are one protein pool (mnp filtering
    # at stage 03 handles divergence; no per-set scoring).  Explicit
    # id=path labels keep the dormant relatedness-weighted multi-set mode.
    labeled = any("=" in it for it in homologs.split(","))
    membership = None
    if len(h_sets) == 1:
        label, hp = next(iter(h_sets.items()))
        check_homolog_headers(hp)
    elif not labeled:
        if not inroot:
            raise InputError(
                "multiple homolog FASTAs require workdir to cache the merged file"
            )
        hp = inroot / "proteins" / "homologs.merged.faa"
        merge_homolog_fastas(h_sets, hp)
        h_sets = {"all": hp}
        log.info("homologs: %d unlabeled file(s) merged into one pool (%s)",
                 len(homologs.split(",")), hp.name)
    else:
        if not inroot:
            raise InputError(
                "multiple homolog FASTAs require workdir to cache the merged file"
            )
        hp = inroot / "proteins" / "homologs.merged.faa"
        _, membership = merge_homolog_fastas(h_sets, hp)

    if busco:
        busco_p = _require(busco, "BUSCO lineage (--busco-lineage)")
    else:
        busco_p = None
    omark_p = _require(omark_db, "OMAmer DB (--omark-db)") if omark_db else None

    samples: list[RnaSample] = parse_rna_entries(rna)
    # a gtf sample may carry several transcript files (ID=a.gtf,b.gtf): merge
    # them into one cached GTF per sample so downstream sees a single file
    for s in samples:
        if s.tag == "gtf" and len(s.files) > 1:
            if not inroot:
                raise InputError("multi-GTF RNA samples require workdir")
            merged = inroot / "rna" / s.name / f"{s.name}.gtf"
            merged.parent.mkdir(parents=True, exist_ok=True)
            if _cache_stale(merged, s.files):
                with open(merged, "w") as out:
                    for f in s.files:
                        out.write(f.read_text())
            manifest.append(("+".join(str(f.resolve()) for f in s.files),
                             str(merged), "generated(rna-gtf-merge)"))
            log.info("rna sample %s: %d GTF(s) merged -> %s",
                     s.name, len(s.files), merged.name)
            s.files = [merged]
    tr = None

    preds = parse_predictors(predictors)
    for name, p in preds.items():
        check_gff_seqids(p, gids, f"predictor {name}")

    anno_p = None
    if annotation:
        anno_p = _require(annotation, "annotation GFF3 (-a)")
        check_gff_seqids(anno_p, gids, "annotation (-a)")

    # ---- materialize every user input under 01.validate/input/ ------------
    if inroot:
        g = _materialize(g, inroot / "genome" / g.name, manifest, "genome")
        preds = {n: _materialize(p, inroot / "denovo" / f"{n}{_fmt_suffix(p)}",
                                 manifest, "denovo")
                 for n, p in preds.items()}
        h_sets = {lab: (p if str(p).startswith(str(inroot))
                        else _materialize(p, inroot / "proteins" / f"{lab}{_fmt_suffix(p)}",
                                          manifest, "proteins"))
                  for lab, p in h_sets.items()}
        # single-set hp IS h_sets' only value (already materialized above);
        # merged pools were generated under inroot and pass through
        hp = next(iter(h_sets.values())) if len(h_sets) == 1 else hp
        for s in samples:
            s.files = [(f if str(f).startswith(str(inroot))
                        else _materialize(f, inroot / "rna" / s.name / f.name,
                                          manifest, "rna"))
                       for f in s.files]
        if busco_p:
            busco_p = _materialize(busco_p, inroot / "busco" / busco_p.name,
                                   manifest, "busco")
        if omark_p:
            omark_p = _materialize(omark_p, inroot / "omark" / omark_p.name,
                                   manifest, "omark")
        if anno_p:
            anno_p = _materialize(anno_p, inroot / "annotation" / anno_p.name,
                                  manifest, "annotation")
        with open(inroot / "manifest.tsv", "w") as mf:
            mf.write("original\tinput_link\tkind\n")
            seen = set()
            for orig, link, kind in manifest:
                if (orig, link) in seen:
                    continue
                seen.add((orig, link))
                mf.write(f"{orig}\t{link}\t{kind}\n")
        log.info("inputs materialized under %s (%d entries)",
                 inroot, len(manifest))

    vi = ValidatedInputs(
        genome=g, genome_ids=gids, homologs=hp,
        homolog_sets=h_sets, homolog_membership=membership,
        predictors=preds,
        annotation=anno_p, rna_samples=samples, transcripts=tr,
        busco=busco_p, omark=omark_p, input_dir=inroot,
    )
    log.info(
        "validation OK: predictors=%d rna_samples=%d transcripts=%s busco=%s",
        len(preds), len(samples), bool(tr), bool(busco_p),
    )
    return vi
