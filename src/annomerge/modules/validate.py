"""Module 1: input validation (fail-fast).

Checks, per CLI design v2.1:
  - genome FASTA exists; seqids extracted
  - -p homolog FASTA is MANDATORY; headers must carry species tags (>id|species=Genus_species)
  - -P predictor list parsing (name=path, comma-separated; unnamed -> filename stem)
  - -r RNA list file parsing (fastq/bam/bam_isoseq/gtf tags), mutually exclusive with --bam/--transcripts
  - --busco / --core mutually exclusive
  - GFF seqids must be a subset of genome seqids (fail-fast)
"""
from __future__ import annotations

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger("annomerge")


class InputError(ValueError):
    """Fatal input validation error (fail-fast)."""


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


@dataclass
class RnaSample:
    files: list[Path]
    tag: str  # fastq | bam | bam_isoseq | gtf

    @property
    def name(self) -> str:
        n = self.files[0].name
        for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam", ".gtf", ".gff"):
            if n.endswith(suf):
                n = n[: -len(suf)]
        for tail in ("_1", "_2", "_R1", "_R2"):
            if n.endswith(tail):
                n = n[: -len(tail)]
        return n


def parse_rna_list(list_file: Path) -> list[RnaSample]:
    samples: list[RnaSample] = []
    with open(list_file) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            tag = parts[-1] if parts[-1] in ("fastq", "bam", "bam_isoseq", "gtf") else "fastq"
            files = parts[:-1] if parts[-1] in ("fastq", "bam", "bam_isoseq", "gtf") else parts
            files = [str(_require(f, f"RNA file (line {ln})")) for f in files]
            if tag == "fastq" and len(files) not in (1, 2):
                raise InputError(f"{list_file}:{ln}: fastq sample needs 1 (SE) or 2 (PE) files")
            if tag in ("bam", "bam_isoseq", "gtf") and len(files) != 1:
                raise InputError(f"{list_file}:{ln}: {tag} sample needs exactly 1 file")
            samples.append(RnaSample(files=[Path(f) for f in files], tag=tag))
    if not samples:
        raise InputError(f"RNA list file is empty: {list_file}")
    return samples


_LABEL_RE = __import__("re").compile(r"^[A-Za-z0-9_.\-]+$")

def parse_homolog_sets(spec: str) -> "dict[str, Path]":
    """Parse '-p id1=xxx.fa,id2=abc.pep' -> {label: path}.

    Same syntax as -P predictors: 'label=path' pairs, comma-separated;
    unlabeled items take their filename stem. Labels must be unique and
    match [A-Za-z0-9_.-] (they become EVM source names / report columns).
    """
    out: "dict[str, Path]" = {}
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            label, p = item.split("=", 1)
            label = label.strip()
        else:
            p = item
            label = Path(p.strip()).stem.split(".")[0]
        if not _LABEL_RE.match(label):
            raise InputError(
                f"invalid homolog set label '{label}': only [A-Za-z0-9_.-] allowed"
            )
        if label in out:
            raise InputError(
                f"duplicate homolog set label '{label}': give each set an explicit "
                "label, e.g. -p scer=scer.pep,sklu=sklu.pep"
            )
        out[label] = _require(p.strip(), f"homolog protein FASTA (set {label})")
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
    core: Path | None = None
    species_map: Path | None = None

    @property
    def has_rna(self) -> bool:
        return bool(self.rna_samples) or self.transcripts is not None

    @property
    def rna_fastqs(self) -> list[RnaSample]:
        return [s for s in self.rna_samples if s.tag == "fastq"]

    @property
    def rna_bams(self) -> list[RnaSample]:
        return [s for s in self.rna_samples if s.tag in ("bam", "bam_isoseq")]

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
    rna_list: str | None = None,
    bam: str | None = None,
    transcripts: str | None = None,
    busco: str | None = None,
    core: str | None = None,
    species_map: str | None = None,
    workdir: str | None = None,
) -> ValidatedInputs:
    """Fail-fast validation shared by merge and evaluate.

    -p accepts 'label=path' items or plain paths, comma-separated
    ('-p scer=scer.pep,sklu=sklu.pep'); plain items take their filename stem.
    Multiple files are header-checked individually and concatenated (cached
    under workdir/input/) so every downstream consumer keeps seeing one
    transparent FASTA; a protein->set membership table is kept alongside.
    """
    log.info("module 1: input validation")

    g = _require(genome, "genome FASTA (-g)")
    gids = set(fasta_seqids(g))
    log.info("genome: %d sequences", len(gids))

    h_sets = parse_homolog_sets(homologs)
    membership = None
    if len(h_sets) == 1:
        label, hp = next(iter(h_sets.items()))
        check_homolog_headers(hp)
    else:
        if not workdir:
            raise InputError(
                "multiple homolog FASTAs require workdir to cache the merged file"
            )
        hp = Path(workdir) / "input" / "homologs.merged.faa"
        _, membership = merge_homolog_fastas(h_sets, hp)

    if busco and core:
        raise InputError("--busco and --core are mutually exclusive")
    busco_p = _require(busco, "BUSCO lineage (--busco)") if busco else None
    core_p = _require(core, "core protein FASTA (--core)") if core else None

    if rna_list and (bam or transcripts):
        raise InputError("-r is mutually exclusive with --bam/--transcripts")

    samples: list[RnaSample] = []
    if rna_list:
        samples = parse_rna_list(_require(rna_list, "RNA list file (-r)"))
    if bam:
        for b in bam.split(","):
            samples.append(RnaSample(files=[_require(b.strip(), "BAM (--bam)")], tag="bam"))
    tr = _require(transcripts, "transcripts GTF (--transcripts)") if transcripts else None

    preds = parse_predictors(predictors)
    for name, p in preds.items():
        check_gff_seqids(p, gids, f"predictor {name}")

    anno_p = None
    if annotation:
        anno_p = _require(annotation, "annotation GFF3 (-a)")
        check_gff_seqids(anno_p, gids, "annotation (-a)")

    sm = _require(species_map, "species map (--species-map)") if species_map else None

    vi = ValidatedInputs(
        genome=g, genome_ids=gids, homologs=hp,
        homolog_sets=h_sets, homolog_membership=membership,
        predictors=preds,
        annotation=anno_p, rna_samples=samples, transcripts=tr,
        busco=busco_p, core=core_p, species_map=sm,
    )
    log.info(
        "validation OK: predictors=%d rna_samples=%d transcripts=%s busco=%s core=%s",
        len(preds), len(samples), bool(tr), bool(busco_p), bool(core_p),
    )
    return vi
