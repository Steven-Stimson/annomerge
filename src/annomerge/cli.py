"""annomerge CLI (typer). Two subcommands: merge / evaluate.

merge:   evidence collection (RNA / homolog proteins / Gold set) -> predictor
         scoring -> EVM weight calibration -> EVM integration -> (--augment)
         isoform resolution / UTR backfill -> sanity repair ->
         confidence tiers -> HTML report + single-file offline browser.
evaluate: the same assessment chain on an existing annotation (read-only).
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import typer
import typer.rich_utils  # noqa: F401  (binds the submodule attr used below)

from . import __version__
from .config import dump_params
from .modules.validate import InputError, validate
from .modules.rna import run_rna_branch
from .modules.homolog import run_homolog_branch
from .modules.goldset import run_goldset
from .utils import (ensure_dir, parse_extra_opts, parse_targeted_opts,
                    set_run_log_dir, setup_logging, tracker)

app = typer.Typer(
    name="annomerge",
    help="""AnnoMerge: EVM weight calibration + evidence auditing for gene annotation.

[blue]Subcommands[/blue]
  merge      Integration pipeline — score de novo predictors, calibrate EVM
             weights, run EVM integration, (optionally) isoform resolution /
             UTR backfill, sanity repair, confidence tiers, HTML report and a
             single-file offline evidence browser.
  evaluate   Read-only assessment of an existing annotation with the same
             evidence chain (never filtered, renamed or coordinate-lifted).

[blue]Global options[/blue]
  -c, --check-dependencies   check all external dependencies and versions
  -v, --version              print version and exit
  -h, --help                 show this help (see <subcommand> --help for options)
""",
    context_settings={"help_option_names": ["-h", "--help"]},
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
)

_ADV_CAUTION = ("[yellow]Caution:[/yellow] only touch these if you fully understand "
                "the consequences")

_TARGETED_NOTE = (" Repeatable. A bare value applies globally; 'NAME=flags' "
                  "targets one RNA sample (explicit ID or filename stem) or "
                  "homolog set label; targeted flags are appended after the "
                  "global ones (last wins). Unknown names are rejected.")

_MNP_HELP = (
    "\b\n"
    "miniprot alignment quality filter (comma-separated tokens).\n"
    "[cyan]Default: 'UnStop,UnShift'[/cyan]; pass [blue]0[/blue] to disable filtering.\n"
    "Switches:\n"
    "  [blue]UnStop[/blue]   drop in-frame-stop alignments\n"
    "  [blue]UnShift[/blue]  drop frameshift alignments\n"
    "Numeric thresholds ([blue]key=value[/blue], 0 = off):\n"
    "  [blue]minIden[/blue]  min identity\n"
    "  [blue]minPos[/blue]   min positive-rate\n"
    "  [blue]minSize[/blue]/[blue]maxSize[/blue]  aligned-CDS length range\n"
    "  [blue]maxRank[/blue]  keep Rank<=N mappings\n"
    "  [blue]topRank[/blue]  keep best N mappings per protein\n"
    "  [blue]R2Q[/blue]/[blue]Q2R[/blue]  lower bounds of the reference-vs-query "
    "CDS length ratio from each side")

_TD_HELP_BASE = (
    "extra flags for TransDecoder.LongOrfs/Predict, e.g. \"-S --single_best_only\". "
    "Predict-only flags (e.g. --single_best_only, --no_refine_starts) are routed "
    "to Predict automatically; all other flags go to both stages. "
    "Protected: -t, -m (use --goldset-min-orf), --retain_blastp_hits. ")

_TD_HELP = _TD_HELP_BASE + _ADV_CAUTION


def _version_callback(value: bool):
    if value:
        typer.echo(f"annomerge {__version__}")
        raise typer.Exit()


def _deps_callback(value: bool):
    if value:
        from .deps import check_all, format_table
        typer.echo(format_table(check_all()))
        n_bad = sum(1 for r in check_all()
                    if not r["ok"] and not r.get("optional"))
        raise typer.Exit(3 if n_bad else 0)


@app.callback()
def _main(
    version: Optional[bool] = typer.Option(None, "--version", "-v",
                                           callback=_version_callback,
                                           is_eager=True, help="report version and exit"),
    check_dependencies: bool = typer.Option(False, "--check-dependencies", "-c",
                                            callback=_deps_callback,
                                            is_eager=True, expose_value=False,
                                            help="verify all external dependencies and exit"),
):
    """AnnoMerge: EVM weight calibration + evidence auditing."""


def _prepare(outdir: str, verbose: bool) -> Path:
    setup_logging(verbose)
    out = ensure_dir(Path(outdir))
    set_run_log_dir(out / "logs" / "run")
    return out


def _norm_mnp(spec: str) -> str:
    """'0' (or none/off) disables the mnp filter entirely."""
    return "" if spec.strip().lower() in ("0", "none", "off") else spec


def _yamlable(v):
    """Make a parsed option value yaml-safe (ToolOpts -> plain dict)."""
    from .utils import ToolOpts
    if isinstance(v, ToolOpts):
        return {"global": v.global_flags, "targeted": v.per_target}
    return v


def _yamlable_key(k: str) -> str:
    """Parsed ToolOpts land under '<name>_resolved' in resolved_params.yaml
    so they are not confused with the raw '<name>_options' CLI values."""
    return k[:-5] + "_resolved" if k.endswith("_opts") else k


def _check_omark_options(spec: str | None) -> str | None:
    """Reject pipeline-managed omark flags; everything else is forwarded."""
    if not spec:
        return None
    import shlex
    protected = {"-f", "--file", "-d", "--database", "-o", "--outputFolder",
                 "-i", "--isoform_file"}
    try:
        toks = shlex.split(spec)
    except ValueError as e:
        typer.secho(f"--omark-options: {e}", fg="red", err=True)
        raise typer.Exit(2)
    for t in toks:
        if t in protected:
            typer.secho(f"--omark-options: {t} is pipeline-managed and "
                        f"rejected (isoform collapsing is auto-generated)",
                        fg="red", err=True)
            raise typer.Exit(2)
    return spec


def _run_evidence_pipeline(
    out: Path, genome: str, homologs: str, rna: Optional[List[str]],
    busco: Optional[str], threads: int,
    goldset_min_orf: int,
    predictors: Optional[str] = None, annotation: Optional[str] = None,
    mnp_filter: str = "",
    miniprot_opts: Optional[List[str]] = None,
    hisat2_opts: Optional[List[str]] = None,
    stringtie_opts: Optional[List[str]] = None,
    td_opts: Optional[List[str]] = None,
    minimap2_opts: Optional[List[str]] = None,
    max_intron: int = 500000,
    min_intron: int = 20,
    partial: int = 1,
):
    tracker.big("01", "validate", "validate inputs (genome, predictions, evidence)")
    vi = validate(
        genome=genome, homologs=homologs, predictors=predictors,
        annotation=annotation, rna=rna, busco=busco,
        workdir=str(out),
    )
    tracker.big_done()
    results: dict = {}

    # module 2: RNA branch
    rna_res = None
    if vi.has_rna:
        tracker.big("02", "rna", "RNA-seq alignment and assembly (hisat2/stringtie/junctions)")
        rna_res = run_rna_branch(vi, out, threads,
                                 hisat2_opts=hisat2_opts,
                                 stringtie_opts=stringtie_opts,
                                 minimap2_opts=minimap2_opts,
                                 max_intron=max_intron,
                                 min_intron=min_intron)
        results["rna"] = rna_res
        tracker.big_done()
    else:
        import logging
        logging.getLogger("annomerge").warning(
            "no RNA evidence; RNA branch skipped (auto-degraded)")

    # module 2: homology branch
    tracker.big("03", "homolog", "homolog protein alignment (miniprot)")
    results["homolog"] = run_homolog_branch(
        vi, out, threads, mnp_filter=mnp_filter, miniprot_opts=miniprot_opts,
        max_intron=max_intron)
    tracker.big_done()

    # module 2b: gold set (built whenever RNA evidence exists)
    if vi.has_rna:
        tracker.big("04", "goldset", "build gold set (gffread + TransDecoder)")
        gtf = rna_res["stringtie_gtf"] if rna_res else (vi.rna_gtfs[0] if vi.rna_gtfs else None)
        results["goldset"] = run_goldset(
            vi, out, threads, stringtie_gtf=gtf,
            referee_faa=None,  # referee protein FASTA extraction: TODO, use full homologs for now
            min_orf_aa=goldset_min_orf,
            td_opts=td_opts,
        )
        tracker.big_done()
    return results, vi


_SCORING_FORMULA = """
[bold]Scoring model[/bold] (module 3 — every de novo predictor gets one P (Prediction score)):
  P = w_r·R_s + w_h·H_s + w_b·B_s + w_o·O_s
  (default w_r=0.45, w_h=0.25, w_b=0.10, w_o=0.20)
Each sub-block maps to exactly ONE evidence source (no double counting):
  R_s  RNA support bundle: junction support + exon precision + single-exon bias
       (internals α/β/γ, see --rna-weights) + GoldMatch (full-length protein
       match vs the RNA-assembled gold set; share --rna-gold-weight)
  H_s  homolog support (structure-aware): fraction of predicted proteins that are
       BOTH full-length matches vs referee proteins AND structurally clean
       (no chimeric/split anomaly from the same miniprot hits). Pure homolog
       evidence; the RNA-derived gold set is mixed in ONLY as a rescue when the
       referee set is smaller than --hs-referee-min proteins.
  B_s  conservation: BUSCO complete recall in proteins mode (needs --busco-lineage)
  O_s  orthology consistency: OMArk Consistent fraction (needs --omark-db)

Sub-blocks without evidence are skipped and the remaining weights are
renormalized.
  --weight-mode adaptive additionally scales the B_s / O_s weights by
  their informative coverage — the fraction of proteins the database
  actually gave a verdict on.
A predictor enters EVM only when its P reaches --min-predictor-score
(a P cannot be computed at all with fewer than 2 usable sub-blocks).
P is then mapped to integer EVM weights (--weight-mapping).
"""

_MERGE_HELP = """[bold]merge[/bold]: score de novo predictors, calibrate EVM
weights and integrate all evidence into a final annotation.

[bold]Examples[/bold]  (repeat a flag or use comma lists for multiple inputs)

  # minimal: one protein pool + one PE short-read sample
  annomerge merge -g genome.fa -P denovo1=a.gff3,denovo2=b.gff3 \\
      -p proteins.pep -r rnaseq_r1.fq,rnaseq_r2.fq --busco-lineage bacteria_odb10

  # multiple protein sets (label=path keeps relatedness-weighted multi-set mode)
  annomerge merge -g genome.fa -P denovo1=a.gff3,denovo3=c.gff3 \\
      -p near=near.pep,far=far.fa -r leaf_r1.fq,leaf_r2.fq

  # multiple unlabeled protein files (merged into one mnp-filtered pool)
  annomerge merge -g genome.fa -P denovo1=a.gff3 -p a.pep,b.pep,c.fa \\
      -r rnaseq_r1.fq,rnaseq_r2.fq

  # multiple RNA samples: repeat -r, one input type per sample (no tag =
  #   inferred: 1-2 fastq/fasta = Illumina SE/PE, 1 bam, >=1 gtf)
  annomerge merge -g genome.fa -P denovo1=a.gff3 -p proteins.pep \
      -r pe_r1.fq,pe_r2.fq -r ccs.fq:isoseq_hifi -r aligned.bam -r asm1.gtf,asm2.gtf

  # long-read / mixed samples via explicit tags
  annomerge merge -g genome.fa -P denovo1=a.gff3 -p proteins.pep \
      -r pe_r1.fq,pe_r2.fq -r ccs.fq:isoseq_hifi -r ont.fq:ont \
      -r "M=il_1.fq,il_2.fq,lr.fq:mix" -r tx.fa:transcripts \
      -r iso.bam:bam_isoseq -r ilmn.bam,lr.bam:bam_mix

  # many RNA samples via a list file (one 'ID=files' sample per line)
  annomerge merge -g genome.fa -P denovo1=a.gff3 -p proteins.pep -r samples.txt

  # --augment: isoform resolution + UTR backfill (needs RNA evidence)
  annomerge merge -g genome.fa -P denovo1=a.gff3 -p proteins.pep \\
      -r rnaseq_r1.fq,rnaseq_r2.fq --augment
"""


class _PanelDescCommand(typer.core.TyperCommand):
    """TyperCommand variant that prints a description paragraph above a named
    rich help panel (typer has no native panel-level description)."""

    PANEL_DESCS: dict = {}

    def format_help(self, ctx, formatter):
        descs = self.PANEL_DESCS
        if not descs:
            return super().format_help(ctx, formatter)
        orig = typer.rich_utils._print_options_panel

        def _with_desc(*args, **kwargs):
            name = kwargs.get("name") or (args[0] if args else "")
            console = kwargs.get("console")
            desc = descs.get(name)
            if desc and console is not None:
                console.print(__import__("rich.panel", fromlist=["Panel"]).Panel(
                    typer.rich_utils._make_rich_text(
                        text=desc, markup_mode=kwargs.get("markup_mode") or "rich"),
                    box=__import__("rich.box", fromlist=["ROUNDED"]).ROUNDED,
                    border_style="dim", padding=(0, 1)))
            return orig(*args, **kwargs)

        typer.rich_utils._print_options_panel = _with_desc
        try:
            super().format_help(ctx, formatter)
        finally:
            typer.rich_utils._print_options_panel = orig


class _MergeHelpCommand(_PanelDescCommand):
    PANEL_DESCS = {
        "Scoring": _SCORING_FORMULA.strip(),
        "Advanced": _ADV_CAUTION,
    }


@app.command(help=_MERGE_HELP, cls=_MergeHelpCommand)
def merge(
    # ---------------------------------------------------------- input files
    genome: Optional[str] = typer.Option(
        None, "--genome", "-g",
        help="[bold green]MANDATORY[/bold green]. Reference genome FASTA.",
        rich_help_panel="Input files"),
    predictors: Optional[List[str]] = typer.Option(
        None, "--predictions", "-P",
        help="[bold green]MANDATORY[/bold green]. De novo gene predictions. "
             "Repeat [blue]-P[/blue] per predictor and/or use a 'name=path' "
             "comma list (unnamed -> filename stem), e.g. "
             "[blue]-P[/blue] denovo1=a.gff3,denovo2=b.gff3 or "
             "[blue]-P[/blue] denovo1=a.gff3 [blue]-P[/blue] denovo2=b.gff3.",
        rich_help_panel="Input files"),
    homologs: Optional[List[str]] = typer.Option(
        None, "--proteins", "-p",
        help="[bold green]MANDATORY[/bold green]. Homolog protein FASTA. "
             "Repeat [blue]-p[/blue] per set and/or use 'label=path' comma lists "
             "(e.g. [blue]-p[/blue] near=near.pep,far=far.fa, or "
             "[blue]-p[/blue] near=near.pep [blue]-p[/blue] far=far.fa; "
             "unlabeled -> filename stem). Files after a label join that "
             "label's set: [blue]-p[/blue] near=n1.pep,n2.pep puts both in set "
             "'near'. Headers: >id|species=Genus_species. "
             "Multiple UNLABELED files merge into one protein pool (mnp-filtered, "
             "unweighted); explicit labels keep the relatedness-weighted multi-set mode.",
        rich_help_panel="Input files"),
    rna: Optional[List[str]] = typer.Option(
        None, "--rna", "-r",
        help="\b\n"
             "RNA evidence; repeat [blue]-r[/blue] per sample. Form:\n"
             "  [blue]-r[/blue] [ID=]file[,file...][:tag]\n"
             "Without a tag the type is inferred: 1-2 fastq/fasta = Illumina\n"
             "(SE/PE), 1 BAM = aligned Illumina, >=1 GTF = assembled\n"
             "transcripts. Explicit tags (long-read & mixed samples):\n"
             "  [blue]isoseq[/blue]/[blue]isoseq_hifi[/blue]  PacBio HiFi (minimap2 splice:hq)\n"
             "  [blue]isoseq_clr[/blue]               PacBio CLR (splice)\n"
             "  [blue]ont[/blue]                      ONT (splice)\n"
             "  [blue]transcripts[/blue]            assembled transcript FASTA\n"
             "                           (splice:hq -u f; excluded from coverage)\n"
             "  [blue]bam_isoseq[/blue]             pre-aligned long-read BAM (stringtie -L)\n"
             "  [blue]mix[/blue]                    one sample sequenced both ways:\n"
             "                           Illumina (1-2 files) then long-read (1),\n"
             "                           assembled with stringtie --mix\n"
             "  [blue]bam_mix[/blue]                pre-aligned pair: Illumina BAM FIRST,\n"
             "                           long-read BAM second\n"
             "A .txt/.list value is a sample sheet, one sample per line:\n"
             "  [ID=]file [file ...] [tag]   (ID needs '=', tag is the last token)\n"
             "Omitting RNA entirely degrades the run to EVM + homolog evidence.",
        rich_help_panel="Input files"),
    busco: Optional[str] = typer.Option(
        None, "--busco-lineage",
        help="Local BUSCO odb lineage directory (offline completeness QC; feeds B_s).",
        rich_help_panel="Input files"),
    omark_db: Optional[str] = typer.Option(
        None, "--omark-db",
        help="OMAmer database (.h5) for OMArk orthology QC (runs after BUSCO).",
        rich_help_panel="Input files"),
    threads: int = typer.Option(
        8, "--threads", "-t", help="Threads (also sets the EVM genome-split granularity).",
        rich_help_panel="Input files"),
    outdir: str = typer.Option(
        "./annomerge_out", "--outdir", "-o", help="Output directory.",
        rich_help_panel="Input files"),
    # ------------------------------------------------------------- scoring
    w_rna: float = typer.Option(
        0.45, "--w-rna",
        help="W_r: weight of the R_s (RNA support) sub-block.",
        rich_help_panel="Scoring"),
    w_homolog: float = typer.Option(
        0.25, "--w-homolog", help="W_h: weight of the H_s (structure-aware homolog "
             "support: full-length AND anomaly-free) sub-block.",
        rich_help_panel="Scoring"),
    w_busco: float = typer.Option(
        0.10, "--w-busco",
        help="W_b: weight of the B_s (BUSCO complete recall) sub-block; "
             "skipped without --busco-lineage.",
        rich_help_panel="Scoring"),
    w_omark: float = typer.Option(
        0.20, "--w-omark",
        help="W_o: weight of the O_s (OMArk consistent fraction) sub-block; "
             "skipped without --omark-db.",
        rich_help_panel="Scoring"),
    rna_weights: str = typer.Option(
        "0.5,0.35,0.15", "--rna-weights",
        help="R_s internals 'alpha,beta,gamma' (junction support, exon coverage, TPM).",
        rich_help_panel="Scoring"),
    rna_gold_weight: float = typer.Option(
        0.25, "--rna-gold-weight",
        help="Share of the gold-set protein match (GoldMatch) inside the R_s "
             "RNA bundle (0-1); the remainder comes from the transcript-level "
             "--rna-weights components. The gold set is RNA-derived, so it "
             "lives in R_s (not H_s) — except in referee-rescue mode.",
        rich_help_panel="Scoring"),
    fulllength_cov: float = typer.Option(
        0.9, "--fulllength-cov",
        help="H_s full-length criterion: a miniprot hit counts as full-length only "
             "when bidirectional protein coverage >= this value.",
        rich_help_panel="Scoring"),
    hs_gold_weight: float = typer.Option(
        0.5, "--hs-gold-weight",
        help="Rescue mode only (referee set smaller than --hs-referee-min): "
             "share of Gold-set hits inside H_s (0-1); the remainder comes from "
             "referee-protein hits. In normal mode the gold set is routed to "
             "R_s instead (see --rna-gold-weight).",
        rich_help_panel="Scoring"),
    hs_referee_min: int = typer.Option(
        1000, "--hs-referee-min",
        help="Minimum referee-set size (proteins) for H_s to rely on the "
             "homolog referee alone. Below this the RNA-derived gold set "
             "rescues H_s (mixed via --hs-gold-weight) instead of feeding R_s.",
        rich_help_panel="Scoring"),
    min_predictor_score: float = typer.Option(
        0.35, "--min-predictor-score",
        help="Minimum P (prediction score) for a predictor to enter EVM integration.",
        rich_help_panel="Scoring"),
    weight_mapping: str = typer.Option(
        "auto", "--weight-mapping",
        help="\b\nP -> EVM integer weight mapping:\n"
             "  [cyan]auto[/cyan]    (default) rank when >=3 predictors — any "
             "number of inputs is spread across the 1..5 weight space; "
             "linear otherwise\n"
             "  [cyan]linear[/cyan]  min-max linear scaling\n"
             "  [cyan]rank[/cyan]    by rank order",
        rich_help_panel="Scoring"),
    weight_mode: str = typer.Option(
        "adaptive", "--weight-mode",
        help="\b\nP aggregation mode:\n"
             "  [cyan]adaptive[/cyan]  (default) scale the B_s / O_s sub-block "
             "weights by their informative coverage — the fraction of "
             "proteins the database actually gave a verdict on\n"
             "  [cyan]fixed[/cyan]     classic fixed weights",
        rich_help_panel="Scoring"),
    ablate: Optional[str] = typer.Option(
        None, "--ablate",
        help="\b\nComma list of sub-blocks to disable — ablation experiments only:\n"
             "  rna, homolog, conserve, omark, structure\n"
             "('structure' only turns off the anomaly screen inside H_s)",
        rich_help_panel="Scoring"),
    # ---------------------------------------------------- EVM integration
    w_gold_evm: int = typer.Option(
        10, "--w-gold-evm", help="EVM evidence weight: TransDecoder Gold set.",
        rich_help_panel="EVM integration"),
    w_transcript_evm: int = typer.Option(
        8, "--w-transcript-evm", help="EVM evidence weight: StringTie transcripts.",
        rich_help_panel="EVM integration"),
    w_protein_evm: int = typer.Option(
        6, "--w-protein-evm",
        help="EVM evidence weight: miniprot homolog proteins (base; per-set weights "
             "apply in labeled multi-set mode).",
        rich_help_panel="EVM integration"),
    # ------------------------------------------------------------- augment
    augment: bool = typer.Option(
        False, "--augment",
        help="Run the augment module after EVM integration: isoform "
             "resolution + UTR backfill + near-duplicate UTR-donor removal. "
             "Requires RNA evidence; auto-falls back to EVM-only without RNA.",
        rich_help_panel="Augment"),
    aug_dedup: str = typer.Option(
        "geometric", "--aug-dedup",
        help="\b\nNear-duplicate UTR-donor removal, two philosophies:\n"
             "  [cyan]geometric[/cyan]  (default) merge donors whose spans sit "
             "within ±--aug-jitter-bp on both ends (assembly-jitter clustering)\n"
             "  [cyan]evidence[/cyan]   EviAnn-style voting — a novel splice "
             "site must pass the splice-PWM gate (--aug-pwm-gate), junction-read "
             "corroboration (--aug-junc-reads) and expression filters "
             "(--aug-tpm-ratio, --aug-min-tpm); no bp threshold\n"
             "Independent of both: --aug-collapse-tol-bp folds "
             "identical-intron-chain isoforms after dedup.",
        rich_help_panel="Augment"),
    aug_jitter_bp: int = typer.Option(
        50, "--aug-jitter-bp",
        help="[cyan]geometric[/cyan] mode only: merge UTR donors closer than this on "
             "both ends (0 = off).",
        rich_help_panel="Augment"),
    aug_pwm_gate: float = typer.Option(
        -3.0, "--aug-pwm-gate",
        help="[cyan]evidence[/cyan] mode only: splice-PWM gate threshold (log2 lift, "
             "calibrated) for novel non-anchor splice sites.",
        rich_help_panel="Augment"),
    aug_junc_reads: int = typer.Option(
        3, "--aug-junc-reads",
        help="[cyan]evidence[/cyan] mode only: minimum junction-read support for a "
             "novel intron to corroborate a donor.",
        rich_help_panel="Augment"),
    aug_tpm_ratio: float = typer.Option(
        0.1, "--aug-tpm-ratio",
        help="[cyan]evidence[/cyan] mode only: keep a donor when TPM >= this ratio x "
             "representative TPM.",
        rich_help_panel="Augment"),
    aug_min_tpm: float = typer.Option(
        0.1, "--aug-min-tpm",
        help="[cyan]evidence[/cyan] mode only: absolute TPM floor for keeping a donor.",
        rich_help_panel="Augment"),
    aug_collapse_tol: int = typer.Option(
        5, "--aug-collapse-tol-bp",
        help="Both modes: isoform collapse tolerance — CDS outer boundaries within "
             "this many bp with an identical intron chain on the same strand count as "
             "one structure (representative = longest total UTR; 0 = exact).",
        rich_help_panel="Augment"),
    # ----------------------------------------------- assessment & report
    viewer: str = typer.Option(
        "lite", "--viewer",
        help="\b\nEvidence browser build:\n"
             "  [cyan]full[/cyan]  embeds the genome sequence (base-level "
             "track, much larger)\n"
             "  [cyan]lite[/cyan]  (default) coordinates only\n"
             "  [cyan]no[/cyan]    skip",
        rich_help_panel="Assessment & report"),
    show_rna_coverage: int = typer.Option(
        1, "--show-rna-coverage",
        help="\b\nRNA coverage tracks embedded in the evidence browser:\n"
             "  [cyan]0[/cyan]  none (smallest browser)\n"
             "  [cyan]1[/cyan]  (default) merged track over all samples only "
             "(a single track)\n"
             "  [cyan]2[/cyan]  merged + one track per sample (largest; "
             "per-sample depth resolution)",
        rich_help_panel="Assessment & report"),
    weight_only: bool = typer.Option(
        False, "--weight-only", help="Stop after the EVM weight table.",
        rich_help_panel="Assessment & report"),
    # ------------------------------------------------------------- advanced
    max_intron: Optional[int] = typer.Option(
        None, "--max-intron",
        help="Max intron length (bp) for hisat2/minimap2/miniprot. Default: "
             "derived from --ploidy as 500000*2/ploidy (diploid = 500000, "
             "after the EviAnn polyploid scaling).",
        rich_help_panel="Advanced"),
    min_intron: int = typer.Option(
        20, "--min-intron",
        help="Min intron length (bp): hisat2 --min-intronlen and EVM "
             "--min_intron_length.",
        rich_help_panel="Advanced"),
    ploidy: int = typer.Option(
        2, "--ploidy",
        help="Genome ploidy; only used to derive the default --max-intron "
             "(500000*2/ploidy). Ignored when --max-intron is given.",
        rich_help_panel="Advanced"),
    partial: int = typer.Option(
        1, "--partial",
        help="\b\nIncomplete-CDS transcripts in the final annotation:\n"
             "  [cyan]0[/cyan]  discard\n"
             "  [cyan]1[/cyan]  (default) keep and tag (partial=true)\n"
             "Keep is the default because --augment isoforms otherwise lose\n"
             "their incomplete-ORF variants (AS structure collapses).",
        rich_help_panel="Advanced"),
    minimap2_options: Optional[List[str]] = typer.Option(
        None, "--minimap2-options",
        help="Extra minimap2 flags for long-read/transcripts/mix alignment, "
             "e.g. \"-x splice:hq\" (overrides the tag-derived preset) or "
             "\"--splice-flank=no\". Protected: -t, -a, -G (use --max-intron), "
             "-u." + _TARGETED_NOTE + " " + _ADV_CAUTION,
        rich_help_panel="Advanced"),
    omark_options: Optional[str] = typer.Option(
        None, "--omark-options",
        help="Extra flags forwarded to the omark assessment step, e.g. "
             "\"-e --taxid 33090 --taxonomic_rank kingdom\" or "
             "\"-e /path/to/taxa.sqlite\". Protected (pipeline-managed): "
             "-f, -d, -o, -i/--isoform_file (isoform collapsing is generated "
             "automatically from each predictor GFF).",
        rich_help_panel="Advanced"),
    hisat2_options: Optional[List[str]] = typer.Option(
        None, "--hisat2-options",
        help="Extra HISAT2 alignment flags, e.g. \"-k 3 --sp "
             "max=15\". Inserted into every per-sample alignment. Protected: "
             "-x, -p, -S, -1, -2, -U, --dta." + _TARGETED_NOTE,
        rich_help_panel="Advanced"),
    stringtie_options: Optional[List[str]] = typer.Option(
        None, "--stringtie-options",
        help="Extra StringTie flags for per-sample assembly, e.g. \"-f 0.1 -c 5\". "
             "Protected: the input BAM, -p, -o, -L (isoseq samples get -L "
             "automatically). The --merge step is pipeline-managed."
             + _TARGETED_NOTE,
        rich_help_panel="Advanced"),
    miniprot_options: Optional[List[str]] = typer.Option(
        None, "--miniprot-options",
        help="Extra miniprot flags, e.g. \"-I -S\". Inserted into both PAF and "
             "GFF runs. Protected: -t, --paf, --gff." + _TARGETED_NOTE,
        rich_help_panel="Advanced"),
    mnp_filter: str = typer.Option(
        "UnStop,UnShift", "--mnp-filter", help=_MNP_HELP,
        rich_help_panel="Advanced"),
    transdecoder_options: Optional[str] = typer.Option(
        None, "--transdecoder-options", help=_TD_HELP_BASE,
        rich_help_panel="Advanced"),
    goldset_min_orf: int = typer.Option(
        50, "--goldset-min-orf",
        help="Gold set: minimum ORF length (aa) for TransDecoder.LongOrfs (-m).",
        rich_help_panel="Advanced"),
    config: Optional[str] = typer.Option(
        None, "--config", help="resolved_params.yaml from a previous run (informational).",
        rich_help_panel="Advanced"),
    verbose: bool = typer.Option(False, "--verbose", help="Verbose logging.",
                                 rich_help_panel="Advanced"),
):

    # repeatable -P/-p values (and comma lists inside them) collapse into the
    # single comma-separated spec the downstream parser has always taken
    predictors = ",".join(predictors) if predictors else None
    homologs = ",".join(homologs) if homologs else None
    missing = [n for n, v in (("--genome/-g", genome), ("--predictions/-P", predictors),
                              ("--proteins/-p", homologs)) if not v]
    if missing:
        typer.secho(f"missing mandatory option(s): {', '.join(missing)} "
                    f"(or use -c/--check-dependencies for a dependency check)",
                    fg="red", err=True)
        raise typer.Exit(2)
    mode = "isoform" if augment else "evm"
    vw = viewer.lower()
    if vw not in ("full", "lite", "no"):
        typer.secho(f"unknown --viewer {viewer!r} (full|lite|no)", fg="red", err=True)
        raise typer.Exit(2)
    wm = weight_mode.lower()
    if wm not in ("fixed", "adaptive"):
        typer.secho(f"unknown --weight-mode {weight_mode!r} (fixed|adaptive)", fg="red", err=True)
        raise typer.Exit(2)
    omark_options = _check_omark_options(omark_options)
    # intron bounds: --max-intron defaults to 500000*2/ploidy (diploid 500k)
    eff_max_intron = max_intron if max_intron else int(500000 * 2 / ploidy)
    try:
        hisat2_opts = parse_targeted_opts(
            hisat2_options, ("-x", "-p", "-S", "-1", "-2", "-U", "--dta"), "hisat2")
        stringtie_opts = parse_targeted_opts(
            stringtie_options, ("-p", "-o", "-L"), "stringtie")
        minimap2_opts = parse_targeted_opts(
            minimap2_options, ("-t", "-a", "-G", "-u"), "minimap2")
        miniprot_opts = parse_targeted_opts(
            miniprot_options, ("-t", "--paf", "--gff"), "miniprot")
        td_opts = parse_extra_opts(
            transdecoder_options, ("-t", "-m", "--retain_blastp_hits"), "TransDecoder")
    except RuntimeError as e:
        typer.secho(str(e), fg="red", err=True)
        raise typer.Exit(2)

    out = _prepare(outdir, verbose)
    from .utils import ToolOpts as _ToolOpts
    params = {(_yamlable_key(k) if isinstance(v, _ToolOpts) else k): _yamlable(v)
              for k, v in locals().items() if k not in ("out", "_ToolOpts")}
    dump_params(params, ensure_dir(out / "01.validate"))
    try:
        results, vi = _run_evidence_pipeline(
            out, genome, homologs, rna, busco, threads,
            goldset_min_orf,
            predictors=predictors, mnp_filter=_norm_mnp(mnp_filter),
            miniprot_opts=miniprot_opts, hisat2_opts=hisat2_opts,
            stringtie_opts=stringtie_opts, td_opts=td_opts,
            minimap2_opts=minimap2_opts, max_intron=eff_max_intron,
            min_intron=min_intron, partial=partial,
        )
    except InputError as e:
        typer.secho(f"input error: {e}", fg="red", err=True)
        raise typer.Exit(2)

    # module 3: predictor scoring + weight mapping
    tracker.big("05", "score", "predictor scoring + EVM weight calibration")
    from .modules.scoring import score_one_predictor, map_weights
    from .modules.validate import parse_predictors

    a, b_, g_ = map(float, rna_weights.split(","))
    cfg = dict(rna_alpha=a, rna_beta=b_, rna_gamma=g_,
               w_rna=w_rna, w_homolog=w_homolog,
               w_busco=w_busco, w_omark=w_omark,
               omark_db=omark_db,
               omark_options=omark_options,
               fulllength_cov=fulllength_cov,
               hs_gold_weight=hs_gold_weight,
               rna_gold_weight=rna_gold_weight,
               hs_referee_min=hs_referee_min,
               weight_mode=weight_mode,
               ablate=(ablate.split(",") if ablate else []))
    preds = parse_predictors(predictors)
    score_dir = ensure_dir(out / "05.score")
    details, p_scores = {}, {}
    # prefetch OMArk for every predictor in parallel (each omamer gets an even
    # share of --threads); the per-predictor compute_os in score_one_predictor
    # then hits the warm *.sum cache instead of running serially. The final
    # annotation's OMArk (report) still runs with the full --threads.
    if cfg.get("omark_db") and "omark" not in set(cfg.get("ablate", [])):
        from concurrent.futures import ThreadPoolExecutor
        from .modules.scoring import (extract_proteins, build_isoform_file,
                                      compute_os)
        import logging as _logging
        _log = _logging.getLogger("annomerge")
        _jobs = list(preds.items())
        _share = max(1, threads // max(1, len(_jobs)))

        def _prefetch_omark(item):
            _name, _gff = item
            _pdir = ensure_dir(score_dir / _name)
            _pf = extract_proteins(vi.genome, _gff, _pdir / "proteins.faa")
            _iso = build_isoform_file(_gff, _pdir / "isoforms.txt")
            try:
                compute_os(_pf, Path(cfg["omark_db"]), _pdir / "omark", _share,
                           isoform_file=_iso,
                           extra_opts=cfg.get("omark_options"))
            except Exception as _e:   # warm-cache best effort; scorer re-checks
                _log.warning("omark prefetch [%s]: %s", _name, _e)

        with ThreadPoolExecutor(max_workers=len(_jobs)) as _pool:
            _futs = [_pool.submit(_prefetch_omark, _item) for _item in _jobs]
            from concurrent.futures import as_completed as _as_completed
            for _i, _f in enumerate(_as_completed(_futs), 1):
                _f.result()
                tracker.parallel("OMArk", _i, len(_futs))
    for name, gff in preds.items():
        d = score_one_predictor(name, gff, vi, score_dir, threads, results, cfg)
        details[name] = d
        p_scores[name] = None if d.get("hard_reject") else d["P"]

    evm_w = map_weights(p_scores, mapping=weight_mapping, min_score=min_predictor_score)
    wtab = score_dir / "evm_weights.tsv"
    with open(wtab, "w") as fh:
        fh.write("predictor\tP\tR_s\tH_s\tB_s\tO_s\teffective_weights\tevm_weight\tnote\n")
        for name in preds:
            d = details[name]
            s = d["subscores"]
            def f(x): return f"{x:.4f}" if isinstance(x, float) else "NA"
            eff = ";".join(f"{k}={v:.3f}" for k, v in
                           (d.get("effective_weights") or {}).items() if v > 0)
            fh.write(f"{name}\t{f(d['P'])}\t{f(s['R_s'])}\t{f(s['H_s'])}\t"
                     f"{f(s['B_s'])}\t{f(s.get('O_s'))}\t"
                     f"{eff or 'NA'}\t"
                     f"{evm_w.get(name) or 'EXCLUDE'}\t"
                     f"{d.get('hard_reject', '')}\n")
    typer.echo(f"\n=== EVM weights ===\n{wtab.read_text()}")
    tracker.big_done()

    # module 4/5: integration (skipped with --weight-only)
    if weight_only:
        typer.secho("weight-only mode: integration skipped", fg="yellow")
        raise typer.Exit(0)
    evm_preds = {name: (preds[name], evm_w[name]) for name in preds
                 if evm_w.get(name) is not None}
    if not evm_preds:
        typer.secho("no predictor passed scoring; integration skipped", fg="red", err=True)
        raise typer.Exit(3)
    # isoform resolution needs RNA evidence (StringTie + gold set); with
    # homolog proteins only the pipeline degrades to EVM-only integration
    # and the viewer shows the EVM result (final = evm.sanity.gff3).
    if not vi.has_rna and augment:
        mode = "evm"
        typer.secho(
            "no RNA evidence: --augment unavailable "
            "(isoform resolution / UTR backfill need RNA); "
            "falling back to EVM-only integration", fg="yellow", err=True)
    gold_dir = out / "04.goldset"
    goldset_dir = gold_dir if (gold_dir / "transcripts.fa.transdecoder.gff3").exists() else None
    stringtie_gtf = (results.get("rna", {}) or {}).get("stringtie_gtf")
    hom_res = results.get("homolog") or {}
    rna_res = results.get("rna") or {}
    tracker.big("06", "merge", "EVM integration + sanity filter + ID normalize")
    from .modules.evm import run_evm
    # homolog-set relatedness weights: per-set EVM protein evidence weights and
    # per-protein weights for confidence assessment (multi-set homolog input)
    set_w = hom_res.get("set_weights") or {}
    membership = hom_res.get("membership") or {}
    multi = len(set_w) > 1
    protein_evm_weights = ({lab: max(1, min(12, round(w_protein_evm * w)))
                            for lab, w in set_w.items()} if multi else None)
    protein_weights = ({p: set_w.get(lab, 1.0) for p, lab in membership.items()}
                       if multi else None)
    evm_res = run_evm(
        vi, out / "06.merge", evm_preds,
        stringtie_gtf=stringtie_gtf,
        miniprot_gff=hom_res.get("gff"),
        goldset_dir=goldset_dir,
        threads=threads, min_intron=min_intron,
        w_gold=w_gold_evm, w_transcript=w_transcript_evm, w_protein=w_protein_evm,
        protein_set_weights=protein_evm_weights,
        protein_membership=membership if multi else None)
    integ_gff3 = evm_res["evm_gff3"]
    integ_d = out / "06.merge"
    typer.secho(f"EVM output: {integ_gff3}", fg="green")
    # module 5b: sanity filter (hard errors removed; notes kept)
    from .modules.sanity import sanity_filter
    sanity_stats = sanity_filter(
        integ_gff3, vi.genome,
        integ_d / "evm.sanity.gff3", integ_d / "sanity_filter.tsv",
        threads=threads, partial=partial)
    typer.secho(f"sanity filter: {sanity_stats or 'skipped (exists)'}", fg="green")
    # module 5c-1: ID normalisation (<seqid>.g<n>.t<m>; Name= keeps old ID)
    from .modules.renamer import normalize_ids
    typer.secho(f"ID normalize (evm): {normalize_ids(integ_d / 'evm.sanity.gff3')}", fg="green")
    tracker.big_done()

    # module 5c: isoform resolution + UTR backfill (isoform mode only)
    final_gff3 = integ_d / "evm.sanity.gff3"
    final_d = integ_d
    if mode == "isoform":
        tracker.big("07", "augment", "isoform resolution + UTR backfill")
        from .modules.augment import run_augment
        aug = run_augment(
            vi, out / "07.augment",
            evm_gff3=integ_d / "evm.sanity.gff3",
            goldset_gff=integ_d / "goldset.genome.gff3",
            stringtie_gtf=stringtie_gtf,
            preds={name: preds[name] for name in evm_preds},
            scores={n: (p_scores.get(n) or 0.0,
                        (details.get(n, {}).get("subscores") or {}).get("R_s") or 0.0)
                    for n in evm_preds},
            jitter_merge_bp=aug_jitter_bp,
            dedup_mode=aug_dedup,
            junctions_bed=rna_res.get("junctions_bed"),
            tpm_ratio=aug_tpm_ratio,
            min_tpm=aug_min_tpm,
            junc_min_reads=aug_junc_reads,
            pwm_gate=aug_pwm_gate,
            collapse_tol_bp=aug_collapse_tol)
        typer.secho(f"augment output: {aug['augment_gff3']}", fg="green")
        aug_sanity = sanity_filter(
            aug["augment_gff3"], vi.genome,
            aug["dir"] / "augment.sanity.gff3", aug["dir"] / "sanity_filter.tsv",
            threads=threads, partial=partial)
        typer.secho(f"augment sanity filter: {aug_sanity or 'skipped (exists)'}", fg="green")
        # module 5c-2: ID normalisation of the final annotation
        typer.secho(f"ID normalize (final): {normalize_ids(aug['dir'] / 'augment.sanity.gff3')}", fg="green")
        tracker.big_done()
        final_gff3 = aug["dir"] / "augment.sanity.gff3"
        final_d = aug["dir"]

    # module 6: coverage-aware evidence tiering (tier, not delete)
    tracker.big("08", "assessment", "confidence tiering + support tagging")
    from .modules.confidence import confidence_annotation
    assess_d = ensure_dir(out / "08.assessment")
    if mode == "evm":
        ensure_dir(out / "07.augment")   # evm mode: create the empty dir only
    conf_stats = confidence_annotation(
        final_gff3, vi.genome,
        rna_res.get("junctions_bed"),
        integ_d / "protein.gff3" if (integ_d / "protein.gff3").exists()
        else hom_res.get("gff"),
        rna_res.get("coverage_bams") or rna_res.get("bams", []),
        final_d / "goldset.genome.gff3" if (final_d / "goldset.genome.gff3").exists()
        else integ_d / "goldset.genome.gff3",
        final_d / "sanity_filter.tsv",
        assess_d / "confidence.tsv",
        protein_weights=protein_weights, threads=threads)
    typer.secho(f"confidence tiers: {conf_stats or 'skipped (exists)'}", fg="green")
    # module 6b: tag the final GFF3 with the evidence verdict (Tier /
    # RNA_support / Protein_support / denovo_support / Flags on every gene
    # and mRNA, viewer wording) + per-transcript CDS numbering (.cds1..N).
    # Always runs (cheap, idempotent) so downstream deliverables stay tagged
    # even when confidence.tsv itself is lazy-skipped.
    from .modules.confidence import annotate_support
    tag_stats = annotate_support(final_gff3, assess_d / "confidence.tsv",
                                 preds or {})
    typer.secho(f"support tags: {tag_stats}", fg="green")
    tracker.big_done()
    # module 7/8: report (run summary + stats + completeness + accuracy + charts)
    import shlex, sys
    tracker.big("09", "result", "report + evidence browser")
    from .modules.report import run_report
    rep = run_report(
        out, final_gff3, vi.genome, final_d,
        busco_lineage=vi.busco, threads=threads,
        omark_db=Path(omark_db) if omark_db else None,
        omark_options=omark_options,
        rna_junctions=rna_res.get("junctions_bed"),
        stringtie_gtf=stringtie_gtf,
        homologs_faa=vi.homologs,
        referee_ids=hom_res.get("referee_ids"),
        run_meta={"mode": mode, "predictors": ",".join(preds),
                  "genome": Path(genome).name, "outdir": str(out),
                  "command": shlex.join(sys.argv),
                  "timestamp": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
                  "threads": threads, "homologs": homologs,
                  "rna": ",".join(rna or []) or "none",
                  "busco": busco or "none",
                  "augment": augment, "viewer": vw,
                  "mnp_filter": _norm_mnp(mnp_filter) or "off",
                  "w_rna": w_rna, "w_homolog": w_homolog,
                  "w_busco": w_busco, "w_omark": w_omark,
                  "min_predictor_score": min_predictor_score,
                  "goldset_min_orf": goldset_min_orf,
                  "weight_mapping": weight_mapping})
    typer.secho(f"report: {rep['html']}", fg="green")
    # evidence browser (igv.js, single-file offline)
    if vw != "no":
        from .modules.viewer import run_viewer
        vw_res = run_viewer(out, vi, results, preds, details,
                            final_gff3, final_d, integ_d, threads,
                            full=(vw == "full"),
                            show_rna_coverage=show_rna_coverage)
        typer.secho(f"evidence browser: {vw_res['html']}", fg="green")
    tracker.big_done()
    # publish layer removed: the module dirs themselves are the numbered
    # deliverables (01.validate ... 09.result)
    typer.secho(f"deliverables: {out}/01.validate ... 09.result", fg="green")


def _annotation_to_gff3(anno: Path, dst: Path) -> Path:
    """Normalize an evaluated annotation to GFF3.

    GFF3 input passes through unchanged: its gene/mRNA structure must be
    preserved (gffread -o silently DROPS gene features, which would empty
    the sanity/confidence model). GTF input (gene_id/transcript_id attrs,
    no ID=/Parent=) is converted by a small pure-Python writer: one gene
    per gene_id (span = union of its features), one mRNA per transcript_id,
    exon/CDS rows re-pointed with Parent=<transcript>.
    """
    if dst.exists():
        return dst
    import logging
    log = logging.getLogger("annomerge")
    is_gtf = False
    with open(anno) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = f[8]
            is_gtf = ("gene_id" in a and "ID=" not in a
                      and "Parent=" not in a)
            break
    if not is_gtf:
        dst.write_bytes(anno.read_bytes())
        log.info("evaluate: annotation already GFF3 -> copied (%s)", dst.name)
        return dst
    # GTF: prefer gffread --keep-genes -F (battle-tested converter; keeps
    # gene records and all attributes), then normalise transcript -> mRNA
    # (the assessment chain parses the mRNA feature type). Pure-Python
    # writer below is the fallback when gffread is unavailable.
    import shutil, subprocess
    gffread = shutil.which("gffread")
    if gffread:
        tmp = dst.with_name(dst.name + ".gffread")
        r = subprocess.run([gffread, "--keep-genes", "-F", "-o", str(tmp), str(anno)],
                           capture_output=True, text=True)
        if r.returncode == 0 and tmp.exists() and tmp.stat().st_size > 0:
            n_g = 0
            with open(tmp) as fh, open(dst, "w") as out:
                out.write("##gff-version 3\n")
                for ln in fh:
                    if ln.startswith("#"):
                        continue
                    f = ln.rstrip("\n").split("\t")
                    if len(f) < 9:
                        continue
                    if f[2] == "transcript":
                        f[2] = "mRNA"
                    if f[2] == "gene":
                        n_g += 1
                    out.write("\t".join(f) + "\n")
            tmp.unlink()
            log.info("evaluate: GTF -> GFF3 via gffread --keep-genes: %d genes", n_g)
            return dst
        log.warning("evaluate: gffread failed (%s); falling back to in-Python GTF conversion",
                    (r.stderr or r.stdout or "").strip()[:200])

    def _attrs(s: str) -> dict:
        av = {}
        for kv in s.rstrip(";").split(";"):
            kv = kv.strip()
            if not kv:
                continue
            if "=" in kv:
                k, v = kv.split("=", 1)
            elif " " in kv:
                k, v = kv.split(None, 1)
            else:
                k, v = kv, ""
            av[k.strip()] = v.strip().strip('"')
        return av

    genes: dict = {}
    with open(anno) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqid, src, ftype = f[0], f[1], f[2]
            start, end = int(f[3]), int(f[4])
            strand, phase, av = f[6], f[7], _attrs(f[8])
            gid = av.get("gene_id")
            if not gid:
                continue
            g = genes.setdefault(gid, {"seqid": seqid, "src": src,
                                       "strand": strand, "start": start,
                                       "end": end, "txs": {}})
            g["start"] = min(g["start"], start)
            g["end"] = max(g["end"], end)
            tid = av.get("transcript_id")
            if not tid:
                continue
            t = g["txs"].setdefault(tid, {"start": start, "end": end,
                                          "strand": strand, "feats": []})
            t["start"] = min(t["start"], start)
            t["end"] = max(t["end"], end)
            t["feats"].append((ftype, start, end, strand, phase))
    n_g = n_t = n_f = 0
    with open(dst, "w") as out:
        out.write("##gff-version 3\n")
        for gid, g in genes.items():
            n_g += 1
            out.write(f'{g["seqid"]}\t{g["src"]}\tgene\t{g["start"]}\t'
                      f'{g["end"]}\t.\t{g["strand"]}\t.\tID={gid}\n')
            for tid, t in g["txs"].items():
                n_t += 1
                out.write(f'{g["seqid"]}\t{g["src"]}\tmRNA\t{t["start"]}\t'
                          f'{t["end"]}\t.\t{t["strand"]}\t.\t'
                          f'ID={tid};Parent={gid}\n')
                for ftype, s, e, strand, phase in t["feats"]:
                    ftype = "exon" if ftype.lower() in ("exon", "cds") \
                        else ftype
                    n_f += 1
                    out.write(f'{g["seqid"]}\t{g["src"]}\t{ftype}\t{s}\t{e}'
                              f'\t.\t{strand}\t{phase}\tParent={tid}\n')
    log.info("evaluate: GTF -> GFF3: %d genes, %d mRNAs, %d sub-features",
             n_g, n_t, n_f)
    return dst


_EVALUATE_HELP = """[bold]evaluate[/bold]: run the merge assessment chain (evidence
collection -> confidence tiers -> per-gene support tags -> HTML report +
evidence browser) on an existing annotation.

The submitted annotation is [bold]read-only[/bold]: it is never filtered,
renamed or coordinate-transformed. The only transformation is a format
normalisation to GFF3 (GFF3 copied verbatim; GTF converted with gffread
--keep-genes -F); assessment tags (Tier / RNA_support / Protein_support /
denovo_support / Flags) are appended to each gene/mRNA. The Gold-set lift
to genome coordinates operates on the StringTie/TransDecoder EVIDENCE
transcripts (never on the submitted annotation) and only feeds the
goldset evidence track. Only -g/-a/-p are mandatory; all evidence inputs
are optional.
"""


@app.command(help=_EVALUATE_HELP)
def evaluate(
    # ---------------------------------------------------------- input files
    genome: Optional[str] = typer.Option(
        None, "--genome", "-g",
        help="[bold green]MANDATORY[/bold green]. Reference genome FASTA.",
        rich_help_panel="Input files"),
    annotation: Optional[str] = typer.Option(
        None, "--annotation", "-a",
        help="[bold green]MANDATORY[/bold green]. Annotation GFF3/GTF to evaluate.",
        rich_help_panel="Input files"),
    homologs: Optional[List[str]] = typer.Option(
        None, "--proteins", "-p",
        help="[bold green]MANDATORY[/bold green]. Homolog protein FASTA "
             "(repeat [blue]-p[/blue] per set, or 'label=path' comma list; "
             "same rules as merge --proteins).",
        rich_help_panel="Input files"),
    predictors: Optional[List[str]] = typer.Option(
        None, "--predictions", "-P",
        help="Optional de novo predictions; repeat [blue]-P[/blue] per predictor "
             "or use a 'name=path' comma list (evidence tracks + denovo_support "
             "attributes only; no scoring in evaluate).",
        rich_help_panel="Input files"),
    rna: Optional[List[str]] = typer.Option(
        None, "--rna", "-r", help="RNA evidence (same syntax as merge --rna).",
        rich_help_panel="Input files"),
    busco: Optional[str] = typer.Option(
        None, "--busco-lineage", help="Local BUSCO odb lineage path (offline).",
        rich_help_panel="Input files"),
    omark_db: Optional[str] = typer.Option(
        None, "--omark-db", help="OMAmer database (.h5) for OMArk orthology QC.",
        rich_help_panel="Input files"),
    omark_options: Optional[str] = typer.Option(
        None, "--omark-options",
        help="Extra flags forwarded to the omark assessment step, e.g. "
             "\"-e --taxid 33090 --taxonomic_rank kingdom\". Protected: "
             "-f, -d, -o, -i/--isoform_file. " + _ADV_CAUTION,
        rich_help_panel="Advanced"),
    mnp_filter: str = typer.Option(
        "UnStop,UnShift", "--mnp-filter", help=_MNP_HELP,
        rich_help_panel="Scoring"),
    goldset_min_orf: int = typer.Option(
        50, "--goldset-min-orf",
        help="Gold set: minimum ORF length (aa) for TransDecoder.LongOrfs (-m).",
        rich_help_panel="Advanced"),
    threads: int = typer.Option(
        8, "--threads", "-t", help="Threads.", rich_help_panel="Input files"),
    outdir: str = typer.Option(
        "./annomerge_eval_out", "--outdir", "-o", help="Output directory.",
        rich_help_panel="Input files"),
    viewer: str = typer.Option(
        "lite", "--viewer",
        help="\b\nEvidence browser build:\n"
             "  [cyan]full[/cyan]  embeds the genome sequence\n"
             "  [cyan]lite[/cyan]  (default) coordinates only\n"
             "  [cyan]no[/cyan]    skip",
        rich_help_panel="Assessment & report"),
    show_rna_coverage: int = typer.Option(
        1, "--show-rna-coverage",
        help="\b\nRNA coverage tracks embedded in the evidence browser:\n"
             "  [cyan]0[/cyan]  none\n"
             "  [cyan]1[/cyan]  (default) merged track over all samples only "
             "(a single track)\n"
             "  [cyan]2[/cyan]  merged + one track per sample",
        rich_help_panel="Assessment & report"),
    max_intron: Optional[int] = typer.Option(
        None, "--max-intron",
        help="Max intron length (bp) for hisat2/minimap2/miniprot. Default: "
             "derived from --ploidy as 500000*2/ploidy (diploid = 500000).",
        rich_help_panel="Advanced"),
    min_intron: int = typer.Option(
        20, "--min-intron",
        help="Min intron length (bp): hisat2 --min-intronlen.",
        rich_help_panel="Advanced"),
    ploidy: int = typer.Option(
        2, "--ploidy",
        help="Genome ploidy; only used to derive the default --max-intron "
             "(500000*2/ploidy). Ignored when --max-intron is given.",
        rich_help_panel="Advanced"),
    minimap2_options: Optional[List[str]] = typer.Option(
        None, "--minimap2-options",
        help="Extra minimap2 flags for long-read/transcripts/mix alignment. "
             "Protected: -t, -a, -G (use --max-intron), -u."
             + _TARGETED_NOTE + " " + _ADV_CAUTION,
        rich_help_panel="Advanced"),
    hisat2_options: Optional[List[str]] = typer.Option(
        None, "--hisat2-options",
        help="Extra HISAT2 alignment flags. Protected: -x, -p, -S, -1, -2, -U, "
             "--dta." + _TARGETED_NOTE + " " + _ADV_CAUTION,
        rich_help_panel="Advanced"),
    stringtie_options: Optional[List[str]] = typer.Option(
        None, "--stringtie-options",
        help="Extra StringTie flags for per-sample assembly. Protected: input "
             "BAM, -p, -o, -L." + _TARGETED_NOTE + " " + _ADV_CAUTION,
        rich_help_panel="Advanced"),
    miniprot_options: Optional[List[str]] = typer.Option(
        None, "--miniprot-options",
        help="Extra miniprot flags, e.g. \"-I -S\". Inserted into both PAF and "
             "GFF runs. Protected: -t, --paf, --gff."
             + _TARGETED_NOTE + " " + _ADV_CAUTION,
        rich_help_panel="Advanced"),
    transdecoder_options: Optional[str] = typer.Option(
        None, "--transdecoder-options", help=_TD_HELP,
        rich_help_panel="Advanced"),
    config: Optional[str] = typer.Option(
        None, "--config", help="resolved_params.yaml from a previous run (informational).",
        rich_help_panel="Advanced"),
    verbose: bool = typer.Option(False, "--verbose", help="Verbose logging.",
                                 rich_help_panel="Advanced"),
):

    # repeatable -P/-p values collapse into one comma-separated spec (as in merge)
    predictors = ",".join(predictors) if predictors else None
    homologs = ",".join(homologs) if homologs else None
    missing = [n for n, v in (("--genome/-g", genome), ("--annotation/-a", annotation),
                              ("--proteins/-p", homologs)) if not v]
    if missing:
        typer.secho(f"missing mandatory option(s): {', '.join(missing)} "
                    f"(or use -c/--check-dependencies for a dependency check)",
                    fg="red", err=True)
        raise typer.Exit(2)
    vw = viewer.lower()
    if vw not in ("full", "lite", "no"):
        typer.secho(f"unknown --viewer {viewer!r} (full|lite|no)", fg="red", err=True)
        raise typer.Exit(2)
    omark_options = _check_omark_options(omark_options)
    eff_max_intron = max_intron if max_intron else int(500000 * 2 / ploidy)
    try:
        hisat2_opts = parse_targeted_opts(
            hisat2_options, ("-x", "-p", "-S", "-1", "-2", "-U", "--dta"), "hisat2")
        stringtie_opts = parse_targeted_opts(
            stringtie_options, ("-p", "-o", "-L"), "stringtie")
        minimap2_opts = parse_targeted_opts(
            minimap2_options, ("-t", "-a", "-G", "-u"), "minimap2")
        miniprot_opts = parse_targeted_opts(
            miniprot_options, ("-t", "--paf", "--gff"), "miniprot")
        td_opts = parse_extra_opts(
            transdecoder_options, ("-t", "-m", "--retain_blastp_hits"), "TransDecoder")
    except RuntimeError as e:
        typer.secho(str(e), fg="red", err=True)
        raise typer.Exit(2)

    out = _prepare(outdir, verbose)
    from .utils import ToolOpts as _ToolOpts
    params = {(_yamlable_key(k) if isinstance(v, _ToolOpts) else k): _yamlable(v)
              for k, v in locals().items() if k not in ("out", "_ToolOpts")}
    dump_params(params, ensure_dir(out / "01.validate"))
    try:
        results, vi = _run_evidence_pipeline(
            out, genome, homologs, rna, busco, threads,
            goldset_min_orf,
            predictors=predictors, annotation=annotation,
            mnp_filter=_norm_mnp(mnp_filter),
            miniprot_opts=miniprot_opts, hisat2_opts=hisat2_opts,
            stringtie_opts=stringtie_opts, td_opts=td_opts,
            minimap2_opts=minimap2_opts, max_intron=eff_max_intron,
            min_intron=min_intron,
        )
    except InputError as e:
        typer.secho(f"input error: {e}", fg="red", err=True)
        raise typer.Exit(2)

    # ---- module 6 (evaluate): normalise the annotation to GFF3 ONLY as a
    # format conversion (GFF3 copied verbatim; GTF -> gffread --keep-genes
    # -F). NO sanity filter / renamer / coordinate lift here: the user's
    # annotation is the assessment target and stays exactly as submitted
    # (original IDs, original models; assessment tags are appended later).
    tracker.big("06", "evaluate", "normalize annotation + lift gold set")
    eval_d = ensure_dir(out / "06.evaluate")
    final_gff3 = _annotation_to_gff3(Path(annotation), eval_d / "annotation.gff3")
    # Gold set lift to genome coordinates: merge does this inside EVM;
    # evaluate has no EVM, so lift directly whenever the Gold set exists.
    # This transforms the StringTie/TransDecoder EVIDENCE transcripts
    # (goldset evidence track) — never the submitted annotation.
    rna_res = results.get("rna") or {}
    gold_dir = out / "04.goldset"
    if (gold_dir / "transcripts.fa.transdecoder.gff3").exists() and rna_res.get("stringtie_gtf"):
        from .modules.evm import lift_goldset
        lift_goldset(gold_dir, rna_res["stringtie_gtf"], eval_d,
                     eval_d / "goldset.genome.gff3", genome_fa=vi.genome)
    tracker.big_done()
    # module 6: coverage-aware evidence tiering (identical to merge;
    # sanity_tsv=None — nothing was filtered)
    tracker.big("08", "assessment", "confidence tiering + support tagging")
    from .modules.confidence import confidence_annotation, annotate_support
    assess_d = ensure_dir(out / "08.assessment")
    hom_res = results.get("homolog") or {}
    set_w = hom_res.get("set_weights") or {}
    membership = hom_res.get("membership") or {}
    protein_weights = ({p: set_w.get(lab, 1.0) for p, lab in membership.items()}
                       if len(set_w) > 1 else None)
    conf_stats = confidence_annotation(
        final_gff3, vi.genome,
        rna_res.get("junctions_bed"),
        hom_res.get("gff"),
        rna_res.get("coverage_bams") or rna_res.get("bams", []),
        eval_d / "goldset.genome.gff3",
        None,
        assess_d / "confidence.tsv",
        protein_weights=protein_weights, threads=threads)
    typer.secho(f"confidence tiers: {conf_stats or 'skipped (exists)'}", fg="green")
    # module 6b: Tier / RNA_support / Protein_support / denovo_support /
    # Flags on every gene and mRNA + CDS numbering (same as merge)
    tag_stats = annotate_support(final_gff3, assess_d / "confidence.tsv",
                                 vi.predictors or {})
    typer.secho(f"support tags: {tag_stats}", fg="green")
    tracker.big_done()
    # module 7/8: report (gene set stats + completeness + charts)
    import shlex, sys
    tracker.big("09", "result", "report + evidence browser")
    from .modules.report import run_report
    rep = run_report(
        out, final_gff3, vi.genome, eval_d,
        busco_lineage=vi.busco, threads=threads,
        omark_db=Path(omark_db) if omark_db else None,
        omark_options=omark_options,
        rna_junctions=rna_res.get("junctions_bed"),
        stringtie_gtf=rna_res.get("stringtie_gtf"),
        homologs_faa=vi.homologs,
        referee_ids=hom_res.get("referee_ids"),
        run_meta={"mode": "evaluate", "annotation": Path(annotation).name,
                  "predictors": ",".join(vi.predictors or {}),
                  "genome": Path(genome).name, "outdir": str(out),
                  "command": shlex.join(sys.argv),
                  "timestamp": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
                  "threads": threads, "homologs": homologs,
                  "rna": ",".join(params.get("rna") or []) or "none",
                  "busco": busco or "none", "omark": omark_db or "none",
                  "augment": False, "viewer": vw,
                  "mnp_filter": _norm_mnp(mnp_filter) or "off"})
    typer.secho(f"report: {rep['html']}", fg="green")
    # evidence browser (igv.js, single-file offline) — same as merge
    if vw != "no":
        from .modules.viewer import run_viewer
        vw_res = run_viewer(out, vi, results, vi.predictors or {}, {},
                            final_gff3, eval_d, eval_d, threads,
                            full=(vw == "full"),
                            show_rna_coverage=show_rna_coverage)
        typer.secho(f"evidence browser: {vw_res['html']}", fg="green")
    tracker.big_done()
    typer.secho(f"deliverables: {out}/01.validate ... 09.result", fg="green")


if __name__ == "__main__":
    app()
