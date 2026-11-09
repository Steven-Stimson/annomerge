"""annomerge CLI (typer). Two subcommands: merge / evaluate.

Current implementation scope: Module 1 (validation) + Module 2 evidence
sub-pipelines (RNA branch, miniprot homology branch, Gold set construction).
Predictor scoring / EVM / augment / confidence annotation are stubs that
log a clear "not implemented" notice.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from . import __version__
from .config import dump_params
from .modules.validate import InputError, validate
from .modules.rna import run_rna_branch
from .modules.homolog import run_homolog_branch
from .modules.goldset import run_goldset
from .utils import ensure_dir, setup_logging

app = typer.Typer(
    name="annomerge",
    help="AnnoMerge: EVM weight calibration + evidence auditing for gene annotation.",
    add_completion=False,
    no_args_is_help=True,
)


def _version_callback(value: bool):
    if value:
        typer.echo(f"annomerge {__version__}")
        raise typer.Exit()


@app.callback()
def _main(
    version: Optional[bool] = typer.Option(None, "--version", callback=_version_callback,
                                           is_eager=True, help="report version and exit"),
):
    """AnnoMerge: EVM weight calibration + evidence auditing."""


def _prepare(outdir: str, verbose: bool) -> Path:
    setup_logging(verbose)
    out = ensure_dir(Path(outdir))
    return out

def _stub(stage: str):
    import logging
    logging.getLogger("annomerge").warning(
        "stage %r is not implemented in this development build; skipped.", stage)


def _run_evidence_pipeline(
    out: Path, genome: str, homologs: str, rna_list: Optional[str],
    bam: Optional[str], transcripts: Optional[str], busco: Optional[str],
    core: Optional[str],    species_map: Optional[str], threads: int,
    referee_min_wi: float, goldset: bool, goldset_min_orf: int,
    predictors: Optional[str] = None, annotation: Optional[str] = None,
):
    vi = validate(
        genome=genome, homologs=homologs, predictors=predictors,
        annotation=annotation, rna_list=rna_list, bam=bam,
        transcripts=transcripts, busco=busco, core=core, species_map=species_map,
        workdir=str(out),
    )
    results: dict = {}

    # module 2: RNA branch
    rna_res = None
    if vi.has_rna:
        rna_res = run_rna_branch(vi, out, threads)
        results["rna"] = rna_res
    else:
        import logging
        logging.getLogger("annomerge").warning("no RNA evidence; RNA branch skipped (auto-degraded)")

    # module 2: homology branch
    results["homolog"] = run_homolog_branch(vi, out, threads, referee_min_wi=referee_min_wi)

    # module 2b: gold set
    if goldset and vi.has_rna:
        gtf = rna_res["stringtie_gtf"] if rna_res else (vi.rna_gtfs[0] if vi.rna_gtfs else None)
        results["goldset"] = run_goldset(
            vi, out, threads, stringtie_gtf=gtf,
            referee_faa=None,  # referee protein FASTA extraction: TODO, use full homologs for now
            min_orf_aa=goldset_min_orf,
        )
    return results, vi


@app.command(help="Predictor scoring + EVM weight calibration (full pipeline by default).")
def merge(
    genome: Optional[str] = typer.Option(None, "-g", help="MANDATORY (unless -c): reference genome FASTA"),
    predictors: Optional[str] = typer.Option(None, "-P", help="MANDATORY (unless -c): de novo predictions, name=path comma list"),
    homologs: Optional[str] = typer.Option(None, "-p", help="MANDATORY (unless -c): homolog protein FASTA, 'label=path' comma list OK (e.g. -p near=near.pep,far=far.fa; unlabeled -> filename stem). headers: >id|species=Genus_species. Multiple sets get per-set miniprot runs + relatedness-weighted evidence"),
    rna_list: Optional[str] = typer.Option(None, "-r", help="RNA list file (fastq/bam/gtf tags)"),
    bam: Optional[str] = typer.Option(None, "--bam", help="sorted+indexed RNA-seq BAM(s), comma list"),
    transcripts: Optional[str] = typer.Option(None, "--transcripts", help="StringTie transcripts GTF/GFF"),
    busco: Optional[str] = typer.Option(None, "--busco", help="local BUSCO odb lineage path (offline)"),
    core: Optional[str] = typer.Option(None, "--core", help="custom core conserved protein FASTA"),
    paf: Optional[str] = typer.Option(None, "--paf", help="precomputed miniprot PAF"),
    species_map: Optional[str] = typer.Option(None, "--species-map", help="queryID -> species TSV"),
    weight_only: bool = typer.Option(False, "--weight-only", help="stop after weight table"),
    threads: int = typer.Option(8, "-t", help="threads"),
    outdir: str = typer.Option("./annomerge_out", "-o", help="output directory"),
    referee_min_wi: float = typer.Option(0.6, "--referee-min-wi", help="referee filter wi=cov*ident"),
    goldset: bool = typer.Option(True, "--goldset/--no-goldset", help="build Gold set (module 2b)"),
    goldset_min_orf: int = typer.Option(50, "--goldset-min-orf", help="min ORF length (aa) in Gold set"),
    w_rna: float = typer.Option(0.45, "--w-rna", help="RNA sub-block weight"),
    w_homolog: float = typer.Option(0.25, "--w-homolog", help="homolog sub-block weight"),
    w_conserve: float = typer.Option(0.15, "--w-conserve", help="conservation sub-block weight"),
    w_structure: float = typer.Option(0.15, "--w-structure", help="structure sub-block weight"),
    rna_weights: str = typer.Option("0.5,0.35,0.15", "--rna-weights", help="R_s internals alpha,beta,gamma"),
    fulllength_cov: float = typer.Option(0.9, "--fulllength-cov", help="bidirectional coverage for full-length match"),
    hs_gold_weight: float = typer.Option(0.5, "--hs-gold-weight", help="Gold set share inside H_s (0-1)"),
    max_anomaly_rate: float = typer.Option(0.08, "--max-anomaly-rate", help="hard-reject anomaly threshold"),
    min_predictor_score: float = typer.Option(0.35, "--min-predictor-score", help="min W_raw to enter EVM"),
    weight_mapping: str = typer.Option("linear", "--weight-mapping", help="linear|rank"),
    ablate: Optional[str] = typer.Option(None, "--ablate", help="comma list: rna,homolog,conserve,structure,me"),
    w_gold_evm: int = typer.Option(10, "--w-gold-evm", help="EVM weight: TransDecoder gold set"),
    w_transcript_evm: int = typer.Option(8, "--w-transcript-evm", help="EVM weight: StringTie transcripts"),
    w_protein_evm: int = typer.Option(6, "--w-protein-evm", help="EVM weight: miniprot homolog proteins"),
    min_intron_evm: int = typer.Option(20, "--min-intron-evm", help="EVM min intron length"),
    evm_parts: int = typer.Option(0, "--evm-parts",
                                  help="EVM genome split: 0=auto (=threads), 1=off (single run), "
                                       "N=split whole chromosomes into N size-balanced parts run in parallel"),
    aug_jitter_bp: int = typer.Option(50, "--aug-jitter-bp",
                                      help="augment: merge UTR donors closer than this on both ends (0=off)"),
    aug_dedup: str = typer.Option("geometric", "--aug-dedup",
                                  help="augment donor dedup: geometric|evidence "
                                       "(evidence = EviAnn-style PWM gate + junction vote + TPM-ratio, no bp threshold)"),
    aug_tpm_ratio: float = typer.Option(0.1, "--aug-tpm-ratio",
                                        help="evidence dedup: keep donor if TPM >= ratio * representative TPM"),
    aug_min_tpm: float = typer.Option(0.1, "--aug-min-tpm",
                                      help="evidence dedup: absolute TPM floor for keeping a donor"),
    aug_junc_reads: int = typer.Option(3, "--aug-junc-reads",
                                       help="evidence dedup: junction-read corroboration min reads"),
    aug_collapse_tol: int = typer.Option(5, "--aug-collapse-tol-bp",
        help="isoform collapse tolerance: CDS outer boundaries within this "
             "many bp (same intron chain, same strand) count as one "
             "structure; representative = longest total UTR (0 = exact)"),
    aug_pwm_gate: float = typer.Option(-3.0, "--aug-pwm-gate",
                                       help="evidence dedup: splice-PWM gate threshold (log2 lift, calibrated)"),
    check_dependencies: bool = typer.Option(False, "--check-dependencies", "-c",
                                            help="verify all external dependencies and exit"),
    merge_mode: str = typer.Option("isoform", "--merge-mode",
                                   help="isoform|evm (isoform: EVM + isoform resolution + UTR backfill + augment; evm: EVM integration only)"),
    omark_db: Optional[str] = typer.Option(None, "--omark-db", help="OMAmer database (.h5) for OMArk orthology QC (after BUSCO)"),
    viewer_off: bool = typer.Option(False, "--no-viewer", help="skip the evidence browser build"),
    viewer_full: bool = typer.Option(False, "--viewer-full",
                                     help="evidence browser: embed genome sequence (base-level track). "
                                          "Default is the lite build (coordinates only, much smaller)"),
    config: Optional[str] = typer.Option(None, "--config", help="resolved_params.yaml from a previous run"),
    verbose: bool = typer.Option(False, "--verbose"),
    version: Optional[bool] = typer.Option(None, "--version", callback=_version_callback, is_eager=True),
):
    if check_dependencies:
        from .deps import check_all, format_table
        table = format_table(check_all())
        typer.echo(table)
        n_bad = sum(1 for r in check_all() if not r["ok"])
        raise typer.Exit(3 if n_bad else 0)
    missing = [n for n, v in (("genome/-g", genome), ("predictors/-P", predictors),
                              ("homologs/-p", homologs)) if not v]
    if missing:
        typer.secho(f"missing mandatory option(s): {', '.join(missing)} "
                    f"(or use -c/--check-dependencies for a dependency check)",
                    fg="red", err=True)
        raise typer.Exit(2)
    out = _prepare(outdir, verbose)
    params = {k: v for k, v in locals().items() if k not in ("out",)}
    dump_params(params, ensure_dir(out / "01.validate"))
    try:
        results, vi = _run_evidence_pipeline(
            out, genome, homologs, rna_list, bam, transcripts, busco, core,
            species_map, threads, referee_min_wi, goldset, goldset_min_orf,
            predictors=predictors,
        )
    except InputError as e:
        typer.secho(f"input error: {e}", fg="red", err=True)
        raise typer.Exit(2)

    # module 3-5: predictor scoring + weight mapping
    from .modules.scoring import score_one_predictor, map_weights
    from .modules.validate import parse_predictors

    a, b_, g_ = map(float, rna_weights.split(","))
    cfg = dict(rna_alpha=a, rna_beta=b_, rna_gamma=g_,
               w_rna=w_rna, w_homolog=w_homolog, w_conserve=w_conserve,
               w_structure=w_structure, fulllength_cov=fulllength_cov,
               hs_gold_weight=hs_gold_weight, max_anomaly_rate=max_anomaly_rate,
               ablate=(ablate.split(",") if ablate else []))
    preds = parse_predictors(predictors)
    score_dir = ensure_dir(out / "05.score")
    details, wraw = {}, {}
    for name, gff in preds.items():
        d = score_one_predictor(name, gff, vi, score_dir, threads, results, cfg)
        details[name] = d
        wraw[name] = None if d.get("hard_reject") else d["W_raw"]

    evm_w = map_weights(wraw, mapping=weight_mapping, min_score=min_predictor_score)
    wtab = score_dir / "evm_weights.tsv"
    with open(wtab, "w") as fh:
        fh.write("predictor\tW_raw\tR_s\tH_s\tB_c\tC_s\tevm_weight\tnote\n")
        for name in preds:
            d = details[name]
            s = d["subscores"]
            def f(x): return f"{x:.4f}" if isinstance(x, float) else "NA"
            fh.write(f"{name}\t{f(d['W_raw'])}\t{f(s['R_s'])}\t{f(s['H_s'])}\t"
                     f"{f(s['B_c'])}\t{f(s['C_s'])}\t{evm_w.get(name) or 'EXCLUDE'}\t"
                     f"{d.get('hard_reject', '')}\n")
    typer.echo(f"\n=== EVM weights ===\n{wtab.read_text()}")

    # module 4/5: integration (skipped with --weight-only)
    if weight_only:
        typer.secho("weight-only mode: integration skipped", fg="yellow")
        raise typer.Exit(0)
    evm_preds = {name: (preds[name], evm_w[name]) for name in preds
                 if evm_w.get(name) is not None}
    if not evm_preds:
        typer.secho("no predictor passed scoring; integration skipped", fg="red", err=True)
        raise typer.Exit(3)
    # merge mode: isoform (default, EVM + isoform resolution + augment) or evm
    n_pred = len(evm_preds)
    mode = merge_mode
    if mode not in ("evm", "isoform"):
        typer.secho(f"unknown --merge-mode {merge_mode!r} (isoform|evm)",
                    fg="red", err=True)
        raise typer.Exit(2)
    gold_dir = out / "04.goldset"
    goldset_dir = gold_dir if (gold_dir / "transcripts.fa.transdecoder.gff3").exists() else None
    stringtie_gtf = (results.get("rna", {}) or {}).get("stringtie_gtf")
    hom_res = results.get("homolog") or {}
    rna_res = results.get("rna") or {}
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
        threads=threads, min_intron=min_intron_evm,
        w_gold=w_gold_evm, w_transcript=w_transcript_evm, w_protein=w_protein_evm,
        evm_parts=(evm_parts if evm_parts > 0 else threads),
        protein_set_weights=protein_evm_weights,
        protein_membership=membership if multi else None)
    integ_gff3 = evm_res["evm_gff3"]
    integ_gff3 = evm_res["evm_gff3"]
    integ_d = out / "06.merge"
    typer.secho(f"EVM output: {integ_gff3}", fg="green")
    # module 5b: sanity filter (hard errors removed; notes kept)
    from .modules.sanity import sanity_filter
    sanity_stats = sanity_filter(
        integ_gff3, vi.genome,
        integ_d / "evm.sanity.gff3", integ_d / "sanity_filter.tsv")
    typer.secho(f"sanity filter: {sanity_stats or 'skipped (exists)'}", fg="green")
    # module 5c-1: ID normalisation (<seqid>.g<n>.t<m>; Name= keeps old ID)
    from .modules.renamer import normalize_ids
    typer.secho(f"ID normalize (evm): {normalize_ids(integ_d / 'evm.sanity.gff3')}", fg="green")

    # module 5c: isoform resolution + UTR backfill (isoform mode only)
    final_gff3 = integ_d / "evm.sanity.gff3"
    final_d = integ_d
    if mode == "isoform":
        from .modules.augment import run_augment
        aug = run_augment(
            vi, out / "07.augment",
            evm_gff3=integ_d / "evm.sanity.gff3",
            goldset_gff=integ_d / "goldset.genome.gff3",
            stringtie_gtf=stringtie_gtf,
            preds={name: preds[name] for name in evm_preds},
            scores={n: (wraw.get(n) or 0.0,
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
            aug["dir"] / "augment.sanity.gff3", aug["dir"] / "sanity_filter.tsv")
        typer.secho(f"augment sanity filter: {aug_sanity or 'skipped (exists)'}", fg="green")
        # module 5c-2: ID normalisation of the final annotation
        typer.secho(f"ID normalize (final): {normalize_ids(aug['dir'] / 'augment.sanity.gff3')}", fg="green")
        final_gff3 = aug["dir"] / "augment.sanity.gff3"
        final_d = aug["dir"]

    # module 6: coverage-aware evidence tiering (tier, not delete)
    from .modules.confidence import confidence_annotation
    assess_d = ensure_dir(out / "08.assessment")
    if mode == "evm":
        ensure_dir(out / "07.augment")   # evm 模式: 仅建空目录
    conf_stats = confidence_annotation(
        final_gff3, vi.genome,
        rna_res.get("junctions_bed"),
        integ_d / "protein.gff3" if (integ_d / "protein.gff3").exists()
        else hom_res.get("gff"),
        rna_res.get("bams", []),
        final_d / "goldset.genome.gff3" if (final_d / "goldset.genome.gff3").exists()
        else integ_d / "goldset.genome.gff3",
        final_d / "sanity_filter.tsv",
        assess_d / "confidence.tsv",
        protein_weights=protein_weights)
    typer.secho(f"confidence tiers: {conf_stats or 'skipped (exists)'}", fg="green")
    # module 7/8: report (run summary + stats + completeness + accuracy + charts)
    import shlex, sys
    from .modules.report import run_report
    rep = run_report(
        out, final_gff3, vi.genome, final_d,
        busco_lineage=vi.busco, threads=threads,
        omark_db=Path(omark_db) if omark_db else None,
        run_meta={"mode": mode, "predictors": ",".join(preds),
                  "genome": Path(genome).name, "outdir": str(out),
                  "command": shlex.join(sys.argv),
                  "timestamp": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
                  "threads": threads, "homologs": homologs,
                  "rna": rna_list or bam or transcripts or "none",
                  "busco": busco or "none",
                  "merge_mode": merge_mode,
                  "w_rna": w_rna, "w_homolog": w_homolog,
                  "w_conserve": w_conserve, "w_structure": w_structure,
                  "min_predictor_score": min_predictor_score,
                  "referee_min_wi": referee_min_wi, "goldset": goldset,
                  "min_intron_evm": min_intron_evm,
                  "weight_mapping": weight_mapping})
    typer.secho(f"report: {rep['html']}", fg="green")
    # evidence browser (igv.js, single-file offline) — generated automatically
    if not viewer_off:
        from .modules.viewer import run_viewer
        vw = run_viewer(out, vi, results, preds, details,
                        final_gff3, final_d, integ_d, threads,
                        full=viewer_full)
        typer.secho(f"evidence browser: {vw['html']}", fg="green")
    # publish layer removed: the module dirs themselves are the numbered
    # deliverables (01.validate ... 09.result)
    typer.secho(f"deliverables: {out}/01.validate ... 09.result", fg="green")


@app.command(help="Evidence auditing / QC of an existing annotation (no EVM).")
def evaluate(
    genome: Optional[str] = typer.Option(None, "-g", help="MANDATORY (unless -c): reference genome FASTA"),
    annotation: Optional[str] = typer.Option(None, "-a", help="MANDATORY (unless -c): annotation GFF3/GTF to evaluate"),
    homologs: Optional[str] = typer.Option(None, "-p", help="MANDATORY (unless -c): homolog protein FASTA"),
    check_dependencies: bool = typer.Option(False, "--check-dependencies", "-c",
                                            help="verify all external dependencies and exit"),
    rna_list: Optional[str] = typer.Option(None, "-r", help="RNA list file"),
    bam: Optional[str] = typer.Option(None, "--bam"),
    transcripts: Optional[str] = typer.Option(None, "--transcripts"),
    busco: Optional[str] = typer.Option(None, "--busco"),
    core: Optional[str] = typer.Option(None, "--core"),
    predictors: Optional[str] = typer.Option(None, "-P", help="optional de novo GFF3s (support rating only)"),
    paf: Optional[str] = typer.Option(None, "--paf"),
    species_map: Optional[str] = typer.Option(None, "--species-map"),
    skip_omark: bool = typer.Option(False, "--skip-omark"),
    threads: int = typer.Option(8, "-t"),
    outdir: str = typer.Option("./annomerge_eval_out", "-o"),
    referee_min_wi: float = typer.Option(0.6, "--referee-min-wi"),
    goldset: bool = typer.Option(True, "--goldset/--no-goldset"),
    goldset_min_orf: int = typer.Option(50, "--goldset-min-orf"),
    config: Optional[str] = typer.Option(None, "--config"),
    verbose: bool = typer.Option(False, "--verbose"),
    version: Optional[bool] = typer.Option(None, "--version", callback=_version_callback, is_eager=True),
):
    if check_dependencies:
        from .deps import check_all, format_table
        table = format_table(check_all())
        typer.echo(table)
        n_bad = sum(1 for r in check_all() if not r["ok"])
        raise typer.Exit(3 if n_bad else 0)
    missing = [n for n, v in (("genome/-g", genome), ("annotation/-a", annotation),
                              ("homologs/-p", homologs)) if not v]
    if missing:
        typer.secho(f"missing mandatory option(s): {', '.join(missing)} "
                    f"(or use -c/--check-dependencies for a dependency check)",
                    fg="red", err=True)
        raise typer.Exit(2)
    out = _prepare(outdir, verbose)
    params = {k: v for k, v in locals().items() if k not in ("out",)}
    dump_params(params, ensure_dir(out / "01.validate"))
    try:
        _run_evidence_pipeline(
            out, genome, homologs, rna_list, bam, transcripts, busco, core,
            species_map, threads, referee_min_wi, goldset, goldset_min_orf,
            annotation=annotation,
        )
    except InputError as e:
        typer.secho(f"input error: {e}", fg="red", err=True)
        raise typer.Exit(2)
    for stage in ("confidence annotation", "report"):
        _stub(stage)


if __name__ == "__main__":
    app()
