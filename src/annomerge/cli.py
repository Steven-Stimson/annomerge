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
    help="""AnnoMerge: EVM weight calibration + evidence auditing for gene annotation.

Subcommands:
  merge     Integration pipeline — score de novo predictors and calibrate
            EVM weights -> EVM integration -> isoform resolution / UTR
            backfill -> sanity repair -> confidence tiers -> HTML report
            + single-file offline evidence browser.
            Inputs: -g genome -P de novo predictions -p homolog proteins
            [-r/--bam/--transcripts RNA] [--busco lineage] [--omark-db db]
  evaluate  Assessment pipeline — run the SAME assessment chain as merge on
            an existing annotation (evidence collection -> confidence tiers
            -> per-gene support tags -> report + evidence browser).
            The annotation is read-only: never filtered, renamed or
            coordinate-lifted.
            Inputs: -g genome -a annotation GFF3/GTF -p homolog proteins
            [-r/--bam/--transcripts RNA] [-P de novo evidence tracks]
            [--busco lineage] [--omark-db db]

Global options:
  -c, --check-dependencies  check all external dependencies and versions
                            (exit 0 = all present, 3 = missing)
  -v, --version             print version and exit
  -h, --help                show this help (use <subcommand> --help for
                            subcommand options)
""",
    context_settings={"help_option_names": ["-h", "--help"]},
    add_completion=False,
    no_args_is_help=True,
)


def _version_callback(value: bool):
    if value:
        typer.echo(f"annomerge {__version__}")
        raise typer.Exit()


def _deps_callback(value: bool):
    if value:
        from .deps import check_all, format_table
        typer.echo(format_table(check_all()))
        n_bad = sum(1 for r in check_all() if not r["ok"])
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
    return out

def _stub(stage: str):
    import logging
    logging.getLogger("annomerge").warning(
        "stage %r is not implemented in this development build; skipped.", stage)


def _run_evidence_pipeline(
    out: Path, genome: str, homologs: str, rna_list: Optional[str],
    bam: Optional[str], transcripts: Optional[str], busco: Optional[str],
    species_map: Optional[str], threads: int,
    referee_min_wi: float, goldset: bool, goldset_min_orf: int,
    predictors: Optional[str] = None, annotation: Optional[str] = None,
    mnp_filter: str = "",
):
    vi = validate(
        genome=genome, homologs=homologs, predictors=predictors,
        annotation=annotation, rna_list=rna_list, bam=bam,
        transcripts=transcripts, busco=busco, species_map=species_map,
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
    results["homolog"] = run_homolog_branch(vi, out, threads, referee_min_wi=referee_min_wi,
                       mnp_filter=mnp_filter)

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
    homologs: Optional[str] = typer.Option(None, "-p", help="MANDATORY (unless -c): homolog protein FASTA, 'label=path' comma list OK (e.g. -p near=near.pep,far=far.fa; unlabeled -> filename stem). headers: >id|species=Genus_species. Multiple UNLABELED files merge into one protein pool (mnp-filtered, unweighted); explicit id=path labels keep the relatedness-weighted multi-set mode"),
    rna_list: Optional[str] = typer.Option(None, "-r", help="RNA list file (fastq/bam/gtf tags)"),
    bam: Optional[str] = typer.Option(None, "--bam", help="sorted+indexed RNA-seq BAM(s), comma list"),
    transcripts: Optional[str] = typer.Option(None, "--transcripts", help="StringTie transcripts GTF/GFF"),
    busco: Optional[str] = typer.Option(None, "--busco", help="local BUSCO odb lineage path (offline)"),
    paf: Optional[str] = typer.Option(None, "--paf", help="precomputed miniprot PAF"),
    mnp_filter: Optional[str] = typer.Option(None, "--mnp-filter",
        help="miniprot homolog quality filter, comma tokens: UnStop,UnShift,Q2R=0.5,topRank=1 ... "
             "(empty = no filter)"),
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
    merge_mode: str = typer.Option("isoform", "--merge-mode",
                                   help="isoform|evm (isoform: EVM + isoform resolution + UTR backfill + augment; evm: EVM integration only)"),
    omark_db: Optional[str] = typer.Option(None, "--omark-db", help="OMAmer database (.h5) for OMArk orthology QC (after BUSCO)"),
    viewer_off: bool = typer.Option(False, "--no-viewer", help="skip the evidence browser build"),
    viewer_full: bool = typer.Option(False, "--viewer-full",
                                     help="evidence browser: embed genome sequence (base-level track). "
                                          "Default is the lite build (coordinates only, much smaller)"),
    config: Optional[str] = typer.Option(None, "--config", help="resolved_params.yaml from a previous run"),
    verbose: bool = typer.Option(False, "--verbose"),
):
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
            out, genome, homologs, rna_list, bam, transcripts, busco,
            species_map, threads, referee_min_wi, goldset, goldset_min_orf,
            predictors=predictors, mnp_filter=mnp_filter or "",
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
    # isoform resolution needs RNA evidence (StringTie + gold set); with
    # homolog proteins only the pipeline degrades to EVM-only integration
    # and the viewer shows the EVM result (final = evm.sanity.gff3).
    if not vi.has_rna and mode == "isoform":
        mode = "evm"
        typer.secho(
            "no RNA evidence: --merge-mode isoform unavailable "
            "(isoform resolution / UTR backfill need RNA); "
            "falling back to EVM-only integration", fg="yellow", err=True)
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
    # module 6b: tag the final GFF3 with the evidence verdict (Tier /
    # RNA_support / Protein_support / denovo_support / Flags on every gene
    # and mRNA, viewer wording) + per-transcript CDS numbering (.cds1..N).
    # Always runs (cheap, idempotent) so downstream deliverables stay tagged
    # even when confidence.tsv itself is lazy-skipped.
    from .modules.confidence import annotate_support
    tag_stats = annotate_support(final_gff3, assess_d / "confidence.tsv",
                                 preds or {})
    typer.secho(f"support tags: {tag_stats}", fg="green")
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
                  "merge_mode": mode,
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


@app.command(help="Evidence auditing / QC of an existing annotation (no EVM integration).")
def evaluate(
    genome: Optional[str] = typer.Option(None, "-g", help="MANDATORY (unless -c): reference genome FASTA"),
    annotation: Optional[str] = typer.Option(None, "-a", help="MANDATORY (unless -c): annotation GFF3/GTF to evaluate"),
    homologs: Optional[str] = typer.Option(None, "-p", help="MANDATORY (unless -c): homolog protein FASTA, 'label=path' comma list OK (same rules as merge)"),
    predictors: Optional[str] = typer.Option(None, "-P", help="optional de novo predictions 'name=path' comma list (evidence tracks + denovo_support attributes only; no scoring in evaluate)"),
    rna_list: Optional[str] = typer.Option(None, "-r", help="RNA list file (fastq/bam/gtf tags)"),
    bam: Optional[str] = typer.Option(None, "--bam", help="sorted+indexed RNA-seq BAM(s), comma list"),
    transcripts: Optional[str] = typer.Option(None, "--transcripts", help="StringTie transcripts GTF/GFF"),
    busco: Optional[str] = typer.Option(None, "--busco", help="local BUSCO odb lineage path (offline)"),
    omark_db: Optional[str] = typer.Option(None, "--omark-db", help="OMAmer database (.h5) for OMArk orthology QC (after BUSCO)"),
    paf: Optional[str] = typer.Option(None, "--paf", help="precomputed miniprot PAF"),
    mnp_filter: Optional[str] = typer.Option(None, "--mnp-filter",
        help="miniprot homolog quality filter, comma tokens: UnStop,UnShift,Q2R=0.5,topRank=1 ... "
             "(empty = no filter; same semantics as merge)"),
    species_map: Optional[str] = typer.Option(None, "--species-map", help="queryID -> species TSV"),
    referee_min_wi: float = typer.Option(0.6, "--referee-min-wi", help="referee filter wi=cov*ident"),
    goldset: bool = typer.Option(True, "--goldset/--no-goldset",
                                 help="build the Gold set from RNA evidence (enables goldset-origin support; needs RNA input)"),
    goldset_min_orf: int = typer.Option(50, "--goldset-min-orf", help="min ORF length (aa) in Gold set"),
    threads: int = typer.Option(8, "-t", help="threads"),
    outdir: str = typer.Option("./annomerge_eval_out", "-o", help="output directory"),
    viewer_off: bool = typer.Option(False, "--no-viewer", help="skip the evidence browser build"),
    viewer_full: bool = typer.Option(False, "--viewer-full",
                                     help="evidence browser: embed genome sequence (base-level track). "
                                          "Default is the lite build (coordinates only, much smaller)"),
    config: Optional[str] = typer.Option(None, "--config", help="resolved_params.yaml from a previous run"),
    verbose: bool = typer.Option(False, "--verbose"),
):
    """Evaluate an existing annotation with the SAME assessment chain as the
    integration flow: evidence collection (RNA / homolog proteins / Gold
    set) -> coverage-aware confidence tiers -> support tags on every gene
    and mRNA -> HTML report (BUSCO + OMArk + stats) -> single-file offline
    evidence browser.

    The submitted annotation is treated as READ-ONLY: it is never filtered,
    renamed or coordinate-transformed (no sanity_filter / renamer here —
    those exist to repair pipeline-generated models, not user data). The
    only transformation is a format normalisation to GFF3 (GFF3 copied
    verbatim; GTF converted with gffread --keep-genes -F), and assessment
    tags (Tier / RNA_support / Protein_support / denovo_support / Flags)
    are APPENDED to each gene/mRNA. The Gold-set lift to genome coordinates
    operates on the StringTie/TransDecoder EVIDENCE transcripts (not on
    the submitted annotation) and only feeds the goldset evidence track.
    RNA inputs, -P predictors, --busco and --omark-db are all optional;
    only -g/-a/-p are mandatory.
    """
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
        results, vi = _run_evidence_pipeline(
            out, genome, homologs, rna_list, bam, transcripts, busco,
            species_map, threads, referee_min_wi, goldset, goldset_min_orf,
            predictors=predictors, annotation=annotation,
            mnp_filter=mnp_filter or "",
        )
    except InputError as e:
        typer.secho(f"input error: {e}", fg="red", err=True)
        raise typer.Exit(2)

    # ---- module 6 (evaluate): normalise the annotation to GFF3 ONLY as a
    # format conversion (GFF3 copied verbatim; GTF -> gffread --keep-genes
    # -F). NO sanity filter / renamer / coordinate lift here: the user's
    # annotation is the assessment target and stays exactly as submitted
    # (original IDs, original models; assessment tags are appended later).
    eval_d = ensure_dir(out / "06.evaluate")
    final_gff3 = _annotation_to_gff3(Path(annotation), eval_d / "annotation.gff3")
    # Gold set lift to genome coordinates: merge does this inside EVM;
    # evaluate has no EVM, so lift directly whenever the Gold set exists.
    # This transforms the StringTie/TransDecoder EVIDENCE transcripts
    # (goldset evidence track) — never the submitted annotation.
    rna = results.get("rna") or {}
    gold_dir = out / "04.goldset"
    if (gold_dir / "transcripts.fa.transdecoder.gff3").exists() and rna.get("stringtie_gtf"):
        from .modules.evm import lift_goldset
        lift_goldset(gold_dir, rna["stringtie_gtf"], eval_d,
                     eval_d / "goldset.genome.gff3", genome_fa=vi.genome)
    # module 6: coverage-aware evidence tiering (identical to merge;
    # sanity_tsv=None — nothing was filtered)
    from .modules.confidence import confidence_annotation, annotate_support
    assess_d = ensure_dir(out / "08.assessment")
    hom_res = results.get("homolog") or {}
    set_w = hom_res.get("set_weights") or {}
    membership = hom_res.get("membership") or {}
    protein_weights = ({p: set_w.get(lab, 1.0) for p, lab in membership.items()}
                       if len(set_w) > 1 else None)
    conf_stats = confidence_annotation(
        final_gff3, vi.genome,
        rna.get("junctions_bed"),
        hom_res.get("gff"),
        rna.get("bams", []),
        eval_d / "goldset.genome.gff3",
        None,
        assess_d / "confidence.tsv",
        protein_weights=protein_weights)
    typer.secho(f"confidence tiers: {conf_stats or 'skipped (exists)'}", fg="green")
    # module 6b: Tier / RNA_support / Protein_support / denovo_support /
    # Flags on every gene and mRNA + CDS numbering (same as merge)
    tag_stats = annotate_support(final_gff3, assess_d / "confidence.tsv",
                                 vi.predictors or {})
    typer.secho(f"support tags: {tag_stats}", fg="green")
    # module 7/8: report (gene set stats + completeness + charts)
    import shlex, sys
    from .modules.report import run_report
    rep = run_report(
        out, final_gff3, vi.genome, eval_d,
        busco_lineage=vi.busco, threads=threads,
        omark_db=Path(omark_db) if omark_db else None,
        run_meta={"mode": "evaluate", "annotation": Path(annotation).name,
                  "predictors": ",".join(vi.predictors or {}),
                  "genome": Path(genome).name, "outdir": str(out),
                  "command": shlex.join(sys.argv),
                  "timestamp": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
                  "threads": threads, "homologs": homologs,
                  "rna": rna_list or bam or transcripts or "none",
                  "busco": busco or "none", "omark": omark_db or "none",
                  "merge_mode": "evaluate", "goldset": goldset,
                  "referee_min_wi": referee_min_wi})
    typer.secho(f"report: {rep['html']}", fg="green")
    # evidence browser (igv.js, single-file offline) — same as merge
    if not viewer_off:
        from .modules.viewer import run_viewer
        vw = run_viewer(out, vi, results, vi.predictors or {}, {},
                        final_gff3, eval_d, eval_d, threads,
                        full=viewer_full)
        typer.secho(f"evidence browser: {vw['html']}", fg="green")
    typer.secho(f"deliverables: {out}/01.validate ... 09.result", fg="green")


if __name__ == "__main__":
    app()
