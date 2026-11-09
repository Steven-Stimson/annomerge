"""Module 4/5: EVidenceModeler integration.

Evidence hierarchy (user-designated): TransDecoder gold set (10) > StringTie
transcripts (8) > homolog proteins via miniprot (6) > de novo predictions
(weights 1-5 from module 3 scoring).

Steps:
  1. prep evidence: lift TransDecoder ORFs to genome coords, convert miniprot
     GFF to EVM PROTEIN (match/match_part), stringtie GTF -> GFF3
  2. write weights file
  3. partition -> write commands -> execute -> recombine (standard EVM flow)
"""
from __future__ import annotations

import logging
import re
import shutil
from pathlib import Path

from ..utils import ensure_dir, lazy, run, which
from .sanity import _read_genome as _read_fasta

log = logging.getLogger("annomerge")

# directory dependencies resolved via the registry (env EVM / TRANSDECODER,
# see deps.py); replaced the previously hard-coded absolute paths
def _evm_utils() -> Path:
    from ..deps import utility
    return Path(utility("evm", "EvmUtils"))


def _td_home() -> Path:
    from ..deps import home
    return home("transdecoder")


# ---------------------------------------------------------------- prep

def rewrite_source(gff_in: Path, source: str, gff_out: Path) -> Path:
    """Copy GFF3 with column 2 replaced by the canonical source name."""
    with open(gff_in) as fh, open(gff_out, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) == 9:
                f[1] = source
                line = "\t".join(f) + "\n"
            out.write(line)
    return gff_out


def stringtie_to_gff3(gtf: Path, out: Path) -> Path:
    if not lazy([out], desc="stringtie gtf->gff3"):
        run([which("gffread"), "-T", str(gtf), "-o", str(out)],
            desc="gffread -T stringtie")
    return out


def lift_goldset(goldset_dir: Path, stringtie_gtf: Path, d: Path,
                 out_gff3: Path, genome_fa: Path | None = None) -> Path:
    """TransDecoder transcript-coord ORF gff3 -> genome coords, in Python.

    Replaces TransDecoder's gtf_to_alignment_gff3.pl + cdna_alignment_orf_to_genome_orf.pl,
    which mishandle StringTie '.'-strand transcripts (single-exon, no junction
    evidence): gtf_to_alignment_gff3.pl maps any non-'+' strand to '-'
    (`($orient eq '+') ? ... : ...`), so '.' transcripts are lifted with the
    wrong strand and shifted CDS coords, and EVM validation then rejects ~95%
    of the gold set on yeast (98% single-exon genes).

    Strand logic here: gffread -w writes '.' transcripts as the plus-strand
    sequence, so ORF genome strand = stringtie_strand XOR orf_strand
    ('.' counts as '+'). CDS segments are mapped through the transcript exon
    chain; phases are recomputed in transcription order. If genome_fa is given,
    every lifted ORF is translated and validation stats are logged.
    """
    if lazy([out_gff3], desc="lift goldset to genome (python)"):
        return out_gff3

    # --- stringtie transcripts: exons + strand
    tx_exons: dict[str, list] = {}
    tx_strand: dict[str, str] = {}
    tx_seqid: dict[str, str] = {}
    for line in open(stringtie_gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        tid = m.group(1)
        tx_exons.setdefault(tid, []).append((int(f[3]), int(f[4])))
        tx_strand[tid] = f[6]
        tx_seqid[tid] = f[0]

    # --- TransDecoder ORFs on cDNA coords
    orf_cds: dict[str, list] = {}   # orf_id -> [(cs, ce)]
    orf_meta: dict[str, tuple] = {} # orf -> (cdna_id, orf_strand, name)
    for line in open(goldset_dir / "transcripts.fa.transdecoder.gff3"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if f[2] == "mRNA":
            oid = re.search(r"ID=([^;]+)", f[8]).group(1)
            nm = re.search(r'Name="([^"]+)"', f[8])
            orf_meta[oid] = (f[0], f[6], nm.group(1) if nm else "")
        elif f[2] == "CDS":
            oid = re.search(r"Parent=([^;]+)", f[8]).group(1)
            orf_cds.setdefault(oid, []).append((int(f[3]), int(f[4])))

    genome = _read_fasta(genome_fa) if genome_fa else None
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    n_ok = n_bad = 0

    with open(out_gff3, "w") as out:
        out.write("##gff-version 3\n")
        for oid, (cdna_id, orf_strand, name) in sorted(orf_meta.items()):
            segs = orf_cds.get(oid)
            exons = tx_exons.get(cdna_id)
            if not segs or not exons:
                continue
            s_T = tx_strand[cdna_id]
            exons = sorted(exons)
            # cDNA position (1-based) -> genome position
            if s_T == "-":
                # gffread -w writes minus-strand transcripts reverse-complemented
                chain = []
                for (s, e) in reversed(exons):
                    chain.extend(range(e, s - 1, -1))  # genome desc
            else:
                chain = []
                for (s, e) in exons:
                    chain.extend(range(s, e + 1))      # genome plus coords
            g_strand = {"+": "-", "-": "+"}.get(orf_strand, orf_strand) \
                if s_T == "-" else orf_strand
            # map CDS segments, splitting at exon boundaries: a cDNA CDS span
            # crossing an intron must become multiple genome CDS segments,
            # otherwise EVM validation rejects the model (no donor/acceptor
            # at the implied intron)
            # exon boundaries in cDNA coords (cumulative exon lengths), in
            # TRANSCRIPTION order — must match chain order: for '-' transcripts
            # chain walks exons genome-descending, so accumulate reversed;
            # a wrong order leaves CDS segments spanning introns (validated:
            # STRG.403.1 minus-strand ORF got one intron-spanning CDS segment)
            bounds = []
            acc = 0
            for (s, e) in (reversed(exons) if s_T == "-" else exons):
                acc += e - s + 1
                bounds.append(acc)  # cDNA position of exon end
            gsegs = []
            for (a, b) in segs:
                # split [a,b] at every exon boundary inside it
                cuts = [a] + [x + 1 for x in bounds if a <= x < b] + [b + 1]
                for u, v in zip(cuts, cuts[1:]):
                    ga, gb = chain[u - 1], chain[v - 2]
                    gsegs.append((min(ga, gb), max(ga, gb)))
            gsegs = sorted(set(gsegs))
            seqid = tx_seqid[cdna_id]
            # keep the FULL transcript exon chain (incl. UTR exons): dropping
            # UTR exons would lose UTR introns and break intron-chain matches
            # against the reference annotation
            ts, te = exons[0][0], exons[-1][1]
            gid = f"GOLD.{oid}"
            # short display name: the ORF/transcript id. TransDecoder's
            # Name carries the full ORF/hit description (long); keep the
            # browser label and popup readable.
            nm_attr = f';Name="{oid}"'
            out.write(f"{seqid}\tgoldset\tgene\t{ts}\t{te}\t.\t{g_strand}\t.\t"
                      f"ID={gid}{nm_attr}\n")
            out.write(f"{seqid}\tgoldset\tmRNA\t{ts}\t{te}\t.\t{g_strand}\t.\t"
                      f"ID={oid};Parent={gid}{nm_attr}\n")
            for (s, e) in exons:
                out.write(f"{seqid}\tgoldset\texon\t{s}\t{e}\t.\t{g_strand}\t.\t"
                          f"ID={oid}.exon{s}_{e};Parent={oid}\n")
            # phases in transcription order: phase of segment i+1 = bases to
            # skip so the codon boundary continues after segment i
            order = gsegs if g_strand == "+" else list(reversed(gsegs))
            phase = 0
            for (s, e) in order:
                out.write(f"{seqid}\tgoldset\tCDS\t{s}\t{e}\t.\t{g_strand}\t"
                          f"{phase}\tID=cds.{oid};Parent={oid}\n")
                phase = (3 - ((e - s + 1 - phase) % 3)) % 3
            # validation: translate
            if genome is not None:
                nuc = "".join(
                    genome[seqid][s - 1:e] if g_strand == "+"
                    else genome[seqid][s - 1:e].translate(comp)[::-1]
                    for s, e in order).upper()
                prot_ok = (nuc[:3] == "ATG" and nuc[-3:] in
                           ("TAA", "TAG", "TGA") and len(nuc) % 3 == 0)
                if prot_ok:
                    n_ok += 1
                else:
                    n_bad += 1
    if genome is not None:
        log.info("lift validation: %d valid (ATG..stop, mod3), %d invalid",
                 n_ok, n_bad)
    return out_gff3


# =================== genome splitting for parallel EVM =====================
def split_genome_balanced(genome_fa: Path, n_parts: int) -> list[list[str]]:
    """Group WHOLE accessions into n_parts bins, balanced by total length
    (first-fit decreasing bin packing): e.g. n_parts=3 with one very long
    chromosome -> [chr1], [chr2, chr3], [chr4, ...all small contigs].
    Chromosomes are never split, so no gene can be cut across parts."""
    lens: list[tuple[str, int]] = []
    name, L = None, 0
    for line in open(genome_fa):
        if line.startswith(">"):
            if name is not None:
                lens.append((name, L))
            name, L = line[1:].split()[0], 0
        else:
            L += len(line.strip())
    if name is not None:
        lens.append((name, L))
    n_parts = max(1, min(n_parts, len(lens)))
    bins: list[list[str]] = [[] for _ in range(n_parts)]
    sizes = [0] * n_parts
    for acc, ln in sorted(lens, key=lambda x: -x[1]):
        i = sizes.index(min(sizes))
        bins[i].append(acc)
        sizes[i] += ln
    parts = [b for b in bins if b]
    log.info("genome split into %d size-balanced parts: %s",
             len(parts), ", ".join(f"[{len(b)} acc]" for b in parts))
    return parts


def _write_part_genome(genome_fa: Path, accessions: set, out: Path) -> Path:
    keep = False
    with open(genome_fa) as fh, open(out, "w") as o:
        for line in fh:
            if line.startswith(">"):
                keep = line[1:].split()[0] in accessions
            if keep:
                o.write(line)
    return out


def _filter_gff_by_accessions(gff_in: Path, accessions: set, out: Path) -> Path:
    with open(gff_in) as fh, open(out, "w") as o:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.split("\t", 1)[0] in accessions:
                o.write(line)
    return out


def _evm_run_one(genome_fa: Path, combined: Path, evidence: dict,
                 wfile: Path, workdir: Path, min_intron: int) -> dict:
    """Run the EVM partition/write/execute/recombine/convert flow inside
    workdir on an arbitrary genome (whole genome or one size-balanced part).
    Returns {accession: evm.out.gff3 path}; resumable via per-accession
    evm.out.gff3 existence (partition dirs are not re-entrant)."""
    d = ensure_dir(workdir)
    plist = d / "partitions.list"
    cmds = d / "commands.list"

    def _accession_dirs() -> list[str]:
        if not plist.exists():
            return []
        seen: dict = {}
        for line in open(plist):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                seen.setdefault(f[1])
        return list(seen)

    evm_done = bool(_accession_dirs()) and all(
        (Path(a) / "evm.out.gff3").exists() for a in _accession_dirs())
    if evm_done:
        log.info("SKIP EVM partition/execute in %s (per-accession gff3 exist)", d)
    else:
        if plist.exists():                    # stale partition from a failed run
            for a in _accession_dirs():
                shutil.rmtree(a, ignore_errors=True)
            plist.unlink()
            cmds.unlink(missing_ok=True)
        # partition_EVM_inputs.pl dies if accession dirs already exist
        for line in open(genome_fa):
            if line.startswith(">"):
                acc_dir = d / line[1:].split()[0]
                if acc_dir.is_dir():
                    shutil.rmtree(acc_dir, ignore_errors=True)
        part_args = [which("perl"), str(_evm_utils() / "partition_EVM_inputs.pl"),
                     "--genome", str(genome_fa),
                     "--gene_predictions", str(combined),
                     "--segmentSize", "100000", "--overlapSize", "10000",
                     "--partition_listing", str(plist)]
        if "stringtie" in evidence:
            part_args += ["--transcript_alignments", str(evidence["stringtie"])]
        if "miniprot" in evidence:
            part_args += ["--protein_alignments", str(evidence["miniprot"])]
        run(part_args, log_file=d / "partition.log", desc="EVM partition", cwd=d)
        wc_args = [which("perl"), str(_evm_utils() / "write_EVM_commands.pl"),
                   "--partitions", str(plist), "--genome", str(genome_fa),
                   "--gene_predictions", str(combined), "--weights", str(wfile),
                   "--output_file_name", "evm.out",
                   "--min_intron_length", str(min_intron)]
        if "stringtie" in evidence:
            wc_args += ["--transcript_alignments", str(evidence["stringtie"])]
        if "miniprot" in evidence:
            wc_args += ["--protein_alignments", str(evidence["miniprot"])]
        run(" ".join(wc_args) + f" > {cmds}", shell=True,
            log_file=d / "write_commands.log", desc="EVM write commands", cwd=d)
        run([which("perl"), str(_evm_utils() / "execute_EVM_commands.pl"), str(cmds)],
            log_file=d / "execute.log", desc="EVM execute (serial)", cwd=d)
        run([which("perl"), str(_evm_utils() / "recombine_EVM_partial_outputs.pl"),
             "--partitions", str(plist), "--output_file_name", "evm.out"],
            log_file=d / "recombine.log", desc="EVM recombine", cwd=d)
        run([which("perl"), str(_evm_utils() / "convert_EVM_outputs_to_GFF3.pl"),
             "--partitions", str(plist), "--output_file_name", "evm.out",
             "--genome", str(genome_fa)],
            log_file=d / "convert.log", desc="EVM convert to GFF3", cwd=d)
    out: dict = {}
    for line in open(plist):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            pg = Path(f[1]) / "evm.out.gff3"
            if pg.exists():
                out[f[1].split("/")[-1]] = pg
    return out


def miniprot_to_protein_gff(miniprot_gff: Path, out_gff3: Path,
                            membership: dict | None = None) -> Path:
    """Convert miniprot mRNA/CDS GFF3 to EVM PROTEIN evidence (match/match_part).
    Protein coordinates per match_part are reconstructed from CDS order and
    phase (approximation, sufficient for EVM weighting).

    membership: protein id -> homolog set label. When provided (multi-set
    homolog input), each evidence line's source becomes 'miniprot_<label>'
    so EVM can weight each set separately via the weights file.
    """
    if lazy([out_gff3], desc="miniprot->EVM protein gff"):
        return out_gff3
    feats: dict[str, dict] = {}
    with open(miniprot_gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                m = re.search(r"ID=([^;]+)", f[8])
                t = re.search(r"Target=(\S+)\s+(\d+)\s+(\d+)", f[8])
                if m and t:
                    feats[m.group(1)] = {
                        "seqid": f[0], "s": int(f[3]), "e": int(f[4]),
                        "strand": f[6], "prot": t.group(1),
                        "plen": int(t.group(3)), "cds": [],
                    }
            elif f[2] == "CDS":
                m = re.search(r"Parent=([^;]+)", f[8])
                if m and m.group(1) in feats:
                    feats[m.group(1)]["cds"].append(
                        (int(f[3]), int(f[4]), int(f[7]) if f[7] != "." else 0))
    n_match = 0
    with open(out_gff3, "w") as out:
        out.write("##gff-version 3\n")
        for mid, r in feats.items():
            cds = sorted(r["cds"])
            if r["strand"] == "-":
                cds.reverse()
            ppos = 1
            parts = []
            for i, (s, e, phase) in enumerate(cds):
                aa = (e - s + 1 - (phase if i == 0 else 0)) // 3
                parts.append((s, e, ppos, ppos + max(aa, 1) - 1))
                ppos += max(aa, 1)
            label = membership.get(r["prot"]) if membership else None
            src = f"miniprot_{label}" if label else "miniprot"
            out.write(f"{r['seqid']}\t{src}\tmatch\t{r['s']}\t{r['e']}\t.\t"
                      f"{r['strand']}\t.\tID={mid}.match;Name={r['prot']};"
                      f"Target={r['prot']} 1 {r['plen']}\n")
            for j, (s, e, ps, pe) in enumerate(parts):
                # EVM requires ID= on EVERY evidence line; Parent= links children
                out.write(f"{r['seqid']}\t{src}\tmatch_part\t{s}\t{e}\t.\t"
                          f"{r['strand']}\t.\tID={mid}.part{j};"
                          f"Parent={mid}.match;Target={r['prot']} {ps} {pe}\n")
            n_match += 1
    log.info("EVM protein evidence: %d matches (%s)", n_match,
             "per-set sources" if membership else "single source")
    return out_gff3


# ---------------------------------------------------------------- weights + run
def run_evm(vi, evm_dir: Path, preds: dict[str, tuple[Path, int]],
            stringtie_gtf: Path, miniprot_gff: Path, goldset_dir: Path,
            threads: int, min_intron: int = 20,
            w_gold: int = 10, w_transcript: int = 8, w_protein: int = 6,
            evm_parts: int = 1,
            protein_set_weights: dict | None = None,
            protein_membership: dict | None = None) -> dict:
    """preds: {name: (gff3_path, evm_weight)} — already filtered by scoring.

    protein_set_weights: {label: weight} — when given with >1 entries, the
    protein evidence carries per-set sources 'miniprot_<label>' (via
    protein_membership: {protein: label}) and the weights file gets one
    PROTEIN line per set, so closer relatives count more inside EVM.
    """
    d = ensure_dir(evm_dir)
    log.info("module 4/5 [evm]: %d predictions + transcripts + proteins + goldset",
             len(preds))

    # 1. prep evidence files
    evidence = {}
    if stringtie_gtf and Path(stringtie_gtf).exists():
        # EVM transcript evidence must be match-style GFF3 with ID=/Target= —
        # gtf_to_alignment_gff3.pl gives exactly that (also reused by lift_goldset)
        aln = d / "stringtie.alignment.gff3"
        if not lazy([aln], desc="gtf->alignment gff3"):
            run(f"{_td_home()}/util/gtf_to_alignment_gff3.pl {stringtie_gtf} > {aln}",
                shell=True, desc="gtf_to_alignment_gff3")
        evidence["stringtie"] = rewrite_source(aln, "stringtie", d / "stringtie.src.gff3")
    if goldset_dir and (goldset_dir / "transcripts.fa.transdecoder.gff3").exists():
        g = lift_goldset(goldset_dir, stringtie_gtf, d, d / "goldset.genome.gff3",
                         genome_fa=vi.genome)
        evidence["goldset"] = rewrite_source(g, "goldset", d / "goldset.src.gff3")
    multi = protein_set_weights and len(protein_set_weights) > 1
    if miniprot_gff and Path(miniprot_gff).exists():
        evidence["miniprot"] = miniprot_to_protein_gff(
            Path(miniprot_gff), d / "protein.gff3",
            membership=protein_membership if multi else None)
    pred_files = {}
    for name, (gff, w) in preds.items():
        pred_files[name] = (rewrite_source(Path(gff), name, d / f"pred.{name}.gff3"), w)

    # 2. weights file
    wfile = d / "weights.txt"
    lines = []
    for name, (_, w) in pred_files.items():
        lines.append(f"ABINITIO_PREDICTION\t{name}\t{w}")
    if "goldset" in evidence:
        lines.append(f"ABINITIO_PREDICTION\tgoldset\t{w_gold}")
    if "stringtie" in evidence:
        lines.append(f"TRANSCRIPT\tstringtie\t{w_transcript}")
    if multi:
        for label, w in protein_set_weights.items():
            lines.append(f"PROTEIN\tminiprot_{label}\t{w}")
    elif "miniprot" in evidence:
        lines.append(f"PROTEIN\tminiprot\t{w_protein}")
    wfile.write_text("\n".join(lines) + "\n")
    log.info("weights:\n%s", wfile.read_text().strip())

    final = d / "evm.out.gff3"
    if lazy([final], desc="EVM run"):
        return {"evm_gff3": final, "weights": wfile}

    # combine all gene predictions (incl. goldset) into one file — EVM takes a single input
    combined = d / "pred.all.gff3"
    if not lazy([combined], desc="combine predictions"):
        with open(combined, "w") as out:
            out.write("##gff-version 3\n")
            for src in [g for _, (g, _) in pred_files.items()] + \
                       [evidence.get("goldset")]:
                if not src:
                    continue
                for line in open(src):
                    if line.startswith("#"):
                        continue
                    out.write(line)

    # ---- EVM execution: whole-genome (evm_parts<=1) or size-balanced parts
    parts = split_genome_balanced(vi.genome, evm_parts) if evm_parts > 1 else [None]
    acc_gffs: dict[str, Path] = {}
    if len(parts) <= 1:
        acc_gffs = _evm_run_one(vi.genome, combined, evidence, wfile, d, min_intron)
    else:
        from concurrent.futures import ThreadPoolExecutor
        jobs = []
        for i, accs in enumerate(parts, 1):
            pd = ensure_dir(d / f"part{i:02d}")
            pgenome = _write_part_genome(vi.genome, set(accs), pd / "genome.part.fa")
            pcombined = _filter_gff_by_accessions(combined, set(accs), pd / "pred.all.gff3")
            pev = {}
            for ename, epath in evidence.items():
                pev[ename] = _filter_gff_by_accessions(
                    epath, set(accs), pd / f"{ename}.part.gff3")
            jobs.append((pgenome, pcombined, pev, pd))
        log.info("EVM split: %d parts, running with %d workers",
                 len(jobs), min(len(jobs), threads))
        with ThreadPoolExecutor(max_workers=min(len(jobs), threads)) as ex:
            results = list(ex.map(lambda j: _evm_run_one(j[0], j[1], j[2], wfile,
                                                         j[3], min_intron), jobs))
        for r in results:
            acc_gffs.update(r)
    # concatenate per-accession gff3s in the ORIGINAL genome order
    if not lazy([final], desc="EVM final concat"):
        acc_order = []
        for line in open(vi.genome):
            if line.startswith(">"):
                acc_order.append(line[1:].split()[0])
        with open(final, "w") as out:
            header_done = False
            for acc in acc_order:
                pg = acc_gffs.get(acc)
                if pg is None or not pg.exists():
                    continue
                for gline in open(pg):
                    if gline.startswith("##gff-version"):
                        if not header_done:
                            out.write(gline)
                            header_done = True
                        continue
                    out.write(gline)
    log.info("module 4/5 [evm] done: %s", final)
    return {"evm_gff3": final, "weights": wfile}
