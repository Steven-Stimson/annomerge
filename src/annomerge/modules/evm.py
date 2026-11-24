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

from ..utils import ensure_dir, lazy, run, tracker, which
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
                 out_gff3: Path, genome_fa: Path | None = None,
                 isoseq_gtf: Path | None = None,
                 isoseq_tol_bp: int = 25) -> Path:
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
    # optional Iso-Seq origin tagging (S2): an ORF is Iso-Seq-origin when
    # its source transcript's full intron chain matches a transcript of the
    # Iso-Seq-only merged assembly (same matching semantics as
    # augment._iso_match).  Written as ;origin=isoseq on the mRNA row;
    # consumed by run_evm (EVM weight split) and the evidence browser
    # (GoldSet track split / gene-list filter).  Short-origin rows carry no
    # attribute so non-Iso-Seq runs stay byte-identical.
    iso_idx = None
    if isoseq_gtf and Path(isoseq_gtf).exists():
        from .augment import _parse_isoseq_gtf, _iso_index
        iso_idx = _iso_index(_parse_isoseq_gtf(isoseq_gtf))

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
            org_attr = ""
            if iso_idx is not None:
                from .augment import _iso_match
                din = tuple((exons[i][1] + 1, exons[i + 1][0] - 1)
                            for i in range(len(exons) - 1))
                if _iso_match(iso_idx, seqid, s_T, din, (ts, te),
                              isoseq_tol_bp) is not None:
                    org_attr = ";origin=isoseq"
            gid = f"GOLD.{oid}"
            # short display name: the ORF/transcript id. TransDecoder's
            # Name carries the full ORF/hit description (long); keep the
            # browser label and popup readable.
            nm_attr = f';Name="{oid}"'
            out.write(f"{seqid}\tgoldset\tgene\t{ts}\t{te}\t.\t{g_strand}\t.\t"
                      f"ID={gid}{nm_attr}\n")
            out.write(f"{seqid}\tgoldset\tmRNA\t{ts}\t{te}\t.\t{g_strand}\t.\t"
                      f"ID={oid};Parent={gid}{nm_attr}{org_attr}\n")
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
        tracker.note("lift validation: %d valid (ATG..stop, mod3), %d invalid",
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
    tracker.note("genome split into %d size-balanced parts: %s",
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
                     "--partition_listing", str(plist),
                     "--partition_dir", str(d)]
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
                # key = partitions.list 第 1 列的原始 accession（带点号）；
                # 目录名（第 2 列末段）被 EVM 把 '.' 换成了 '_'，不能直接当 key
                out[f[0]] = pg
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
    tracker.note("EVM protein evidence: %d matches (%s)", n_match,
             "per-set sources" if membership else "single source")
    return out_gff3


def _split_goldset_by_origin(lifted: Path, d: Path) -> dict:
    """Split a lifted goldset GFF3 into per-platform prediction sources by
    the ;origin=isoseq marker written by lift_goldset.  Returns
    {"goldset": path|None, "isoseq_goldset": path|None} — whole ORF blocks
    (gene/mRNA/exon/CDS) are routed together; the source column becomes the
    EVM prediction-source name so the weights file can score them apart."""
    out_paths = {"goldset": d / "goldset.short.src.gff3",
                 "isoseq_goldset": d / "goldset.isoseq.src.gff3"}
    bufs: dict = {"goldset": [], "isoseq_goldset": []}
    cur = "goldset"
    for line in open(lifted):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9:
            if f[2] == "mRNA":
                cur = "isoseq_goldset" if "origin=isoseq" in f[8] else "goldset"
            f[1] = cur
            line = "\t".join(f) + "\n"
        bufs[cur].append(line)
    res = {}
    for key, lines in bufs.items():
        if lines:
            with open(out_paths[key], "w") as fh:
                fh.write("##gff-version 3\n")
                fh.writelines(lines)
            res[key] = out_paths[key]
    return res


# ---------------------------------------------------------------- weights + run
def run_evm(vi, evm_dir: Path, preds: dict[str, tuple[Path, int]],
            stringtie_gtf: Path, miniprot_gff: Path, goldset_dir: Path,
            threads: int, min_intron: int = 20,
            w_gold: int = 10, w_transcript: int = 8, w_protein: int = 6,
            protein_set_weights: dict | None = None,
            protein_membership: dict | None = None,
            short_gtf: Path | None = None,
            isoseq_gtf: Path | None = None,
            isoseq_tol_bp: int = 25,
            w_isoseq_transcript: int = 10,
            w_isoseq_gold: int = 12) -> dict:
    """preds: {name: (gff3_path, evm_weight)} — already filtered by scoring.

    protein_set_weights: {label: weight} — when given with >1 entries, the
    protein evidence carries per-set sources 'miniprot_<label>' (via
    protein_membership: {protein: label}) and the weights file gets one
    PROTEIN line per set, so closer relatives count more inside EVM.
    """
    d = ensure_dir(evm_dir)
    tracker.note("evm: %d predictions + transcripts + proteins + goldset",
             len(preds))

    # 1. prep evidence files
    evidence = {}
    _iso_ok = bool(isoseq_gtf) and Path(isoseq_gtf).exists()
    if _iso_ok:
        # platform-split transcript evidence (S1): the Iso-Seq-only merged
        # assembly becomes its own TRANSCRIPT source with a higher weight —
        # full-length molecules outrank short-read assemblies.  Both sources
        # live in ONE --transcript_alignments file (EVM scores by the source
        # column, weights.txt gets one row per source).  A structure
        # assembled by both platforms appears twice (once per source) — the
        # EVM consensus then naturally favours dual-platform models.
        _parts = []
        if short_gtf and Path(short_gtf).exists():
            aln_s = d / "stringtie.short.alignment.gff3"
            if not lazy([aln_s], desc="gtf->alignment gff3 (short)"):
                run(f"{_td_home()}/util/gtf_to_alignment_gff3.pl {short_gtf} > {aln_s}",
                    shell=True, desc="gtf_to_alignment_gff3 short")
            _parts.append(rewrite_source(aln_s, "stringtie",
                                         d / "stringtie.short.src.gff3"))
        aln_i = d / "isoseq.alignment.gff3"
        if not lazy([aln_i], desc="gtf->alignment gff3 (isoseq)"):
            run(f"{_td_home()}/util/gtf_to_alignment_gff3.pl {isoseq_gtf} > {aln_i}",
                shell=True, desc="gtf_to_alignment_gff3 isoseq")
        _parts.append(rewrite_source(aln_i, "isoseq", d / "isoseq.src.gff3"))
        _tx_src = d / "transcript.src.gff3"
        if not lazy([_tx_src], desc="transcript evidence concat"):
            with open(_tx_src, "w") as out:
                out.write("##gff-version 3\n")
                for sp in _parts:
                    for line in open(sp):
                        if not line.startswith("#"):
                            out.write(line)
        evidence["stringtie"] = _tx_src
    elif stringtie_gtf and Path(stringtie_gtf).exists():
        # EVM transcript evidence must be match-style GFF3 with ID=/Target= —
        # gtf_to_alignment_gff3.pl gives exactly that (also reused by lift_goldset)
        aln = d / "stringtie.alignment.gff3"
        if not lazy([aln], desc="gtf->alignment gff3"):
            run(f"{_td_home()}/util/gtf_to_alignment_gff3.pl {stringtie_gtf} > {aln}",
                shell=True, desc="gtf_to_alignment_gff3")
        evidence["stringtie"] = rewrite_source(aln, "stringtie", d / "stringtie.src.gff3")
    goldset_files: dict = {}
    if goldset_dir and (goldset_dir / "transcripts.fa.transdecoder.gff3").exists():
        g = lift_goldset(goldset_dir, stringtie_gtf, d, d / "goldset.genome.gff3",
                         genome_fa=vi.genome,
                         isoseq_gtf=isoseq_gtf if _iso_ok else None,
                         isoseq_tol_bp=isoseq_tol_bp)
        if _iso_ok:
            # per-platform prediction sources (S2): Iso-Seq-origin ORFs get
            # their own ABINITIO_PREDICTION source with a higher weight
            goldset_files = _split_goldset_by_origin(g, d)
        else:
            goldset_files = {"goldset": rewrite_source(
                g, "goldset", d / "goldset.src.gff3")}
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
    if goldset_files.get("goldset"):
        lines.append(f"ABINITIO_PREDICTION\tgoldset\t{w_gold}")
    if goldset_files.get("isoseq_goldset"):
        lines.append(f"ABINITIO_PREDICTION\tisoseq_goldset\t{w_isoseq_gold}")
    if "stringtie" in evidence:
        lines.append(f"TRANSCRIPT\tstringtie\t{w_transcript}")
    if _iso_ok:
        lines.append(f"TRANSCRIPT\tisoseq\t{w_isoseq_transcript}")
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
            for src in ([g for _, (g, _) in pred_files.items()]
                        + [goldset_files.get("goldset"),
                           goldset_files.get("isoseq_goldset")]):
                if not src:
                    continue
                for line in open(src):
                    if line.startswith("#"):
                        continue
                    out.write(line)

    # ---- EVM execution: whole-genome (single thread) or size-balanced
    # parts, one part per thread (genome split granularity matches the
    # parallelism the user asked for with --threads)
    parts = split_genome_balanced(vi.genome, threads) if threads > 1 else [None]
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
        tracker.note("EVM split: %d parts, running with %d workers",
                     len(jobs), min(len(jobs), threads))
        with ThreadPoolExecutor(max_workers=min(len(jobs), threads)) as ex:
            futs = [ex.submit(_evm_run_one, j[0], j[1], j[2], wfile,
                              j[3], min_intron) for j in jobs]
            from concurrent.futures import as_completed
            for i, fut in enumerate(as_completed(futs), 1):
                acc_gffs.update(fut.result())
                tracker.parallel("EVM", i, len(futs))
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
    tracker.note("evm done: %s", final)
    return {"evm_gff3": final, "weights": wfile}
