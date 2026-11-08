"""Module 5c: gene / transcript ID normalisation.

Rewrites the IDs of a GFF3 in place to a compact, coordinate-based scheme:

  gene        <seqid>.g<n>        (n = 1.. per seqid, by gene start coordinate)
  transcript  <seqid>.g<n>.t<m>   (m = 1.. within the gene, by mRNA start)

Exon / CDS / UTR features keep their ordinal suffix but inherit the new
transcript prefix (``<new_tid>.exon3``, ``<new_tid>.cds``, ...).  The
previous ID is preserved in ``Name=`` when the feature has no Name yet, so
the provenance of every feature stays traceable.

Idempotent: a file whose gene IDs already follow the scheme is left
untouched (no mtime bump), so cached downstream steps are not invalidated.

Runs twice in the pipeline: on ``evm.sanity.gff3`` (module 5b output) and
on the post-augment final annotation (isoform merge mode only; in evm mode
the final file IS the renamed ``evm.sanity.gff3``).
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

log = logging.getLogger("annomerge")

_GID_NUM = re.compile(r"^\.g(\d+)$")


def _attr_get(attrs: str, key: str) -> str | None:
    for kv in attrs.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            if k.strip() == key:
                return v.strip()
    return None


def _attr_set(attrs: str, key: str, value: str) -> str:
    """Replace key's value, or append the pair if absent."""
    out, found = [], False
    for kv in attrs.split(";"):
        if "=" in kv and kv.split("=", 1)[0].strip() == key:
            out.append(f"{key}={value}")
            found = True
        else:
            out.append(kv)
    if not found:
        out.append(f"{key}={value}")
    return ";".join(out)


def _already_normalized(genes: dict) -> bool:
    """True when every gene ID is <its own seqid>.g<n>."""
    for seqid, glist in genes.items():
        for _, gid in glist:
            if not gid.startswith(seqid + ".g"):
                return False
            if not _GID_NUM.match(gid[len(seqid):]):
                return False
    return True


def normalize_ids(gff3: Path) -> dict:
    """Rename IDs in place. Returns a small stats dict."""
    genes: dict = {}          # seqid -> [(start, old_gid)]
    txs: dict = {}            # old_gid -> [(start, old_tid)]
    tid2gene: dict = {}
    lines = open(gff3).read().splitlines()

    for ln in lines:
        if ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 9:
            continue
        ftype, attrs = f[2], f[8]
        gid = _attr_get(attrs, "ID") if ftype == "gene" else None
        if gid:
            genes.setdefault(f[0], []).append((int(f[3]), gid))
            continue
        if ftype in ("mRNA", "transcript"):
            tid = _attr_get(attrs, "ID")
            parent = (_attr_get(attrs, "Parent") or "").split(",")[0]
            if tid:
                txs.setdefault(parent, []).append((int(f[3]), tid))
                tid2gene[tid] = parent

    if not genes:
        return {"skipped": "no gene features"}
    if _already_normalized(genes):
        return {"skipped": "already normalized"}

    # assign new IDs
    gid_map: dict = {}        # old_gid -> new_gid
    tid_map: dict = {}        # old_tid -> new_tid
    n_genes = 0
    for seqid, glist in genes.items():
        for i, (_s, old_gid) in enumerate(sorted(glist), 1):
            new_gid = f"{seqid}.g{i}"
            gid_map[old_gid] = new_gid
            n_genes += 1
    for old_gid, tlist in txs.items():
        new_gid = gid_map.get(old_gid)
        if new_gid is None:
            continue
        for j, (_s, old_tid) in enumerate(sorted(tlist), 1):
            tid_map[old_tid] = f"{new_gid}.t{j}"

    def _remap_parent(attrs: str) -> str:
        parent = _attr_get(attrs, "Parent")
        if not parent:
            return attrs
        new_parents = []
        for p in parent.split(","):
            p = p.strip()
            new_parents.append(tid_map.get(p, gid_map.get(p, p)))
        return _attr_set(attrs, "Parent", ",".join(new_parents))

    out_lines = []
    for ln in lines:
        if ln.startswith("#"):
            out_lines.append(ln)
            continue
        f = ln.split("\t")
        if len(f) < 9:
            out_lines.append(ln)
            continue
        ftype, attrs = f[2], f[8]
        old_id = _attr_get(attrs, "ID")
        if ftype == "gene" and old_id in gid_map:
            new_id = gid_map[old_id]
            attrs = _attr_set(attrs, "ID", new_id)
            if _attr_get(attrs, "Name") is None:
                attrs = _attr_set(attrs, "Name", old_id)
        elif ftype in ("mRNA", "transcript") and old_id in tid_map:
            new_id = tid_map[old_id]
            attrs = _attr_set(attrs, "ID", new_id)
            if _attr_get(attrs, "Name") is None:
                attrs = _attr_set(attrs, "Name", old_id)
            attrs = _remap_parent(attrs)
        else:
            # exon / CDS / UTR / anything else carrying a Parent
            parent = _attr_get(attrs, "Parent")
            attrs = _remap_parent(attrs)
            if old_id and parent:
                # swap the old transcript prefix for the new one, keep suffix;
                # longest prefix first so t1 does not shadow t10
                for p_old, p_new in sorted(tid_map.items(),
                                           key=lambda kv: -len(kv[0])):
                    if old_id.startswith(p_old):
                        attrs = _attr_set(attrs, "ID",
                                          p_new + old_id[len(p_old):])
                        break
        f[8] = attrs
        out_lines.append("\t".join(f))

    gff3.write_text("\n".join(out_lines) + "\n")
    stats = {"genes_renamed": n_genes, "transcripts_renamed": len(tid_map)}
    log.info("ID normalize %s: %s", gff3.name, stats)
    return stats
