"""Synthetic end-to-end viewer test: S5 assembly tracks, S7 goldset
four-state, S10 dual RNA columns, S13 panel order, S12 popup patch."""
import sys, json, re, tempfile
from pathlib import Path
sys.path.insert(0, "/work/shidishen/Annomergev02/annomerge/src")
from annomerge.modules import viewer

d = Path(tempfile.mkdtemp(prefix="vw_test_"))
data = d / "data"
data.mkdir()

# genome
(data / "genome.fa").write_text(">c1\n" + "ATGGTA" * 5000 + "\n")
(data / "genome.fa.fai").write_text("c1\t30000\t4\t30000\t30001\n")

# final.gff3: gene g1 with two transcripts (t1 gold-matches isoseq-origin
# goldset tx, t2 gold-matches short-origin), gene g2 non-gold
(data / "final.gff3").write_text("""##gff-version 3
c1\tAnnoMerge\tgene\t100\t2000\t.\t+\t.\tID=g1
c1\tAnnoMerge\tmRNA\t100\t2000\t.\t+\t.\tID=t1;Parent=g1
c1\tAnnoMerge\tCDS\t100\t900\t.\t+\t0\tID=cd1a;Parent=t1
c1\tAnnoMerge\tCDS\t1100\t2000\t.\t+\t0\tID=cd1b;Parent=t1
c1\tAnnoMerge\tmRNA\t100\t1500\t.\t+\t.\tID=t2;Parent=g1
c1\tAnnoMerge\tCDS\t100\t1500\t.\t+\t0\tID=cd2;Parent=t2
c1\tAnnoMerge\tgene\t5000\t6000\t.\t-\t.\tID=g2
c1\tAnnoMerge\tmRNA\t5000\t6000\t.\t-\t.\tID=t3;Parent=g2
c1\tAnnoMerge\tCDS\t5000\t6000\t.\t-\t0\tID=cd3;Parent=t3
""")

# goldset.gff3: tx_a (origin=isoseq) == t1's CDS structure; tx_b (no
# origin attr = short) == t2's
(data / "goldset.gff3").write_text("""##gff-version 3
c1\tgoldset\tgene\t100\t2000\t.\t+\t.\tID=gg1
c1\tgoldset\tmRNA\t100\t2000\t.\t+\t.\tID=tx_a;Parent=gg1;origin=isoseq
c1\tgoldset\tCDS\t100\t900\t.\t+\t0\tID=ga1;Parent=tx_a
c1\tgoldset\tCDS\t1100\t2000\t.\t+\t0\tID=ga2;Parent=tx_a
c1\tgoldset\tgene\t100\t1500\t.\t+\t.\tID=gg2
c1\tgoldset\tmRNA\t100\t1500\t.\t+\t.\tID=tx_b;Parent=gg2
c1\tgoldset\tCDS\t100\t1500\t.\t+\t0\tID=gb1;Parent=tx_b
c1\tgoldset\tgene\t100\t1500\t.\t+\t.\tID=gg3
c1\tgoldset\tmRNA\t100\t1500\t.\t+\t.\tID=tx_c;Parent=gg3;origin=isoseq
c1\tgoldset\tCDS\t100\t1500\t.\t+\t0\tID=gc1;Parent=tx_c
c1\tgoldset\tfive_prime_UTR\t100\t199\t.\t+\t.\tParent=tx_a
c1\tgoldset\tthree_prime_UTR\t1400\t1500\t.\t+\t.\tParent=tx_b
""")

# confidence.tsv with 20 columns (has_long run)
hdr = ("gene_id\ttranscript_id\tseqid\tstrand\tstart\tend\ttier\t"
       "goldset_origin\trna_support_frac\trna_basis\tlocus_cov_frac\t"
       "locus_mean_depth\trna_reached\tprotein_support_frac\tprotein_basis\t"
       "protein_reached\tsanity_notes\tintrinsic_flags\t"
       "rna_support_short_frac\trna_support_long_frac")
rows = [
    "g1\tt1\tc1\t+\t100\t2000\tstrong\tyes\t1.000\tS:0/1;L:1/1\t0.9\t30.0\tyes\t0.500\t1/2\tyes\t.\t.\t0.000\t1.000",
    "g1\tt2\tc1\t+\t100\t1500\tmoderate\tyes\t0.900\tS:0.900;L:NA\t0.8\t25.0\tyes\t0.000\t0/1\tno\t.\tsingle_exon\t0.900\tNA",
    "g2\tt3\tc1\t-\t5000\t6000\tweak\tno\t0.000\tcoverage\t0.1\t2.0\tno\t0.000\t0/1\tno\t.\tsingle_exon\t0.000\tNA",
]
(data / "confidence.tsv").write_text(hdr + "\n" + "\n".join(rows) + "\n")

# denovo predictor
(data / "denovo_augustus.gff3").write_text("""##gff-version 3
c1\taugustus\tgene\t100\t2000\t.\t+\t.\tID=dg1
c1\taugustus\tmRNA\t100\t2000\t.\t+\t.\tID=dt1;Parent=dg1
c1\taugustus\tCDS\t100\t900\t.\t+\t0\tID=dcd1;Parent=dt1
c1\taugustus\tCDS\t1100\t2000\t.\t+\t0\tID=dcd2;Parent=dt1
""")

# assembly gtfs (merged short, merged isoseq, one per-sample)
(data / "stringtie.gtf").write_text(
    'c1\tStringTie\ttranscript\t100\t2000\t1000\t+\t.\tgene_id "STG.1"; transcript_id "STG.1.1";\n'
    'c1\tStringTie\texon\t100\t2000\t1000\t+\t.\tgene_id "STG.1"; transcript_id "STG.1.1";\n')
(data / "isoseq.stringtie.gtf").write_text(
    'c1\tStringTie\ttranscript\t100\t2000\t1000\t+\t.\tgene_id "ISG.1"; transcript_id "ISG.1.1";\n'
    'c1\tStringTie\texon\t100\t2000\t1000\t+\t.\tgene_id "ISG.1"; transcript_id "ISG.1.1";\n')
(data / "assembly_s1.gtf").write_text(
    'c1\tStringTie\ttranscript\t100\t2000\t1000\t+\t.\tgene_id "S1G.1"; transcript_id "S1G.1.1";\n'
    'c1\tStringTie\texon\t100\t2000\t1000\t+\t.\tgene_id "S1G.1"; transcript_id "S1G.1.1";\n')

# tracks.json with assembly metadata (as run_viewer writes it)
tracks = {"denovo": [{"pred": "augustus", "w": 0.9, "file": "denovo_augustus.gff3", "color": "#1f77b4"}],
          "protein": [], "coverage": {},
          "assembly": [
              {"name": "Transcript Assembly (StringTie Short-Reads)", "platform": "short", "file": "stringtie.gtf"},
              {"name": "Transcript Assembly (StringTie Long-Reads)", "platform": "isoseq", "file": "isoseq.stringtie.gtf"},
              {"name": "Assembly (StringTie Short-Reads): s1", "platform": "short", "file": "assembly_s1.gtf"}]}
(data / "tracks.json").write_text(json.dumps(tracks))

# genes.json via build_genes_json (S7/S10 path)
n = viewer.build_genes_json(data / "confidence.tsv", data / "genes.json",
                            final_gff3=data / "final.gff3",
                            denovo_gff3={"augustus": data / "denovo_augustus.gff3"},
                            goldset_gff3=data / "goldset.gff3")
genes = json.loads((data / "genes.json").read_text())
assert n == 2, n
g1 = [g for g in genes if g["gene"] == "g1"][0]
g2 = [g for g in genes if g["gene"] == "g2"][0]
assert g1["goldset"] == "both", g1["goldset"]
assert g2["goldset"] == "no", g2["goldset"]
iso = {t["tid"]: t for t in g1["isoforms"]}
assert iso["t1"]["gold"] == "isoseq", iso["t1"]["gold"]
assert iso["t2"]["gold"] == "both", iso["t2"]["gold"]
assert iso["t1"]["rna_long"] == "1.000" and iso["t1"]["rna_short"] == "0.000"
assert g1["rna_long"] == "1.000" and g1["rna_short"] == "0.900"
print("genes.json S7/S10 OK")

# split the display goldset into per-platform tracks (S-new)
sp = viewer._split_goldset_display(data / "goldset.gff3", data)
assert sp["isoseq"] == "goldset.isoseq.gff3" and sp["short"] == "goldset.short.gff3"
iso_txt = (data / "goldset.isoseq.gff3").read_text()
short_txt = (data / "goldset.short.gff3").read_text()
assert "tx_a" in iso_txt and "tx_c" in iso_txt and "tx_b" not in iso_txt
assert "tx_b" in short_txt and "tx_a" not in short_txt and "tx_c" not in short_txt
# UTR rows appended at file end must follow their PARENT mRNA's platform
assert "Parent=tx_a" in iso_txt and "Parent=tx_a" not in short_txt
assert "Parent=tx_b" in short_txt and "Parent=tx_b" not in iso_txt
# no orphan features: every Parent= reference resolves within the same file
import re as _re
for name, txt in (("isoseq", iso_txt), ("short", short_txt)):
    ids = set(_re.findall(r"ID=([^;\n]+)", txt))
    for parent in set(_re.findall(r"Parent=([^;\n]+)", txt)):
        assert parent in ids, (name, parent)
print("goldset split + UTR routing OK")

# render the offline browser
out_html = viewer.render_offline_html(data, d / "browser.html",
                                      ref_name="test.fa",
                                      fallback_locus="c1:1-200000",
                                      coverage={}, full=False, theme="npg")
html = out_html.read_text()
print("html size: %.1f MB" % (len(html) / 1048576))

# S5: assembly tracks present with platform colors
assert "Transcript Assembly (StringTie Long-Reads)" in html
assert "Transcript Assembly (StringTie Short-Reads)" in html
assert "Assembly (StringTie Short-Reads): s1" in html
th = viewer.theme_colors("npg")
assert f'color: {json.dumps(th["isoseq"])}' in html, "isoseq color missing"
assert "stringtie.gtf" not in [m for m in re.findall(r'url: "(data/[^"]+)"', html)], "leftover data refs"
print("S5 assembly tracks OK")

# S12: popup ID patch inlined
assert '(()=>{const _m=String(this.attributeString).match(/(?:ID=|transcript_id "|gene_id ")([^;"]+)/)' in html
print("S12 popup patch OK")

# S13: panel order — Overview < Structure < Gene list indicators <
# (Orthology conditional absent here) < De novo < Flags < Transcripts
order = [html.find(x) for x in (
    '<div class="dt-title">Overview</div>',
    '<div class="dt-title">Structure</div>',
    '<div class="dt-title">Gene list indicators</div>',
    'De novo evidence（CDS coverage）',
    '<div class="dt-title">Flags</div>',
    '<div class="dt-title">Transcripts (')]
assert all(x > 0 for x in order), order
assert order == sorted(order), order
assert "(from final GFF3)" not in html
print("S13 panel order OK")

# S7: gold filter panel present; four-state labels + no-GoldSet option
assert "flt-gold-panel" in html and "long-read + short-read" in html
assert "no GoldSet evidence" in html
print("S7 filter OK")

# per-platform GoldSet tracks (long-read above short-read)
i_lr = html.find("GoldSet (TransDecoder Long-Reads)")
i_sr = html.find("GoldSet (TransDecoder Short-Reads)")
assert 0 < i_lr < i_sr
print("goldset tracks OK")

# no leftover unresolved tokens
for tok in ("__ASSEMBLY_TRACKS__", "__PROTEIN_TRACKS__", "__COVERAGE_TRACK__",
            "__STRINGTIE_COLOR__", "__GOLDSET_COLOR__"):
    assert tok not in html, tok
leftover = [l for l in html.splitlines() if 'url: "data/' in l or 'fetch("data/' in l]
assert not leftover, leftover[:3]
print("tokens/refs OK")
print("VIEWER TEST PASS ->", out_html)
