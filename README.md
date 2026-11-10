# AnnoMerge

AnnoMerge: automated EVM weight calibration and evidence auditing for gene annotation in non-model organisms.

> Current status: both entry points run end to end. `merge` scores de novo
> predictors, calibrates EVM weights and integrates evidence (EVM + isoform
> resolution + augment) into a final annotation; `evaluate` audits an
> existing annotation through the identical assessment chain (evidence ->
> sanity -> confidence tiers -> support tags -> HTML report + offline
> evidence browser).

## Design documents

- 技术路线: ../AnnoMerge_技术路线.html
- CLI 设计 (v2.1): ../AnnoMerge_CLI设计_v2.1.md

## Installation

```bash
# external tools + python deps (pod environment)
bash scripts/setup_env.sh
source env.sh
```

## Usage
```bash
# merge (integration): scoring -> EVM -> augment -> confidence tiers -> report + browser
annomerge merge -g genome.fa -P helixer=helixer.gff,annevo=annevo.gff \
    -p homologs.faa -r rnaseq.list --busco odb10_lineage/ -t 16 -o result/

# evaluate (audit an existing annotation with the SAME assessment chain;
# only -g/-a/-p are mandatory, RNA / predictors / BUSCO / OMArk are optional)
annomerge evaluate -g genome.fa -a annotation.gff3 -p homologs.faa \
    -r rnaseq.list --busco odb10_lineage/ -t 16 -o result_eval/
```

## Repository layout

```
src/annomerge/        python package (typer CLI + modules)
  modules/validate.py input validation (fail-fast)
  modules/rna.py      RNA branch: HISAT2 -> samtools -> StringTie -> junctions.bed
  modules/homolog.py  miniprot -> PAF parsing -> wi = cov x ident referee filter
  modules/goldset.py  Gold set: gffread -> TransDecoder -> DIAMOND validation
scripts/              environment / utility scripts
tests/                pytest suite
env.sh                PATH / PYTHONPATH for the pod environment
```

## License

TBD
