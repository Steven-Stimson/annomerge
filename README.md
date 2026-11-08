# AnnoMerge

AnnoMerge: automated EVM weight calibration and evidence auditing for gene annotation in non-model organisms.

> Current status: **early development**. This repository currently implements and tests the
> evidence-preparation layer (Module 1 input validation, Module 2 evidence sub-pipelines:
> RNA-seq alignment/assembly, miniprot homology alignment with referee-set filtering,
> and the Gold set construction via StringTie + TransDecoder + homology validation).
> Predictor scoring, EVM integration, augment and confidence annotation modules are under development.

## Design documents

- 技术路线: ../AnnoMerge_技术路线.html
- CLI 设计 (v2.1): ../AnnoMerge_CLI设计_v2.1.md

## Installation

```bash
# external tools + python deps (pod environment)
bash scripts/setup_env.sh
source env.sh
```

## Usage (current scope)

```bash
# evidence preparation on the S. cerevisiae test dataset
annomerge evaluate -g genome.fa -a annotation.gff \
    -r rnaseq.list -p homologs.faa --skip-omark -t 16 -o result/
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
