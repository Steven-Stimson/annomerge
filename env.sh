# AnnoMerge pod environment: PATH for external tools, PYTHONPATH for python deps
# usage: source env.sh
export PATH=/home/shidishen/software/miniprot:/home/shidishen/software/TransDecoder:\
/public/tools/hisat2/2.2.1:/public/tools/htslib/1.24/bin:/public/tools/stringtie/3.0.3:/public/tools/samtools/1.21/bin:\
/public/tools/gffcompare/0.12.10:/public/tools/gffread/0.12.9:\
/public/tools/hmmer/3.4:/public/tools/diamond/2.2.1:/public/tools/busco/5.7.1:$PATH
export PYTHONPATH=/work/shidishen/AnnoMerge/annomerge/src:\
/public/tools/python-packages/py3.12-bio:/public/tools/python-packages/py3.12-proteinlm:\
/home/shidishen/software/pylibs/annomerge

# directory-based dependencies (resolved by src/annomerge/deps.py)
export EVM=/home/shidishen/software/EVidenceModeler
export TRANSDECODER=/home/shidishen/software/TransDecoder
