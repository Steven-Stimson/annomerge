# AnnoMerge v02 改造 — 变更说明（交付）

本目录 `annomerge/` 为工作副本，本地 git 记录全部改动（11 个 commit，`git log --oneline` 可查）。
远端 origin 指向 github.com/Steven-Stimson/annomerge，需在本机（有 GitHub 凭据处）执行
`git push origin master` 完成推送。

## 改造内容总览

### 1. RNA 输入：12 种 tag 数据模型（validate.py）

`-r` 支持内联 tag（`file:isoseq`）与 list 行尾 tag（`name: f1.fq,f2.fq :fastq`）。
tag 一览（`annomerge merge --help` 的 `-r` 帮助有完整表）：

| tag | 数据 | 比对器 | 组装 |
|---|---|---|---|
| fastq（默认） | Illumina SE/PE | hisat2 | stringtie 默认 |
| fasta | Illumina fasta | hisat2 -f | stringtie 默认 |
| isoseq | PacBio Iso-Seq | minimap2 -x splice:hq | stringtie -L |
| isoseq_hifi / isoseq_clr / ont | 长读细分 | minimap2 -x splice:hq / splice / splice | stringtie -L |
| transcripts | 已组装转录本 | minimap2 -x splice:hq -u f | stringtie 默认 -c 1 -j 1（不进 coverage） |
| mix | Illumina + 长读同一样本 | hisat2 + minimap2 双轨 | stringtie --mix |
| bam / bam_isoseq / bam_mix | 预比对 BAM | — | 直通 / --mix |
| gtf | 预组装转录本 | — | 直通 |

**行为变更**：单个 fastq 不带 tag 时，旧代码当作 Iso-Seq 长读（hisat2 直接报错），
现在按 Illumina SE 处理（hisat2）。带 `:isoseq` 可显式还原旧意图。

解析规则：未知 `:tag` 明确报错；list 行中裸 token 一律视为文件、样本 ID 必须带 `=`；
`bam_mix` 约定 Illumina BAM 在前、长读 BAM 在后。

### 2. minimap2 成为新依赖（deps.py）

长读 / transcripts / mix 样本需要 minimap2。解析顺序：`MINIMAP2`/`MINIMAP2_HOME`
环境变量 → PATH → conda 布局。Pod 内已装官方静态二进制 v2.30-r1287
（/home/shidishen/software/minimap2/minimap2）。

### 3. 新增 CLI 参数（merge / evaluate）

- `--max-intron INT`（Advanced）：统一 intron 长度上限，默认 500000*2/ploidy，
  落点 hisat2 --max-intronlen、minimap2 -G、miniprot -G、EVM max_intron。
- `--ploidy INT`（Advanced，默认 2）：参与 max-intron 默认推导。
- `--min-intron INT`（Advanced，默认 20）：落 hisat2 --min-intronlen 与 EVM
  --min_intron_length。**注意 StringTie 没有 intron 长度参数，此参数不落 StringTie。**
- `--minimap2-options STR`（Advanced）：保护 -t/-a/-G/-u（管道管理，用户不可覆盖）。
- `--partial {0,1}`（merge 独有，默认 0）：不完整 CDS（no_stop / non_atg_start）处理。
  0=丢弃（记 tsv 的 removed/partial_cds(<原因>)）；1=保留并在 mRNA 打 `partial=true`。
  extended_stop 修复成功不算 partial。evaluate 无此参数（read-only）。
- 删除 `--evm-options`（原白名单仅 min_intron 一个 key，已提升为 --min-intron）。

### 4. coverage 语义

confidence 分层与浏览器 coverage 轨道只统计真实测序深度的 BAM：
排除 transcripts（天然 1x 假象）与 mix/bam_mix 的长读侧。

### 5. 本步（e2e）新增修复

- **diamond 死锁修复**（scoring.py `_sanitize_faa`）：diamond 2.2 对含 `.` 残基的
  蛋白序列会在 "Loading query sequences" 处 futex 死锁（ath TAIR10 注释实测
  195/48133 条含点，来自 gffread -y 对模糊/trans-spliced CDS 的翻译）。
  extract_proteins 现在幂等清洗：小写残基视为合法（判定前先 upper()）；`.`、`*` 与 U/B/Z/J 等非标准字母 → X。
- **samtools sort 内存封顶**：rna.py 两处 `samtools sort` 加 `-m 256M`
  （默认 768MB/线程在小内存容器会把 hisat2/minimap2 挤到 OOM）。

## 验证

- 解析单测 tag_test/test_parse.py 20 项 PASS；--partial 单测 test_partial.py 6 项 PASS。
- toy 基因组 e2e（test_rna_e2e.py）12 项 PASS：11 个 tag 路径全覆盖。
- ath 真实数据：merge_multipool 全量（21953 基因，P=0.8321）+ ath_eval evaluate
  （isoseq :isoseq 路径，全链路到 report + 89MB 证据浏览器）；上游产物 md5 回归一致。

## 已知限制 / 后续建议

1. **evaluate 模式对 NCBI 风格注释需预清洗**：trans-splicing / RNA-editing 基因
   （`?` 链向、end<start、mRNA 缺失的孤儿 exon/CDS）会使 gffread rc=1 终止。
   目前未在管道内做预清洗（evaluate 设计为 read-only），建议后续加一个可选的
   `--preclean` 或在文档中给出清洗脚本（ath_eval/prep_fixture2.sh 可作参考）。
2. merge 模式 evaluate 侧的 BUSCO/OMArk 输入（longest_isoforms.faa）由 Python
   翻译（未知密码子→X），不经过 gffread，天然无 `.` 问题；但交付物 annomerge.pep
   在 evaluate + NCBI 注释下可能含 `.`（仅展示用途，不喂工具）。
3. miniprot 对大规模 homolog 集（~10 万蛋白）在 16c 下约需 1 小时级，属正常。

## 提交清单（git log --oneline，新→旧）

```
5d9745c step 9: review follow-ups + delivery notes
4a140a2 scoring: sanitize gffread proteins for diamond (step 8 e2e finding)
ea6a7e0 rna: docstring — per-unit stringtie opts are pipeline-managed
70d5d14 step 7: partial tsv fix + toy e2e tests
4839bd1 step 6: --partial semantics in sanity/cli
ede75ef step 5: cli merge/evaluate params + -r help
021ac1c step 4: miniprot -G max_intron passthrough
474de56 step 3: rna.py minimap2/mix/transcripts rewrite
67245f5 step 2: deps register minimap2
98e0c1e step 1: RnaSample tag model in validate.py
```
（以实际 git log 为准；以上为本次改造的全部 commit）
