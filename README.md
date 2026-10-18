# WGS-nextflow

Illumina 150 bp paired-end WGS workflow scaffold for FASTQ data. It performs adapter/QC trimming, BWA-MEM2 alignment, GATK HaplotypeCaller gVCF generation (SNV/INDEL), CNVkit CNV calling, Manta SV calling, and an optional Arriba fusion step.

## 运行环境

推荐在 Linux x86_64 云主机或 HPC 节点运行：

- Nextflow 23.10 或更高版本
- Docker 20.10+（或将 `nextflow.config` 中的 Docker 配置改为 Singularity/Apptainer）
- Java 11 或更高版本（Nextflow 依赖）
- 至少 16 CPU、64 GB RAM 和足够的本地临时空间；全基因组样本建议 32 CPU、128 GB RAM
- 能访问 Docker Hub、Quay.io 和 GitHub；首次运行会下载工具容器

安装 Nextflow 的示例：

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
docker --version
```

生产环境建议固定容器 digest，并通过机构镜像仓库缓存容器。不要关闭 TLS、镜像签名或校验。

## 参考文件

`--reference` 指向 GRCh38 或 GRCh37 的参考 FASTA。运行前需要在同一目录准备：

```bash
bwa-mem2 index GRCh38.fa
samtools faidx GRCh38.fa
gatk CreateSequenceDictionary -R GRCh38.fa
```

CNVkit 还需要适配所选参考版本的正常样本 reference `.cnn`；Arriba 需要对应版本的基因组、注释和 blacklist 文件。参考版本必须和注释数据库、样本分析方案保持一致。

## Run

The input naming convention is `data/{sample}_R1.fastq.gz` and `data/{sample}_R2.fastq.gz`.

```bash
nextflow run main.nf -profile docker \
  --reads 'data/*_{R1,R2}.fastq.gz' \
  --reference /refs/GRCh38.fa \
  --outdir results
```

如果使用 Apptainer，可在配置中将 `docker.enabled` 改为 `false`、`singularity.enabled` 改为 `true`，并设置机构允许的容器缓存目录。重跑时可使用 `-resume`：

```bash
nextflow run main.nf -profile docker -resume \
  --reads 'data/*_{R1,R2}.fastq.gz' \
  --reference /refs/GRCh38.fa --outdir results
```

FASTQ 文件名必须成对出现，例如 `sample01_R1.fastq.gz` 和 `sample01_R2.fastq.gz`。运行日志、任务报告和 DAG 会写入 `results/pipeline_info/`。

Reference FASTA indexes (`.fai`, BWA-MEM2 indexes, and a sequence dictionary) and tool databases must be prepared before running. Containers should be pinned to validated versions in a production deployment.

Fusion detection with Arriba is meaningful for RNA-seq input; WGS-only FASTQ cannot reliably support gene-fusion calling. For WGS structural rearrangements, use the Manta VCF and annotate/filter with a caller such as GRIDSS2 or SvABA as a separate validation stage.

The current file is an extensible scaffold. Before clinical use, add cohort-aware joint genotyping, base-quality recalibration/known sites, contamination and sex checks, matched-normal tumor calling where applicable, CNV reference cohorts, SV quality filters, and VEP/ANNOVAR annotation.
