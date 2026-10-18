# WGS-nextflow

Illumina 150 bp paired-end WGS workflow scaffold for FASTQ data. It performs adapter/QC trimming, BWA-MEM2 alignment, GATK HaplotypeCaller gVCF generation (SNV/INDEL), CNVkit CNV calling, Manta SV calling, and an optional Arriba fusion step.

## Run

The input naming convention is `data/{sample}_R1.fastq.gz` and `data/{sample}_R2.fastq.gz`.

```bash
nextflow run main.nf -profile docker \
  --reads 'data/*_{R1,R2}.fastq.gz' \
  --reference /refs/GRCh38.fa \
  --outdir results
```

Reference FASTA indexes (`.fai`, BWA-MEM2 indexes, and a sequence dictionary) and tool databases must be prepared before running. Containers should be pinned to validated versions in a production deployment.

Fusion detection with Arriba is meaningful for RNA-seq input; WGS-only FASTQ cannot reliably support gene-fusion calling. For WGS structural rearrangements, use the Manta VCF and annotate/filter with a caller such as GRIDSS2 or SvABA as a separate validation stage.

The current file is an extensible scaffold. Before clinical use, add cohort-aware joint genotyping, base-quality recalibration/known sites, contamination and sex checks, matched-normal tumor calling where applicable, CNV reference cohorts, SV quality filters, and VEP/ANNOVAR annotation.
