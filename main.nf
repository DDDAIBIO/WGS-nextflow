nextflow.enable.dsl=2

params.reads = params.reads ?: 'data/*_{R1,R2}.fastq.gz'
params.outdir = params.outdir ?: 'results'

process FASTP {
  tag "$sample_id"
  label 'high_memory'
  container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_3'
  publishDir "${params.outdir}/fastq", mode: params.publish_mode
  input:
  tuple val(sample_id), path(reads)
  output:
  tuple val(sample_id), path("${sample_id}.R{1,2}.fastq.gz"), emit: reads
  path "${sample_id}.fastp.html", emit: report
  script:
  """
  fastp -i ${reads[0]} -I ${reads[1]} \\
    -o ${sample_id}.R1.fastq.gz -O ${sample_id}.R2.fastq.gz \\
    --detect_adapter_for_pe --thread ${task.cpus} \\
    -h ${sample_id}.fastp.html -j ${sample_id}.fastp.json
  """
}

process ALIGN {
  tag "$sample_id"
  label 'high_memory'
  container 'quay.io/biocontainers/bwa-mem2:2.2.1--he70b90d_5'
  publishDir "${params.outdir}/bam", mode: params.publish_mode
  input:
  tuple val(sample_id), path(reads)
  path reference
  output:
  tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
  script:
  """
  bwa-mem2 mem -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} | \\
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
  samtools index ${sample_id}.sorted.bam
  """
}

process GATK_VARIANTS {
  tag "$sample_id"
  label 'high_memory'
  container 'broadinstitute/gatk:4.5.0.0'
  publishDir "${params.outdir}/snv_indel", mode: params.publish_mode
  input:
  tuple val(sample_id), path(bam), path(bai)
  path reference
  output:
  tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf
  script:
  """
  gatk HaplotypeCaller -R ${reference} -I ${bam} -O ${sample_id}.g.vcf.gz \\
    -ERC GVCF --native-pair-hmm-threads ${task.cpus}
  tabix -p vcf ${sample_id}.g.vcf.gz
  """
}

process CNVKIT {
  tag "$sample_id"
  label 'high_memory'
  container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'
  publishDir "${params.outdir}/cnv", mode: params.publish_mode
  input:
  tuple val(sample_id), path(bam), path(bai)
  path reference
  output:
  path "${sample_id}.cnv.*"
  script:
  """
  cnvkit.py batch ${bam} --method hybrid --fasta ${reference} \\
    --output-reference ${sample_id}.reference.cnn --output-dir . --processes ${task.cpus}
  cnvkit.py call ${sample_id}.cnr -o ${sample_id}.cnv.call.cns
  """
}

process MANTA {
  tag "$sample_id"
  publishDir "${params.outdir}/sv", mode: params.publish_mode
  container 'quay.io/biocontainers/manta:1.6.0--py39hdfd78af_4'
  input:
  tuple val(sample_id), path(bam), path(bai)
  path reference
  output:
  tuple val(sample_id), path("results/variants/diploidSV.vcf.gz"), emit: sv
  script:
  """
  configManta.py --bam ${bam} --referenceFasta ${reference} --runDir manta
  manta/runWorkflow.py -m local -j ${task.cpus}
  """
}

process ARRIBA {
  tag "$sample_id"
  errorStrategy 'ignore'
  container 'quay.io/biocontainers/arriba:2.4.0--hdfd78af_0'
  publishDir "${params.outdir}/fusion", mode: params.publish_mode
  input:
  tuple val(sample_id), path(reads)
  output:
  path "${sample_id}.arriba.fusions.tsv", emit: fusions
  script:
  """
  arriba -x ${reads[0]} -o ${sample_id}.arriba.fusions.tsv -O ${sample_id}.arriba.discarded.tsv \\
    -a ${params.arriba_annotation ?: 'blacklist.tsv.gz'} \\
    -b ${params.arriba_blacklist ?: 'blacklist_hg38.tsv.gz'} \\
    -g ${params.arriba_genome ?: 'assembly.fa'}
  """
}

workflow {
  if (!params.reference) error 'Provide --reference reference.fa (and its indexes in the same directory).'
  reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true, flat: true)
  FASTP(reads_ch)
  aligned = ALIGN(FASTP.out.reads, file(params.reference))
  GATK_VARIANTS(aligned, file(params.reference))
  CNVKIT(aligned, file(params.reference))
  MANTA(aligned, file(params.reference))
  ARRIBA(FASTP.out.reads)
}
