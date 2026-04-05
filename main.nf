#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Minimal pipeline: run FastQC on input reads

process FASTQC {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastqc/${sample_id}", mode: 'copy'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"),  emit: zip

    script:
    """
    fastqc ${reads} --threads ${task.cpus}
    """
}

workflow {
    reads_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, [file(row.fastq_r1), file(row.fastq_r2)]) }

    FASTQC(reads_ch)
}