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