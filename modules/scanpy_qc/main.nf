process SCANPY_QC {
    tag "${sample_id}"
    publishDir "${params.outdir}/scanpy_qc/${sample_id}", mode: 'copy'

    container 'murtiabhishek/scanpy-qc:1.0.1'

    cpus   4
    memory '32 GB'

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), path("${sample_id}_clean.h5ad"), emit: h5ad

    script:
    """
    python /usr/local/bin/scanpy_qc.py \
        --input ${h5ad} \
        --output ${sample_id}_clean.h5ad \
        --sample-id ${sample_id} \
        --min-genes ${params.min_genes} \
        --max-genes ${params.max_genes} \
        --max-mito ${params.max_mito} \
        --min-counts ${params.min_counts}
    """
}