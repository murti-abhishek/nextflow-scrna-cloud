process INGEST_SOMA {
    tag "${sample_id}"

    container 'murtiabhishek/soma-ingest:1.0.2'

    cpus   4
    memory { 32.GB * task.attempt }
    time   '2h'

    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'finish' }
    maxRetries 3

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), val("${params.atlas_uri}/${sample_id}"), emit: soma_uri

    script:
    """
    python /usr/local/bin/ingest_soma.py \
        --input ${h5ad} \
        --uri ${params.atlas_uri} \
        --sample-id ${sample_id} \
        --tissue ${params.tissue} \
        --organism "${params.organism}" \
        --assay "${params.assay}"
    """
}