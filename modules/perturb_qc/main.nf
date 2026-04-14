process PERTURB_QC {
    tag "${sample_id}"
    publishDir "${params.outdir}/perturb_qc/${sample_id}", mode: 'copy'

    container 'murtiabhishek/perturb-qc:1.0.0'

    cpus   8
    memory '32 GB'

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), path("${sample_id}_clean.h5ad"), emit: h5ad

    script:
    """
    python /usr/local/bin/perturb_qc.py \
        --input ${h5ad} \
        --output ${sample_id}_clean.h5ad \
        --sample-id ${sample_id} \
        --n-top-perturbs ${params.n_top_perturbs} \
        --min-cells-per-perturb ${params.min_cells_per_perturb} \
        --min-genes ${params.min_genes} \
        --max-genes ${params.max_genes} \
        --max-mito ${params.max_mito} \
        --min-counts ${params.min_counts}
    """
}