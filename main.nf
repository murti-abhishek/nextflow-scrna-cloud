#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC     } from './modules/fastqc/main'
include { SCANPY_QC  } from './modules/scanpy_qc/main'
include { INGEST_SOMA } from './modules/ingest_soma/main'

workflow {
    input_ch = channel.fromPath(params.input)
        .splitCsv(header: true)

    fastq_ch = input_ch
        .filter { row -> row.containsKey('fastq_r1') }
        .map { row -> tuple(row.sample_id, [file(row.fastq_r1), file(row.fastq_r2)]) }

    h5ad_ch = input_ch
        .filter { row -> row.containsKey('h5ad') }
        .map { row -> tuple(row.sample_id, file(row.h5ad)) }

    FASTQC(fastq_ch)
    SCANPY_QC(h5ad_ch)
    INGEST_SOMA(SCANPY_QC.out.h5ad)
}