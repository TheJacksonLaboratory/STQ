#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STAGE } from String.format('./subworkflows/local/stage%s', params.stage)

workflow {

    Channel
    .from(file(params.input))
    .splitCsv(header:true, sep:',')
    .set{ sample_ids }

    STAGE ( sample_ids )

}
