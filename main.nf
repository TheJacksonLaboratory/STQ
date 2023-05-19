#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAIN } from './workflows/full'
include { WTILES } from './workflows/sample_tiles'

workflow {

    Channel
    .from(file(params.input))
    .splitCsv(header:true, sep:',')
    .set{ samples }

    MAIN ( samples )

}
