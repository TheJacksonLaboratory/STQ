#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAIN } from './workflows/full'

workflow {

    Channel
    .from(file(params.input))
    .splitCsv(header:true, sep:',')
    .set{ slides }

    MAIN ( slides )

}
