#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARBITRARY_GRID } from './workflows/arbitrary_grid'
include { ONE_REFERENCE } from './workflows/one_reference'
include { TWO_REFERENCES } from './workflows/two_references'

workflow {

    Channel
    .from(file(params.input))
    .splitCsv(header:true, sep:',')
    .set{ samples }

    // ARBITRARY_GRID ( samples )

    // ONE_REFERENCE ( samples )

    TWO_REFERENCES ( samples )

}
