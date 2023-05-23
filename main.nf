#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARB } from './workflows/arbitrary_grid'
include { ONE } from './workflows/one_reference'
include { TWO } from './workflows/two_references'

workflow {

    Channel
    .from(file(params.input))
    .splitCsv(header:true, sep:',')
    .map( { it -> [ (it.sample), it ] } )
    .set{ samples }

    if ( params.workflow == "arbitrary_grid" ) {
        ARB ( samples )
    }

    if ( params.workflow == "one_reference" ) {
        ONE ( samples )
    }

    if ( params.workflow == "two_references" ) {
        TWO ( samples )
    }

}
