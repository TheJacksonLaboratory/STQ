#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARB } from './workflows/arbitrary_grid'
include { ONE } from './workflows/one_reference'
include { TWO } from './workflows/two_references'
include { XINDEX } from './subworkflows/xenome_index'

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

    if ( params.workflow == "xenome_indices" ) {
        if ( params.deconvolution_tool == "xenome" ) {
            if ( !file("${params.deconvolution_indices_path}/${params.deconvolution_indices_name}-both.kmers.low-bits.lwr").exists() ) {
                XINDEX ( )
            }
        }
        else if ( params.deconvolution_tool == "xengsort" ) {
            if ( !file("${params.deconvolution_indices_path}/${params.deconvolution_indices_name}-xind.hash").exists() ) {
                XINDEX ( )
            }        
        }
    }

}
