#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMPLE_TILES } from './workflows/sample_tiles'
include { MAIN } from './workflows/full'
        
workflow {

    Channel
    .from(file(params.input))
    .splitCsv(header:true, sep:',')
    .set{ slides }

    //SAMPLE_TILES ( slides )
    
    MAIN ( slides )

}
