
include { XINDEX } from '../subworkflows/xenome_index'
include { SEQ } from '../subworkflows/sequencing'
include { IMG } from '../subworkflows/imaging'

workflow TWO {

    take:
        samples

    main:
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

        SEQ ( samples )

        if ( params.do_img_subworkflow ) {
            IMG ( samples
                  .join(SEQ.out) )
        }

}
