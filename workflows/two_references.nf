
include { XINDEX } from '../subworkflows/xenome_index'
include { SEQ } from '../subworkflows/sequencing'
include { IMG } from '../subworkflows/imaging'

workflow TWO {

    take:
        samples

    main:
        // Need to accomodate XENGSORT structure
        if ( !file("${params.deconvolution_indices_path}/${params.deconvolution_indices_name}-both.kmers.low-bits.lwr").exists() ) {
            XINDEX ( )
        }

        SEQ ( samples )

        if ( params.do_img_subworkflow ) {
            IMG ( samples
                  .join(SEQ.out) )
        }

}
