
include { XINDEX } from '../subworkflows/xenome_index'
include { SEQ } from '../subworkflows/sequencing'
include { IMG } from '../subworkflows/imaging'

workflow TWO {

    take:
        samples

    main:
        if ( !file("${params.xenome_indices_path}/${params.xenome_indices_name}-both.kmers.low-bits.lwr").exists() ) {
            XINDEX ( )
        }

        SEQ ( samples )

        if ( params.do_img_subworkflow ) {
            IMG ( samples
                  .join(SEQ.out) )
        }

}
