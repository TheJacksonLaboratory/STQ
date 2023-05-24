
include { SEQ } from '../subworkflows/sequencing_single'
include { IMG } from '../subworkflows/imaging'

workflow ONE {

    take:
        samples

    main:
        SEQ ( samples )

        if ( params.do_img_subworkflow ) {
            IMG ( samples
                  .join(SEQ.out) )
        }

}
