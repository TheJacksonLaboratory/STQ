
include { SEQ } from '../subworkflows/sequencing'
include { IMG } from '../subworkflows/imaging'

workflow ONE {

    take:
        samples

    main:
        SEQ ( samples )

        IMG ( samples
              .join(SEQ.out))

}
