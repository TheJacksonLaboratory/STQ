
include { IMAGING } from '../subworkflows/imaging'

workflow ARBITRARY_GRID {

    take:
        samples

    main:
        IMAGING ( samples )

}
