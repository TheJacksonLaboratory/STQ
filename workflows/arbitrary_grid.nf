
include { IMAGING } from './workflows/imaging'

workflow ARBITRARY_GRID {

    take:
        samples

    main:
        IMAGING ( samples )

}
