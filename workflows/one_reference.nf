
include { SEQUENCING } from './workflows/sequencing'
include { IMAGING } from './workflows/imaging'

workflow ONE_REFERENCE {

    take:
        samples

    main:
        SEQUENCING ( samples )

        IMAGING ( SEQUENCING.out.grid )

}
