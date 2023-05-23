
include { SEQUENCING } from '../subworkflows/sequencing'
include { IMAGING } from '../subworkflows/imaging'

workflow ONE_REFERENCE {

    take:
        samples

    main:
        SEQUENCING ( samples )

        IMAGING ( SEQUENCING.out.grid )

}
