
include { XENOME_INDEX } from './workflows/xenome_index'
include { SEQUENCING } from './workflows/sequencing'
include { IMAGING } from './workflows/imaging'

workflow TWO_REFERENCES {

    take:
        samples

    main:
        if ( path("${params.xenome_indices_path}").exists() ) {
            indices_path = params.xenome_indices_path
        }
        else {
            XENOME_INDEX ( )
            indices_path = XENOME_INDEX.out
        }

        println indices_path

        // SEQUENCING ( samples, indices_path )

        // IMAGING ( SEQUENCING.out.grid )

}
