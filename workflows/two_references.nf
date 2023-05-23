
include { XENOME_INDEX } from '../subworkflows/xenome_index'
include { SEQUENCING } from '../subworkflows/sequencing'
include { IMAGING } from '../subworkflows/imaging'

workflow TWO_REFERENCES {

    take:
        samples

    main:
        if ( !file("${params.xenome_indices_path}/${params.xenome_indices_name}-both.kmers.low-bits.lwr").exists() ) {
            // XENOME_INDEX ( )
            println DDD
        }

        SEQUENCING ( samples )
        
        SEQUENCING.out.view()

        // IMAGING ( samples
        //           .join(SEQUENCING.out) )

}
