
include { XENOME_GENERATE_INDEX;
        } from '../modules/local/deconvolution'
        
workflow XENOME_INDEX {

    main:
        XENOME_GENERATE_INDEX ( params.xenome_reference_host, params.xenome_reference_graft, params.xenome_kmer_size )

    emit:
       XENOME_GENERATE_INDEX.out.indices_path

}
