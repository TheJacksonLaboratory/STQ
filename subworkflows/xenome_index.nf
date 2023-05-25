
include { GUNZIP_FASTA as UNPACK_HOST;
          GUNZIP_FASTA as UNPACK_GRAFT;
        } from '../modules/local/gunzip'

include { XENOME_GENERATE_INDEX;
        } from '../modules/local/deconvolution'
        
workflow XINDEX {

    main:
        if ( file(params.xenome_reference_host).getExtension() == "gz" ) {
            UNPACK_HOST ( params.xenome_reference_host )
            reference_host = UNPACK_HOST.out
        }
        else {
            reference_host = params.xenome_reference_host
        }

        if ( file(params.xenome_reference_graft).getExtension() == "gz" ) {
            UNPACK_GRAFT ( params.xenome_reference_graft )
            reference_graft = UNPACK_GRAFT.out
        }
        else {
            reference_graft = params.xenome_reference_graft
        }

        XENOME_GENERATE_INDEX ( reference_host, reference_graft, params.xenome_kmer_size )
    
    emit:
       XENOME_GENERATE_INDEX.out.indices_path

}
