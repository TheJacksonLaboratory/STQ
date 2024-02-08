
include { GUNZIP_FASTA as UNPACK_HOST;
          GUNZIP_FASTA as UNPACK_GRAFT;
        } from '../modules/local/gunzip'

include { XENOME_GENERATE_INDEX;
          XENGSORT_GENERATE_INDEX;
        } from '../modules/local/deconvolution'
        
workflow XINDEX {

    main:
        if ( file(params.deconvolution_reference_host).getExtension() == "gz" ) {
            UNPACK_HOST ( params.deconvolution_reference_host )
            reference_host = UNPACK_HOST.out
        }
        else {
            reference_host = params.deconvolution_reference_host
        }

        if ( file(params.deconvolution_reference_graft).getExtension() == "gz" ) {
            UNPACK_GRAFT ( params.deconvolution_reference_graft )
            reference_graft = UNPACK_GRAFT.out
        }
        else {
            reference_graft = params.deconvolution_reference_graft
        }

        if ( params.deconvolution_tool == "xenome" ) {
            XENOME_GENERATE_INDEX ( reference_host, reference_graft, params.deconvolution_kmer_size )
            
            output = XENOME_GENERATE_INDEX.out.indices_path
        }
        else if ( params.deconvolution_tool == "xengsort" ) {
            XENGSORT_GENERATE_INDEX ( reference_host, reference_graft, params.deconvolution_kmer_size )
            
            output = XENGSORT_GENERATE_INDEX.out.indices_path
        }
    
    emit:
       output

}
