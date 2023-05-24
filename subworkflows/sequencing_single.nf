
include { LOAD_SAMPLE_INFO;
        } from '../modules/local/load'
        
include { GUNZIP as UNPACK_FASTQ;
        } from '../modules/local/gunzip'
        
include { SPACERANGER;
          RETURN_SPACERANGER_ALIGNMENT_SINGLE;
        } from '../modules/local/spaceranger'
     
include { GET_REFERENCE_PILEUP;
          GET_PILEUP_OF_BAM;
          GET_SNV_FROM_PILEUP;
        } from '../modules/local/bafextract'
        
include { CELLSORT_BAM;
          SPLICING_QUANTIFICATION;
        } from '../modules/local/velocyto'
        
include { RETURN_MTX;
        } from '../modules/local/merge'        


workflow SEQ {

    take:
        samples

    main:   
        fastqs = samples.map{[it[0], (it[1].fastq)]}
        images = samples.map{[it[0], (it[1].image)]}
    
        LOAD_SAMPLE_INFO ( samples
                           .join(fastqs)
                           .join(images) )

        UNPACK_FASTQ ( LOAD_SAMPLE_INFO.out.fastq )
        
        SPACERANGER ( UNPACK_FASTQ.out
                      .join(LOAD_SAMPLE_INFO.out.image),
                      file("${params.reference_genome}") ) 
        
        RETURN_SPACERANGER_ALIGNMENT_SINGLE ( SPACERANGER.out.metrics
                                              .join(SPACERANGER.out.spatial) )
                                              
        RETURN_MTX ( SPACERANGER.out.mtx )
        
        
        if ( params.do_snv_extract ) {
        
            GET_REFERENCE_PILEUP ( file("${params.reference_genome}") )
        
        
            GET_PILEUP_OF_BAM ( SPACERANGER.out.bam,
                                GET_REFERENCE_PILEUP.out )
        
        
            GET_SNV_FROM_PILEUP ( GET_PILEUP_OF_BAM.out,
                                  GET_REFERENCE_PILEUP.out,
                                  "baf" )

        }
        
        if ( params.do_splicing_quantification ) {
            CELLSORT_BAM ( SPACERANGER.out.bam )
        
            SPLICING_QUANTIFICATION ( CELLSORT_BAM.out
                                      .join(SPACERANGER.out.bam)
                                      .join(SPACERANGER.out.mtx),
                                      file("${params.reference_genome}"),
                                      "baf" )

        }

    emit:
        SPACERANGER.out.spatial
}
