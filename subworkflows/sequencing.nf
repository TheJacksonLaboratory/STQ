
include { LOAD_SAMPLE_INFO } from '../modules/local/load'
include { GUNZIP as UNPACK_FASTQ } from '../modules/local/gunzip'
include { DECONVOLUTION;
          SORT_FASTQ as SORT_FASTQ_MOUSE;
          SORT_FASTQ as SORT_FASTQ_HUMAN;
        } from '../modules/local/deconvolution'
        
include { SPACERANGER as SPACERANGER_MOUSE;
          SPACERANGER as SPACERANGER_HUMAN;
          RETURN_SPACERANGER_ALIGNMENT;
        } from '../modules/local/spaceranger'
     
include { GET_REFERENCE_PILEUP as GET_REFERENCE_PILEUP_MOUSE;
          GET_REFERENCE_PILEUP as GET_REFERENCE_PILEUP_HUMAN;
          GET_PILEUP_OF_BAM as GET_PILEUP_OF_BAM_MOUSE;
          GET_PILEUP_OF_BAM as GET_PILEUP_OF_BAM_HUMAN;
          GET_SNV_FROM_PILEUP as GET_SNV_FROM_PILEUP_MOUSE;
          GET_SNV_FROM_PILEUP as GET_SNV_FROM_PILEUP_HUMAN;
        } from '../modules/local/bafextract'
        
include { CELLSORT_BAM as CELLSORT_BAM_MOUSE;
          CELLSORT_BAM as CELLSORT_BAM_HUMAN;
          SPLICING_QUANTIFICATION as SPLICING_QUANTIFICATION_MOUSE;
          SPLICING_QUANTIFICATION as SPLICING_QUANTIFICATION_HUMAN;
        } from '../modules/local/velocyto'
        
include { MERGE_MTX;
          RETURN_SEPARATE_MTX;
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
        
        DECONVOLUTION ( UNPACK_FASTQ.out,
                        file("${params.xenome_indices_path}"),
                        params.xenome_indices_name )
        
        
        SORT_FASTQ_MOUSE ( DECONVOLUTION.out.mouse )     
         
        SORT_FASTQ_HUMAN ( DECONVOLUTION.out.human )


        SPACERANGER_MOUSE ( SORT_FASTQ_MOUSE.out
                            .join(LOAD_SAMPLE_INFO.out.image),
                            file("${params.mouse_reference_genome}") )
        
        SPACERANGER_HUMAN ( SORT_FASTQ_HUMAN.out
                            .join(LOAD_SAMPLE_INFO.out.image),
                            file("${params.human_reference_genome}") ) 
        
        RETURN_SPACERANGER_ALIGNMENT ( SPACERANGER_MOUSE.out.metrics
                                       .join(SPACERANGER_HUMAN.out.metrics)
                                       .join(SPACERANGER_HUMAN.out.spatial) )
        
        if ( params.do_snv_extract ) {
            GET_REFERENCE_PILEUP_MOUSE ( file("${params.mouse_reference_genome}") )
        
            GET_REFERENCE_PILEUP_HUMAN ( file("${params.human_reference_genome}") )
        
        
            GET_PILEUP_OF_BAM_MOUSE ( SPACERANGER_MOUSE.out.bam,
                                      GET_REFERENCE_PILEUP_MOUSE.out )
        
            GET_PILEUP_OF_BAM_HUMAN ( SPACERANGER_HUMAN.out.bam,
                                      GET_REFERENCE_PILEUP_HUMAN.out )
        
        
            GET_SNV_FROM_PILEUP_MOUSE ( GET_PILEUP_OF_BAM_MOUSE.out,
                                        GET_REFERENCE_PILEUP_MOUSE.out,
                                        "mouse" )
        
            GET_SNV_FROM_PILEUP_HUMAN ( GET_PILEUP_OF_BAM_HUMAN.out,
                                        GET_REFERENCE_PILEUP_HUMAN.out,
                                        "human" )
        }
        
        if ( params.do_splicing_quantification ) {
            CELLSORT_BAM_MOUSE ( SPACERANGER_MOUSE.out.bam )
        
            CELLSORT_BAM_HUMAN ( SPACERANGER_HUMAN.out.bam )
        
        
            SPLICING_QUANTIFICATION_MOUSE ( CELLSORT_BAM_MOUSE.out
                                            .join(SPACERANGER_MOUSE.out.bam)
                                            .join(SPACERANGER_MOUSE.out.mtx),
                                            file("${params.mouse_reference_genome}"),
                                            "mouse" )
       
            SPLICING_QUANTIFICATION_HUMAN ( CELLSORT_BAM_HUMAN.out
                                            .join(SPACERANGER_HUMAN.out.bam)
                                            .join(SPACERANGER_HUMAN.out.mtx),
                                            file("${params.human_reference_genome}"),
                                            "human" )
        }

        if ( params.do_merge_mtx ) {
            MERGE_MTX ( SPACERANGER_MOUSE.out.mtx
                        .join(SPACERANGER_HUMAN.out.mtx) )
        }
        else {
            RETURN_SEPARATE_MTX ( SPACERANGER_MOUSE.out.mtx
                                  .join(SPACERANGER_HUMAN.out.mtx) )
        }

    emit:
        SPACERANGER_HUMAN.out.spatial
}
