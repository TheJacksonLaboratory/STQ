
include { LOAD_SAMPLE_INFO;
          GET_IMAGE_SIZE;
          EXTRACT_ROI;
          COLOR_NORMALIZATION;
          STAIN_NORMALIZATION;
          CONVERT_TO_TILED_TIFF;
          GET_PIXEL_MASK;
          TILE_WSI;
          GET_TILE_MASK;
          GET_TISSUE_MASK;
          SELECT_SAVE_TILES;
          GET_INCEPTION_FEATURES_TILES;
          GET_INCEPTION_FEATURES;
        } from '../modules/local/tasks'

include { CHECK_FOCUS;
        } from '../modules/local/focus'
                
include { SUPERPIXELATION;
          EXPORT_DOWN_IMAGE_FOR_CONTOURS;
          CALCULATE_CELLS_OD;
          ASSIGN_NUCLEI_TO_SUPERPIXELS;
          EXPORT_SUPERPIXELATION_CONTOURS;
        } from '../modules/local/superpixel'

include { GET_NUCLEI_MASK_FROM_HOVERNET_JSON;  
          INFER_HOVERNET_TILES;
          GET_NUCLEI_TYPE_COUNTS;
          INFER_HOVERNET;
          INFER_STARDIST;
          COMPRESS_JSON_FILE;
          COMPUTE_SEGMENTATION_DATA;
          GENERATE_PERSPOT_SEGMENTATION_DATA;
        } from '../modules/local/hovernet'
        
include { MERGE_IMAGING_DATA;
          CONVERT_CSV_TO_ANNDATA;
        } from '../modules/local/merge'
        
workflow IMG {

    take:
        samples

    main:
        //if (params.hovernet_segmentation && params.do_superpixels) {
        //    error "NotImplementedError: HoVer-Net nuclear segmentation and superpixels subworkflow are incompatible."
        //}
    
        images = samples.map{[it[0], (it[1].image)]}
        
        LOAD_SAMPLE_INFO ( samples
                           .join(images) )
         
        GET_IMAGE_SIZE ( LOAD_SAMPLE_INFO.out.main )
        
        EXTRACT_ROI ( LOAD_SAMPLE_INFO.out.main
                      .join(GET_IMAGE_SIZE.out) )

        if ( params.stain_normalization ) {
            if ( params.macenko_normalization ) {
                STAIN_NORMALIZATION ( EXTRACT_ROI.out.image
                                      .join(GET_IMAGE_SIZE.out) )
                
                normimage = STAIN_NORMALIZATION.out
                }
            else {
                COLOR_NORMALIZATION ( EXTRACT_ROI.out.image
                                      .join(GET_IMAGE_SIZE.out) )
                
                normimage = COLOR_NORMALIZATION.out
                }
            
            CONVERT_TO_TILED_TIFF ( normimage )
            }
        else
            CONVERT_TO_TILED_TIFF ( EXTRACT_ROI.out.image )
        
        if ( params.check_focus ) {
            CHECK_FOCUS ( CONVERT_TO_TILED_TIFF.out.full
                          .join(CONVERT_TO_TILED_TIFF.out.size) )
        }
        
        if ( params.do_superpixels ) {
            SUPERPIXELATION ( CONVERT_TO_TILED_TIFF.out.full
                              .join(CONVERT_TO_TILED_TIFF.out.size) )
            
            if ( params.export_superpixels_contours ) {
                EXPORT_DOWN_IMAGE_FOR_CONTOURS ( CONVERT_TO_TILED_TIFF.out.full
                                  .join(CONVERT_TO_TILED_TIFF.out.size) )
            
                EXPORT_SUPERPIXELATION_CONTOURS ( SUPERPIXELATION.out.main
                                                  .join(CONVERT_TO_TILED_TIFF.out.size) )
                }
            }


        GET_PIXEL_MASK ( CONVERT_TO_TILED_TIFF.out.thumb
                         .join(CONVERT_TO_TILED_TIFF.out.size) )
        
        TILE_WSI ( CONVERT_TO_TILED_TIFF.out.full
                  .join(LOAD_SAMPLE_INFO.out.grid)
                  .join(CONVERT_TO_TILED_TIFF.out.size) )
        
        GET_TILE_MASK ( CONVERT_TO_TILED_TIFF.out.thumb
                        .join(GET_PIXEL_MASK.out)
                        .join(TILE_WSI.out.grid) 
                        .join(CONVERT_TO_TILED_TIFF.out.size))       
        
        
        if ( params.extract_tile_features ) {
            GET_INCEPTION_FEATURES ( CONVERT_TO_TILED_TIFF.out.full
                                     .join(GET_TILE_MASK.out.mask)
                                     .join(TILE_WSI.out.grid)
                                     .join(LOAD_SAMPLE_INFO.out.grid)
                                     .join(CONVERT_TO_TILED_TIFF.out.size) )
        }      

        
        if ( params.do_nuclear_sementation ) {
        
            GET_TISSUE_MASK ( TILE_WSI.out.grid
                      .join(GET_TILE_MASK.out.mask)
                      .join(CONVERT_TO_TILED_TIFF.out.size) )
                      
            if ( params.hovernet_segmentation ) {
                INFER_HOVERNET ( CONVERT_TO_TILED_TIFF.out.full
                                 .join(GET_TISSUE_MASK.out)
                                 .join(CONVERT_TO_TILED_TIFF.out.size) )
                
                jsonout = INFER_HOVERNET.out.json
                                 
                GET_NUCLEI_MASK_FROM_HOVERNET_JSON ( CONVERT_TO_TILED_TIFF.out.full
                                                     .join(jsonout)
                                                     .join(CONVERT_TO_TILED_TIFF.out.size) )
                
                segmaskout = GET_NUCLEI_MASK_FROM_HOVERNET_JSON.out
            }
            else {
                INFER_STARDIST ( CONVERT_TO_TILED_TIFF.out.full
                                 .join(GET_TISSUE_MASK.out)
                                 .join(CONVERT_TO_TILED_TIFF.out.size) )
            
                jsonout = INFER_STARDIST.out.json
                segmaskout = INFER_STARDIST.out.mask
            }

            COMPRESS_JSON_FILE ( jsonout )
        
            COMPUTE_SEGMENTATION_DATA ( jsonout
                                        .join(CONVERT_TO_TILED_TIFF.out.size) )
            
            if ( params.do_superpixels ) {
                CALCULATE_CELLS_OD ( CONVERT_TO_TILED_TIFF.out.full
                                     .join(segmaskout)
                                     .join(COMPUTE_SEGMENTATION_DATA.out)
                                     .join(CONVERT_TO_TILED_TIFF.out.size) )
                                     
                ASSIGN_NUCLEI_TO_SUPERPIXELS ( SUPERPIXELATION.out.main
                                               .join(CALCULATE_CELLS_OD.out)
                                               .join(CONVERT_TO_TILED_TIFF.out.size) )
                }

            GENERATE_PERSPOT_SEGMENTATION_DATA ( TILE_WSI.out.grid
                                             .join(COMPUTE_SEGMENTATION_DATA.out)
                                             .join(CONVERT_TO_TILED_TIFF.out.size) )

            if ( params.extract_tile_features ) {
                MERGE_IMAGING_DATA ( GET_INCEPTION_FEATURES.out
                                     .join(GENERATE_PERSPOT_SEGMENTATION_DATA.out.data)
                                     .join(CONVERT_TO_TILED_TIFF.out.size) )

                if ( params.do_imaging_anndata ) {
                    CONVERT_CSV_TO_ANNDATA ( MERGE_IMAGING_DATA.out )
                }
            }
        }
        else {
            if ( params.extract_tile_features ) {
                if ( params.do_imaging_anndata ) {
                    CONVERT_CSV_TO_ANNDATA ( GET_INCEPTION_FEATURES.out )
                }
            }        
        }
}
