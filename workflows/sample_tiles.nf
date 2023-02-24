
include { LOAD_SAMPLE_INFO;
          EXTRACT_ROI;
          STAIN_NORMALIZATION;
          CONVERT_TO_TILED_TIFF;
          CREATE_THUMBNAIL_TIFF;
          GET_PIXEL_MASK;
          TILE_WSI;
          GET_TILE_MASK;
          SELECT_SAVE_TILES;
          GET_INCEPTION_FEATURES_TILES;
          GET_INCEPTION_FEATURES;
        } from '../modules/local/tasks'

include { GET_HOVERNET_MASK;  
          CHECK_MASK;
          INFER_HOVERNET_TILES;
          GET_NUCLEI_TYPE_COUNTS;
          INFER_HOVERNET;
          INFER_STARDIST;
          COMPRESS_JSON_FILE;
          COMPUTE_SEGMENTATION_DATA;
          GENERATE_PERSPOT_SEGMENTATION_DATA;
        } from '../modules/local/hovernet'
        
include { MERGE_IMAGING_DATA
        } from '../modules/local/merge'
        
workflow WTILES {

    take:
        dataset

    main:
        LOAD_SAMPLE_INFO ( dataset )
        
        EXTRACT_ROI ( LOAD_SAMPLE_INFO.out.main )        

        if ( params.stain_normalization ) {
            STAIN_NORMALIZATION ( EXTRACT_ROI.out )
            
            CONVERT_TO_TILED_TIFF ( STAIN_NORMALIZATION.out )
            }
        else
            CONVERT_TO_TILED_TIFF ( EXTRACT_ROI.main )
        
        
        CREATE_THUMBNAIL_TIFF ( CONVERT_TO_TILED_TIFF.out )
        
        GET_PIXEL_MASK ( CREATE_THUMBNAIL_TIFF.out )
        
        TILE_WSI ( CONVERT_TO_TILED_TIFF.out
                  .join(LOAD_SAMPLE_INFO.out.grid) )
        
        GET_TILE_MASK ( CREATE_THUMBNAIL_TIFF.out
                        .join(GET_PIXEL_MASK.out)
                        .join(TILE_WSI.out.grid) )
                        
                        
        // Tilitng sub-worflow for a small number of tiles        
        
        SELECT_SAVE_TILES ( CONVERT_TO_TILED_TIFF.out
                            .join(TILE_WSI.out.grid)
                            .join(GET_TILE_MASK.out.mask) )
        
        GET_INCEPTION_FEATURES_TILES ( SELECT_SAVE_TILES.out.tiles )
                                         
        INFER_HOVERNET_TILES ( SELECT_SAVE_TILES.out.tiles )
        
        GET_NUCLEI_TYPE_COUNTS ( INFER_HOVERNET_TILES.out.json )
        
        
        // Feature extraction for all tiles
        
        GET_INCEPTION_FEATURES ( CONVERT_TO_TILED_TIFF.out
                                 .join(GET_TILE_MASK.out.mask)
                                 .join(TILE_WSI.out.grid)
                                 .join(LOAD_SAMPLE_INFO.out.grid) )
        
        
        if ( params.hovernet_segmentation ) {
            GET_HOVERNET_MASK ( CONVERT_TO_TILED_TIFF.out )

            GET_TISSUE_MASK ( TILE_WSI.out.grid
                              .join(GET_TILE_MASK.out.mask) )

            INFER_HOVERNET ( CONVERT_TO_TILED_TIFF.out
                             .join(GET_TISSUE_MASK.out) )
            
            jsonout = INFER_HOVERNET.out.json
        }
        else {
            INFER_STARDIST ( CONVERT_TO_TILED_TIFF.out )
            
            jsonout = INFER_STARDIST.out.json
        }

        COMPRESS_JSON_FILE ( jsonout )
        
        COMPUTE_SEGMENTATION_DATA ( jsonout )

        GENERATE_PERSPOT_SEGMENTATION_DATA ( TILE_WSI.out.grid
                                         .join(COMPUTE_SEGMENTATION_DATA.out) )
        
        MERGE_IMAGING_DATA ( GET_INCEPTION_FEATURES.out
                             .join(GENERATE_PERSPOT_SEGMENTATION_DATA.out.data) )        

}
