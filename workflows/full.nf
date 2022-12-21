
include { LOAD_SAMPLE_INFO;
          STAIN_NORMALIZATION;
          CONVERT_TO_TILED_TIFF;
          CREATE_THUMBNAIL_TIFF;
          GET_PIXEL_MASK;
          TILE_WSI;
          GET_TILE_MASK;
          GET_INCEPTION_FEATURES;
        } from '../modules/local/tasks'

include { GET_HOVERNET_MASK;  
          CHECK_MASK;
          INFER_HOVERNET;
          INFER_STARDIST;
          COMPRESS_JSON_FILE;
          COMPUTE_SEGMENTATION_DATA;
          GENERATE_PERSPOT_SEGMENTATION_DATA;
        } from '../modules/local/hovernet'
        
include { MERGE_IMAGING_DATA
        } from '../modules/local/merge'
        
workflow MAIN {

    take:
        dataset

    main:
        LOAD_SAMPLE_INFO ( dataset )
        
        if ( params.stain_normalization ) {
            STAIN_NORMALIZATION ( LOAD_SAMPLE_INFO.out.main )
            
            CONVERT_TO_TILED_TIFF ( STAIN_NORMALIZATION.out )
            }
        else
            CONVERT_TO_TILED_TIFF ( LOAD_SAMPLE_INFO.out.main )
        
        
        CREATE_THUMBNAIL_TIFF ( CONVERT_TO_TILED_TIFF.out )
        
        GET_PIXEL_MASK ( CREATE_THUMBNAIL_TIFF.out )
        
        TILE_WSI ( CONVERT_TO_TILED_TIFF.out
                  .join(LOAD_SAMPLE_INFO.out.grid) )
        
        GET_TILE_MASK ( CREATE_THUMBNAIL_TIFF.out
                        .join(GET_PIXEL_MASK.out)
                        .join(TILE_WSI.out.grid) )
        
        GET_INCEPTION_FEATURES ( CONVERT_TO_TILED_TIFF.out
                                 .join(GET_TILE_MASK.out.mask)
                                 .join(TILE_WSI.out.grid)
                                 .join(LOAD_SAMPLE_INFO.out.grid) )
        
        
        if ( params.hovernet_segmentation ) {
            GET_HOVERNET_MASK ( CONVERT_TO_TILED_TIFF.out )
        
            CHECK_MASK ( GET_HOVERNET_MASK.out
                         .join(CREATE_THUMBNAIL_TIFF.out) )
                                     
            INFER_HOVERNET ( CONVERT_TO_TILED_TIFF.out
                             .join(CHECK_MASK.out) )
                             
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
