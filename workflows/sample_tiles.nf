
include { LOAD_SAMPLE_INFO;
          STAIN_NORMALIZATION;
          CONVERT_TO_TILED_TIFF;
          CREATE_THUMBNAIL_TIFF;
          GET_PIXEL_MASK;
          TILE_WSI;
          GET_TILE_MASK;
          SELECT_SAVE_TILES;
          GET_INCEPTION_FEATURES_TILES;
        } from '../modules/local/tasks'

include { GET_HOVERNET_MASK;  
          CHECK_MASK;
          INFER_HOVERNET_TILES;
          GET_NUCLEI_TYPE_COUNTS;
        } from '../modules/local/hovernet'
        
workflow SAMPLE_TILES {

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
        
        SELECT_SAVE_TILES ( CONVERT_TO_TILED_TIFF.out
                            .join(TILE_WSI.out.grid)
                            .join(GET_TILE_MASK.out.mask) )
        
        GET_INCEPTION_FEATURES_TILES ( SELECT_SAVE_TILES.out.tiles )
                                         
        INFER_HOVERNET_TILES ( SELECT_SAVE_TILES.out.tiles )
        
        GET_NUCLEI_TYPE_COUNTS ( INFER_HOVERNET_TILES.out.json )

}
