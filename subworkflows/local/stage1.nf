
include { LOAD_SAMPLE_INFO;
          CONVERT_TO_TILED_TIFF;
          CREATE_THUMBNAIL_TIFF;
          GET_PIXEL_MASK;
          TILE_WSI;
          GET_TILE_MASK;
          GET_INCEPTION_FEATURES;
        } from '../../modules/local/tasks'

workflow STAGE {

    take:
        dataset

    main:
        LOAD_SAMPLE_INFO ( dataset )
    
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

}