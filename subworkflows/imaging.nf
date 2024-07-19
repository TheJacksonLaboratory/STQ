
include { LOAD_SAMPLE_INFO;
          GET_IMAGE_SIZE;
          EXTRACT_ROI;
          COLOR_NORMALIZATION;
          STAIN_NORMALIZATION;
          CONVERT_TO_TILED_TIFF;
          RESIZE_IMAGE;
          GET_THUMB;
          GET_PIXEL_MASK;
          TILE_WSI;
          GET_TILE_MASK;
          GET_TISSUE_MASK;
          SELECT_SAVE_TILES;
          GET_INCEPTION_FEATURES_TILES;
          GET_INCEPTION_FEATURES;
          GET_CTRANSPATH_FEATURES;
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
          INFER_PREP_HOVERNET;
          INFER_STARDIST;
          COMPRESS_JSON_FILE;
          COMPUTE_SEGMENTATION_DATA;
          GENERATE_PERSPOT_SEGMENTATION_DATA;
        } from '../modules/local/hovernet'

include { CONVERT_TO_PYRAMIDAL_OME;
          EXTRACT_IMAGE_METADATA;
        } from '../modules/local/ome'
        
include { CONVERT_SEGMENTATION_DATA;
          CONVERT_CSV_TO_ANNDATA;
        } from '../modules/local/merge'

include { DIMRED_CLUSTER;
          DIMRED_CLUSTER_MORPH;
        } from '../modules/local/postprocessing'
        
workflow IMG {

    take:
        samples

    main:    
        images = samples.map{[it[0], (it[1].image)]}
        
        LOAD_SAMPLE_INFO ( samples
                           .join(images) )
         
        GET_IMAGE_SIZE ( LOAD_SAMPLE_INFO.out.main )
        
        if ( params.short_workflow ) {
            GET_THUMB ( LOAD_SAMPLE_INFO.out.image )

            convertedimage = LOAD_SAMPLE_INFO.out.image
            thumbimage = GET_THUMB.out
            imagesize = GET_IMAGE_SIZE.out
        }
        else {
            if ( params.export_image_metadata ) {
                EXTRACT_IMAGE_METADATA ( LOAD_SAMPLE_INFO.out.main
                                         .join(GET_IMAGE_SIZE.out) )
            }
            
            EXTRACT_ROI ( LOAD_SAMPLE_INFO.out.main
                          .join(GET_IMAGE_SIZE.out) )
                          
            RESIZE_IMAGE ( EXTRACT_ROI.out.image )
            
            imageroi = RESIZE_IMAGE.out.full
            imagesize = RESIZE_IMAGE.out.size
    
            if ( params.stain_normalization ) {
                if ( params.macenko_normalization ) {
                    STAIN_NORMALIZATION ( imageroi
                                          .join(imagesize) )
                    
                    normimage = STAIN_NORMALIZATION.out
                    }
                else {
                    COLOR_NORMALIZATION ( imageroi
                                          .join(imagesize) )
                    
                    normimage = COLOR_NORMALIZATION.out
                    }
                
                CONVERT_TO_TILED_TIFF ( normimage )
                }
            else
                CONVERT_TO_TILED_TIFF ( imageroi )
            
            convertedimage = CONVERT_TO_TILED_TIFF.out.full
            thumbimage = CONVERT_TO_TILED_TIFF.out.thumb
            
            if ( params.export_image ) {
                CONVERT_TO_PYRAMIDAL_OME ( convertedimage )
            }
        }


        if ( params.check_focus ) {
            CHECK_FOCUS ( convertedimage
                          .join(imagesize) )
        }
        
        if ( params.do_superpixels ) {
            SUPERPIXELATION ( convertedimage
                              .join(imagesize) )
            
            if ( params.export_superpixels_contours ) {
                EXPORT_DOWN_IMAGE_FOR_CONTOURS ( convertedimage
                                  .join(imagesize) )
            
                EXPORT_SUPERPIXELATION_CONTOURS ( SUPERPIXELATION.out.main
                                                  .join(imagesize) )
                }
            }


        GET_PIXEL_MASK ( thumbimage
                         .join(imagesize) )
        
        TILE_WSI ( convertedimage
                  .join(LOAD_SAMPLE_INFO.out.grid)
                  .join(imagesize)
                  .join(LOAD_SAMPLE_INFO.out.mpp) )
        
        GET_TILE_MASK ( thumbimage
                        .join(GET_PIXEL_MASK.out)
                        .join(TILE_WSI.out.grid) 
                        .join(imagesize))    
                        
                        
        // Tilitng sub-workflow for a small number of tiles
        if ( params.sample_tiles_subworkflow ) {
            SELECT_SAVE_TILES ( convertedimage
                                .join(TILE_WSI.out.grid)
                                .join(GET_TILE_MASK.out.mask) )
                                             
            INFER_HOVERNET_TILES ( SELECT_SAVE_TILES.out.tiles )
            
            GET_NUCLEI_TYPE_COUNTS ( INFER_HOVERNET_TILES.out.json )
        }
      
        
        if ( params.extract_tile_features ) {        
            if (params.extract_inception_features) {
                GET_INCEPTION_FEATURES ( convertedimage
                                         .join(GET_TILE_MASK.out.mask)
                                         .join(TILE_WSI.out.grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                features_out = GET_INCEPTION_FEATURES.out
            }

            if (params.extract_transpath_features) {                   
                GET_CTRANSPATH_FEATURES ( convertedimage
                                         .join(GET_TILE_MASK.out.mask)
                                         .join(TILE_WSI.out.grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                if (params.extract_inception_features) {
                    features_out = features_out.concat( GET_CTRANSPATH_FEATURES.out )
                }
                else {
                    features_out = GET_CTRANSPATH_FEATURES.out
                }
            }
            
            if ( params.do_imaging_anndata ) {
                CONVERT_CSV_TO_ANNDATA ( features_out
                .filter{ it[2]== params.expansion_factor_for_clustering }
                .filter{ it[3] == params.suffix_for_clustering } )
            }
        }

        if ( params.do_nuclear_segmentation ) {
        
            GET_TISSUE_MASK ( TILE_WSI.out.grid
                      .join(GET_TILE_MASK.out.mask)
                      .join(imagesize) )
                      
            if ( params.hovernet_segmentation ) {
                INFER_PREP_HOVERNET ( convertedimage
                                   .join(GET_TISSUE_MASK.out)
                                   .join(imagesize) )

                INFER_HOVERNET ( convertedimage
                                 .join(GET_TISSUE_MASK.out)
                                 .join(imagesize)
                                 .join(INFER_PREP_HOVERNET.out) )
                
                jsonout = INFER_HOVERNET.out.json
                                 
                GET_NUCLEI_MASK_FROM_HOVERNET_JSON ( convertedimage
                                                     .join(jsonout)
                                                     .join(imagesize) )
                
                segmaskout = GET_NUCLEI_MASK_FROM_HOVERNET_JSON.out
            }
            else {
                INFER_STARDIST ( convertedimage
                                 .join(GET_TISSUE_MASK.out)
                                 .join(imagesize) )
            
                jsonout = INFER_STARDIST.out.json
                segmaskout = INFER_STARDIST.out.mask
            }

            COMPRESS_JSON_FILE ( jsonout )
        
            COMPUTE_SEGMENTATION_DATA ( jsonout
                                        .join(imagesize) )
            
            if ( params.do_superpixels ) {
                CALCULATE_CELLS_OD ( convertedimage
                                     .join(segmaskout)
                                     .join(COMPUTE_SEGMENTATION_DATA.out)
                                     .join(imagesize) )
                                     
                ASSIGN_NUCLEI_TO_SUPERPIXELS ( SUPERPIXELATION.out.main
                                               .join(CALCULATE_CELLS_OD.out)
                                               .join(imagesize) )
                }

            GENERATE_PERSPOT_SEGMENTATION_DATA ( TILE_WSI.out.grid
                                             .join(COMPUTE_SEGMENTATION_DATA.out)
                                             .join(imagesize) )

            if ( params.do_clustering ) {
                if ( params.do_imaging_anndata ) {
                    features_selected_out = CONVERT_CSV_TO_ANNDATA.out
                    .filter{ it[2]== params.expansion_factor_for_clustering }
                    .filter{ it[3] == params.suffix_for_clustering }

                    DIMRED_CLUSTER_MORPH ( TILE_WSI.out.grid
                                         .join(thumbimage)
                                         .join(GENERATE_PERSPOT_SEGMENTATION_DATA.out.data)
                                         .join(features_selected_out) )
                }
            }
        }
        else {
            if ( params.do_clustering ) {
                if ( params.do_imaging_anndata ) {
                    features_selected_out = CONVERT_CSV_TO_ANNDATA.out
                    .filter{ it[2]== params.expansion_factor_for_clustering }
                    .filter{ it[3] == params.suffix_for_clustering }

                    DIMRED_CLUSTER ( TILE_WSI.out.grid
                                     .join(thumbimage)
                                     .join(features_selected_out) )
                }
            }
        }
}
