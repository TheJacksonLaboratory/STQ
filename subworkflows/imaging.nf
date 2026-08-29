
include { LOAD_SAMPLE_INFO;
          GET_IMAGE_SIZE;
          EXTRACT_ROI;
          STAIN_NORMALIZATION;
          CONVERT_TO_TILED_TIFF;
          RESIZE_IMAGE;
          GET_THUMB;
          MAKE_TINY_THUMB;
          GET_PIXEL_MASK;
          TILE_WSI;
          GET_TILE_MASK;
          GET_TISSUE_MASK;
          SELECT_SAVE_TILES;
          SELECT_SAVE_TILES_RAW;
          ASSEMBLE_TILES_CELLS;
          GET_INCEPTION_FEATURES_TILES;
        } from '../modules/local/tasks'

include { GET_UNI_FEATURES;
          GET_UNI2_FEATURES;
          GET_HOPTIMUS0_FEATURES;
          GET_CTRANSPATH_FEATURES;
          GET_MOCOV3_FEATURES;
          GET_CONCH_FEATURES;
          GET_INCEPTION_FEATURES;
          GET_VIRCHOW_FEATURES;
          GET_VIRCHOW2_FEATURES;
          GET_PHIKON_FEATURES;
          GET_PHIKON2_FEATURES;
          GET_HOPTIMUS1_FEATURES;
          GET_H0MINI_FEATURES;
          GET_MIDNIGHT12K_FEATURES;
          GET_GIGAPATH_FEATURES;
          GET_GIGAPATHF_FEATURES;
        } from '../modules/local/features'

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
          INFER_HOVERNET_TILES_RAW;
          GET_NUCLEI_TYPE_COUNTS;
          INFER_HOVERNET;
          INFER_PREP_HOVERNET;
          INFER_STARDIST;
          COMPRESS_JSON_FILE;
          COMPUTE_SEGMENTATION_DATA;
          GENERATE_PERSPOT_SEGMENTATION_DATA;
        } from '../modules/local/hovernet'

include { TIFFFILE_OMETIFF;
          EXTRACT_IMAGE_METADATA;
        } from '../modules/local/ome'
        
include { CONVERT_SEGMENTATION_DATA;
          CONVERT_CSV_TO_ANNDATA;
        } from '../modules/local/merge'

include { CREATE_CELLSZARRZIP;
        } from '../modules/local/zarr'

include { DIMRED_CLUSTER;
          DIMRED_CLUSTER_MORPH;
          MAKE_CLUSTER_VECTORS;
          MAKE_SAMPLER_VECTOR;
        } from '../modules/local/postprocessing'
        
workflow IMG {

    take:
        samples

    main:    
        images = samples.map{[it[0], (it[1].image)]}
        
        LOAD_SAMPLE_INFO ( samples
                           .join(images) )

        mpp = LOAD_SAMPLE_INFO.out.mpp
         
        GET_IMAGE_SIZE ( LOAD_SAMPLE_INFO.out.main )
        
        if ( params.short_workflow ) {
            if ( params.reuse_previous_run ) {
                thumbimage = LOAD_SAMPLE_INFO.out.image.map { sample_id, image ->
                    tuple(sample_id, image.getParent().resolve("thumbnail.tiff"))
                }
            }
            else {
                GET_THUMB ( LOAD_SAMPLE_INFO.out.image )
                thumbimage = GET_THUMB.out
            }

            convertedimage = LOAD_SAMPLE_INFO.out.image
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
                SELECT_SAVE_TILES_RAW ( imageroi )
                                                
                INFER_HOVERNET_TILES_RAW ( SELECT_SAVE_TILES_RAW.out.tiles )

                ASSEMBLE_TILES_CELLS ( SELECT_SAVE_TILES_RAW.out.tiles
                                        .join(INFER_HOVERNET_TILES_RAW.out.json) )

                STAIN_NORMALIZATION ( imageroi
                                        .join(ASSEMBLE_TILES_CELLS.out)
                                        .join(imagesize) )
                
                normimage = STAIN_NORMALIZATION.out

                CONVERT_TO_TILED_TIFF ( normimage )
                }
            else
                CONVERT_TO_TILED_TIFF ( imageroi )
            
            convertedimage = CONVERT_TO_TILED_TIFF.out.full
            thumbimage = CONVERT_TO_TILED_TIFF.out.thumb
            
            if ( params.export_image ) {
                TIFFFILE_OMETIFF ( convertedimage
                                   .join(imagesize) )
            }
        }
        
        if ( !params.reuse_previous_run ) {
            MAKE_TINY_THUMB ( thumbimage )
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


        if ( params.reuse_previous_run ) {
            grid = LOAD_SAMPLE_INFO.out.image.map { sample_id, image ->
                tuple(sample_id, image.getParent().resolve("grid/grid.csv"), image.getParent().resolve("grid/grid.json"))
            }
            mask = LOAD_SAMPLE_INFO.out.image.map { sample_id, image ->
                tuple(sample_id, image.getParent().resolve("mask/tile_mask.csv"))
            }
        }
        else {
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

            grid = TILE_WSI.out.grid
            mask = GET_TILE_MASK.out.mask
        }
                        
                        
        // Tilitng sub-workflow for a small number of tiles
        if ( params.sample_tiles_subworkflow ) {
            SELECT_SAVE_TILES ( convertedimage
                                .join(grid)
                                .join(mask) )
                                             
            INFER_HOVERNET_TILES ( SELECT_SAVE_TILES.out.tiles )
            
            GET_NUCLEI_TYPE_COUNTS ( INFER_HOVERNET_TILES.out.json )
        }
      
        
        if ( params.extract_tile_features ) {

            if (params.extract_transpath_features) {                   
                GET_CTRANSPATH_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                ctranspath_features_out = GET_CTRANSPATH_FEATURES.out
            }

            if (params.extract_mocov3_features) {                   
                GET_MOCOV3_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                mocov3_features_out = GET_MOCOV3_FEATURES.out
            }

            if (params.extract_inception_features) {
                GET_INCEPTION_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                inception_features_out = GET_INCEPTION_FEATURES.out
            }

            if (params.extract_uni_features) {                   
                GET_UNI_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                uni_features_out = GET_UNI_FEATURES.out
            }

            if (params.extract_uni2_features) {                   
                GET_UNI2_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                uni2_features_out = GET_UNI2_FEATURES.out
            }

            if (params.extract_conch_features) {                   
                GET_CONCH_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                conch_features_out = GET_CONCH_FEATURES.out
            }           

            if (params.extract_hoptimus0_features) {                   
                GET_HOPTIMUS0_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                hoptimus0_features_out = GET_HOPTIMUS0_FEATURES.out
            }

            if (params.extract_virchow_features) {                   
                GET_VIRCHOW_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                virchow_features_out = GET_VIRCHOW_FEATURES.out
            }
            
            if (params.extract_virchow2_features) {                   
                GET_VIRCHOW2_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                virchow2_features_out = GET_VIRCHOW2_FEATURES.out
            }

            if (params.extract_phikon_features) {                   
                GET_PHIKON_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                phikon_features_out = GET_PHIKON_FEATURES.out
            }

            if (params.extract_phikon2_features) {                   
                GET_PHIKON2_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                phikon2_features_out = GET_PHIKON2_FEATURES.out
            }

            if (params.extract_hoptimus1_features) {                   
                GET_HOPTIMUS1_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                hoptimus1_features_out = GET_HOPTIMUS1_FEATURES.out
            }

            if (params.extract_h0mini_features) {                   
                GET_H0MINI_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                h0mini_features_out = GET_H0MINI_FEATURES.out
            }

            if (params.extract_midnight12k_features) {                   
                GET_MIDNIGHT12K_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                midnight12k_features_out = GET_MIDNIGHT12K_FEATURES.out
            }

            if (params.extract_gigapath_features) {                   
                GET_GIGAPATH_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                gigapath_features_out = GET_GIGAPATH_FEATURES.out
            }

            if (params.extract_gigapathf_features) {                   
                GET_GIGAPATHF_FEATURES ( convertedimage
                                         .join(mask)
                                         .join(grid)
                                         .join(LOAD_SAMPLE_INFO.out.grid)
                                         .join(imagesize)
                                         .combine(Channel.fromList(params.expansion_factor)) )
                
                gigapathf_features_out = GET_GIGAPATHF_FEATURES.out
            }


            // Default features
            features_out = channel.empty()

            if (params.extract_transpath_features) {
                features_out = features_out.concat( ctranspath_features_out )
            }
            if (params.extract_mocov3_features) {
                features_out = features_out.concat( mocov3_features_out )
            }
            if (params.extract_inception_features) {
                features_out = features_out.concat( inception_features_out )
            }
            if (params.extract_uni_features) {
                features_out = features_out.concat( uni_features_out )
            }
            if (params.extract_uni2_features) {
                features_out = features_out.concat( uni2_features_out )
            }
            if (params.extract_conch_features) {
                features_out = features_out.concat( conch_features_out )
            }
            if (params.extract_hoptimus0_features) {
                features_out = features_out.concat( hoptimus0_features_out )
            }
            if (params.extract_virchow_features) {
                features_out = features_out.concat( virchow_features_out )
            }
            if (params.extract_virchow2_features) {
                features_out = features_out.concat( virchow2_features_out )
            }
            if (params.extract_phikon_features) {
                features_out = features_out.concat( phikon_features_out )
            }
            if (params.extract_phikon2_features) {
                features_out = features_out.concat( phikon2_features_out )
            }
            if (params.extract_hoptimus1_features) {
                features_out = features_out.concat( hoptimus1_features_out )
            }
            if (params.extract_h0mini_features) {
                features_out = features_out.concat( h0mini_features_out )
            }
            if (params.extract_midnight12k_features) {
                features_out = features_out.concat( midnight12k_features_out )
            }
            if (params.extract_gigapath_features) {
                features_out = features_out.concat( gigapath_features_out )
            }
            if (params.extract_gigapathf_features) {
                features_out = features_out.concat( gigapathf_features_out )
            }

            MAKE_SAMPLER_VECTOR ( samples
                                 .combine(features_out, by: 0) )

            if ( params.do_clustering ) {
                if ( params.do_imaging_anndata ) {
                    CONVERT_CSV_TO_ANNDATA ( features_out
                    .filter{ it[2]== params.expansion_factor_for_clustering }
                    .filter{ it[3] == params.suffix_for_clustering } )
                }
            }
        }

        if ( params.do_nuclear_segmentation ) {
        
            GET_TISSUE_MASK ( grid
                      .join(mask)
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

            CREATE_CELLSZARRZIP ( mpp
                                .join(segmaskout) )

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

            GENERATE_PERSPOT_SEGMENTATION_DATA ( grid
                                             .join(COMPUTE_SEGMENTATION_DATA.out)
                                             .join(imagesize) )

            if ( params.do_clustering ) {
                if ( params.do_imaging_anndata ) {
                    features_selected_out = CONVERT_CSV_TO_ANNDATA.out
                    .filter{ it[2]== params.expansion_factor_for_clustering }
                    .filter{ it[3] == params.suffix_for_clustering }

                    DIMRED_CLUSTER_MORPH ( grid
                                         .join(thumbimage)
                                         .join(GENERATE_PERSPOT_SEGMENTATION_DATA.out.data)
                                         .join(features_selected_out) )

                    MAKE_CLUSTER_VECTORS ( thumbimage
                                         .join(GENERATE_PERSPOT_SEGMENTATION_DATA.out.data)
                                         .join(features_selected_out)
                                         .join(DIMRED_CLUSTER_MORPH.out.clusters) )
                }
            }
        }
        else {
            if ( params.do_clustering ) {
                if ( params.do_imaging_anndata ) {
                    features_selected_out = CONVERT_CSV_TO_ANNDATA.out
                    .filter{ it[2]== params.expansion_factor_for_clustering }
                    .filter{ it[3] == params.suffix_for_clustering }

                    DIMRED_CLUSTER ( grid
                                     .join(thumbimage)
                                     .join(features_selected_out) )
                }
            }
        }
}
