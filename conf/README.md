
# Description of the pipeline parameters

To find out the default value for each parameter, see `conf/analysis.config`.

### Sequencing analysis parameters

+ **`do_merge_mtx`** Merge graft and host MTX (gene by spot) matrices into one MTX matrix

+ **`do_splicing_quantification`** Run splicing quantification with velocyto. The pipeline also sorts by cell barcodes the BAM file produced by Space Ranger.

+ **`do_snv_extract`** Run the BAF extraction sub-workflow to get bulk-level SNV.

+ **`reference_genome`** Path to the reference genome to use for Space Ranger reads alignment in one-reference analysis route. See https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build for Space Ranger requirements of the reference genomes.

+ **`mouse_reference_genome`** Path to the mouse reference genome for Space Ranger reads alignment in two-reference analysis route.

+ **`human_reference_genome`** Path to the human reference genome for Space Ranger reads alignment in two-reference analysis route.

+ **`deconvolution_reference_graft`** Path to a graft (e.g., human) reference genome (e.g., *.fa, *.fna, *.fa.gz, *.fna.gz) to build xenome or xengsort indices. If the indices supplied in `nextflow.config` already exits, then this parameter is ignored.

+ **`deconvolution_reference_host`** Path to a host (e.g., mouse) reference genome (e.g., *.fa, *.fna, *.fa.gz, *.fna.gz) to build xenome or xengsort indices. If the indices supplied in `nextflow.config` already exits, then this parameter is ignored.

+ **`deconvolution_kmer_size`** K-mer size for building xenome or xengsort indices. See https://github.com/data61/gossamer/blob/master/docs/xenome.md for a detailed description.

+ **`deconvolution_indices_path`** Path to save deconvolution indices.

+ **`deconvolution_indices_name`** Name of the indices.

+ **`xengsort_n`** Xengsort-specific parameter. See https://gitlab.com/genomeinformatics/xengsort for details.


##### See https://github.com/akdess/BAFExtract for the description of the following filtering parameters:

+ **`bafextract_minimum_mapping_quality`** 

+ **`bafextract_minimum_base_quality`** 

+ **`bafextract_min_coverage_per_SNV`** 

+ **`bafextract_min_MAF_covg_per_SNV`** 

+ **`bafextract_min_MAF`** 


### Imaging analysis parameters

+ **`do_img_subworkflow`** Run the imaging sub-workflow to generate imaging and nuclear morphometric features for each spot on the grid.

+ **`short_workflow`** Run short imaging workflow instead of the full imaging workflow. See config for details.

+ **`do_imaging_anndata`** Create an AnnData object (e.g., for use with Scanpy) from the *.csv.gz data file with imaging and nuclear morphometric features

+ **`do_nuclear_sementation`** Perform nuclear segmentation (use either HoVer-Net or StarDist to segment nuclei) of the entire WSI.

+ **`target_mpp`** desired image resolution for scaling the images. Note that specific DL and ML models require full-resolution images, and the supplied pre-trained models are designed for images with a resolution of around 0.25 (mpp). In case a low-magnification image is supplied (e.g., mpp is 0.5) while target_mpp is 0.25, the image is upsampled and will have doubled dimensions.

+ **`tiled_tiff_tile_size`** The TIFF WSI is internally stored in blocks (for memory management). The tile size determines the block size. This parameter is not the size of tiles used for feature extraction or segmentation aggregation. The grid parameter `grid_spot_diamter` (in micrometers) and resolution parameter `target_mpp` define the scaled image tile size.

+ **`thumbnail_downsample_factor`** A factor used to reduce the WSI dimensions to create a low-resolution slide representation.

+ **`check_focus`** Run DeepFocus module to assess focus (blurryness) of the whole slide image.

+ **`deepfocus_model_path`** Path to DeepFocus checkpoint to use.



+ **`stain_normalization`** Whether to do any stain or color normalization.

+ **`stainnet`** Path to checkpoint for stain normalization model.

+ **`macenko_normalization`** If true, then use Macenko stain normalization. If false, use StainNet color normalization. This parameter is ignored if `stain_normalization` is false.

+ **`stain_reference_image`** Reference image (or a small patch, e.g., 2000 by 2000 pixels) to use with Macenko stain normalization.

+ **`stain_patch_size`** Macenco stain normalization patch size.


+ **`mask_background_cutoff`** Parameter for detecting image background with HoVer-Net.

+ **`pixel_mask_threshold_low`** Parameter for detecting tissue pixels on the low-resolution image.

+ **`pixel_mask_threshold_high`** Parameter for detecting tissue pixels on the low-resolution image.

+ **`fraction_for_mask`** Fraction of pixels in tissue required to call tile in tissue.


+ **`use_provided_grid`** Whether to use the grid provided in the input sample sheet. If false and no Space Ranger alignment is done, then a new grid of tiles is generated based on the grid parameters.

+ **`grid_type`** Type of the grid of tiles to generate. it can be hex, square, or random.

+ **`grid_spot_diamter`** Diameter of the spot (dimension of a tile) in micrometers.

+ **`grid_spot_horizontal_spacing`** Horizontal center-to-center distance between adjacent spots (or tiles).

+ **`grid_aspect_correction`** Factor to correct Visium slide aspect ratio.


+ **`overlap_scale_factor`** Imaging features extraction parameter. If the factor is 1, then features are extracted from the tile of the ST spot dimension.
 

+ **`hovernet_segmentation`** Do HoVer-Net segmetation. If false do StarDist segmentation. 

+ **`nuclei_segmentation_dir`** name of directory to save segmentation.

+ **`hovernet_batch_size`** Parameter of HoVer-Net segmetation. This parameter is ignored when segmentation is done with StarDist.

+ **`hovernet_num_inference_workers`** Parameter of HoVer-Net segmetation. This parameter is ignored when segmentation is done with StarDist.

+ **`hovernet_chunk_size`** Parameter of HoVer-Net segmetation. This parameter is ignored when segmentation is done with StarDist.

+ **`hovernet_tile_size`** Parameter of HoVer-Net segmetation. This parameter is ignored when segmentation is done with StarDist.

+ **`stardist_model`** Path to checkpoint of stardist model.

+ **`stardist_block_size`** Size of the image block to run segmentation. Blocks are merged internally at the end of segmentation.

+ **`stardist_expand_size`** Size of cytoplasm arouhd nucleus in pixels.



+ **`hovernet_spot_assignment_factor`** Used for either HoVer-Net or StarDist segmentation postprocessing. Scaling factor of the boundary limiting the inclusion of nuclei to an ST spot. A value of 1 means the boundary size equals ST spot size.

+ **`hovernet_spot_assignment_shape`** Used for either HoVer-Net or StarDist segmentation postprocessing. The shape of the boundary, either square or disk.

+ **`hovernet_min_cell_type_prob`** Used for either HoVer-Net or StarDist segmentation postprocessing. This filtering parameteris used to remove nuclei assigned with low confidence.


+ **`extract_tile_features`** Extract (generate) imaging features for all tiles.

+ **`extract_inception_features`** If `extract_tile_features` then do Inception V3 features.

+ **`extract_transpath_features`** If `extract_tile_features` then do TransPath features.

+ **`extract_uni_features`** If `extract_tile_features` then do UNI features.

+ **`extract_conch_features`** If `extract_tile_features` then do CONCH features.

+ **`transpath_features_model`** One of 'CTransPath' or 'MoCoV3'.

+ **`use_conch_normalizer`** Use specialized CONCH normalizer, instead of the standard normalizer used with UNI and CTransPath.

+ **`uni_model_checkpoint`** Path to downloaded CONCH checkpoint. Download requires registration https://huggingface.co/MahmoodLab/UNI/blob/main/pytorch_model.bin.

+ **`conch_model_checkpoint`** Path to downloaded CONCH checkpoint. Download requires registration https://huggingface.co/MahmoodLab/CONCH/blob/main/pytorch_model.bin.



+ **`do_superpixels`** Do superpixel segmentation using SNIC algorithm.

+ **`export_superpixels_contours`** If true, export superpixel contours in JSON format.

+ **`superpixel_compactness`** Superpixel compactness parameter, see details of SNIC algorithm.

+ **`pixels_per_segment`** Number of pixels per superpixel segment, i.e., superpixel size.

+ **`superpixel_patch_size`** Superpixel patch size. Warning: patches boundaries are kept flat.

+ **`superpixel_downsampling_factor`** Superpixel downsampling factor for the input image downsampling .

+ **`od_block_size`** Block size for OD calculation.

+ **`expand_nuclei_distance`** Distance in pixels to expand the nuclei mask.




+ **`export_image`** Export the resized and normalized image in OME-TIFF format.

+ **`export_image_metadata`** Export input image metadata in OME-XML format.

+ **`compression`** Compression library to use with OME-TIFF, e.g., 'LZW'.


+ **`downsample_expanded_tile`** Downsample expanded tile.

+ **`expansion_factor`** Tile is read from expanded area around the tile center.

+ **`subtiling`** If true, split tile into subtiles, then extract features and compute average across the subtiles.

+ **`subcoords_factor`** Factor that defines the size of subtiles.

+ **`subcoords_list`** Centers of the subtiles within a tile.



+ **`do_clustering`** Do dimensionality reduction and clustering. Generate spatial and UMAP plots of imaging feature clusters as well as nucler morphometric features and classification results.

+ **`expansion_factor_for_clustering`** Features of the specified expansion factor are used for clustering.

+ **`suffix_for_clustering`** Features of this type are used for clustering.

+ **`plot_dpi`** DPI (dots per inch) of the figures.



+ **`hovernet_device_mode`** GPU or CPU device for use with HoVer-Net.

+ **`ctranspath_device_mode`** GPU or CPU device for use with TransPath inference models.



+ **`sample_tiles_subworkflow`** Run a subworkflow where a small number of tiles is saved, along with the HoVer-Net classification data.

+ **`tiles_per_slide`** Number of randomly selected tiles to use in the sampling tiles subworkflow.



+ **`do_segmentation_anndata`**   DEPRECATED parameter, will be removed in future.
