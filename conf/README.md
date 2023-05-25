
# Description of the pipeline parameters

To find out the default value for each parameter, see `conf/analysis.config`.

### Sequencing analysis parameters

+ **`do_merge_mtx`** Merge graft and host MTX (gene by spot) matrices into one MTX matrix

+ **`do_splicing_quantification`** Run splicing quantification with velocyto. Also sorts by cell barcodes the BAM file produced by Space Ranger.

+ **`do_snv_extract`** Run the BAF extraction sub-workflow to get bulk level SNV.

+ **`reference_genome`** Path to the reference genome to use for Space Ranger reads alignment in one-reference analysis route. See https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build for Space Ranger requirements of the reference genomes.

+ **`mouse_reference_genome`** Path to the mouse reference genome to use for Space Ranger reads alignment in two-reference analysis route.

+ **`human_reference_genome`** Path to the human reference genome to use for Space Ranger reads alignment in two-reference analysis route.

+ **`xenome_reference_graft`** Path to a graft (e.g., human) reference genome (e.g., *.fa, *.fna, *.fa.gz, *.fna.gz) to build xenome indices. If the indices supplied in `nextflow.config` already exits then this parameter is ignored.

+ **`xenome_reference_host`** Path to a host (e.g., mouse) reference genome (e.g., *.fa, *.fna, *.fa.gz, *.fna.gz) to build xenome indices. If the indices supplied in `nextflow.config` already exits then this parameter is ignored.

+ **`xenome_kmer_size`** K-mer size for building xenome indices. See https://github.com/data61/gossamer/blob/master/docs/xenome.md for detailed description.


##### See https://github.com/akdess/BAFExtract for the description of following filtering parameters:

+ **`bafextract_minimum_mapping_quality`** 

+ **`bafextract_minimum_base_quality`** 

+ **`bafextract_min_coverage_per_SNV`** 

+ **`bafextract_min_MAF_covg_per_SNV`** 

+ **`bafextract_min_MAF`** 


### Imaging analysis parameters

+ **`do_img_subworkflow`** Run the imaging sub-workflow to generate imaging and nuclear morphometric features for each spot on the grid.

+ **`do_imaging_anndata`** Create AnnData object (e.g., for use with Scanpy) from the *.csv.gz data file with imaging and nuclear morphometric features

+ **`do_nuclear_sementation`** Perform nuclear segmentation (use either HoVer-Net or StarDist to segment nuclei) of the entire WSI.

+ **`target_mpp`** desired image resolution for scaling the images. Note specific DL and ML model require full-resolution image, and the supplied pre-trained models are designed for images with resolution around 0.25 (mpp). In case a low-magnification image is supplied (e.g. mpp is 0.5) while target_mpp is 0.25, the image is upsampled and will have doubled dimensions.

+ **`tiled_tiff_tile_size`** The TIFF WSI is internally stored in blocks (for memory management). The tile size determines the block size. This is not the size of tiles used for feature extraction or segmentation aggregation. The grid parameter `grid_spot_diamter` (in mircometers) and resolution parameter `target_mpp` define the scaled image tile size.

+ **`thumbnail_downsample_factor`** A factor used to reduce the WSI dimensions to create a low resolution slide representation.

+ **`stain_normalization`** Whether to do any stain or color normalization.

+ **`macenko_normalization`** If true then use Macenko stain normalization. If false use StainNet color normalization. This parameter is ignored if `stain_normalization` if false.

+ **`stain_reference_image`** Reference image (or a small patch, e.g., 2000 by 2000 pixels) to use with Macenko stain normalization.

+ **`stain_patch_size`** Macenco stain normalization patch size.


+ **`mask_background_cutoff`** Parameter for detecting image background with HoVer-Net.

+ **`pixel_mask_threshold_low`** Parameter for detecting tissue pixels on the low resolution image.

+ **`pixel_mask_threshold_high`** Parameter for detecting tissue pixels on the low resolution image.

+ **`fraction_for_mask`** Fraction of pixels in tissue required to call tile in tissue.


+ **`use_provided_grid`** Whether to use the grid provided in the input samplesheet. If false and no Space Ranger alignment is done then a new grid of tiles si generated based on the grid parameters.

+ **`grid_type`** Type of the grid of tiles to generate. Can be hex, square, or random.

+ **`grid_spot_diamter`** Diameter of the spot (dimension of a tile) in micrometers.

+ **`grid_spot_horizontal_spacing`** Horizontal center-to-center distance between adjacent spots (or tiles).

+ **`grid_aspect_correction`** Factor to correct Visium slide aspect ratio.


+ **`overlap_scale_factor`** Imaging features extraction parameter. If the factor is 1 then features are extracted from the tile of the ST spot dimension.
 

+ **`hovernet_segmentation`** Do HoVer-Net segmetation. If false do StarDist segmentation. 

+ **`hovernet_batch_size`** Parameter of HoVer-Net segmetation. This parameter is ignored when segmentation is done with StarDist.

+ **`hovernet_num_inference_workers`** Parameter of HoVer-Net segmetation. This parameter is ignored when segmentation is done with StarDist.

+ **`hovernet_chunk_size`** Parameter of HoVer-Net segmetation. This parameter is ignored when segmentation is done with StarDist.

+ **`hovernet_tile_size`** Parameter of HoVer-Net segmetation. This parameter is ignored when segmentation is done with StarDist.



+ **`hovernet_spot_assignment_factor`** Used for either HoVer-Net or StarDist segmentation postprocessing. Scaling factor of the boundary limiting inclusion of nuclei to an ST spot. Value of 1 means the boundary size equals ST spot size.

+ **`hovernet_spot_assignment_shape`** Used for either HoVer-Net or StarDist segmentation postprocessing. Shape of the boundary, either square or disk.

+ **`hovernet_min_cell_type_prob`** Used for either HoVer-Net or StarDist segmentation postprocessing. Filtering parameter to remove nuclei assigned with low confidence.
