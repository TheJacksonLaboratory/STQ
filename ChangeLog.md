
v0.4.0 (pre-release)
+ Added a large set of new histopathology foundation models as tile-feature extractors: UNI2, Virchow, Virchow2, H-optimus-0, H-optimus-1, H0-mini, Phikon, Phikon-v2, prov-gigapath, prov-gigapath-flash, and Kaiko-Midnight12k (see docs/foundation-models.txt for the full list with license and origin).
+ Any of the feature extractors can run on either GPU (faster) or CPU, configurable.
+ Refactored the Inception v3 feature extractor from a Keras/TensorFlow implementation to a PyTorch/timm implementation, unifying it with the other foundation models under a single feature-extraction script (bin/run-extract.py) and container.
+ Optimized model-loading code paths shared across all foundation-model feature extractors to reduce redundant checkpoint loading overhead.
+ Added a new singularity container (assets/container-singularity-prov-gigapath-flash.def) dedicated to the prov-gigapath-flash model, exposed as params.container_gigapathf.
+ Improved H&E stain normalization robustness by estimating the reference stain matrix from an assembled canvas of representative nuclei/cell tiles (ASSEMBLE_TILES_CELLS), instead of relying on a single representative patch.
+ Added a reference stain matrix to use if the user does not provide a reference image.
+ Began refactoring pipeline configuration files toward compliance with Nextflow v26 conventions; this migration is in progress and not yet complete across all config files.
+ Updated documentation (README.md, docs/foundation-models.txt) to describe the expanded set of supported foundation models and the HuggingFace gated-access requirement for most non-Inception checkpoints.

v0.3.0
+ Updated imaging workflow output directory structure.
+ Added export of pipeline parameters, metadata items, and additional image outputs.
+ Added multiscale feature extraction inspired by SAMPLER work (PMID: 37577691).
+ Added experimental option of sub-tiling for use with feature extraction.
+ Added slide focus checking.
+ Added CTransPath, MoCoV3, UNI, CONCH as an extractors of imaging features.
+ Added subworkflow for sampling of tiles.
+ Optimized HoVer-Net segmentation steps. Added option for GPU-based segmentation.
+ Added postprocessing step to visualize outputs.
+ Added the AOI (automatic object identification) util for preparing ROI JSON for STQ.

v0.2.0
+ This version was referenced in the publication (PMID: 38626768).
+ Refactored and optimized codebase.
+ Added Xengsort read classification option.
+ Added documentation details to improve user experience.

v0.1.0
+ Initial release
