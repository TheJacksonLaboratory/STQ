# Nextflow Pipeline for Visium and H&E Data Processing

- [Overview](#overview)
- [Motivation](#Motivation)
- [Documentation](#Documentation)
- [Output](#Output)
- [Running the piepline](#running-the-piepline)
- [Tools used in the pipeline](#tools-used-in-the-pipeline)
  * fastq-tools
  * xenome
  * spaceranger
  * velocyto
  * bafextract
  * Inception v3
  * HoVer-Net
  * Stardist
- [Nextflow pipeline data flow](#nextflow-pipeline-data-flow)
- [Nextflow pipeline resources](#nextflow-pipeline-resources)
- [Glossary of Terms](#glossary-of-terms)

## Routes of analysis

<p>
    <img src="docs/route-map.png" width="800"/>
</p>

## Overview

This repository contains the source code of the nextflow implementation of the 10x Visium Spatial Gene Expression data and full-resolution H&E-stained Whole Slide Images (WSI) processing developed at [The Jackson Laboratory](https://www.jax.org/ "The Jackson Laboratory"). The overview of the pipeline is shown above. The primary input consists of compressed FASTQ files, reference FASTA files, and a full-resolution image of the 10x Visium Slide sample. Additional required inputs include either pre-built Xenome indices or host and graft genome assemblies, mouse and human reference transcriptomes for reads mapping, DL pre-trained model weights, and singularity containers with software tools.

## Motivation

Most of the steps implemented in our pipeline are computationally expensive and must be carried out on high-performance computer (HPC) systems. The most computationally intensive pipeline steps include RNA-seq reads mapping, full-resolution image alignment, preprocessing for RNA-velocity calculation, and preprocessing for RNA-based CNV inference, deep learning imaging features, and nuclear morphometrics data extraction. The pipeline generates a standardized set of files that can be used for downstream analysis using R-based Seurat of Python-based Scanpy or any other available environments.

## Documentation

The description of the pipeline components, parameters, analysis routes, required resources, and configuration guide are provided in this repository. The documentation files are [README.md](README.md), [conf/README.md](conf/README.md), and [workflows/README.md](workflows/README.md).

## Running the piepline

![Flow](docs/flow.gif)

<details closed><summary>Click to see all steps</summary><p>

<img src="docs/flow-static.png" width="1200"/>

</p></details>


##### Prerequisites
+ HPC environment with sufficient CPU and RAM and temporary storage resources

Processing 1 sample requires approximately 100+ CPU hours of computing time. Some of the processes need 1 CPU others need 4 CPUs or 8 CPUs as specified in the `nexflow.config` file. The temporary storage requires roughly 250 GB per sample for the pipeline to run. For example, if 32 samples are processed simultaneously, about 8TB of storage will be used until the pipeline completes.

+ Nextflow <img src="https://www.nextflow.io/img/nextflow2014_no-bg.png" height="30"/>

      https://www.nextflow.io/docs/latest/getstarted.html#installation

+ Singularity <img src="https://docs.sylabs.io/guides/3.0/user-guide/_static/logo.png" height="30"/>

      https://docs.sylabs.io/guides/3.0/user-guide/installation.html


+ The pipeline source code <img src="https://github.githubassets.com/images/modules/logos_page/Octocat.png" height="30"/> (this repository)

      mkdir my-piepline-run
      cd my-piepline-run
      git clone https://github.com/TheJacksonLaboratory/STQ.git
      cd STQ

+ Singularity software containers used in this pipeline

The singularity containers used in our pipeline can be downloaded or built with the definition `*.def` files and recipes contained in the directory `assets`.

> ***Note*** 10x Genomics requires that any software containers with Space Ranger are not shared publicly. We provide an example of a definition file for building a Space Ranger container with singularity: [assets/container-singularity-spaceranger.def](conf/container-singularity-spaceranger.def) that pulls a standard `debian:buster-slim` container from docker and installs all necessary Linux libraries. After that, a copy of Space Ranger is downloaded and installed from the 10x Genomics download portal. To obtain a download link for a specific version of Space Ranger user must navigate to https://www.10xgenomics.com/support/software/space-ranger/downloads, register, review, and accept any required user agreements from 10x, and copy the download link. Next, paste the link to a copy of the `def` file. Finally, build a container with any desired resource, for example, https://cloud.sylabs.io/builder.

+ Reference genomes for mouse and human (see tool `spaceranger count` below)

+ Graft and host reference FASTA files.





##### Samplesheet

The pipeline is designed to take one or more samples and process them in parallel. The sample sheet must have a `csv` format and must contain a header, as shown below:

| sample | fastq | image | grid | roifile | mpp |
|--------|-------|-------|------|---------|-----|
| sample_A   | path/to/sample1/fatqs/   | path/to/fullres/image1.tiff    |  |  | 0.22075 |
| sample_B   | path/to/sample2/fatqs/   | path/to/fullres/image2.tiff    |  |  | 0.22075 |
| sample_C   | path/to/sample3/fatqs/   | path/to/fullres/image3.tiff    |  |  | 0.22075 |
| ...    | ...   | ...   | ...  | ...     | ... |

> ***Tip*** Column headers must be named as shown above but not necessarily preserve the order, except that the first column must contain the sample identifiers. Even if a certain column is not necessary for a certain pipeline run, that column must be present in the CSV with an empty value, e.g., if roifile is left blank, then the entire WSI image is used in the analysis.

> ***Tip*** Sample identifiers in column 1 must be unique. The output directory will contain sub-directories named by the sample identifiers. The pipeline output reports are also based on these identifiers.

> ***Note*** The fastq files are expected to be compressed `fastq.gz` paired-end reads. Two files are expected in the fastq directory.


The column "grid" points to a path of directories where the sample spatial image alignment is stored. The value in the column "grid" can be empty if the WSI is not from the 10x Visium Spatial Gene Expression Slide or Space Ranger image alignment was not done.

The column "roifile" is used to specify the path to a JSON file that contains information on how to crop the WSI for processing in the imaging sub-workflow. Note that the sequencing sub-workflow image must contain the slide fiducials for Space Ranger to align the image and grid. When using a grid from the sequencing sub-workflow roifile value is ignored.

The column "mpp" is used to specify WSI resolution in microns per pixel. The recommended resolution of images that can be used with the imaging sub-workflow is approximately 0.25 microns per pixel. Images are downsampled or upsampled to the pre-specified in the configuration target mpp value.


##### Configure the pipeline

Edit file `nextflow.config` to specify paths to singularity containers, reference genome sequences, and deconvolution indices. If necessary, adjust any of the resources allowed to be allocated to the pipeline processes.

Check the defaults or edit the file `conf/analysis.config` to adjust any customizable parameters. See file `conf/README.md` for the description of each parameter.

For JAX users, preparation of the pipeline can be done by editing the file `submit.sb` and modifying the following lines:
+ workflow="two_references" ### "two_references" "one_reference" "arbitrary_grid" "xenome_indices"
+ samplesheet="/path/to/samplesheet.csv"
+ workdir="/flashscratch/[jaxuser]/some_work_run"
+ outdir="/path/to/results_my_analysis"

> ***Tip*** Use unique paths for "workdir" and "outdir" for each pipeline run to make the cleanup of temporary files easier. This strategy also allows the separation of pipeline run statistics and output samples.

> ***Note*** `workdir` must point to a device with a large storage volume and fast I/O access capability, e.g. flashscratch or fastscratch for JAX users.

##### Run the pipeline

###### Interactive mode (this mode use is discouraged)

For JAX users, start an interactive session and run the script:

        srun -p compute -q batch -t 6:00:00 --cpus-per-task=1 --mem=2G -J ijob --pty /bin/bash
        cd my-piepline-run/STQ/
        ./submit.sb

> ***Note*** The exemplified command above requests 6 hours of wall-time; this can be adjusted to a specific necessary run time.

> ***Tip*** Submitting via srun interactive session will show interactively updated progress of the pipeline run. We do not recommend using srun for our pipeline since any network interruption will cause the ssh connection to drop and the pipeline to fail. In such an unfortunate case user can rerun the pipeline. However, any previously unfinished processes will restart, while all finished processes resume from cached data.

###### Detached mode (preferred)

For JAX users use any HPC node (login node use is acceptable too) and submit the script to slurm:

       cd my-piepline-run/STQ/
       sbatch submit.sb

In this detached mode, the above command submits the pipeline to the HPC slurm system, which creates a low resource but long wall time job (duration specified in the submit.sb file), which manages the nextflow pipeline run. Nextflow manages all the pipeline processes and monitors the execution progress. Users can periodically monitor the contents of the `slurm-*.out` file to see progress:

       tail -n 50 slurm-<your_job_id>.out


## Output

Example output directory structure:

<pre>
.
├── <b>pipeline_info</b>
│   ├── <i>execution_report_2024-02-09_10-57-36.html</i>
│   ├── <i>execution_timeline_2024-02-09_10-57-36.html</i>
│   ├── <i>execution_trace_2024-02-09_10-57-36.txt</i>
│   └── <i>pipeline_dag_2024-02-09_10-57-36.svg</i>
│
└── <b>WM4237_TE_S1_ST</b>
    ├── <b>raw_feature_bc_matrix</b>
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── <b>spatial</b>
    │   ├── aligned_fiducials.jpg
    │   ├── detected_tissue_image.jpg
    │   ├── scalefactors_json.json
    │   ├── tissue_hires_image.png
    │   ├── tissue_lowres_image.png
    │   └── tissue_positions_list.csv
    │
    ├── xenome.summary.txt
    ├── <b>human</b>
    │   ├── extracted.baf
    │   ├── metrics_summary.csv
    │   ├── velocyto.loom
    │   └── web_summary.html
    ├── <b>mouse</b>
    │   ├── extracted.baf
    │   ├── metrics_summary.csv
    │   ├── velocyto.loom
    │   └── web_summary.html
    │
    ├── <b>grid</b>
    │   ├── grid.csv
    │   └── grid.json
    ├── <b>mask</b>
    │   ├── pixel_mask.csv
    │   ├── tile_mask.csv
    │   ├── tile_mask.png
    │   └── tissue_mask.png
    ├── <b>features</b>
    │   └── inception_features.tsv.gz
    ├── <b>nucseg</b>
    │   ├── outfile.json.gz
    │   ├── per_nucleus_data.csv.gz
    │   ├── per_nucleus_data.csv.gz.csv
    │   └── per_spot_data.csv
    |
    ├── data.csv.gz
    └── thumbnail.tiff
</pre>





## Tools used in the pipeline

1. **`fastq-tools`** (https://github.com/dcjones/fastq-tools)

2. **`xenome classify`** (https://github.com/data61/gossamer)

> Conway T, Wazny J, Bromage A, et al. Xenome--a tool for classifying reads from xenograft samples. Bioinformatics (Oxford, England). 2012 Jun;28(12):i172-8. DOI: 10.1093/bioinformatics/bts236. PMID: 22689758; PMCID: PMC3371868.

`xenome classify` is designed to classify xenograft-derived RNA-seq reads to deconvolve the graft (human) from the host (mouse) reads. Xenome defines classes of reads: definitely human, probably human, definitely mouse, probably mouse, both, ambiguous, neither. In xenome classes definitely human and probably human are combined into human; classes definitely mouse and probably mouse are combined into mouse. We discard reads classified as both, ambiguous, or neither. The statistics of Xenome reads deconvolution is generated in the sample output directory in file `xenome.summary.txt`.

`xenome classify` requires indices generated by `xenome index` as an input. The indices used in our analysis were built with `-K 35`, `-H Custom_Genomes/R84-REL1505/NOD_ShiLtJ/NOD_ShiLtJ.fa`, and `-G GCA_009914755.4_T2T-CHM13v2.0_genomic.fna`.

3. **`spaceranger count`** (https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count) is a pipeline developed by 10x Genomics based on a Martian pipeline (https://martian-lang.org/). 

> Martian pipeline is implemented in GOYACC (https://pkg.go.dev/golang.org/x/tools/cmd/goyacc), i.e., YACC (Yet Another Compiler-Compiler) for GO (https://go.dev/). We build `spaceranger count` pipeline into our pipeline as a nextflow process. 

`space ranger count` is designed to process one capture area of a Visium Spatial Gene Expression Slide for fresh frozen (FF) and formalin fixed paraffin embedded FFPE tissue samples. We use automatic capture area image alignment, where fiducials and the grid are detected by SpaceRanger. Note the image alignment is independent of the reference transcriptome or any sequencing parameters. The image alignment result is generated in the sample output folder under `spatial` and contains six files:
+ `aligned_fiducials.jpg`
+ `detected_tissue_image.jpg`
+ `scalefactors_json.json`
+ `tissue_hires_image.png`
+ `tissue_lowres_image.png`
+ `tissue_positions_list.csv`

Space Ranger requires the mouse, human, or combined reference genomes as input. These can be downloaded from the 10x Genomics website:

    wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
    wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
    wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-and-mm10-2020-A.tar.gz

and unpacked to create two directories, `refdata-gex-GRCh38-2020-A` and `refdata-gex-mm10-2020-A`:

    tar -xf refdata-gex-GRCh38-2020-A.tar.gz
    tar -xf refdata-gex-mm10-2020-A.tar.gz
    tar -xf refdata-gex-GRCh38-and-mm10-2020-A.tar.gz

> ***Note*** Each unpacked reference genome takes about 20 GB of storage space.

The Xenome-classified mouse reads alignment to the reference transcriptome is done using mouse reference, and the Xenome-classified human reads alignment to the reference transcriptome is done using human reference. The output mouse and human gene count matrices in the MTX format are merged into one and generated as output in the sample folder under `raw_feature_bc_matrix`. The "raw" means that all 4992 ST spots are present in the matrix. Summaries of the human and mouse read quantification are generated in the `human` and `mouse` directories of the sample output directory.
 
4. **`velocyto run10x`** (http://velocyto.org/velocyto.py/index.html)

Velocyto is run for all 4992 spots using the human or mouse genome reference and barcode-sorted BAM file. The latter is generated from a possorted BAM file which is the output from `spaceranger count`. The output is a `human/velocyto.loom` and `mouse/velocyto.loom` in the sample output directory. The `*.loom` files can be loaded with `velocyto.py` or `scvelo` (https://scvelo.readthedocs.io/) or other compatible toolsets for analysis of RNA-velocity.

> To learn more about RNA velocity, see: La Manno, G., Soldatov, R., Zeisel, A. et al. RNA velocity of single cells. Nature 560, 494-498 (2018). https://doi.org/10.1038/s41586-018-0414-6

> ***Note*** Analysis of the spliced/unspliced RNA counts generated in our pipeline is outside the scope of this pipeline and is carried out in the downstream steps.


5. **`BAFextract`** (https://github.com/akdess/CaSpER)

BAFextract is a tool designed by the authors of CaSpER and is intended to extract B-allele frequencies (BAF) from BAM files in order to estimate CNV events in the downstream steps with CaSpER.

> To learn about CaSpER, see: Serin Harmanci, A., Harmanci, A.O. & Zhou, X. CaSpER identifies and visualizes CNV events by integrative analysis of single-cell or bulk RNA-sequencing data. Nat Commun 11, 89 (2020). https://doi.org/10.1038/s41467-019-13779-x

BAFextract includes DNA scaffolds and mitochondrion (MT) DNA if those are present in the species genome reference. The output is a file `extracted.baf` generated in the human and mouse sub-directories of the sample output directory.


> ***Note*** CaSpER RNA-based CNV-inference is outside the scope of this pipeline and is carried out in the downstream steps.


6. **`Inception v3`** (https://keras.io/api/applications/inceptionv3/)

Inception v3 is a convolutional neural network model developed and trained by Google on more than a million images from the ImageNet database. In our pipeline, the Inception v3 step is designed to take a whole slide image and grid information in the same format as generated by the SpaceRanger pipeline for a 10x Visium ST image. The 2048 imaging features (or mones) are generated for each image tile defined by the provided tiling grid. 

<p>
    <img src="docs/mones-per-tile.png" width="450"/>
</p>

The default setting is to generate a grid that has the same geometry as 10x Visium Slide, except that tiles are set to be square instead of the round shape of ST spots. Note, in this default setting, tiles are not covering the image entirely. Users can change grid parameters in the `analysis.config` file.


<p>
    <img src="docs/example ST wsi.png" width="600"/>
</p>

7. **`HoVer-Net`** (https://github.com/vqdang/hover_net)

HoVer-Net nuclear segmentation and classification inference step was originally designed for GPU devices, but in our pipeline, is optimized for CPU. This pipeline step is computationally demanding.

> *Publication:* Simon Graham, Quoc Dang Vu, Shan E Ahmed Raza, Ayesha Azam, Yee Wah Tsang, Jin Tae Kwak, and Nasir Rajpoot, HoVer-Net: Simultaneous Segmentation and Classification of Nuclei in Multi-Tissue Histology Images, Medical Image Analysis, Volume 58, December 2019, 101563.

The inference calculation is carried out on the masked region (the white region is the tissue in the example below). The whole slide image is split into chunks; each chunk is split into tiles for segmentation. Ultimately HoVer-Net inference step produces a list of nuclear boundaries and assigned labels by using a pre-trained convolutional neural network (CNN) model on the PankNuke dataset, a dataset containing nearly 200,000 nuclei from 19 tissues from diverse cancer types that has undergone clinical annotation in a semi-automated manner and quality control by clinicians.

> *Publication:* Jevgenij Gamper, Navid Alemi Koohbanani, Ksenija Benes, Simon Graham, Mostafa Jahanifar, Seyyed Ali, Khurram, Ayesha Azam, Katherine Hewitt, and Nasir Rajpoot, PanNuke Dataset Extension, Insights and Baselines, arXiv:2003.10778v7 [eess.IV] 22 Apr 2020.

In-tissue mask:
<p>
    <img src="docs/hovernet-tissue-mask.png" width="300"/>
</p>

> ***Note*** The PanNuke dataset and the derived weights for HoVer-Net is licensed under Attribution-NonCommercial-ShareAlike 4.0 International.


8. **`StarDist`** (https://github.com/stardist/stardist)

StarDist is a method for the detection of star-convex-shaped objects from images. We use a pre-trained StarDist 2D model generated by the developers of StarDist.


## Nextflow pipeline data flow

Execution of the pipeline and data flow through processes enable generating the Directed Acyclic Graph (DAG). The inputs and outputs are interconnected with the logic encoded by or nextflow pipeline. The DAGs represent a snapshot (or analysis trace) based on user-specified or default pipeline parameters. Below are the DAGs for the three analysis routes: two references, one reference, and an arbitrary grid:

![Dag Arb](docs/dag-two.svg)

<details closed><summary>Directed Acyclic Graph for one-reference workflow</summary><p>

![Dag Two](docs/dag-one.svg)

</p></details>


<details closed><summary>Directed Acyclic Graph for arbitrary grid workflow</summary><p>

![Dag Two](docs/dag-arb.svg)

</p></details>


> ***Note*** DAG is generated during each pipeline run and is determined by the configuration and input parameters of the pipeline. For example, if specific optional processes are toggled "off," those processes and relevant graph connections will not show in the DAG.

## Nextflow pipeline resources

Our pipeline enables detailed tracing of resource usage. The run HTML report in the `pipeline_info` directory contains the time and resources spent by each process and each sample. A DAG and detailed raw execution trace are also generated in `pipeline_info` upon a successful pipeline run.

Example of CPU usage by various processes in the pipeline. A screenshot from the HTML report shows that SPACERANGER_HUMAN was the most CPU-intensive step in this example run:

<p>
    <img src="docs/Example_CPU_usage.png" width="1000"/>
</p>

## Temporary files

The Nextflow pipeline is designed to cache each step in the work directory, path to which is specified when the pipeline is invoked. The temporary files can take a large amount of storage space. After the pipeline runs are completed for all necesary analysis routes, the user has to remove the temporary files.


#### Glossary of Terms:
+ STQ - Spatial Transcriptomics Quantification
+ WSI - Whole Slide Image
+ HPC - High-Performance Computer
+ VIPS - open source toolset for image processing (IP), based on libvips. V stands for "vf", the first version of VIPS GUI in SunView
+ TF - TensorFlow
+ HoVer-Net - Horizontal and vertical distances prediction network
+ FF - Flash Frozen
+ FFPE - Formalin-Fixed Paraffin-Embedded
+ DAG - Directed Acyclic Graph
+ JSON - JavaScript Object Notation
+ CPU - Central Processing Unit
+ GPU - Graphical Processing Unit
+ CNV - Copy Number Variation
+ BAM - Binary sequence Alignment Map
+ mpp - microns per pixel
+ CNV - Copy Number Variation
+ SNV - Single Nucleotide Variant
+ BAF - B-allele frequency
+ YACC - Yet Another Compiler-Compiler
+ ST - Spatial Transcriptomics
+ Full-resolution image - capture area image supplied to SpaceRanger, e.g., 30000 x 30000 pixels
+ Low-resolution image - downsampled full-resolution capture area image, e.g., 600 x 600 pixels
+ High-resolution image - downsampled full-resolution capture area image, e.g., 2000 x 2000 pixels
