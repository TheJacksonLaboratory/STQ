
process LOAD_SAMPLE_INFO {

    tag "$sample_id"

    input:
    tuple val(sample_id), val(meta), path(srgrid), path(image)
    
    output:
    tuple val(sample_id), file(image), file("roifile.json"), env(mpp), emit: main
    tuple val(sample_id), file("tissue_positions_list.csv"), file("scalefactors_json.json"), emit: grid
    tuple val(sample_id), env(mpp), emit: mpp
    tuple val(sample_id), file(image), emit: image
    
    script:
    """
    mpp=${meta.mpp}

    if [ ! "${meta.roifile}" = "" ];
    then
        cp "${meta.roifile}" "roifile.json"
    else
        echo '{"0": {"location": 0, "size": 1}, "1": {"location": 0, "size": 1}}' > "roifile.json"
    fi
    
    if [ "${meta.grid}" = "" -o $params.use_provided_grid == false ];
    then
        if [ ! "${srgrid}" = "" ];
        then
            echo "Using a spaceranger grid"
            echo '{"0": {"location": 0, "size": 1}, "1": {"location": 0, "size": 1}}' > "roifile.json"
        else
            echo "Proceeding with a new grid"
            echo "" >  "tissue_positions_list.csv"
            echo "" >  "scalefactors_json.json"
        fi
    else
        echo "Using a grid from samplesheet"
        cp "${meta.grid}/scalefactors_json.json" .

        if [ -s "${meta.grid}/tissue_positions_list.csv" ]
        then
            cp "${meta.grid}"/tissue_positions_list.csv .
        else
            echo "Removed spots header"
            tail -n +2 "${meta.grid}"/tissue_positions.csv > ./tissue_positions_list.csv
        fi
    
    fi
    """
}


process GET_IMAGE_SIZE {

    tag "$sample_id"
    label 'process_estimate_size'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }

    input:
    tuple val(sample_id), path(fileslide), path(roifile), val(mpp)
    
    output:
    tuple val(sample_id), env(size)
    
    script:
    """  
    python -u ${projectDir}/bin/extractROI.py --fileslide="${fileslide}" --roifile="${roifile}" --sizefile="size.txt" --wholeside
    size=`cat size.txt`
    """
}


process EXTRACT_ROI {

    tag "$sample_id"
    label 'process_extract'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 18.GB }

    input:
    tuple val(sample_id), path(fileslide), path(roifile), val(mpp), val(size)
    
    output:
    tuple val(sample_id), file("outfile.tiff"), val(mpp), emit: image
    
    script:
    """
    echo $size
    python -u ${projectDir}/bin/extractROI.py --fileslide="${fileslide}" --roifile="${roifile}" --outfile="outfile.tiff" --extract=True
    """
}


process COLOR_NORMALIZATION {

    tag "$sample_id"
    label 'color_normalization_process'
    errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 30.GB }

    input:
    tuple val(sample_id), path("outfile.tiff"), val(size)
    
    output:
    tuple val(sample_id), file("output_images/outfile.tiff")

    script:    
    """
    [ ! -d "output_images" ] && mkdir "output_images"
    
    python -u "${projectDir}/bin/StainNetNorm.py" \
    --source_dir "." \
    --save_dir "output_images/" \
    --model_path ${params.stainnet}
    """
}


process STAIN_NORMALIZATION {

    tag "$sample_id"
    label 'stain_normalization_process'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 12.GB }

    input:
    tuple val(sample_id), path("outfile.tiff"), path("canvas.tif"), val(size)
    
    output:
    tuple val(sample_id), file("output_images/outfile.tiff")

    script:    
    """
    [ ! -d "output_images" ] && mkdir "output_images"
    
    python -u "${projectDir}/bin/StainToolsNorm.py" \
    --referenceImagePath "${params.stain_reference_image}" \
    --inputImagePath "outfile.tiff" \
    --outputImageName "output_images/outfile.tiff" \
    --canvasImagePath "canvas.tif" \
    --s ${params.stain_patch_size}
    """
}


process RESIZE_IMAGE {

    tag "$sample_id"
    label 'vips_process'
    errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
    
    input:
    tuple val(sample_id), path(image), val(mpp)
    
    output:
    tuple val(sample_id), file("tempfile.tiff"), emit: full
    tuple val(sample_id), env(size), emit: size

    script:    
    """
    f=`echo ${mpp} / ${params.target_mpp} | bc -l`
    vips resize ${image} tempfile.tiff \$f

    w=`vipsheader -f width tempfile.tiff`
    h=`vipsheader -f height tempfile.tiff`
    
    size=`echo "\$w * \$h / 1000000" | bc -l`
    size=`echo "\$size/1" | bc`
    """
}


process CONVERT_TO_TILED_TIFF {

    tag "$sample_id"
    label 'vips_process'
    errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'thumbnail.tiff', mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image)
    
    output:
    tuple val(sample_id), file("converted/outfile.tiff"), emit: full
    tuple val(sample_id), file("thumbnail.tiff"), emit: thumb

    script:    
    """
    [ ! -d "converted" ] && mkdir "converted"
    
    vips resize ${image} thumbnail.tiff ${params.thumbnail_downsample_factor}
    vips tiffsave ${image} converted/outfile.tiff --compression none --tile --tile-width ${params.tiled_tiff_tile_size} --tile-height ${params.tiled_tiff_tile_size} --bigtiff
    """
}


process GET_THUMB {

    tag "$sample_id"
    label 'process_extract'
    errorStrategy  { task.attempt <= 0  ? 'retry' : 'finish' }
    memory { 64.GB }
    publishDir "${params.outdir}/${sample_id}", pattern: 'thumbnail.tiff', mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image)
    
    output:
    tuple val(sample_id), file("thumbnail.tiff")

    script:    
    """
    #!/usr/bin/env python

    import tifffile
    from tifffile import TiffFile
    import numpy as np
    import cv2

    with TiffFile("${image}") as imgh:
        dims0 = imgh.pages[0].tags[256].value, imgh.pages[0].tags[257].value
        num_levels = len(imgh.series[0].levels)
        print(dims0, num_levels)

    # Try reading the smallest level of the image and resize it to thumbnail
    # If that doesn't work, read level 0, downsample it and resize it to thumbnail
    try:
        img = tifffile.imread("${image}", level=num_levels-1)

        # Adjust order of dimensions from OME-TIFF input
        if img.shape[0] == 3:
            img = np.moveaxis(img, 0, -1)
    except Exception as exception:
        print('Reading last level failed:', exception)
        img = tifffile.imread("${image}")

        assert img.shape[2] == 3
        f = 10
        img = img[::f, ::f, :]

    f = ${params.thumbnail_downsample_factor}
    target_dims = tuple((np.array(dims0) * f).astype(int))
    thumb = cv2.resize(img, dsize=target_dims, interpolation=cv2.INTER_CUBIC)
    print(thumb.shape)

    tifffile.imwrite("thumbnail.tiff", thumb, bigtiff=False)
    """
}


process MAKE_TINY_THUMB {

    tag "$sample_id"
    label 'process_extract'
    errorStrategy  { task.attempt <= 0  ? 'retry' : 'finish' }
    memory { 2.GB }
    publishDir "${params.outdir}/${sample_id}", pattern: 'thumbnail.jpeg', mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image)
    
    output:
    tuple val(sample_id), file("thumbnail.jpeg")

    script:    
    """
    #!/usr/bin/env python
    
    import tifffile
    img = tifffile.imread("${image}")
    tifffile.imwrite("thumbnail.jpeg", img, compression=('jpeg', 85), bigtiff=False)
    """
}


process GET_PIXEL_MASK {

    tag "$sample_id"
    label 'python_process_low'
    errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), path(image), val(size)
    
    output:
    tuple val(sample_id), file("mask/pixel_mask.csv")
    
    script:
    """
    #!/usr/bin/env python
    
    import sys
    sys.path.append("${projectDir}/lib")
    from wsiMask import getInTissuePixelMask

    getInTissuePixelMask(low_res_image = "${image}",
                         low = ${params.pixel_mask_threshold_low},
                         kernel_size=${params.pixel_mask_kernel_size},
                         savepath = 'mask/', sname = 'pixel_mask')
    """
}


process GET_TISSUE_MASK {

    tag "$sample_id"
    label 'python_process_low'
    errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 6.GB * task.attempt }
    
    input:
    tuple val(sample_id), path(meta_grid_csv), path(meta_grid_json), path(tile_mask), val(size)
    
    output:
    tuple val(sample_id), file('mask/tissue_mask.png')
    
    script:
    """
    #!/usr/bin/env python
    
    import os
    import sys
    sys.path.append("${projectDir}/lib")
    from wsiMask import makeTissueMaskFromTileMask

    if not os.path.exists("mask/"):
        os.makedirs("mask/")
    
    makeTissueMaskFromTileMask("${meta_grid_csv}",
                               "${meta_grid_json}",
                               "${tile_mask}",
                               savePath='mask/tissue_mask.png')
    """
}


process TILE_WSI {

    tag "$sample_id"
    label 'python_process_low'
    errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 6.GB * task.attempt }
    
    input:
    tuple val(sample_id), path(image), path(meta_grid_csv), path(meta_grid_json), val(size), val(mpp)
    
    output:
    tuple val(sample_id), file("grid/grid.csv"), file("grid/grid.json"), emit: grid
    tuple val(sample_id), file("grid/grid.png"), emit: plot, optional: true
    
    script:
    """
    #!/usr/bin/env python        
       
    import sys
    sys.path.append("${projectDir}/lib")
    from wsiGrid import getGrid, plotGrid
    import os
    import openslide   
    import PIL.Image
    PIL.Image.MAX_IMAGE_PIXELS = None
    with openslide.open_slide("${image}") as slide:
        slide_dimensions = slide.dimensions
    
    savepath = 'grid/'
    
    # Spaceranger output is given
    if os.path.getsize("${meta_grid_csv}") > 10:
        import os
        if not os.path.exists(savepath):
            os.makedirs(savepath)
            
        f = ${mpp} / ${params.target_mpp}
        ${params.thumbnail_downsample_factor}
          
        import pandas as pd
        grid = pd.read_csv("${meta_grid_csv}", index_col=0, header=None) 
        grid.index.name = 'id'
        grid.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
        grid['pxl_row_in_fullres'] = (grid['pxl_row_in_fullres'] * f).astype(int)
        grid['pxl_col_in_fullres'] = (grid['pxl_col_in_fullres'] * f).astype(int)
        grid.to_csv(savepath + 'grid.csv', header=False)
        
        import json
        with open("${meta_grid_json}") as tempfile:
            info_dict_sr = json.load(tempfile)
        print(info_dict_sr)
        
        info_dict = dict()
        info_dict['tissue_lowres_scalef'] = float(${params.thumbnail_downsample_factor})
        info_dict['spot_diameter_fullres'] = info_dict_sr['spot_diameter_fullres'] * f
        
        tile_size = info_dict['spot_diameter_fullres']
        info_dict['x'], info_dict['y'] = slide_dimensions
        
        with open(savepath + 'grid.json', 'w') as outfile:        
            outfile.write(json.dumps(info_dict))
    else:
        grid, tile_size = getGrid(*slide_dimensions,
                                  grid_type="${params.grid_type}",
                                  spot_diamter=${params.grid_spot_diamter},
                                  spot_horizontal_spacing=${params.grid_spot_horizontal_spacing},
                                  resolution=1./${params.target_mpp},
                                  aspect_correction=${params.grid_aspect_correction},
                                  savepath=savepath, sname='grid')
    if False:
        plotGrid(grid, *slide_dimensions,
                 size=tile_size, show_spot_labels=False,
                 show=False, object_shape='square',
                 savepath=savepath, sname='grid')
    """
}


process GET_TILE_MASK {

    tag "$sample_id"
    label 'python_process_low'
    errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 6.GB * task.attempt }
    
    input:
    tuple val(sample_id), path(low_res_image), path(pixel_mask_csv), path(grid_csv), path(grid_json), val(size)
    
    output:
    tuple val(sample_id), file("mask/tile_mask.csv"), emit: mask
    tuple val(sample_id), file("mask/tile_mask.png"), emit: plot
    
    script:
    """
    #!/usr/bin/env python
    
    import sys
    sys.path.append("${projectDir}/lib")
    from wsiMask import getInTissueTileMask
            
    getInTissueTileMask(pixel_mask_csv = "${pixel_mask_csv}", grid_csv = "${grid_csv}",
                        grid_json = "${grid_json}", low_res_image = "${low_res_image}", 
                        fraction = ${params.fraction_for_mask},
                        plot_mask = True, show = False,
                        savepath = 'mask/', sname = 'tile_mask')
    """
}


process GET_INCEPTION_FEATURES {

    tag "$sample_id"
    label 'process_inception'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    memory { 36.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 14.GB }
    //publishDir "${params.outdir}/${sample_id}", pattern: 'features/*.tsv.gz', mode: 'copy', overwrite: true
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.tsv.gz', saveAs: { filename -> "${params.subtiling}-${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/inception_features.tsv.gz"), val(expansion_factor), val("inception"), optional: true
    
    script:
    """
    # If the grid from SpaceRanger, then don't owerwrite the tile mask with my mask
    filesize=`wc -c <"${meta_grid_csv}"`
    if [ \$filesize -ge 10 ];
    then
        vtilemask=None
    else
        vtilemask="${tile_mask}"
    fi
    
    [ ! -d "features" ] && mkdir "features"

    python -u ${projectDir}/bin/run-inception-v3.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/inception_features" \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile}
    """   
}


process GET_CTRANSPATH_FEATURES {

    tag "$sample_id"
    label 'process_ctranspath'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    //publishDir "${params.outdir}/${sample_id}", pattern: 'features/*.tsv.gz', mode: 'copy', overwrite: true
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.tsv.gz', saveAs: { filename -> "${params.subtiling}-${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/ctranspath_features.tsv.gz"), val(expansion_factor), val("ctranspath"), optional: true
    
    script:
    """
    # If the grid from SpaceRanger, then don't owerwrite the tile mask with my mask
    filesize=`wc -c <"${meta_grid_csv}"`
    if [ \$filesize -ge 10 ];
    then
        vtilemask=None
    else
        vtilemask="${tile_mask}"
    fi
    
    [ ! -d "features" ] && mkdir "features"

    python -u ${projectDir}/bin/run-ctranspath.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/ctranspath_features" \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --subtiling=${params.subtiling} \
    --subcoords-factor=${params.subcoords_factor} \
    --subcoords-list="${params.subcoords_list}" \
    --model="CTransPath" \
    --cuda-visible-devices=""
    """   
}


process GET_MOCOV3_FEATURES {

    tag "$sample_id"
    label 'process_mocov3'
    errorStrategy  { task.attempt <= 0  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.tsv.gz', saveAs: { filename -> "${params.subtiling}-${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/mocov3_features.tsv.gz"), val(expansion_factor), val("mocov3"), optional: true
    
    script:
    """
    CUDEV="\$CUDA_VISIBLE_DEVICES"
    
    # If the grid from SpaceRanger, then don't owerwrite the tile mask with my mask
    filesize=`wc -c <"${meta_grid_csv}"`
    if [ \$filesize -ge 10 ];
    then
        vtilemask=None
    else
        vtilemask="${tile_mask}"
    fi
    
    [ ! -d "features" ] && mkdir "features"

    python -u ${projectDir}/bin/run-ctranspath.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/mocov3_features" \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --subtiling=${params.subtiling} \
    --subcoords-factor=${params.subcoords_factor} \
    --subcoords-list="${params.subcoords_list}" \
    --model="MoCoV3" \
    --cuda-visible-devices="\$CUDEV"
    """   
}


process GET_UNI_FEATURES {

    tag "$sample_id"
    label 'process_uni'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.tsv.gz', saveAs: { filename -> "${params.subtiling}-${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/uni_features.tsv.gz"), val(expansion_factor), val("uni"), optional: true
    
    script:
    """    
    # If the grid from SpaceRanger, then don't owerwrite the tile mask with my mask
    filesize=`wc -c <"${meta_grid_csv}"`
    if [ \$filesize -ge 10 ];
    then
        vtilemask=None
    else
        vtilemask="${tile_mask}"
    fi
    
    [ ! -d "features" ] && mkdir "features"

    python -u ${projectDir}/bin/run-uni.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/uni_features" \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --subtiling=${params.subtiling} \
    --subcoords-factor=${params.subcoords_factor} \
    --subcoords-list="${params.subcoords_list}" \
    --model-checkpoint-path=${params.uni_model_checkpoint} \
    """   
}


process GET_CONCH_FEATURES {

    tag "$sample_id"
    label 'process_conch'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    memory { 60.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.tsv.gz', saveAs: { filename -> "${params.subtiling}-${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/conch_features.tsv.gz"), val(expansion_factor), val("conch"), optional: true
    
    script:
    """    
    # If the grid from SpaceRanger, then don't owerwrite the tile mask with my mask
    filesize=`wc -c <"${meta_grid_csv}"`
    if [ \$filesize -ge 10 ];
    then
        vtilemask=None
    else
        vtilemask="${tile_mask}"
    fi
    
    [ ! -d "features" ] && mkdir "features"

    python -u ${projectDir}/bin/run-conch.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/conch_features" \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --subtiling=${params.subtiling} \
    --subcoords-factor=${params.subcoords_factor} \
    --subcoords-list="${params.subcoords_list}" \
    --model-checkpoint-path=${params.conch_model_checkpoint} \
    --use-conch-normalizer=${params.use_conch_normalizer}
    """   
}


process SELECT_SAVE_TILES_RAW {

    tag "$sample_id"
    label 'python_low_process'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 2.GB }
    
    input:
    tuple val(sample_id), path(image)
    
    output:
    tuple val(sample_id), file("tiles/*.tif"), emit: tiles
    
    script:
    """
    #!/usr/bin/env python
    
    import os
    import numpy as np
    import zarr
    import tifffile

    if not os.path.exists('tiles/'):
        os.makedirs('tiles/')

    with tifffile.TiffFile("${image}") as tif:
        dims = tif.pages[0].tags[256].value, tif.pages[0].tags[257].value
    print(dims)

    tile_size = 224 # Pixels
    s = tile_size // 2
    N_tiles = 100

    np.random.seed(0)
    cx = np.random.randint(s, dims[0] - s, N_tiles)
    cy = np.random.randint(s, dims[1] - s, N_tiles)

    with tifffile.imread("${image}", aszarr=True) as store:
        z = zarr.open(store, mode='r')
        if z.shape[0] <= 4:
            z = np.moveaxis(z, 0, -1)
        print(z.shape)

        for i in range(N_tiles):
            print(i, cx[i], cy[i])
            tile = z[cy[i]-s:cy[i]+s, cx[i]-s:cx[i]+s, :]
            tifffile.imwrite(f'tiles/tile-{i}.tif', tile)
    """
}


process ASSEMBLE_TILES_CELLS {

    tag "$sample_id"
    label 'python_low_process'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'canvas.tif', mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 4.GB }
    
    input:
    tuple val(sample_id), path("tiles/"), path("json/")
    
    output:
    tuple val(sample_id), file("canvas.tif")
    
    script:
    """
    #!/usr/bin/env python
    
    import os
    import json
    import tifffile
    import numpy as np

    patches = []
    N = len([f for f in os.listdir("json/") if f.endswith('.json')])
    hcsize = 24 # pixels

    for i in range(N):
        try:
            with open(f"json/tile-{i}.json") as f:
                d = json.load(f)
            centroids = np.array([c['centroid'] for i, c in d['nuc'].items()]).astype(int)

            tile = tifffile.imread(f"tiles/tile-{i}.tif")
            msize = tile.shape[0]
            # Adjust centroids such that they are at least hcsize pixels from the edge of the tile
            centroids = np.clip(centroids, hcsize, msize - hcsize)

            # Extract patches around centroids
            for x, y in centroids:
                patches.append(tile[y-hcsize:y+hcsize, x-hcsize:x+hcsize])
        except Exception as e:
            print(f"Error processing tile {i}: {e}")
    print(f"Extracted {len(patches)} patches.")

    # Collage patches into a single image
    S = int(np.floor(np.sqrt(len(patches))))
    canvas = np.zeros((2*hcsize*S, 2*hcsize*S, patches[0].shape[2]), dtype=patches[0].dtype)
    for i, patch in enumerate(patches[:S*S]):
        x = (i % S) * 2 * hcsize
        y = (i // S) * 2 * hcsize
        canvas[y:y+2*hcsize, x:x+2*hcsize] = patch

    # Save the canvas as a TIFF file
    tifffile.imwrite("canvas.tif", canvas)
    """
}


process SELECT_SAVE_TILES {

    tag "$sample_id"
    label 'python_process_low'
    errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'tiles/*.csv', mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 2.GB }
    
    input:
    tuple val(sample_id), path(image), path(grid_csv), path(grid_json), path(tile_mask)
    
    output:
    tuple val(sample_id), file("tiles/*.tif"), emit: tiles
    
    script:
    """
    #!/usr/bin/env python
    
    import os
    import pandas as pd
    import numpy as np
    import json
    import openslide
    import PIL
    import PIL.Image
    PIL.Image.MAX_IMAGE_PIXELS = None
    
    if not os.path.exists('tiles/'):
        os.makedirs('tiles/')
     
    se_mask = pd.read_csv("${tile_mask}", index_col=1, header=None)[0].xs(1)
    np.random.seed(0)
    sel_tiles = se_mask.sample(min(${params.tiles_per_slide}, se_mask.shape[0]))
    sel_tiles.to_csv('tiles/tiles.csv', index=False, header=False)
    
    df_grid = pd.read_csv("${grid_csv}", index_col=0, header=None)[[4, 5]].loc[sel_tiles.values]
    
    with open("${grid_json}", 'r') as tempfile:
        s = int(json.loads(tempfile.read())['spot_diameter_fullres'])
    
    slide = openslide.open_slide("${image}")  
    for id in df_grid.index:
        cy, cx = df_grid.loc[id]
        print(id, s, cx, cy)
        slide.read_region((int(cx - s / 2), int(cy - s / 2)), 0, (int(s), int(s))).convert('RGB').save('tiles/%s.tif' % id)
    """
}


process GET_INCEPTION_FEATURES_TILES {

    tag "$sample_id"
    label 'process_inception'
    errorStrategy  { task.attempt <= 3  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'tiles/*.csv.gz', mode: 'copy', overwrite: params.overwrite_files_on_publish
    memory { 4.GB }
    
    input:
    tuple val(sample_id), path("tiles/")
    
    output:
    tuple val(sample_id), file("tiles/features.csv.gz")
    
    script:
    """
    python -u ${projectDir}/bin/run-inception-v3-tiles.py \
    --input-path="tiles/" \
    --output-path="tiles/features.csv.gz"
    """   
}
