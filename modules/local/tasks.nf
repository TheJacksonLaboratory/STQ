
process LOAD_SAMPLE_INFO {

    tag "$sample_id"

    input:
    tuple val(sample_id), val(meta), path(srgrid), path(image)
    
    output:
    tuple val(sample_id), file(image), file("roifile.json"), env(mpp), emit: main
    tuple val(sample_id), file("tissue_positions_list.csv"), file("scalefactors_json.json"), emit: grid
    
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
        cp "${meta.grid}"/scalefactors_json.json .
        cp "${meta.grid}"/tissue_positions_list.csv .
    fi
    """
}


process GET_IMAGE_SIZE {

    tag "$sample_id"
    label 'process_estimate_size'
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }

    input:
    tuple val(sample_id), path(fileslide), path(roifile), val(mpp)
    
    output:
    tuple val(sample_id), env(size)
    
    script:
    """  
    python -u ${projectDir}/bin/extractROI.py --fileslide="${fileslide}" --roifile="${roifile}" --sizefile="size.txt"
    size=`cat size.txt`
    """
}


process EXTRACT_ROI {

    tag "$sample_id"
    label 'process_extract'
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 12.GB }

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
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 30.GB }

    input:
    tuple val(sample_id), path("outfile.tiff"), val(mpp), val(size)
    
    output:
    tuple val(sample_id), file("output_images/outfile.tiff"), val(mpp)

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
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 12.GB }

    input:
    tuple val(sample_id), path("outfile.tiff"), val(mpp), val(size)
    
    output:
    tuple val(sample_id), file("output_images/outfile.tiff"), val(mpp)

    script:    
    """
    [ ! -d "output_images" ] && mkdir "output_images"
    
    python -u "${projectDir}/bin/StainToolsNorm.py" \
    --referenceImagePath "${params.stain_reference_image}" \
    --inputImagePath "outfile.tiff" \
    --outputImageName "output_images/outfile.tiff" \
    --s ${params.stain_patch_size}
    """
}


process CONVERT_TO_TILED_TIFF {

    tag "$sample_id"
    label 'vips_process'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'thumbnail.tiff', mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(image), val(mpp)
    
    output:
    tuple val(sample_id), file("converted/outfile.tiff"), emit: full
    tuple val(sample_id), file("thumbnail.tiff"), emit: thumb
    tuple val(sample_id), env(size), emit: size
    tuple val(sample_id), val(mpp), emit: mpp

    script:    
    """
    f=`echo ${mpp} / ${params.target_mpp} | bc -l`
    echo vips resize ${image} tempfile.tiff \$f
    vips resize ${image} tempfile.tiff \$f
    [ ! -d "converted" ] && mkdir "converted"
    
    vips resize tempfile.tiff thumbnail.tiff ${params.thumbnail_downsample_factor}
    vips tiffsave tempfile.tiff converted/outfile.tiff --compression none --tile --tile-width ${params.tiled_tiff_tile_size} --tile-height ${params.tiled_tiff_tile_size} --bigtiff
    rm tempfile.tiff
    
    w=`vipsheader -f width converted/outfile.tiff`
    h=`vipsheader -f height converted/outfile.tiff`
    
    size=`echo "\$w * \$h / 1000000" | bc -l`
    size=`echo "\$size/1" | bc`
    """
}


process GET_PIXEL_MASK {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
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
                         high = ${params.pixel_mask_threshold_high},
                         savepath = 'mask/', sname = 'pixel_mask')
    """
}


process GET_TISSUE_MASK {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
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
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), path(image), path(meta_grid_csv), path(meta_grid_json), val(size)
    
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
            
        import pandas as pd
        grid = pd.read_csv("${meta_grid_csv}", index_col=0, header=None) 
        grid.index.name = 'id'
        grid.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
        grid.to_csv(savepath + 'grid.csv', header=False)
        
        import json
        with open("${meta_grid_json}") as f:
            info_dict = json.load(f)
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
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
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
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 14.GB }
        
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size)
    
    output:
    tuple val(sample_id), file("inception/inception_features.tsv.gz"), optional: true
    
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
    
    [ ! -d "inception" ] && mkdir "inception"

    python -u ${projectDir}/bin/run-inception-v3.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="inception/inception_features" \
    --overlap-scale-factor=${params.overlap_scale_factor} 
    """   
}


process SELECT_SAVE_TILES {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'tiles/*.csv', mode: 'copy', overwrite: true
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
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'tiles/*.csv.gz', mode: 'copy', overwrite: true
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
