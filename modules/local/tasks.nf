
process LOAD_SAMPLE_INFO {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(image), val(meta_grid)
    
    output:
    tuple val(sample_id), file(image), emit: main
    tuple val(sample_id), file("tissue_positions_list.csv"), file("scalefactors_json.json"), emit: grid
    
    script:
    """
    if [ ! "${meta_grid}" = "" ];
    then
        cp "${meta_grid}"/scalefactors_json.json .
        cp "${meta_grid}"/tissue_positions_list.csv .
    else
        echo "" >  "tissue_positions_list.csv"
        echo "" >  "scalefactors_json.json"
    fi
    """
}


process CONVERT_TO_TILED_TIFF {

    tag "$sample_id"
    label 'vips_process'

    input:
    tuple val(sample_id), path(image)
    
    output:
    tuple val(sample_id), file("outfile.tiff")

    script:
    """
    vips tiffsave ${image} ./outfile.tiff --compression none --tile --tile-width ${params.tiled_tiff_tile_size} --tile-height ${params.tiled_tiff_tile_size}
    """
}


process CREATE_THUMBNAIL_TIFF {

    tag "$sample_id"
    label 'vips_process'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(image)
    
    output:
    tuple val(sample_id), file("thumbnail.tiff")

    script:
    """
    vips resize ${image} thumbnail.tiff ${params.thumbnail_downsample_factor}
    """
}


process GET_PIXEL_MASK {

    tag "$sample_id"
    label 'python_process_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(image)
    
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


process TILE_WSI {

    tag "$sample_id"
    label 'python_process_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(image), path(meta_grid_csv), path(meta_grid_json)
    
    output:
    tuple val(sample_id), file("grid/grid.csv"), file("grid/grid.json"), emit: grid
    tuple val(sample_id), file("grid/grid.png"), emit: plot
    
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
                                  savepath=savepath, sname='grid')
    
    plotGrid(grid, *slide_dimensions,
             size=tile_size, show_spot_labels=False,
             show=False, object_shape='square',
             savepath=savepath, sname='grid')
    """
}


process GET_TILE_MASK {

    tag "$sample_id"
    label 'python_process_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(low_res_image), path(pixel_mask_csv), path(grid_csv), path(grid_json)
    
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
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
        
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json)
    
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
    
    mkdir inception

    python -u ${projectDir}/bin/run-inception-v3.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="inception/inception_features" \
    --overlap-scale-factor=${params.overlap_scale_factor} 
    """   
}


process INFER_HOVERNET {

    tag "$sample_id"
    label 'process_hovernet'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(image)
    
    output:
    tuple val(sample_id), file("hovernet/json/outfile.json"), emit: json
    tuple val(sample_id), file("hovernet/thumb/outfile.png"), file("hovernet/mask/outfile.png"), emit: mask
    
    script:
    """
    mkdir data
    cp "${image}" ./data/. 
    
    python -c 'import os; import tifftools; from tifftools.constants import Tag; \
    tifftools.tiff_set(os.getcwd() + "/data/outfile.tiff", overwrite=True, setlist=[(Tag[269], "40")]);'
      
    python /hover_net/run_infer.py \
    --gpu=0 \
    --model_mode=fast \
    --nr_types=6 \
    --type_info_path=/hover_net/type_info.json \
    --model_path=/hovernet_fast_pannuke_type_tf2pytorch.tar \
    wsi \
    --input_dir="./data/" \
    --output_dir=hovernet/ \
    --proc_mag=${params.magnification} \
    --save_thumb \
    --save_mask \
    > /dev/null 2>&1
    """ 
}


process COMPUTE_HOVERNET_DATA {

    tag "$sample_id"
    label 'python_process_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(hovernet_json)
    
    output:
    tuple val(sample_id), file("hovernet/per_nucleus_data.csv.gz")
    
    script:
    """
    #!/usr/bin/env python
    
    import sys
    sys.path.append("${projectDir}/lib")
    from hovernetConv import loadNuclei
    
    loadNuclei("${hovernet_json}", savepath='hovernet/', sname='per_nucleus_data')   
    """
}


process GENERATE_PERSPOT_HOVERNET_DATA {

    tag "$sample_id"
    label 'python_process_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(grid_csv), path(grid_json), path(hovernet_csv)
    
    output:
    tuple val(sample_id), file("hovernet/per_nucleus_data.csv.gz.csv"), emit: assignment
    tuple val(sample_id), file("hovernet/per_spot_data.csv"), emit: data
    
    script:
    """
    #!/usr/bin/env python
    
    import sys
    sys.path.append("${projectDir}/lib")
    from hovernetConv import assignNuceiToSpots, calculateAggregateValues

    df = assignNuceiToSpots(grid_file_path="${grid_csv}",
                       scalefactors_json_file="${grid_json}",
                       hovernet_data_file_path="${hovernet_csv}",
                       spot_shape="${params.hovernet_spot_assignment_shape}",
                       factor=${params.hovernet_spot_assignment_factor},
                       savepath='hovernet/')
                       
    calculateAggregateValues(df, min_cell_type_prob = ${params.hovernet_min_cell_type_prob},
                             savepath = 'hovernet/', sname = 'per_spot_data')
    """
}


process MERGE_IMAGING_DATA {

    tag "$sample_id"
    label 'python_process_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(inception_csv), path(hovernet_csv)
    
    output:
    tuple val(sample_id), file("data.csv.gz")
    
    script:
    """
    #!/usr/bin/env python
    
    import pandas as pd

    df_inception = pd.read_csv("${inception_csv}", index_col=0, header=0).sort_index() 
    df_hovernet = pd.read_csv("${hovernet_csv}", index_col=0, header=0).sort_index().reindex(df_inception.index)
    df = pd.concat([df_inception, df_hovernet], axis=1, sort=False)
    df.index.name = 'id'    
    df.to_csv("data.csv.gz")
    """
}
