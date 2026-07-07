process GET_UNI2SUB_FEATURES {

    tag "$sample_id"
    label 'process_uni2'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.parquet', saveAs: { filename -> "${params.subtiling}-${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), path(segmask), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/uni2_features_cell_coords.parquet"), val(expansion_factor), val("uni2"), optional: true

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

    python -u ${projectDir}/bin/run-uni2-sub.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/uni2_features" \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --subtiling=${params.subtiling} \
    --subcoords-factor=${params.subcoords_factor} \
    --subcoords-list="${params.subcoords_list}" \
    --segmentation="${segmask}" \
    --upscale-factor=${params.subfeatures_upscale_factor} \
    --model-checkpoint-path=${params.uni2_model_checkpoint} \
    """
}



process GET_CTRANSPATHSUB_FEATURES {

    tag "$sample_id"
    label 'process_ctranspath'
    errorStrategy  { task.attempt <= 2  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.parquet', saveAs: { filename -> "${params.subtiling}-${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), path(segmask), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/ctranspath_features_cell_coords.parquet"), val(expansion_factor), val("ctranspath"), optional: true
    
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

    python -u ${projectDir}/bin/run-ctranspath-sub.py \
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
    --segmentation="${segmask}" \
    --upscale-factor=${params.subfeatures_upscale_factor} \
    --model="CTransPath"
    """   
}