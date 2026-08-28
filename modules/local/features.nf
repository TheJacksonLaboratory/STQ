process GET_UNI_FEATURES {

    tag "$sample_id"
    label 'process_uni'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/uni_features.csv.gz"), val(expansion_factor), val("uni"), optional: true
    
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
    
    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/uni_features" \
    --model-name="uni" \
    --model-suffix="uni" \
    --batch-size=${params.uni_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.uni_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_UNI2_FEATURES {

    tag "$sample_id"
    label 'process_uni2'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/uni2_features.csv.gz"), val(expansion_factor), val("uni2"), optional: true
    
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
    
    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/uni2_features" \
    --model-name="uni2" \
    --model-suffix="uni2" \
    --batch-size=${params.uni2_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.uni2_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_HOPTIMUS0_FEATURES {

    tag "$sample_id"
    label 'process_hoptimus0'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/hoptimus0_features.csv.gz"), val(expansion_factor), val("hoptimus0"), optional: true
    
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
    
    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/hoptimus0_features" \
    --model-name="hoptimus0" \
    --model-suffix="hoptimus0" \
    --batch-size=${params.hoptimus0_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.hoptimus0_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_CTRANSPATH_FEATURES {

    tag "$sample_id"
    label 'process_ctranspath'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/ctranspath_features.csv.gz"), val(expansion_factor), val("ctranspath"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/ctranspath_features" \
    --model-name="ctranspath" \
    --model-suffix="ctranspath" \
    --batch-size=${params.ctranspath_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.ctranspath_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_MOCOV3_FEATURES {

    tag "$sample_id"
    label 'process_mocov3'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/mocov3_features.csv.gz"), val(expansion_factor), val("mocov3"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/mocov3_features" \
    --model-name="mocov3" \
    --model-suffix="mocov3" \
    --batch-size=${params.mocov3_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.mocov3_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_CONCH_FEATURES {

    tag "$sample_id"
    label 'process_conch'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 60.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/conch_features.csv.gz"), val(expansion_factor), val("conch"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/conch_features" \
    --model-name="conch" \
    --model-suffix="conch" \
    --batch-size=${params.conch_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.conch_model_checkpoint} \
    --use-conch-normalizer=${params.use_conch_normalizer} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_INCEPTION_FEATURES {

    tag "$sample_id"
    label 'process_inception'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/inception_features.csv.gz"), val(expansion_factor), val("inception"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/inception_features" \
    --model-name="inception" \
    --model-suffix="inception" \
    --batch-size=${params.inception_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.inception_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_VIRCHOW_FEATURES {

    tag "$sample_id"
    label 'process_virchow'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/virchow_features.csv.gz"), val(expansion_factor), val("virchow"), optional: true
    
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
    
    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/virchow_features" \
    --model-name="virchow" \
    --model-suffix="virchow" \
    --batch-size=${params.virchow_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.virchow_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_VIRCHOW2_FEATURES {

    tag "$sample_id"
    label 'process_virchow2'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish
    
    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)
    
    output:
    tuple val(sample_id), file("features/virchow2_features.csv.gz"), val(expansion_factor), val("virchow2"), optional: true
    
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
    
    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/virchow2_features" \
    --model-name="virchow2" \
    --model-suffix="virchow2" \
    --batch-size=${params.virchow2_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.virchow2_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_PHIKON_FEATURES {

    tag "$sample_id"
    label 'process_phikon'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/phikon_features.csv.gz"), val(expansion_factor), val("phikon"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/phikon_features" \
    --model-name="phikon" \
    --model-suffix="phikon" \
    --batch-size=${params.phikon_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.phikon_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_PHIKON2_FEATURES {

    tag "$sample_id"
    label 'process_phikon2'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/phikon2_features.csv.gz"), val(expansion_factor), val("phikon2"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/phikon2_features" \
    --model-name="phikon2" \
    --model-suffix="phikon2" \
    --batch-size=${params.phikon2_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.phikonv2_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_HOPTIMUS1_FEATURES {

    tag "$sample_id"
    label 'process_hoptimus1'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/hoptimus1_features.csv.gz"), val(expansion_factor), val("hoptimus1"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/hoptimus1_features" \
    --model-name="hoptimus1" \
    --model-suffix="hoptimus1" \
    --batch-size=${params.hoptimus1_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.hoptimus1_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_H0MINI_FEATURES {

    tag "$sample_id"
    label 'process_h0mini'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/h0mini_features.csv.gz"), val(expansion_factor), val("h0mini"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/h0mini_features" \
    --model-name="h0mini" \
    --model-suffix="h0mini" \
    --batch-size=${params.h0mini_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.h0mini_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_MIDNIGHT12K_FEATURES {

    tag "$sample_id"
    label 'process_midnight12k'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/midnight12k_features.csv.gz"), val(expansion_factor), val("midnight12k"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/midnight12k_features" \
    --model-name="kaikomidnight12k" \
    --model-suffix="midnight12k" \
    --batch-size=${params.midnight12k_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.midnight_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_GIGAPATH_FEATURES {

    tag "$sample_id"
    label 'process_gigapath'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/gigapath_features.csv.gz"), val(expansion_factor), val("gigapath"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/gigapath_features" \
    --model-name="provgigapath" \
    --model-suffix="gigapath" \
    --batch-size=${params.gigapath_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.gigapath_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}

process GET_GIGAPATHF_FEATURES {

    tag "$sample_id"
    label 'process_gigapathf'
    errorStrategy  { task.attempt <= 1  ? 'retry' : 'finish' }
    memory { 56.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 28.GB }
    publishDir "${params.outdir}/${sample_id}/features", pattern: 'features/*.csv.gz', saveAs: { filename -> "${expansion_factor}-${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
    tuple val(sample_id), path(image), path(tile_mask), path(grid_csv), path(grid_json), path(meta_grid_csv), path(meta_grid_json), val(size), val(expansion_factor)

    output:
    tuple val(sample_id), file("features/gigapathf_features.csv.gz"), val(expansion_factor), val("gigapathf"), optional: true

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

    CUDEV="\${CUDA_VISIBLE_DEVICES:-}"

    python -u ${projectDir}/bin/run-extract.py \
    --wsi-file="${image}" \
    --positions-list-file="${grid_csv}" \
    --tile-mask="\${vtilemask}" \
    --scalefactors-json-file="${grid_json}" \
    --output-path="features/gigapathf_features" \
    --model-name="provgigapathflash" \
    --model-suffix="gigapathf" \
    --batch-size=${params.gigapathf_batch_size} \
    --num-workers=${task.cpus} \
    --expansion-factor=${expansion_factor} \
    --downsample-expanded=${params.downsample_expanded_tile} \
    --model-checkpoint-path=${params.gigapathf_model_checkpoint} \
    --cuda-visible-devices="\${CUDEV}"
    """   
}
