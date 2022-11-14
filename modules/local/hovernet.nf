
process GET_HOVERNET_MASK {

    tag "$sample_id"
    label 'process_hovernet_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(image)
    
    output:
    tuple val(sample_id), file("hovernet/mask/outfile.png"), emit: mask
    
    script:
    """   
    python /hover_net/run_infer.py \
    --gpu="" \
    --device_mode="cpu" \
    --cpu_count=1 \
    --save_mask_and_exit \
    --model_mode=fast \
    --nr_inference_workers=1 \
    --nr_post_proc_workers=1 \
    --nr_types=6 \
    --type_info_path=/hover_net/type_info.json \
    --model_path=/hovernet_fast_pannuke_type_tf2pytorch.tar \
    --batch_size=1 \
    wsi \
    --input_dir="./data/${image}" \
    --output_dir=hovernet/ \
    --proc_mag=${params.magnification} \
    --chunk_shape=${params.hovernet_chunk_size} \
    --tile_shape=${params.hovernet_tile_size} \
    --save_mask
    """ 
}


process INFER_HOVERNET {

    tag "$sample_id"
    label 'process_hovernet'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(image), path(mask)
    
    output:
    tuple val(sample_id), file("hovernet/json/outfile.json"), emit: json
    
    script:
    """      
    python /hover_net/run_infer.py \
    --gpu="" \
    --device_mode="cpu" \
    --cpu_count=${task.cpus} \
    --model_mode=fast \
    --nr_inference_workers=${params.hovernet_num_inference_workers} \
    --nr_post_proc_workers=${task.cpus} \
    --nr_types=6 \
    --type_info_path=/hover_net/type_info.json \
    --model_path=/hovernet_fast_pannuke_type_tf2pytorch.tar \
    --batch_size=${params.hovernet_batch_size} \
    wsi \
    --input_dir="./data/${image}" \
    --output_dir=hovernet/ \
    --input_mask_dir=${mask} \
    --proc_mag=${params.magnification} \
    --chunk_shape=${params.hovernet_chunk_size} \
    --tile_shape=${params.hovernet_tile_size}
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
