
process MERGE_IMAGING_DATA {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    memory { (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), path(inception_csv), path(hovernet_csv), val(size)
    
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


process MERGE_MTX {

    tag "$sample_id"
    label "python_low_process"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path('mtx_mouse/*'), path('mtx_human/*')

    output:
    tuple val(sample_id), file("raw_feature_bc_matrix/*")
    
    script:
        
    """    
    #!/usr/bin/env python
    
    import sys
    sys.path.append("${projectDir}/bin")   
    from mtx_tools import read_mtx_combine_and_write_mtx as combine
    from mtx_tools import read_sc_from_mtx as read
    
    combine(read('mtx_mouse/'), read('mtx_human/'), saveDataDir='raw_feature_bc_matrix/')
    """
}