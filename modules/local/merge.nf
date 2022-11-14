
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
