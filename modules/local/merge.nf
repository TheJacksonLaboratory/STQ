
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
