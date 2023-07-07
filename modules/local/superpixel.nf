
process SUPERPIXELATION {

    tag "$sample_id"
    label 'process_inception'
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 12.GB }
    publishDir "${params.outdir}/${sample_id}/superpixels", pattern: "superpixelation_*.png", mode: 'copy', overwrite: true
    cpus 1
    
    input:
    tuple val(sample_id), path(image), val(size)
    
    output:
    tuple val(sample_id), file("segmentation.npy"), emit: main
    tuple val(sample_id), file("superpixelation_*.png"), emit: images optional true

    script:    
    """
    python -u "${projectDir}/bin/superpixelation.py" \
    --inputImagePath "${image}" \
    --segmentationSavePath "segmentation.npy" \
    --pixelsPerSegment ${params.pixels_per_segment} \
    --compactness ${params.superpixel_compactness} \
    --s ${params.superpixel_patch_size} \
    --downsamplingFactor ${params.superpixel_downsampling_factor}
    """
}


process EXPORT_DOWN_IMAGE_FOR_CONTOURS {

    tag "$sample_id"
    label 'process_inception'
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 5.GB }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    cpus 1
    
    input:
    tuple val(sample_id), path(image), val(size)
    
    output:
    tuple val(sample_id), file("im_down.tiff")

    script:    
    """
    #!/usr/bin/env python
    import tifffile
    import numpy as np
    
    # Convert image to numpy array to remove OME TIFF metadata
    f = ${params.superpixel_downsampling_factor}
    img = np.array(tifffile.imread("${image}"))[::f, ::f, :3]
    print(img.shape)

    print('Saving downsampled image', flush=True)
    tifffile.imwrite("im_down.tiff", img, bigtiff=True)
    """
}


process EXPORT_SUPERPIXELATION_CONTOURS {

    tag "$sample_id"
    label 'process_inception'
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    cpus 1
    
    input:
    tuple val(sample_id), path(superpixelation), val(size)
    
    output:
    tuple val(sample_id), file("superpixelation.json.gz")

    script:    
    """
    #!/usr/bin/env python
    import numpy as np
    
    import sys
    sys.path.append("${projectDir}/lib")
    from superpixels import get_countours_from_mask, save_contours
    
    with open("${superpixelation}", 'rb') as tempfile:
        superpixelation = np.load(tempfile)
    print('Superpixelation mask shape:', superpixelation.shape)
    
    print('Computing contours', flush=True)
    contours = get_countours_from_mask(superpixelation)
    
    print('Saving contours', flush=True)
    save_contours(contours, filename='superpixelation.json.gz')
    """
}


process CALCULATE_CELLS_OD {

    tag "$sample_id"
    label 'process_inception'
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 10.GB }
    cpus 1
    
    input:
    tuple val(sample_id), path(img), path(nuclei), path(nuc_seg_measures), val(size)
    
    output:
    tuple val(sample_id), path("qp.csv")
    
    script:
    """
    #!/usr/bin/env python
    import numpy as np
    import tifffile
    import pandas as pd
    from tqdm import tqdm
    
    import sys
    sys.path.append("${projectDir}/lib")
    from hovernetConv import calculate_H_E_OD_quantities
    
    # Load nuclei mask
    with open("${nuclei}", 'rb') as tempfile:
        nuclei = np.load(tempfile)
    print('Nuclei segmetation mask shape:', nuclei.shape, flush=True)
        
    # Load nuc_seg_measures
    df_nuc_seg_measures = pd.read_csv("${nuc_seg_measures}", index_col=0)
    df_nuc_seg_measures.index = df_nuc_seg_measures.index.astype(str)
    print(df_nuc_seg_measures, flush=True)
    
    # Load image
    print('Loading full resolution image', flush=True)
    img = np.array(tifffile.imread("${img}"))[:, :, :3]
    dims = img.shape[0], img.shape[1]
    print(dims, flush=True)
    
    # Prepare image patches coordinates
    s = ${params.stardist_block_size}
    r = [np.append(s*np.array(range(0, int(np.floor(dims[i]/s))+1)), [dims[i]]) for i in range(2)]
    coords = [(i,j) for i in range(len(r[0])-1) for j in range(len(r[1])-1)]
    print(coords, flush=True)
    
    # Calculate HE OD quantities by patches
    dfs = []
    for ipatch, (i, j) in enumerate(tqdm(coords)):
        df_OD_temp = calculate_H_E_OD_quantities(img[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1], :],
                                                 nuclei[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1]],
                                                 (r[0][i], r[1][j]),
                                                 df_nuc_seg_measures,
                                                 expand_nuclei_distance=${params.expand_nuclei_distance})
        dfs.append(df_OD_temp)
    
    # Merge patches data, average cells fragments due to patching
    df_OD = pd.concat(dfs)
    df_OD = df_OD.groupby(level=0).mean()
    print(df_OD, flush=True)
    
    # Save data
    df_OD.to_csv("qp.csv")
    """
}


process ASSIGN_NUCLEI_TO_SUPERPIXELS {

    tag "$sample_id"
    label 'process_inception'
    maxRetries 1
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    cpus 1
    
    input:
    tuple val(sample_id), path(superpixelation), path(qp_od), val(size)
    
    output:
    tuple val(sample_id), path("od_per_cell.csv.gz")
    
    script:
    """
    #!/usr/bin/env python
    import numpy as np
    import pandas as pd
    
    # Load pre-calculated nuclear and cytoplasmic quantities
    df_qp_od = pd.read_csv("${qp_od}", index_col=0)
    print(df_qp_od)
    
    # Load superpixelation mask
    with open("${superpixelation}", 'rb') as tempfile:
        superpixelation = np.load(tempfile)
    print('Superpixelation mask shape:', superpixelation.shape)
    
    # Assign each cell to a superpixel; superpixelation is in downsampled coordinates
    df_qp_od['xd'] = (df_qp_od['x'] / ${params.superpixel_downsampling_factor}).astype(int)
    df_qp_od['yd'] = (df_qp_od['y'] / ${params.superpixel_downsampling_factor}).astype(int)
    df_qp_od['uspx'] = df_qp_od.apply(lambda se: superpixelation[int(se['yd']), int(se['xd'])], axis=1)
    df_qp_od['ipatch'] = df_qp_od['uspx'].apply(lambda v: int((v - v % 1000) / 1000))
    
    # Sort dataframe index
    df_qp_od.index = df_qp_od.index.astype(int)
    df_qp_od = df_qp_od.sort_index()
    df_qp_od.index = df_qp_od.index.astype(str)
    print(df_qp_od)
    
    # Normalize quantities within each ipatch identifier
    for col in df_qp_od.columns[~df_qp_od.columns.isin(['x', 'y', 'xd', 'yd', 'uspx', 'ipatch'])]:
        df_qp_od[col] = df_qp_od[col] * 0.18 / df_qp_od[col].quantile(0.5)
        
    se_sizes = pd.Series(superpixelation.ravel()).value_counts().sort_index()*(${params.superpixel_downsampling_factor}**2)
    df_qp_od['spx_size'] = (se_sizes.loc[df_qp_od['uspx'].values]).values
    
    # Export data
    df_qp_od.to_csv('od_per_cell.csv.gz')
    """    
}
