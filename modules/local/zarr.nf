
process CREATE_CELLSZARRZIP {

    tag "$sample_id"
    label 'zarr'
    maxRetries 0
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", pattern: 'cells.zarr.zip', mode: 'copy', overwrite: params.overwrite_files_on_publish
    publishDir "${params.outdir}/${sample_id}", pattern: 'cells-index.pkl', mode: 'copy', overwrite: params.overwrite_files_on_publish

    input:
        tuple val(sample_id), val(mpp), path(mfile)

    output:
        tuple val(sample_id), file("cells.zarr.zip")
        tuple val(sample_id), file("cells-index.pkl")

    script:
    """
    #!/usr/bin/env python
    import pickle
    import pandas as pd
    import numpy as np
    from skimage.segmentation import expand_labels
    import sys
    sys.path.append("${projectDir}/lib")
    from parquetcellsbuilder import parquet_to_cells_zarr
    from contours import get_contours

    # Make a parquet file from the mfile (mask.npy)
    pfile = 'boundaries.parquet'
    nuclear_expansion = int(float(${params.nuclear_expansion}) / float(${params.target_mpp}))  # in pixels
    mask = np.load("${mfile}")
    print(f"Loaded mask with shape {mask.shape} and unique labels: {np.unique(mask)}")
    mask = expand_labels(mask, distance=nuclear_expansion)
    print(f"Expanded mask with nuclear expansion of {nuclear_expansion} pixels.")
    contours, centroids = get_contours(mask)
    cell_ids = sorted(contours.keys())
    print(f"Extracted contours for {len(contours)} cells.")

    df = pd.DataFrame({
        'cell_id': [f"C{cid:08d}" for cid in cell_ids for _ in range(len(contours[cid]))],
        'vertex_x': [v[0] for cid in cell_ids for v in contours[cid]],
        'vertex_y': [v[1] for cid in cell_ids for v in contours[cid]],
    })
    df.to_parquet(pfile, index=False)
    print(f"Created DataFrame with {len(df)} rows for cell boundaries.")

    parquet_to_cells_zarr(
        parquet_path  = pfile,
        output_path   = "cells.zarr.zip",
        boundary_set  = 1,
        scalefactor   = ${params.target_mpp})
    print("Converted parquet to cells.zarr.zip")

    index = pd.Index(pd.read_parquet(pfile)['cell_id'])
    index = index[~index.duplicated(keep='first')]
    with open("cells-index.pkl", 'wb') as f:
        pickle.dump(index, f)
    print("Created cells-index.pkl")
    """
}
