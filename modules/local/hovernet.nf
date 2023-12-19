
process GET_HOVERNET_MASK {

    tag "$sample_id"
    label 'process_hovernet_low'
    errorStrategy 'finish'
    
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
    --input_dir="./${image}" \
    --output_dir=hovernet/ \
    --slide_mag=40 \
    --proc_mag=40 \
    --chunk_shape=${params.hovernet_chunk_size} \
    --tile_shape=${params.hovernet_tile_size} \
    --save_mask
    """ 
}


process CHECK_MASK {

    tag "$sample_id"
    label 'process_hovernet_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'finish'
    
    input:
    tuple val(sample_id), path(mask), path(thumbnail)
    
    output:
    tuple val(sample_id), file("hovernet/mask/outfile.png")
    
    script:
    """
    #!/usr/bin/env python
    
    import sys
    sys.path.append("${projectDir}/lib")
    from hovernetConv import checkMask
    
    checkMask("${thumbnail}", "${mask}", "hovernet/mask/", bc=${params.mask_background_cutoff}) 
    """
}


process INFER_HOVERNET {

    tag "$sample_id"
    label 'process_hovernet'
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 2.GB }
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    
    input:
    tuple val(sample_id), path(image), path(mask), val(size)
    
    output:
    tuple val(sample_id), file("hovernet/outfile.json"), emit: json
    
    script:
    """
    [ ! -d "mask" ] && mkdir "mask"
    cp ${mask} mask/outfile.png
     
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
    --input_dir="./${image}" \
    --output_dir=hovernet/ \
    --input_mask_dir=mask/ \
    --slide_mag=40 \
    --proc_mag=40 \
    --chunk_shape=${params.hovernet_chunk_size} \
    --tile_shape=${params.hovernet_tile_size}
    """ 
}


process GET_NUCLEI_MASK_FROM_HOVERNET_JSON {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 0
    debug true
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 6.GB }
    
    input:
    tuple val(sample_id), path(image), file(json), val(size)
    
    output:
    tuple val(sample_id), file("nuclei.npy")
    
    script:
    """
    #!/usr/bin/env python
    import sys
    import json
    import cv2
    import tifffile
    import numpy as np

    sys.path.append("${projectDir}/lib")
    from hovernetConv import close_contour

    imgshape = tifffile.TiffFile("${image}").pages[0].shape
    print(imgshape)
    
    nuclei = np.zeros(imgshape, dtype=np.int32)

    with open("${json}", 'r') as tempfile:
        data = json.loads(tempfile.read())
    
    #for id in sorted(list(data['nuc'].keys())):
    #    print(id, end='\t')
    #    c = np.array([data['nuc'][id]['contour']], dtype=np.int32)
    #    vmin = tuple(c.min(axis=1)[0])
    #    vmax = tuple(c.max(axis=1)[0])
    #    c[0] -= vmin[0]
    #    c[1] -= vmin[1]    
    #    temp = np.zeros((vmax[0]+1, vmax[1]+1), dtype=np.int32)
    #    cv2.fillPoly(temp, c, 1)
    #    wh = np.where(temp!=0)
    #    nuclei[wh[0] + vmin[0], wh[1] + vmin[1]] = int(id)

    with open('nuclei.npy', 'wb') as tempfile:
        np.save(tempfile, nuclei)
    """
}


process INFER_HOVERNET_TILES {

    tag "$sample_id"
    label 'process_hovernet'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}/tiles", pattern: 'temp/overlay/*.png', saveAs: { filename -> "${filename.split("/")[filename.split("/").length - 1]}" }, mode: 'copy', overwrite: true
    memory { 8.GB }
    
    input:
    tuple val(sample_id), path("tiles/")
    
    output:
    tuple val(sample_id), file("temp/json/*.json"), emit: json
    tuple val(sample_id), file("temp/overlay/*.png"), emit: png
    
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
    tile \
    --input_dir="tiles/" \
    --output_dir="temp/" \
    --mem_usage=0.2
    
    rm -R temp/mat/
    """ 
}


process GET_NUCLEI_TYPE_COUNTS {

    tag "$sample_id"
    label 'process_hovernet_low'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path("json/")
    
    output:
    tuple val(sample_id), file("tiles/classes.csv.gz")
    
    script:
    """
    #!/usr/bin/env python
    
    import os
    import json
    import pandas as pd
    
    fnames = [fname for fname in os.listdir("json/") if fname[-len('.json'):]=='.json']
    num_images = len(fnames)
    print("Number of json files:", num_images)
    
    ses = []
    sen = []
    for fname in fnames:
        with open("json/" + fname, 'r') as tempfile:
            s = json.loads(tempfile.read())
        if len(s['nuc'].keys()) > 0:
            df = pd.DataFrame([(i, s['nuc'][i]['type'], s['nuc'][i]['type_prob']) for i in s['nuc'].keys()])
            se = df.loc[df[2]>=0.75][1].value_counts()
            ses.append(se)
            sen.append(fname[:-len('.json')])
    
    ses = pd.concat(ses, axis=1)
    ses.columns = sen
    ses = ses.fillna(0).sort_index()
    print(ses)
    
    if not os.path.exists("tiles/"):
        os.makedirs("tiles/")
    
    ses.T.to_csv('tiles/classes.csv.gz')
    """

}


process INFER_STARDIST {

    tag "$sample_id"
    label 'process_stardist'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    memory { 6.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 16.GB }
    
    input:
    tuple val(sample_id), path(image), path(mask), val(size)
    
    output:
    tuple val(sample_id), file("outfile.json"), emit: json
    tuple val(sample_id), file("nuclei.npy"), emit: mask
    
    script:
    """
    #!/usr/bin/env python
    import sys
    import json
    import shutil
    sys.path.append("${projectDir}/lib")
    from hovernetConv import makeJSONoutput
    import tifffile
    import numpy as np
    from csbdeep.utils import normalize
    from stardist.data import test_image_nuclei_2d
    from stardist.models import StarDist2D
    import matplotlib.pyplot as plt
    
    shutil.copytree("${params.stardist_model}", "custom_model/")
    model = StarDist2D(None, name="custom_model/")
    
    img = tifffile.imread("${image}")[..., :3]
    print(img.shape)
    
    tissuemask = plt.imread("${mask}")
    mask_reduction_factor = int(img.shape[0] / tissuemask.shape[0])
    print('Input tissue mask:', tissuemask.shape, mask_reduction_factor)
    
    #img[tissuemask==0] = 0
    
    from csbdeep.data import Normalizer, normalize_mi_ma
    class MyNormalizer(Normalizer):
        def __init__(self, mi, ma):
                self.mi, self.ma = mi, ma
        def before(self, x, axes):
            return normalize_mi_ma(x, self.mi, self.ma, dtype=np.float32)
        def after(*args, **kwargs):
            assert False
        @property
        def do_after(self):
            return False
    
    normalizer = MyNormalizer(0, 255)
    
    nuclei, details = model.predict_instances_big(img, axes='YXC', block_size=int($params.stardist_block_size / np.power(2, $task.attempt)), min_overlap=128, context=128, normalizer=normalizer, n_tiles=(4,4,1))
    
    data = makeJSONoutput(details)
    
    with open('outfile.json', 'w') as f:
        f.write(json.dumps(data))

    # Save nuclei mask
    with open('nuclei.npy', 'wb') as tempfile:
        np.save(tempfile, nuclei)
    """ 
}


process COMPRESS_JSON_FILE {

    tag "$sample_id"
    label 'vips_process'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(json_file)
    
    output:
    tuple val(sample_id), file("hovernet/outfile.json.gz")
    
    script:
    """
    [ ! -d "hovernet" ] && mkdir "hovernet"
    gzip -c ${json_file} > hovernet/outfile.json.gz
    """ 
}


process COMPUTE_SEGMENTATION_DATA {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), path(hovernet_json), val(size)
    
    output:
    tuple val(sample_id), file("hovernet/per_nucleus_data.csv.gz")
    
    script:
    """
    #!/usr/bin/env python
    
    import sys
    sys.path.append("${projectDir}/lib")
    from hovernetConv import loadNuclei
    print()
    
    loadNuclei("${hovernet_json}",
              savepath='hovernet/',
              sname='per_nucleus_data',
              original_mpp=${params.target_mpp})   
    """
}


process GENERATE_PERSPOT_SEGMENTATION_DATA {

    tag "$sample_id"
    label 'python_process_low'
    maxRetries 3
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true
    memory { 3.GB + (Float.valueOf(size) / 1000.0).round(2) * params.memory_scale_factor * 3.GB }
    
    input:
    tuple val(sample_id), path(grid_csv), path(grid_json), path(hovernet_csv), val(size)
    
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
