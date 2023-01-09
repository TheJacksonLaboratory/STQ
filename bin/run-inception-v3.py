import os
import argparse
import pandas as pd
import numpy as np
import cv2
import tensorflow as tf
import datatable as dt
import openslide
from pathlib import Path
import itertools
import PIL
import json
from scipy.ndimage import gaussian_filter
from tqdm import tqdm

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute Inception V3 features on tiles that cover a _single_ spatial transcriptomic spot')
    parser.add_argument('--wsi-file', dest='wsi_file', action='store',
                        required=True,
                        help="""The path to the whole slide image (WSI) in a format readable by openslide (e.g., svs or ndpi).""")
    parser.add_argument('--positions-list-file', dest='positions_list_file', action='store',
                        required=True,
                        help="""The positions_list.csv file output by spaceranger that has one row per spot and columns indicating whether the spot is within the tissue and its x and y coordinates in pixels.""")
    parser.add_argument('--scalefactors-json-file', dest='scalefactors_json_file', action='store',
                        required=True,
                        help="""The scalefactors_json.json file output by spaceranger that defines the spot diameter in spaceranger's full resolution (i.e., the resolution of the file input to spaceranger, which may or may not be wsi_file).""")
    parser.add_argument('--output-path', dest='output_path', action='store',
                        required=True,
                        help="""Name of _CSV_ file in which to store the feature matrix (rows are tiles, cols are features). 
                                The file will be compressed if it is named *.gz""")
    parser.add_argument('--tile-mask', dest='tile_mask', default=None, action='store', required=False)
    parser.add_argument('--sigma', dest='sigma', action='store', default=0, required=False,
                        help="""Gaussian smoothing parameter for background subtraction""")
    parser.add_argument('--overlap-scale-factor', dest='overlap', action='store',
                        required=True,
                        help="""Overlap factor, 1 means no overlap""")
    args = parser.parse_args()
    overlap = float(args.overlap)
    sigma = float(args.sigma)
    wsi_file = args.wsi_file
    positions_list_file = args.positions_list_file
    scalefactors_json_file = args.scalefactors_json_file
    output_path = args.output_path
    # Read in the spaceranger positions list file
    pos = pd.read_csv(positions_list_file, header=None)
    pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    
    if args.tile_mask != 'None':
        print('Received tile mask %s' % args.tile_mask)
        mask = pd.read_csv(args.tile_mask, index_col=0, header=None)
        pos['in_tissue'] = mask.reindex(pos['barcode'].values).values

    # Read the spot diameter at spaceranger's "full resolution" from the scalefactors_json file
    # output by spaceranger, i.e., in the resolution of the file passed to spaceranger, which may not
    # be the same resolution of wsi_file.
    with open(scalefactors_json_file) as f:
        scalefactors_tbl = json.load(f)
    spot_diameter_fullres = scalefactors_tbl['spot_diameter_fullres']
    
    
    # scale_factor = ratio of resolution of 'wsi_file' to resolution of "fullres" image input to spaceranger.
    # scale_factor = 4
    # NB: ideally, this code would accept the full resolution image along with the wsi_file and compare their sizes.
    # You would do that with something like (wait ... probably the full resolution image is a png/jpg/etc not openable by openslide)
    # full_resolution_slide = openslide.open_slide(full_resolution_file)
    # base_magnification = float(full_resolution_slide.properties[openslide.PROPERTY_NAME_OBJECTIVE_POWER])
    scale_factor = 1
    # Define the spot diameter in the resolution of the wsi_file
    spot_diameter_wsi = round(spot_diameter_fullres * scale_factor)
    # Translate the pixel coordinates from full resolution to the resolution of the wsi
    pos['pxl_row_in_wsi'] = pos.pxl_row_in_fullres * scale_factor
    pos['pxl_col_in_wsi'] = pos.pxl_col_in_fullres * scale_factor
    # Create the inception v3 model
    num_dimensions = 3
    num_rows = num_cols = round(spot_diameter_wsi * overlap)
    base_modeli = tf.keras.applications.inception_v3.InceptionV3(include_top=False, weights='imagenet',
                                                                 input_shape=(num_rows, num_cols, num_dimensions),
                                                                 classes=2)
    xi = base_modeli.output
    xi = tf.keras.layers.GlobalAveragePooling2D(data_format=None)(xi)
    model = tf.keras.models.Model(inputs=base_modeli.input, outputs=xi)
    
    num_images = len(pos)
    batch_size = int(10**8 / (num_cols * num_rows))
    num_batches = int(np.ceil(num_images / batch_size))
    
    print('Reading and pocessing tiles:', num_images)
    print('Batch size:', batch_size)
    print('Number of batches:', num_batches)
    
    slide = openslide.open_slide(wsi_file)
    
    w = num_cols
    h = num_rows
    lvl = 0
    features = []
    for ibatch in tqdm(range(num_batches)):
        images = []
        for indx in range(batch_size):
            try:
                cy = pos.loc[indx + ibatch*batch_size, 'pxl_row_in_wsi']
                cx = pos.loc[indx + ibatch*batch_size, 'pxl_col_in_wsi']
                if pos.loc[indx + ibatch*batch_size, 'in_tissue']:
                    images.append(np.array(slide.read_region((int(cx - w / 2), int(cy - h / 2)), lvl, (int(w), int(h))).convert('RGB')))
            except:
                pass
 
        if len(images)>0:
            features.append(model.predict(tf.keras.applications.inception_v3.preprocess_input(np.stack(images)), verbose=0))
    features = np.vstack(features)
    
    # Convert the dictionary of features to a dataframe and name its columns featXXX
    df_features = pd.DataFrame(features)
    df_features.columns = ['feat' + str(i) for i in range(df_features.shape[1])]
    df_features.index = pos.loc[pos['in_tissue']==1].index
    
    # Append the spot position information to each row
    tbl = pd.concat([pos.loc[pos['in_tissue']==1], df_features], axis=1)
    print(tbl)   
    
    # Output the features with spot information
    DT = dt.Frame(tbl)
    ## This will automatically compress if the file suffix is .gz   
    
    overlapAppendName = '' if overlap==1.0 else '_overlap_%s' % overlap
    sigmaAppendName = '' if sigma==0.0 else '_corrected_%s' % args.sigma

    DT.to_csv(path=output_path + '%s%s.tsv.gz' % (sigmaAppendName, overlapAppendName), header=True, compression="auto") 
    print('Successfully wrote ' + output_path)
    
exit(0)
