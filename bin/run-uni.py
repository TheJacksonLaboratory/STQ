import argparse
import pandas as pd
import numpy as np
import openslide
import PIL
import json
from tqdm import tqdm

import timm
import torch
from torchvision import transforms
import torch.nn as nn
import openslide

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

def normalizer(img, mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225), size=224):
    func = transforms.Compose([transforms.Resize(size),
                              transforms.ToTensor(),
                              transforms.Normalize(mean=mean, std=std)])
    return func(img)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute features of each tile')
    parser.add_argument('--wsi-file', dest='wsi_file', action='store',
                        required=True,
                        help="""The path to the whole slide image (WSI) in a format readable by openslide (e.g., svs or ndpi).""")
    parser.add_argument('--model-checkpoint-path', dest='modelPath', action='store',
                        required=True,
                        help="""Path to bin checkpoint.""")
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
    parser.add_argument('--downsample-expanded', dest='downsample', action='store', default=True,
                        required=False,
                        help="""If expansion factor is greater than 1 then downsample the tiles back to the input size""")
    parser.add_argument('--expansion-factor', dest='expansion', action='store',
                        required=True,
                        help="""Expansion factor, 1 means no expansion""")
    parser.add_argument('--subtiling', dest='subtiling', action='store',
                        required=True,
                        help="""Do subtiling""")
    parser.add_argument('--subcoords-factor', dest='subcoordsf', action='store',
                        required=True,
                        help="""Factor for subtiling subtiling""")                        
    parser.add_argument('--subcoords-list', dest='subcoords', action='store',
                        required=True,
                        help="""Subtiling coordinates""")    

    args = parser.parse_args()
    expansion = float(args.expansion)
    downsample = args.downsample=='true'
    subtiling = args.subtiling=='true'

    subcoordsf = int(args.subcoordsf)
    subcoords = json.loads(args.subcoords)

    if expansion == 1.0:
        print('Expansion factor is 1, requested downsampling:', downsample)
        downsample = False
    else:
        if downsample:
            expansion = np.ceil(expansion)
            print('Expansion factor rounded to next interger:', expansion)
            print('Tiles will be expanded and then downsampled')
        else:
            print('Expansion without downsampling is requested')
    
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
    
    if downsample:
        num_rows = num_cols = round(spot_diameter_wsi)
    else:
        num_rows = num_cols = round(spot_diameter_wsi * expansion)

    # Load pre-trained UNI model
    model = timm.create_model("vit_large_patch16_224", img_size=224, patch_size=16, init_values=1e-5, num_classes=0, dynamic_img_size=True)
    model.load_state_dict(torch.load(args.modelPath, map_location="cpu"), strict=True)
    model.eval()

    num_images = len(pos)
    batch_size = int(10**8 / (float(args.expansion) * float(args.expansion) * num_cols * num_rows))
    if subtiling:
        batch_size = int(batch_size / 5)
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
                    if downsample:
                        ew = round(w * expansion)
                        eh = round(h * expansion)
                    else:
                        ew = w
                        eh = h
                    
                    img = np.array(slide.read_region((int(cx - ew / 2), int(cy - eh / 2)), lvl, (int(ew), int(eh))).convert('RGB'))

                    if subtiling:
                        a = int(np.floor(img.shape[0]/subcoordsf))
                        b = int(np.floor(img.shape[1]/subcoordsf))
                        for i, j in subcoords:
                            subimg = img[a*(i-1): a*(i+1), b*(i-1): b*(i+1), :]
                            images.append(subimg)
                    else:
                        # The downsampling is done to save memory
                        if downsample:
                            img = img[::int(expansion), ::int(expansion), :]
                            assert (img.shape[0], img.shape[1])==(w, h), 'Wrong tile dimensions after downsampling!'
                        
                        images.append(img)

            except Exception as exception:
                #print(exception)
                pass    
        print('Number of tiles:', len(images)) 

        if len(images)>0:
            images = torch.cat([normalizer(PIL.Image.fromarray(image))[None, :, :, :] for image in images], 0)
            with torch.inference_mode():
                temp_features = model(images).cpu().numpy()

                # Average the subtiles, e.g., every 5 subtiles
                if subtiling:
                    df_temp = pd.DataFrame(temp_features)
                    temp_features = df_temp.groupby(np.arange(len(df_temp.index))//len(subcoords)).mean().values

                features.append(temp_features)
        
    features = np.vstack(features)
    
    # Convert the dictionary of features to a dataframe and name its columns featXXX
    df_features = pd.DataFrame(features)
    df_features.columns = [f'feat_uni_' + str(i) for i in range(df_features.shape[1])]
    df_features.index = pos.loc[pos['in_tissue']==1].index
    
    # Append the spot position information to each row
    tbl = pd.concat([pos.loc[pos['in_tissue']==1], df_features], axis=1)
    print(tbl)   
    
    # Output the features with spot information
    ## This will automatically compress if the file suffix is .gz

    tbl.to_csv(output_path + '.tsv.gz', index=False)
    print('Successfully wrote ' + output_path)
    
exit(0)
