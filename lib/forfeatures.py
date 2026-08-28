"""Shared helpers reusable across the different feature-extraction backends
(e.g., features.py for torch-based models, featuresinception.py for the
TensorFlow-based InceptionV3 model). Nothing in this module depends on torch
or tensorflow, so it is safe to import from either backend/container.
"""
import os
import numpy as np
import openslide
import pandas as pd

def load_positions(positions_list_file, tile_mask):
    """Load spaceranger positions list and translate pixel coords into wsi resolution."""
    pos = pd.read_csv(positions_list_file, header=None)
    pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col',
                   'pxl_row_in_fullres', 'pxl_col_in_fullres']

    if tile_mask != 'None' and tile_mask is not None:
        print('Received tile mask %s' % tile_mask)
        mask = pd.read_csv(tile_mask, index_col=0, header=None)
        pos['in_tissue'] = mask.reindex(pos['barcode'].values).values

    # scale_factor = ratio of resolution of 'wsi_file' to resolution of "fullres" image input to spaceranger.
    scale_factor = 1
    pos['pxl_row_in_wsi'] = pos.pxl_row_in_fullres * scale_factor
    pos['pxl_col_in_wsi'] = pos.pxl_col_in_fullres * scale_factor
    return pos, scale_factor

def resolve_expansion_and_downsample(expansion, downsample):
    if expansion == 1.0:
        print('Expansion factor is 1, requested downsampling:', downsample)
        downsample = False
    elif downsample:
        expansion = np.ceil(expansion)
        print('Expansion factor rounded to next interger:', expansion)
        print('Tiles will be expanded and then downsampled')
    else:
        print('Expansion without downsampling is requested')
    return expansion, downsample

def read_tile(slide, pos, row_idx, w, h, expansion, downsample, lvl=0):
    """Read and (optionally) downsample a single tile centered on a spot. Returns np array or None."""
    cy = pos.loc[row_idx, 'pxl_row_in_wsi']
    cx = pos.loc[row_idx, 'pxl_col_in_wsi']

    if not pos.loc[row_idx, 'in_tissue']:
        return None

    if downsample:
        ew, eh = round(w * expansion), round(h * expansion)
    else:
        ew, eh = w, h

    img = np.array(slide.read_region(
        (int(cx - ew / 2), int(cy - eh / 2)), lvl, (int(ew), int(eh))).convert('RGB'))

    if downsample:
        img = img[::int(expansion), ::int(expansion), :]
        assert (img.shape[0], img.shape[1]) == (w, h), 'Wrong tile dimensions after downsampling!'

    return img

def _try_read_tile(task):
    """Wrapper so a bad/out-of-range row never kills the thread pool; mirrors the
    original per-tile try/except behavior (silently skip, return None)."""
    slide, pos, row_idx, w, h, expansion, downsample = task
    try:
        return read_tile(slide, pos, row_idx, w, h, expansion, downsample)
    except Exception:
        return None

def save_features(features, pos, modelSuffix, output_path):
    df_features = pd.DataFrame(features)
    df_features.columns = [f'feat_{modelSuffix}_' + str(i) for i in range(df_features.shape[1])]
    df_features.index = pos.loc[pos['in_tissue'] == 1].index

    df_features = pd.concat([pos.loc[pos['in_tissue'] == 1], df_features], axis=1)
    print(df_features)

    output_path_dir = os.path.dirname(output_path)
    if not os.path.exists(output_path_dir):
        os.makedirs(output_path_dir)

    df_features.to_csv(output_path + '.csv.gz', index=False)
    print('Successfully wrote ' + output_path)
    return
