
"""Written by S.Domanskyi, 2022

Module designed to generate a mask for a given grid of centers of tiles from Whole Slide Image (WSI).
Generate updated grid containing mask values, plot mask and low resolution image.

Examples of usage below:

# Load the module
import lib.wsiMask as wsiMask

# List module components
dir(wsiMask)

# Get help on all module functions
help(wsiMask)

## Generate in tissue mask
wsiMask.getInTissueMask(grid_csv='grid_sample.csv',
                        grid_json='grid_sample.json,
                        low_res_image='image_sample.tiff',
                        show=True, savepath='', sname='sample');
"""

import os
import json
import pandas as pd
import numpy as np

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.patches import Circle, Rectangle
from matplotlib.collections import PatchCollection

def plotMask(df, width: int = None, height: int = None, size: float = None, 
             image = None, figdim = 10, object_shape: str = 'spot', spot_alpha: float = 0.4,
             savepath: str = None, sname: str = '', show: bool = True):

    """Plot mask as square tiles or disks/spots
    
    Parameters:
    df: grid produced by function getGrid of wsiGrid module

    width: full resolution image width

    height: full resolution image height

    size: value produced by function getGrid

    image: low resolution image 3D array

    figdim: image scale, the bigger the value, the large mask image will be

    object_shape: ['spot', 'square'] shape of patch to plot as mask

    spot_alpha: transparency of the patches

    savepath: directory to save data files

    sname: identifier for saving data files

    show: display the image, needs interactive backend

    Output:
    None
    """
    
    figdim *= max(width, height) / 30000
    
    fig, ax = plt.subplots(figsize=(figdim, figdim))
    ax.imshow(image, origin='lower', extent=(0, width, 0, height))
    
    v = df[['pxl_row_in_fullres', 'pxl_col_in_fullres']].loc[df['in_tissue']==1].values
    
    if object_shape == 'spot':
        ax.add_collection(PatchCollection([Circle((x1, y1), size/2) for y1, x1 in v], alpha=spot_alpha, color='k', edgecolor=None, linewidth=0))
    elif object_shape == 'square':
        ax.add_collection(PatchCollection([Rectangle((x1-size/2, y1-size/2), size, size) for y1, x1 in v], alpha=spot_alpha, color='k', edgecolor=None, linewidth=0))
    else:
        raise NotImplementedError
        
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_aspect('equal')
    ax.axis('off')
    ax.plot([0, 0, width, width, 0], [0, height, height, 0, 0], linewidth=0.5, c='k', clip_on=False)
    
    fig.tight_layout()
    
    if not savepath is None:
        fig.savefig(savepath + 'mask_%s.png' % sname, dpi=150, facecolor='w')

    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return


def getInTissueMask(grid_csv: str, grid_json: str, low_res_image: str,
                    low: float = 100, high: float = 200, plot_mask: bool = True,
                    fraction: float = 0.1, savepath: str = '', sname: str = '', show: bool = False):

    """Plot mask as square tiles or disks/spots
    
    Parameters:
    grid_csv: csv file produced by function getGrid of wsiGrid module

    grid_json: json file produced by function getGrid of wsiGrid module

    low_res_image: path to file with low resolution image

    plot_mask: whether to make a plot with mask and low resolution image

    fraction: fraction of low resolution pixels in tissue to call patch in_tissue

    savepath: directory to save data files

    sname: identifier for saving data files

    show: display the image, needs interactive backend

    Output:
    df_grid: grid of centers with updated mask column
    """


    
    img_RGB_high_res = plt.imread(low_res_image)
    v = img_RGB_high_res[:, :, :3].copy()
    v = v.mean(axis=2)
    vc = v.copy()
    v[vc < low] = 0
    v[vc > high] = 0
    v[(vc >= low) & (vc <= high)] = 1
    df = pd.DataFrame(v.T)
    
    with open(grid_json) as f:
        info_dict = json.load(f)
    slide_fullres_width = info_dict['x']
    slide_fullres_height = info_dict['y']
    spot_diameter_fullres = info_dict['spot_diameter_fullres']
    
    scale_factor = 0.5 * (img_RGB_high_res.shape[0] / slide_fullres_height) + 0.5 * (img_RGB_high_res.shape[1] / slide_fullres_width)
    
    df_grid = pd.read_csv(grid_csv, header=None, index_col=0)
    df_grid.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    df_grid.index.name = 'barcode'
    
    for tile in df_grid.index[:]:
        tile_x = int(df_grid.loc[tile]['pxl_col_in_fullres'] * scale_factor)
        tile_y = int(df_grid.loc[tile]['pxl_row_in_fullres'] * scale_factor)
        tile_half_size = int(spot_diameter_fullres * scale_factor / 2)
        in_tissue = int(df.iloc[tile_x - tile_half_size : tile_x + tile_half_size,
                                tile_y - tile_half_size : tile_y + tile_half_size].mean().mean() >= fraction)
        df_grid.loc[tile, 'in_tissue'] = in_tissue
        
    if plot_mask:
        plotMask(df_grid, width=slide_fullres_width, height=slide_fullres_height,
                 size=spot_diameter_fullres, image=img_RGB_high_res[:, :, :3], 
                 figdim=10, object_shape='square', savepath=savepath, sname=sname, show=show)
        
    if not savepath is None:
        df_grid.to_csv(savepath + 'grid_with_mask_%s.csv' % sname, header=False)
        
    return df_grid
