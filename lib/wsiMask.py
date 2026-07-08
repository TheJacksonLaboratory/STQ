
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

import cv2
from skimage.transform import resize
from scipy.ndimage import binary_fill_holes
from skimage.draw import disk

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
        if not os.path.exists(savepath):
            os.makedirs(savepath)
            
        fig.savefig(savepath + '%s.png' % sname, dpi=150, facecolor='w')

    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return
        
def getInTissuePixelMask(low_res_image: str, low: float=None, max_color: int=255, sh: int=200,
                        kernel_size: int=31, q: float=99, downsample: float=4.0,
                        savepath: str = None, sname: str = '', verbose=False):

    """Geet mask of pixels in tissue using Otsu's method and morphological operations
    
    Parameters:
    low_res_image: path to file with low resolution image or a numpy array
    
    low: low threshold

    max_color: maximum color value for the mask (default 255)

    sh: size of the border to add around the image for morphological operations (default 200)

    kernel_size: size of the kernel for Gaussian blur (default 31)

    q: percentile for contrast enhancement (default 99)

    downsample: factor to downsample the image for processing (default 4.0)
    
    savepath: directory to save data files

    sname: identifier for saving data files

    Output:
    Pixel in tissue mask
    """

    if isinstance(low_res_image, str):
        img_normalized = plt.imread(low_res_image)
    else:
        img_normalized = low_res_image
    fshape = img_normalized.shape
    fimage = img_normalized.copy()

    q99 = np.percentile(img_normalized, q)
    img_normalized = np.clip(img_normalized.astype(int) + 255 - q99, 0, 255).astype(np.uint8)

    img_normalized = cv2.resize(img_normalized, (0, 0), fx=1./downsample, fy=1./downsample)
    if verbose:
        print(f"Resized image from {fshape} to {img_normalized.shape} for processing.", flush=True)

    # Convert the image to grayscale using the luminosity
    img_normalized = (img_normalized.astype(np.float32) + 1) / 256
    luminance = 0.2989 * img_normalized[:, :, 0] + \
                0.5870 * img_normalized[:, :, 1] + \
                0.1140 * img_normalized[:, :, 2]

    # Apply optical density transformation
    mask = (-np.log(luminance) * 255).astype(np.uint8)

    if not low is None:
        if verbose:
            print('Using provided low threshold:', low, flush=True)
    else:
        # Use Otsu's method to find the optimal threshold
        maskb = cv2.GaussianBlur(mask, (kernel_size,kernel_size), 0)
        low = 0.5 * cv2.threshold(maskb, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)[0]
        if verbose:
            print('Using low threshold from Otsu:', low, flush=True)

    mask[mask < low] = 0
    mask[mask > 0] = max_color

    extended_mask = np.zeros((mask.shape[0]+2*sh, mask.shape[1]+2*sh), np.uint8)
    extended_mask[sh:-sh, sh:-sh] = mask
    mask = extended_mask

    mask = cv2.dilate(mask, None, iterations=2)
    mask = cv2.erode(mask, None, iterations=2)
    mask = cv2.dilate(mask, None, iterations=2)
    mask = cv2.erode(mask, None, iterations=2)

    mask = cv2.GaussianBlur(mask, (kernel_size,kernel_size), 0)
    th, im_th = cv2.threshold(mask, 0, max_color, cv2.THRESH_BINARY+cv2.THRESH_OTSU)

    h, w = im_th.shape[:2]
    mask = np.zeros((h+2, w+2), np.uint8)

    cv2.floodFill(im_th, mask, (0,0), max_color)
    mask = cv2.bitwise_not(mask)[1:-1, 1:-1]
    mask[mask<max_color] = 0

    mask = cv2.dilate(mask, None, iterations=2)
    mask = mask[sh:-sh, sh:-sh]

    if mask.shape[0] != fshape[0] or mask.shape[1] != fshape[1]:
        mask = cv2.resize(mask, (fshape[1], fshape[0]), interpolation=cv2.INTER_NEAREST)
        if verbose:
            print(f"Resized mask back to original image size: {mask.shape}", flush=True)

    df = pd.DataFrame(mask>0).astype(int).T
    
    if not savepath is None:
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        
        filepath = os.path.join(savepath, f"{sname}.csv")
        df.to_csv(filepath, header=False, index=False)

    return df

def getInTissueTileMask(pixel_mask_csv: str, grid_csv: str, grid_json: str, low_res_image: str, plot_mask: bool = True,
                    fraction: float = 0.1, savepath: str = None, sname: str = '', show: bool = False):

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
 
    with open(grid_json) as f:
        info_dict = json.load(f)
    slide_fullres_width = info_dict['x']
    slide_fullres_height = info_dict['y']
    spot_diameter_fullres = info_dict['spot_diameter_fullres']
    
    img_RGB_high_res = plt.imread(low_res_image)[:, :, :3]
    
    scale_factor = 0.5 * (img_RGB_high_res.shape[0] / slide_fullres_height) + 0.5 * (img_RGB_high_res.shape[1] / slide_fullres_width)
    
    df_grid = pd.read_csv(grid_csv, header=None, index_col=0)
    df_grid.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    df_grid.index.name = 'barcode'
    
    df_pixel_mask = pd.read_csv(pixel_mask_csv, index_col=None, header=None)
    
    for tile in df_grid.index[:]:
        tile_x = int(df_grid.loc[tile]['pxl_col_in_fullres'] * scale_factor)
        tile_y = int(df_grid.loc[tile]['pxl_row_in_fullres'] * scale_factor)
        tile_half_size = int(spot_diameter_fullres * scale_factor / 2)
        in_tissue = int(df_pixel_mask.iloc[tile_x - tile_half_size : tile_x + tile_half_size,
                                tile_y - tile_half_size : tile_y + tile_half_size].mean().mean() >= fraction)
        df_grid.loc[tile, 'in_tissue'] = in_tissue
        
    if plot_mask:
        plotMask(df_grid, width=slide_fullres_width, height=slide_fullres_height,
                 size=spot_diameter_fullres, image=img_RGB_high_res[:, :, :3], 
                 figdim=10, object_shape='square', savepath=savepath, sname=sname, show=show)
        
    if not savepath is None:
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        
        df_grid['in_tissue'].to_csv(savepath + '%s.csv' % sname, header=False)
        
    return df_grid

def makeTissueMaskFromTileMask(gridFile, gridInfoFile, tileMaskFile, squarePatch=False, upSizeFactor=1.5, 
                               downSizeChunkPx=1000, kernelSize=20, savePath='tissue_mask.png'):

    with open(gridInfoFile, 'r') as tempfile:
        info = json.loads(tempfile.read())
    s, x, y = int(info['spot_diameter_fullres']), info['x'], info['y']
    print(s, x, y)

    se_mask = pd.read_csv(tileMaskFile, index_col=0, header=None)[1].rename(None)

    df_grid = pd.read_csv(gridFile, index_col=0, header=None)[[4, 5]]
    df_grid.columns = ['x', 'y']
    df_grid.index.name = None

    df_grid = df_grid.loc[se_mask[se_mask==1].index.values]

    downsampleFactor = int(np.ceil(max(x, y) / downSizeChunkPx))
    
    m = np.zeros((x, y), dtype=np.int8)[::downsampleFactor, ::downsampleFactor]
    print(m.shape)

    maxxd = int(x / downsampleFactor) 
    maxyd = int(y / downsampleFactor) 

    for ty, tx in df_grid.values:
        xd = int(tx / downsampleFactor)
        yd = int(ty / downsampleFactor)
        radius = int(s * upSizeFactor / downsampleFactor)
        if squarePatch:
            m[xd-radius: xd+radius, yd-radius: yd+radius] = 1
        else:
            cc, rr = disk((xd, yd), radius)
            cc[cc<0] = 0
            cc[cc>maxxd] = maxxd-1
            rr[rr<0] = 0
            rr[rr>maxyd] = maxyd-1
            try:
                m[cc, rr] = 1
            except:
                pass

    m = m.T * 255
    m = resize(m, (int(np.round(y/downsampleFactor, 0)), int(np.round(x/downsampleFactor, 0))), order=3)
    m[m>=m.max()/2] = 255
    m[m<m.max()/2] = 0

    kernel = np.ones((kernelSize, kernelSize), dtype=np.uint8)
    m = cv2.dilate(m, kernel, iterations=1)
    m = binary_fill_holes(m).astype(np.uint8)
    m = cv2.erode(m, kernel, iterations=1)

    print(m.T.shape)
    
    cv2.imwrite(savePath, m * 255)
    
    return savePath
