# Prepared by Domanskyi
# The superpixelation is done by patches
# Each patch plot of superpixels is generated, small superpixels' identifiers are not shown

import os
import argparse
import tifffile
import numpy as np
from tqdm import tqdm
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from skimage.segmentation import mark_boundaries

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

import skimage
from pysnic.algorithms.snic import snic

def infere_spx(im_down_patch, target_number_of_segments='auto', pixels_per_segment=10000, compactness=1):
    
    lab_image = skimage.color.rgb2lab(im_down_patch).tolist()

    if target_number_of_segments == 'auto':
        target_number_of_segments = int(im_down_patch.shape[0] * im_down_patch.shape[1] / pixels_per_segment)
    
    segmentation, _, centroids = snic(lab_image, target_number_of_segments, compactness, update_func=None)
    segmentation = np.array(segmentation)

    return segmentation
    
def plot_all_spx_nf(im_down, seg, seg_id='', fontcolor='k', fontsize=8, fontweight='demibold',
                 figsize=(10, 10), boundaries_color=(1, 0, 0), min_size=500,
                 pe=path_effects.Stroke(linewidth=2, foreground='w')):

    if im_down.shape[0] > im_down.shape[1]:
        figsize = figsize[0] * im_down.shape[1] / im_down.shape[0], figsize[1]
    else:
        figsize = figsize[0], figsize[1] * im_down.shape[0] / im_down.shape[1]
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(mark_boundaries(im_down, seg, color=boundaries_color))
    
    for s in np.unique(seg.ravel()):
        wh = np.array(np.where(seg==s))
        if len(wh[0]) >= min_size:
            m = wh.mean(axis=1)
            params = dict(va='center', ha='center', color=fontcolor, fontsize=fontsize, fontweight=fontweight)
            ltext = ax.text(m[1], m[0], s, **params)
            ltext.set_path_effects([pe, path_effects.Normal()])
              
    ax.set_aspect('equal')  
    ax.axis('off')
    fig.tight_layout()
    
    plt.savefig(f'superpixelation_{seg_id}.png', facecolor='w', dpi=100, pad_inches=0.01)    
    plt.close(fig)
    
    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--inputImagePath", type=str, required=True, help="input image name")
    parser.add_argument("--segmentationSavePath", type=str, required=True, help="output file name")
    parser.add_argument('--s', type=int, default=4096, help='patch size')
    parser.add_argument('--compactness', required=True, type=float)
    parser.add_argument('--pixelsPerSegment', required=True, type=int)
    parser.add_argument('--downsamplingFactor', required=True, type=int)
    args = parser.parse_args()
    
    print('s:', args.s)

    # If the image is in 40x, the downsampling_factor 4 will bring the resolution to 10x
    img = np.array(tifffile.imread(args.inputImagePath))[::args.downsamplingFactor, ::args.downsamplingFactor, :3]
    print(img.shape)
    
    ## Save downsampled image
    #print('Saving downsampled image')
    #tifffile.imwrite(args.outputImagePath, img, bigtiff=True)
    #print('Done')
    
    dims = img.shape[0], img.shape[1]
    
    # Prepare image patches' coordinates
    r = [np.append(args.s*np.array(range(0, int(np.floor(dims[i]/args.s))+1)), [dims[i]]) for i in range(2)]
    coords = [(i,j) for i in range(len(r[0])-1) for j in range(len(r[1])-1)]
    print(coords)

    segmentation = np.zeros(dims, dtype=np.int32)
    for ipatch, (i, j) in enumerate(tqdm(coords)):
        try:
            seg_patch = infere_spx(img[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1], :],
                                   pixels_per_segment=args.pixelsPerSegment,
                                   compactness=args.compactness)
        
            # There are less than 1000 superpixels in each patch, less than 1000 patches
            # Make each superpixel id unique
            seg_patch += ipatch * 10000
            
            segmentation[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1]] = seg_patch
        except Exception as exception:
            print('Superpixel ERROR:', exception)
            
        try:
            plot_all_spx_nf(img[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1], :], 
                         segmentation[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1]],
                         seg_id=ipatch * 10000)
        except Exception as exception:
            print('Superpixel plot ERROR:', exception)
            
    # Save segmentation mask
    print('\nSegmentation:', segmentation.shape)
    with open(args.segmentationSavePath, 'wb') as tempfile:
        np.save(tempfile, segmentation)

    exit(0)
