
"""Written by S.Domanskyi, 2022

Module designed to generate a grid of centers of tiles from Whole Slide Image (WSI)
Examples of usage below:

# Load the module
import lib.wsiGrid as wsiGrid

# List module components
dir(wsiGrid)

# Get help on all module functions
help(wsiGrid)

# Generate and plot Visium-like hexagonal grid
slide_dimensions = 2900, 3060
grid, tile_size = wsiGrid.getGrid(*slide_dimensions, grid_type='hex')
wsiGrid.plotGrid(grid, *slide_dimensions, size=tile_size)

# Generate and plot Slide-seq-like random grid, save plot to current working directory
slide_dimensions = 2900, 3060
grid, tile_size = wsiGrid.getGrid(*slide_dimensions, grid_type='hex')
grid, tile_size = wsiGrid.perturbGrid(*slide_dimensions, grid, tile_size, delta=0.1)
wsiGrid.plotGrid(grid, *slide_dimensions, size=tile_size, savepath='', show=False)

# Generate and plot square grid, specify magnification and spot diameter
slide_dimensions = 3000, 1100
grid, tile_size = wsiGrid.getGrid(*slide_dimensions, grid_type='square', magnification=20, spot_diamter=55)
wsiGrid.plotGrid(grid, *slide_dimensions, size=tile_size)

# enerate and plot large Visium-like hexagonal grid, save to file to current working directory
slide_dimensions = 29000, 30600
grid, tile_size = wsiGrid.getGrid(*slide_dimensions, savepath='')
wsiGrid.plotGrid(grid, *slide_dimensions, size=tile_size, show_spot_labels=False, savepath='')
"""

import json
import pandas as pd
import numpy as np

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.patches import Circle, Rectangle
from matplotlib.collections import PatchCollection

def getGrid(x: int, y: int, grid_type: str = 'hex', factor: float = 64/39, magnification: float = 40.0,
               resolution: float = 294/55, spot_diamter: float = 55, spot_horizontal_spacing: float = 85,
               aspect_correction: float = 0.95, savepath: str = None, sname: str = ''):
    
    """Generate grid of tile centers

    Parameters:
    x: full resolution image width

    y: full resolution image height

    grid_type: ['hex', 'square']

    factor: Visium Spatial Gene Expression hex grid factor

    magnification: image magnification

    resolution: pixels per micron of sample at 40x magnification

    spot_diamter: spot diameter in microns

    spot_horizontal_spacing: spot horizontal center-to-center distance in microns

    aspect_correction: Visium capture area is not square, even though officially it is 6.5x6.5mm

    savepath: directory to save data files

    sname: identifier for saving data files

    Output:
    grid: pandas.DataFrame

    tile_size_pixels: tile size
    """
    
    tile_size_pixels = resolution * spot_diamter * magnification / 40.0
    tile_horizontal_spacing_pixels = resolution * spot_horizontal_spacing * magnification / 40.0
    if grid_type=='hex':
        tile_vertical_spacing_pixels = 0.5 * factor * tile_horizontal_spacing_pixels / aspect_correction
    elif grid_type=='square':
        tile_vertical_spacing_pixels = tile_horizontal_spacing_pixels
    else:
        raise NotImplementedError
      
    if not savepath is None:
        info_dict = {'grid_type': grid_type, 'factor': factor, 'magnification': magnification,
                     'resolution': resolution, 'aspect_correction': aspect_correction,
                     'spot_diamter': spot_diamter, 'spot_horizontal_spacing': spot_horizontal_spacing,
                     'spot_diameter_fullres': tile_size_pixels, 'x': x, 'y': y}
        with open(savepath + 'grid_%s.json' % sname, 'w') as outfile:
            outfile.write(json.dumps(info_dict))
    
    nx = int(np.ceil(x / tile_horizontal_spacing_pixels))
    ny = int(np.ceil(y / tile_vertical_spacing_pixels))
    
    _grid = [['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']]
    for i in range(nx):
        for j in range(ny):
            temp_x = tile_size_pixels/2. + i * tile_horizontal_spacing_pixels
            if grid_type=='hex':
                if j % 2 == 1:
                    temp_x += 0.5 * tile_horizontal_spacing_pixels
            temp_x = int(temp_x)
            temp_y = int(tile_size_pixels/2. + j * tile_vertical_spacing_pixels)
            if (temp_x + tile_size_pixels/2. <= x) and (temp_y + tile_size_pixels/2. <= y):
                _grid.append([1, j, i, temp_y, temp_x])
                
    _grid = pd.DataFrame(_grid).T.set_index(0).T
    _grid.index = 'tile-' + _grid.index.astype(str).str.pad(8, fillchar='0')
    _grid.index.name = 'barcode'
    
    if not savepath is None:
        _grid.to_csv(savepath + 'grid_%s.csv' % sname, header=False)
    
    return _grid, tile_size_pixels

def perturbGrid(x, y, grid, tile_size_pixels, n_iterations: int = 5, delta: float = 0.5, seed: int = None, dmax: float = 7., verbose: int = 1):
    
    """Random perturbations of the grid produce grid limilar to Slide-Seq ST technology

    Parameters:
    x: full resolution image width

    y: full resolution image height

    grid: produced by function getGrid

    tile_size_pixels: produced by function getGrid

    n_iterations: number of iterations to randomly perturb the drid

    delta: fraction of the tile size that the tile can be displaced along x or y at most in one move

    seed: random seed to have reproducible perturbation

    dmax: max neighbor distance, number of tile sizes away from the spot center

    verbose: set to 0 to suppress print output

    Output:
    grid: pandas.DataFrame

    tile_size_pixels: tile size
    """
    
    if not seed is None:
        np.random.seed(seed)
        
    x_col_ind = np.where(grid.columns=='pxl_col_in_fullres')[0][0]
    y_col_ind = np.where(grid.columns=='pxl_row_in_fullres')[0][0]
    
    if verbose >= 1:
        print('\tComuting neigbors of each spot with dmax: %s' % dmax)
    se = pd.Series(index=range(len(grid)), dtype='object')
    v = grid[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values
    for _i, (_y, _x) in enumerate(v):
        nn = set(np.where((((v.T[0] - _y)**2 + (v.T[1] - _x)**2)**0.5) < (dmax * tile_size_pixels))[0])
        se[_i] = nn.difference({_i})
        
    for iter in range(n_iterations):
        if verbose >= 1:
            print('Iteration: %s' % iter)
        for _i in range(len(grid)):
            _x = grid.iloc[_i, x_col_ind]
            _y = grid.iloc[_i, y_col_ind]
            
            _p = (np.random.rand(2) - 0.5) * delta * tile_size_pixels
            
            nviolations = 0
            
            nn = np.array(list(se[_i]))
            if len(nn) > 0:
                vx = grid.iloc[nn, x_col_ind].values
                vy = grid.iloc[nn, y_col_ind].values
                d = ((vy - _y + _p[0])**2 + (vx - _x + _p[1])**2)**0.5
                nviolations += len(set(np.where(d < tile_size_pixels)[0]))

            if ((_y - _p[0]) > (y - 0.5 * tile_size_pixels)) or ((_y - _p[0]) < (0.5 * tile_size_pixels)):
                nviolations += 1

            if ((_x - _p[1]) > (x - 0.5 * tile_size_pixels)) or ((_x - _p[1]) < (0.5 * tile_size_pixels)):
                nviolations += 1

            if nviolations==0:
                grid.iloc[_i, y_col_ind] -= _p[0]
                grid.iloc[_i, x_col_ind] -= _p[1]
                
    return grid, tile_size_pixels

def plotGrid(grid: pd.DataFrame, x: int, y: int, f: float = 300, object_shape: str = 'spot', size: int = 294, show_spot_labels: bool = True, show: bool = True, savepath: str = None, sname: str = '', verbose: int = 1):

    
    """Plot grid of tiles or spots

    Parameters:
    x: full resolution image width

    y: full resolution image height

    f: figure scaling factor

    object_shape: ['spot', 'square', compatible <matplotlib.patches> object]

    size: tile height and width, or spot diameter

    show_spot_labels: display spot labels

    savepath: directory to save data files

    sname: identifier for saving data files

    verbose: set to 0 to suppress print output

    Output:
    None
    """

    fig, ax = plt.subplots(figsize=(x/f, y/f))
    ax.set_xlim(0, x)
    ax.set_ylim(0, y)
    ax.axis('off')
    ax.set_aspect('equal')
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.plot([0, 0, x, x, 0], [y, 0, 0, y, y], color='k')
    ax.scatter(0, 0, marker='+', s=1000, c='crimson', clip_on=False)
    v = grid[['pxl_row_in_fullres', 'pxl_col_in_fullres']].values

    if object_shape == 'spot':
        ax.add_collection(PatchCollection([Circle((x1, y1), size/2) for y1, x1 in v], alpha=0.9, color='gray'))
    elif object_shape == 'square':
        ax.add_collection(PatchCollection([Rectangle((x1-size/2, y1-size/2), size, size) for y1, x1 in v], alpha=0.9, color='gray'))
    else:
        raise NotImplementedError
    
    if show_spot_labels:
        for i, (y1, x1) in enumerate(v):
            ax.text(x1, y1, i, va='center', ha='center')

    fig.tight_layout()
    
    if not savepath is None:
        fig.savefig(savepath + 'grid_%s.png' % sname, facecolor='w')

    if show:
        plt.show()
    else:
        plt.close(fig)

    return