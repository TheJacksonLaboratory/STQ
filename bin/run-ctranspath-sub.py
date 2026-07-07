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

try:
    from timm.models.layers.helpers import to_2tuple  # older timm versions
except ModuleNotFoundError:
    from timm.layers import to_2tuple  # timm >= ~0.9

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

## here we normalize the image by mean RBG channels, std dev RGB channels, and resizes the pixels
## using the PyTorch package. it makes a function with parameters
## ToTensor renders the image and scales pixel values
def normalizer(img, mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225), size=224):
    func = transforms.Compose([transforms.Resize(size),
                              transforms.ToTensor(),
                              transforms.Normalize(mean=mean, std=std)])
    return func(img)


## CTransPath's custom convolutional patch embedding (replaces the standard
## linear patch projection used by plain ViTs). Copied from run-ctranspath.py
## so this script has no dependency on the /TransPath/ repo.
class ConvStem(nn.Module):
    def __init__(self, img_size=224, patch_size=4, in_chans=3, embed_dim=768, norm_layer=None, flatten=True):
        super().__init__()

        assert patch_size == 4
        assert embed_dim % 8 == 0

        img_size = to_2tuple(img_size)
        patch_size = to_2tuple(patch_size)
        self.img_size = img_size
        self.patch_size = patch_size
        self.grid_size = (img_size[0] // patch_size[0], img_size[1] // patch_size[1])
        self.num_patches = self.grid_size[0] * self.grid_size[1]
        self.flatten = flatten

        stem = []
        input_dim, output_dim = 3, embed_dim // 8
        for l in range(2):
            stem.append(nn.Conv2d(input_dim, output_dim, kernel_size=3, stride=2, padding=1, bias=False))
            stem.append(nn.BatchNorm2d(output_dim))
            stem.append(nn.ReLU(inplace=True))
            input_dim = output_dim
            output_dim *= 2
        stem.append(nn.Conv2d(input_dim, embed_dim, kernel_size=1))
        self.proj = nn.Sequential(*stem)
        self.norm = norm_layer(embed_dim) if norm_layer else nn.Identity()

    def forward(self, x):
        B, C, H, W = x.shape
        assert H == self.img_size[0] and W == self.img_size[1], \
            f"Input image size ({H}*{W}) doesn't match model ({self.img_size[0]}*{self.img_size[1]})."
        x = self.proj(x)
        if self.flatten:
            x = x.flatten(2).transpose(1, 2)  # BCHW -> BNC
        x = self.norm(x)
        return x


## this is a command line argument parser using the argparse library
## these are all the things we ask to input in the CLI, destination, where it's stores, if it's required
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute per-cell spatial features of each tile using CTransPath')
    parser.add_argument('--wsi-file', dest='wsi_file', action='store',
                        required=True,
                        help="""The path to the whole slide image (WSI) in a format readable by openslide (e.g., svs or ndpi).""")
    parser.add_argument('--model-checkpoint-path', dest='modelPath', action='store',
                        required=True,
                        help="""Path to CTransPath checkpoint (e.g. ctranspath.pth). Expected to contain a 'model' key.""")
    parser.add_argument('--positions-list-file', dest='positions_list_file', action='store',
                        required=True,
                        help="""The positions_list.csv file output by spaceranger that has one row per spot and columns indicating whether the spot is within the tissue and its x and y coordinates in pixels.""")
    parser.add_argument('--scalefactors-json-file', dest='scalefactors_json_file', action='store',
                        required=True,
                        help="""The scalefactors_json.json file output by spaceranger that defines the spot diameter in spaceranger's full resolution (i.e., the resolution of the file input to spaceranger, which may or may not be wsi_file).""")
    parser.add_argument('--output-path', dest='output_path', action='store',
                        required=True,
                        help="""Prefix for the output _cell_coords.parquet file.""")
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
    parser.add_argument('--segmentation', dest='segmentation', action='store',
                        required=True,
                        help="""Path to a .npy cell-segmentation mask, same resolution as wsi_file.""")
    parser.add_argument('--upscale-factor', dest='upscale_factor', action='store',
                        required=True,
                        help="""Factor for upscaling the cell patches""")

    ## this extracts the inputted data and puts it in args
    args = parser.parse_args()
    ## here we convert args.expansion into a float, and rename
    expansion = float(args.expansion)
    ## make sure the downsample and subtiling are true before storing
    downsample = args.downsample=='true'
    subtiling = args.subtiling=='true'

    ## make the subcordsf an integer and store it
    subcoordsf = int(args.subcoordsf)
    ## convert the json string into a python argument
    subcoords = json.loads(args.subcoords)

    ## whether or not the expansion factor is 1 tells us if we
    ## can downsample the tiles
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

    ## rename the wsi_file, positions_list_file, output_path, and scale_factors_json_file from args
    wsi_file = args.wsi_file
    positions_list_file = args.positions_list_file
    scalefactors_json_file = args.scalefactors_json_file
    output_path = args.output_path
    # Read in the spaceranger positions list file
    pos = pd.read_csv(positions_list_file, header=None)
    pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']

    ## if there is a tile_mask entered, print so
    ## read the mask, as it's a csv file
    ## align the mask values with the barcode values, then add that data to the in_tissue column
    if args.tile_mask != 'None':
        print('Received tile mask %s' % args.tile_mask)
        mask = pd.read_csv(args.tile_mask, index_col=0, header=None)
        pos['in_tissue'] = mask.reindex(pos['barcode'].values).values

    # Read the spot diameter at spaceranger's "full resolution" from the scalefactors_json file
    with open(scalefactors_json_file) as f:
        scalefactors_tbl = json.load(f)
    spot_diameter_fullres = scalefactors_tbl['spot_diameter_fullres']

    scale_factor = 1
    # Define the spot diameter in the resolution of the wsi_file
    spot_diameter_wsi = round(spot_diameter_fullres * scale_factor)
    # Translate the pixel coordinates from full resolution to the resolution of the wsi
    pos['pxl_row_in_wsi'] = pos.pxl_row_in_fullres * scale_factor
    pos['pxl_col_in_wsi'] = pos.pxl_col_in_fullres * scale_factor

    ## adjust the num rows and columns based on if downsampling or not
    if downsample:
        num_rows = num_cols = round(spot_diameter_wsi)
    else:
        num_rows = num_cols = round(spot_diameter_wsi * expansion)

    # Load pre-trained CTransPath model (Swin-Tiny backbone with a custom
    # convolutional patch embedding). model.head is replaced with Identity
    # so the classifier weights aren't required, but note that model(x) is
    # NOT used for feature extraction below — see the comment at the
    # extraction step for why.
    model = timm.create_model('swin_tiny_patch4_window7_224', embed_layer=ConvStem, pretrained=False)
    model.head = nn.Identity()
    model.load_state_dict(torch.load('/TransPath/ctranspath.pth', map_location="cpu")['model'], strict=True)
    model.eval()

    ## define the number of images as length of position file
    ## define the batch size with how much memory we have
    ## but if subtiling, then divide the batch size by 5
    ## also determine the number of batches
    num_images = len(pos)
    # NB: CTransPath (Swin-Tiny) is far smaller than UNI2 (ViT-giant), so this
    # numerator is set much higher than in run-uni2-sub.py. Tune for your GPU.
    batch_size = int(10**8 / (float(args.expansion) * float(args.expansion) * num_cols * num_rows))
    if subtiling:
        batch_size = int(batch_size / 5)
    num_batches = int(np.ceil(num_images / batch_size))

    print('Reading and pocessing tiles:', num_images)
    print('Batch size:', batch_size)
    print('Number of batches:', num_batches)

    ## open the slide with openslide library
    slide = openslide.open_slide(wsi_file)
    mfull = np.load(args.segmentation)

    border = 48
    A = 7  # CTransPath's Swin-Tiny backbone outputs a native 7x7 spatial
           # token grid at 224x224 input, one token per mosaic cell, so no
           # pooling/aggregation step is needed (unlike UNI2's 16x16 -> 8x8).
    upscale_factor = float(args.upscale_factor)

    ## set zoom level to highest resolution
    w = num_cols
    h = num_rows
    lvl = 0
    features = []
    ## loop through the batches, add progress bar
    cell_coords = {}
    mosaic_ids = []
    for ibatch in tqdm(range(num_batches)):
        images = []
        for indx in range(batch_size):
            if indx + ibatch*batch_size >= pos.shape[0]:
                continue
            try:
                ## cy and cx are the pixel coordinates
                cy = pos.loc[indx + ibatch*batch_size, 'pxl_row_in_wsi']
                cx = pos.loc[indx + ibatch*batch_size, 'pxl_col_in_wsi']

                ## if it's in the tissue, then continue
                if pos.loc[indx + ibatch*batch_size, 'in_tissue']:
                    ## compute expanded extraction size
                    if downsample:
                        ew = round(w * expansion)
                        eh = round(h * expansion)
                    else:
                        ew = w
                        eh = h
                    # TODO: temporarily increase ew and eh, find a more efficient way
                    bew = ew + border * 2
                    beh = eh + border * 2
                    bimg = np.array(slide.read_region((int(cx - bew / 2), int(cy - beh / 2)), lvl, (int(bew), int(beh))).convert('RGB'))
                    m = mfull[int(cy - eh / 2):int(cy + eh / 2), int(cx - ew / 2):int(cx + ew / 2)]

                    ## Make cell mosaic
                    ucells = np.unique(m)
                    cell_centers = {}
                    for cell in ucells:
                        cell_mask = (m == cell)
                        if np.sum(cell_mask) == 0:
                            continue
                        y_coords, x_coords = np.where(cell_mask)
                        center_y = int(np.mean(y_coords))
                        center_x = int(np.mean(x_coords))
                        cell_centers[cell] = (center_y, center_x)

                    mosaic = np.zeros((ew, eh, 3), dtype=np.uint8)
                    max_cells = A**2
                    sel_cells = dict()
                    hs = int(0.5 * ew / A)
                    cell_size = 2 * hs
                    img_h, img_w = bimg.shape[:2]

                    for i in range(A):
                        for j in range(A):
                            center_y = int((i + 0.5) * cell_size)
                            center_x = int((j + 0.5) * cell_size)
                            closest_cell = None
                            closest_distance = float('inf')
                            for cell, (cyc, cxc) in cell_centers.items():
                                #skip cells too close to image borders
                                if (cyc - hs + border < 0 or cyc + hs + border > img_h or
                                        cxc - hs + border < 0 or cxc + hs + border > img_w):
                                    continue
                                distance = np.sqrt(((cyc - center_y)) ** 2 + (cxc - center_x) ** 2)
                                if distance < closest_distance:
                                    closest_distance = distance
                                    closest_cell = cell
                            if closest_cell is not None:
                                cell_center_y, cell_center_x = cell_centers[closest_cell]
                                ## adjusted to account for boundary
                                patch = bimg[cell_center_y - hs + border: cell_center_y + hs + border,
                                             cell_center_x - hs + border: cell_center_x + hs + border]

                                # upscale the patch
                                upscaled_size = int(cell_size * upscale_factor)
                                upscaled = np.array(PIL.Image.fromarray(patch).resize((upscaled_size, upscaled_size), PIL.Image.BILINEAR))

                                # trim back to cell_size from the center
                                trim_start = (upscaled_size - cell_size) // 2
                                trimmed = upscaled[trim_start:trim_start + cell_size, trim_start:trim_start + cell_size]

                                sel_cells[closest_cell] = trimmed
                                try:
                                    mosaic[center_y - hs:center_y + hs, center_x - hs:center_x + hs] = sel_cells[closest_cell]
                                except Exception as e:
                                    print(f'Error filling mosaic for cell {closest_cell}: {e}')

                                mosaic_num = indx + (ibatch * batch_size)
                                cell_coords[(mosaic_num, i, j)] = {"x": int(cell_center_x + cx - ew / 2),
                                                                   "y": int(cell_center_y + cy - eh / 2)}
                    images.append(mosaic)
                    mosaic_ids.append(indx + (ibatch * batch_size))

            except Exception as exception:
                print(exception)
                pass
        print('Number of tiles:', len(images))

        if len(images) > 0:
            timages = torch.cat([normalizer(PIL.Image.fromarray(image))[None, :, :, :] for image in images], 0)
            print('Extracting features')
            with torch.no_grad():
                # IMPORTANT: do NOT call model(x), and don't rely on
                # forward_features() either — behavior of both varies across
                # timm versions and can silently return an already-pooled
                # (B, 768) vector instead of the spatial grid (confirmed to
                # happen in this environment). Calling the submodules
                # directly is version-stable and matches validated testing:
                # patch_embed -> layers -> norm gives the pre-pool spatial
                # patch features. No cls/register tokens to strip here
                # (CTransPath has none), and A already equals the native
                # grid size, so no pooling/aggregation is needed.
                patch_features = model.norm(model.layers(model.patch_embed(timages))).cpu().numpy()  # (B, 49, 768)
                print('Patch features shape:', patch_features.shape)

                Bn, N, C = patch_features.shape
                h_ = w_ = int(round(N ** 0.5))
                assert h_ * w_ == N, f'Unexpected token count {N}, expected a perfect square'
                spatial_features = patch_features.reshape(Bn, h_, w_, C).transpose(0, 3, 1, 2)  # -> (B, 768, 7, 7) = (B, C, H, W)

            features.append(spatial_features)
        print('Processed batch', ibatch + 1, 'of', num_batches)

    features = np.vstack(features)
    mosaic_ids = np.array(mosaic_ids)

    print(features.shape)  # (X, 768, 7, 7)

    mid_to_imid = {mid: imid for imid, mid in enumerate(mosaic_ids)}

    # collect only lightweight scalars/tuples in loop order — no feature slicing here
    records = []
    for imid, mid in enumerate(mosaic_ids):
        for i in range(A):
            for j in range(A):
                key = (mid, i, j)
                coords = cell_coords.get(key)
                if coords is not None:
                    records.append((coords['x'], coords['y'], imid, i, j))

    # build a small dataframe of indices, then dedup keeping the LAST occurrence
    # (this replicates dict[(x,y)] = ... overwrite semantics)
    idx_df = pd.DataFrame(records, columns=['x', 'y', 'imid', 'i', 'j'])
    idx_df = idx_df.drop_duplicates(subset=['x', 'y'], keep='last')

    # single vectorized fancy-index pull for all surviving cells
    feat_arr = features[idx_df['imid'].to_numpy(),
                        :,
                        idx_df['i'].to_numpy(),
                        idx_df['j'].to_numpy()]  # shape (N_unique, 768)

    df = pd.DataFrame(feat_arr, columns=[f'feat_{k}' for k in range(feat_arr.shape[1])])
    df.insert(0, 'y', idx_df['y'].to_numpy())
    df.insert(0, 'x', idx_df['x'].to_numpy())
    df = df.set_index(['x', 'y']).sort_index().reset_index()

    df.to_parquet(args.output_path + '_cell_coords.parquet', index=False)
    print('Successfully wrote ' + args.output_path + '_cell_coords.parquet')

exit(0)