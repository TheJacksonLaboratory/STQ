import os
import argparse
import pandas as pd
import numpy as np
import openslide
import itertools
import PIL
import json
from tqdm import tqdm
import time

import torch
import torchvision
import torch.nn as nn
from torchvision import transforms
import openslide

from functools import partial

import sys
sys.path.append('/TransPath/')
import moco.builder_infence
import vits

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

models = {'CTransPath': {'path': '/TransPath/ctranspath.pth',    'size': 224},
          'MoCoV3':     {'path': '/TransPath/vit_small.pth.tar', 'size': 224},
          'TransPath':  {'path': '/TransPath/checkpoint.pth',    'size': 256}}

def normalizer(img, mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225), size=224):
    func = transforms.Compose([transforms.Resize(size),
                              transforms.ToTensor(),
                              transforms.Normalize(mean=mean, std=std)])
    return func(img)

import timm
from timm.models.layers.helpers import to_2tuple

class MoCo(nn.Module):
    """
    Build a MoCo model with a base encoder, a momentum encoder, and two MLPs
    https://arxiv.org/abs/1911.05722
    """
    def __init__(self, base_encoder, dim=256, mlp_dim=4096, T=1.0):
        """
        dim: feature dimension (default: 256)
        mlp_dim: hidden dimension in MLPs (default: 4096)
        T: softmax temperature (default: 1.0)
        """
        super(MoCo, self).__init__()

        self.T = T

        # build encoders
        self.base_encoder = base_encoder(num_classes=mlp_dim)
        self.momentum_encoder = base_encoder(num_classes=mlp_dim)

        self._build_projector_and_predictor_mlps(dim, mlp_dim)

        for param_b, param_m in zip(self.base_encoder.parameters(), self.momentum_encoder.parameters()):
            param_m.data.copy_(param_b.data)  # initialize
            param_m.requires_grad = False  # not update by gradient

    def _build_mlp(self, num_layers, input_dim, mlp_dim, output_dim, last_bn=True):
        mlp = []
        for l in range(num_layers):
            dim1 = input_dim if l == 0 else mlp_dim
            dim2 = output_dim if l == num_layers - 1 else mlp_dim

            mlp.append(nn.Linear(dim1, dim2, bias=False))

            if l < num_layers - 1:
                mlp.append(nn.BatchNorm1d(dim2))
                mlp.append(nn.ReLU(inplace=True))
            elif last_bn:
                # follow SimCLR's design: https://github.com/google-research/simclr/blob/master/model_util.py#L157
                # for simplicity, we further removed gamma in BN
                mlp.append(nn.BatchNorm1d(dim2, affine=False))

        return nn.Sequential(*mlp)

    def _build_projector_and_predictor_mlps(self, dim, mlp_dim):
        pass

    @torch.no_grad()
    def _update_momentum_encoder(self, m):
        """Momentum update of the momentum encoder"""
        for param_b, param_m in zip(self.base_encoder.parameters(), self.momentum_encoder.parameters()):
            param_m.data = param_m.data * m + param_b.data * (1. - m)

    def contrastive_loss(self, q, k):
        # normalize
        q = nn.functional.normalize(q, dim=1)
        k = nn.functional.normalize(k, dim=1)
        # gather all targets
        k = concat_all_gather(k)
        # Einstein sum is more intuitive
        logits = torch.einsum('nc,mc->nm', [q, k]) / self.T
        N = logits.shape[0]  # batch size per GPU
        labels = (torch.arange(N, dtype=torch.long) + N * torch.distributed.get_rank()).cuda()
        return nn.CrossEntropyLoss()(logits, labels) * (2 * self.T)

    def forward(self, x1):

        return self.base_encoder(x1)

    def forward_train(self, x1, x2, m):
        """
        Input:
            x1: first views of images
            x2: second views of images
            m: moco momentum
        Output:
            loss
        """

        # compute features
        q1 = self.predictor(self.base_encoder(x1))
        q2 = self.predictor(self.base_encoder(x2))

        with torch.no_grad():  # no gradient
            self._update_momentum_encoder(m)  # update the momentum encoder

            # compute momentum features as targets
            k1 = self.momentum_encoder(x1)
            k2 = self.momentum_encoder(x2)

        return self.contrastive_loss(q1, k2) + self.contrastive_loss(q2, k1)


class MoCo_ResNet(MoCo):
    def _build_projector_and_predictor_mlps(self, dim, mlp_dim):
        hidden_dim = self.base_encoder.fc.weight.shape[1]
        del self.base_encoder.fc, self.momentum_encoder.fc # remove original fc layer

        # projectors
        self.base_encoder.fc = self._build_mlp(2, hidden_dim, mlp_dim, dim)
        self.momentum_encoder.fc = self._build_mlp(2, hidden_dim, mlp_dim, dim)

        # predictor
        self.predictor = self._build_mlp(2, dim, mlp_dim, dim, False)


class MoCo_ViT(MoCo):
    def _build_projector_and_predictor_mlps(self, dim, mlp_dim):
        hidden_dim = self.base_encoder.head.weight.shape[1]
        del self.base_encoder.head, self.momentum_encoder.head # remove original fc layer

        # projectors
        self.base_encoder.head = self._build_mlp(3, hidden_dim, mlp_dim, dim)
        self.momentum_encoder.head = self._build_mlp(3, hidden_dim, mlp_dim, dim)

        # predictor
        self.predictor = self._build_mlp(2, dim, mlp_dim, dim)


# utils
@torch.no_grad()
def concat_all_gather(tensor):
    """
    Performs all_gather operation on the provided tensors.
    *** Warning ***: torch.distributed.all_gather has no gradient.
    """
    tensors_gather = [torch.ones_like(tensor)
        for _ in range(torch.distributed.get_world_size())]
    torch.distributed.all_gather(tensors_gather, tensor, async_op=False)

    output = torch.cat(tensors_gather, dim=0)
    return output

print(torch.__version__)
print(torchvision.__version__)
print(openslide.__version__)
print(timm.__version__)

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute features of each tile')
    parser.add_argument('--wsi-file', dest='wsi_file', action='store',
                        required=True,
                        help="""The path to the whole slide image (WSI) in a format readable by openslide (e.g., svs or ndpi).""")
    parser.add_argument('--model', dest='model', action='store',
                        required=True,
                        help="""Name of the pre-trained model: CTransPath, MoCoV3, or TransPath.""")
    parser.add_argument('--cuda-visible-devices', dest='cuda_visible_devices', action='store', default="", required=False,
                        help="""List of GPUs to use.""")
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

    if args.model == 'CTransPath':
        model = timm.create_model('swin_tiny_patch4_window7_224', embed_layer=ConvStem, pretrained=False)
        model.head = nn.Identity()
        model.load_state_dict(torch.load(models[args.model]['path'])['model'], strict=True)

    elif args.model == 'MoCoV3':
        os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda_visible_devices  
        model = MoCo_ViT(partial(vits.__dict__['vit_small'], stop_grad_conv1=True))
        model = nn.DataParallel(model).cuda()
        model.load_state_dict(torch.load(models['MoCoV3']['path'])['state_dict'], strict=True)

    elif args.model == 'TransPath':
        raise NotImplementedError("Forward call to get positional embedding fails.")
        os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda_visible_devices
        from byol_pytorch.byol_pytorch_get_feature import BYOL
        model = BYOL(image_size=256, hidden_layer='to_latent')
        pretext_model = torch.load(models[args.model]['path'])
        model = nn.DataParallel(model).cuda()
        model.load_state_dict(pretext_model, strict=True)
        model.module.online_encoder.net.head = nn.Identity()
        # featrues = model(imgs, return_embedding=True)[1]
        
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
            images = torch.cat([normalizer(PIL.Image.fromarray(image), size=models[args.model]['size'])[None, :, :, :] for image in images], 0)
            with torch.no_grad():
                temp_features = model(images).cpu().numpy()

                # Average the subtiles, e.g., every 5 subtiles
                if subtiling:
                    df_temp = pd.DataFrame(temp_features)
                    temp_features = df_temp.groupby(np.arange(len(df_temp.index))//len(subcoords)).mean().values

                features.append(temp_features)
        
    features = np.vstack(features)
    
    # Convert the dictionary of features to a dataframe and name its columns featXXX
    df_features = pd.DataFrame(features)
    df_features.columns = [f'feat_{args.model}_' + str(i) for i in range(df_features.shape[1])]
    df_features.index = pos.loc[pos['in_tissue']==1].index
    
    # Append the spot position information to each row
    tbl = pd.concat([pos.loc[pos['in_tissue']==1], df_features], axis=1)
    print(tbl)   
    
    # Output the features with spot information
    ## This will automatically compress if the file suffix is .gz

    tbl.to_csv(output_path + '.tsv.gz', index=False)
    print('Successfully wrote ' + output_path)
    
exit(0)
