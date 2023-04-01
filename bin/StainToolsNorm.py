# Prepared by Domanskyi
# The normalization is done by patches
# Some stitching lines may be visible

import os
import argparse
import tifffile
import staintools
import numpy as np
from tqdm import tqdm
from PIL import Image

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

import spams

def get_concentrations(I, stain_matrix, regularizer=0.01):
    OD = convert_RGB_to_OD(I).reshape((-1, 3))
    return spams.lasso(X=OD.T, D=stain_matrix.T, mode=2, lambda1=regularizer, pos=True).toarray().T
    
def convert_RGB_to_OD(I):
    mask = (I == 0)
    I[mask] = 1
    return np.maximum(-1 * np.log(I / 255), 1e-6)
    
def convert_OD_to_RGB(OD):
    assert OD.min() >= 0, "Negative optical density."
    OD = np.maximum(OD, 1e-6)
    return (255 * np.exp(-1 * OD)).astype(np.uint8)

class StainNormalizer(staintools.StainNormalizer):

    def __init__(self, method):
        super().__init__(method)

    def estimate(self, I):
        stain_matrix_source = self.extractor.get_stain_matrix(I)
        print(stain_matrix_source.dtype)
        source_concentrations = get_concentrations(I, stain_matrix_source)
        maxC_source = np.percentile(source_concentrations, 99, axis=0).reshape((1, 2))
        return stain_matrix_source, maxC_source

    def transform(self, I, stain_matrix_source, maxC_source):
        source_concentrations = get_concentrations(I, np.array(stain_matrix_source))        
        source_concentrations *= self.maxC_target / maxC_source
        tmp = 255 * np.exp(-1 * np.dot(source_concentrations, self.stain_matrix_target))
        return tmp.reshape(I.shape).astype(np.uint8)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--referenceImagePath", type=str, required=True, help="path to reference or target image")
    parser.add_argument("--inputImagePath", type=str, required=True, help="input image name")
    parser.add_argument("--outputImageName", type=str, required=True, help="output image name")
    parser.add_argument('--s', type=int, default=4096, help='patch size')
    parser.add_argument('--low', type=int, default=100, help='low threshold')
    parser.add_argument('--high', type=int, default=200, help='high threshold')
    parser.add_argument('--qfraction', type=float, default=0.75, help='quantile of fraction for tissue')
    args = parser.parse_args()
    
    print('s:', args.s)

    target = tifffile.imread(args.referenceImagePath)[:,:,:3][:,:,::-1]
    normalizer = StainNormalizer(method='macenko')
    normalizer.fit(target)
    
    img = tifffile.imread(args.inputImagePath)[:,:,:3][:,:,::-1]
    dims = img.shape[0], img.shape[1]
    r = [np.append(args.s*np.array(range(0, int(np.floor(dims[i]/args.s))+1)), [dims[i]]) for i in range(2)]
    coords = [(i,j) for i in range(len(r[0])-1) for j in range(len(r[1])-1)]
    print(coords)

    # Determine representative patch
    coordsf = []
    fractions = []
    for i, j in tqdm(coords):
        try:
            # Get in_tissue flags for patch
            v = img[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1], :].mean(axis=2)
            vc = v.copy()
            v[vc < args.low] = 0
            v[vc > args.high] = 0
            v[(vc >= args.low) & (vc <= args.high)] = 1
            f = v.ravel().mean()
            print(i, j, f)
            if f==f:
                coordsf.append((i, j, f))
                fractions.append(f)
        except Exception as exception:
            print('Exception:', exception)
        
    fcutoff = np.quantile(fractions, args.qfraction)
    fcutoff = min(fcutoff, 0.5)
    print('fcutoff:', fcutoff)
    
    ms = []
    cs = []
    for i, j, f in tqdm(coordsf):
        # Get in_tissue flags for patch
        in_tissue = f >= fcutoff
        print(i, j, in_tissue, f)
        if in_tissue:
            try:
                m, c = normalizer.estimate(img[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1], :])
                ms.append(m)
                cs.append(c)
            except Exception as exception:
                print('Exception:', exception)

    ms = np.dstack(ms)
    msm = np.median(ms, axis=2)
    cs = np.vstack(cs)
    csm = np.median(cs, axis=0)
    print(ms.shape, cs.shape)
    closest = np.argsort(np.sqrt((np.array([(ms - msm[:, :, None])[:, :, i].ravel() for i in range(ms.shape[2])])**2).sum(axis=1)))[0]
    m, c = ms[:, :, closest], cs[closest, :]
    print(m)
    print(c)

    # Normalize all patches
    for i, j in tqdm(coords):
        try:
            img[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1], :] = normalizer.transform(img[r[0][i]:r[0][i+1], r[1][j]:r[1][j+1], :], m, c)[:,:,::-1]
        except Exception as exception:
            print('Exception:', exception)
    
    tifffile.imwrite(args.outputImageName, img, bigtiff=True)
    
    exit(0)
