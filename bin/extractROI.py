import os
import openslide
import json
import tifffile
from tifffile import TiffFile
import numpy as np
import argparse
import cv2

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--fileslide", type=str, required=True, help="")
    parser.add_argument("--roifile", type=str, required=True, help="")
    parser.add_argument('--wholeside', default=False, action=argparse.BooleanOptionalAction, help="")
    parser.add_argument('--sizefile', type=str, default="size.txt", help="")
    parser.add_argument('--outfile', type=str, default="outfile.tiff", help="")
    parser.add_argument('--extract', type=str, default="False", help="")
    args = parser.parse_args()

    fileslide = args.fileslide.replace("\\", "")
    
    try:
        slide = openslide.open_slide(fileslide)
        dims0 = slide.dimensions
    except Exception as exception:
        print(exception)
        # If the slide is too large openslide may fail to read
        with TiffFile(fileslide) as imgh:
            dims0 = imgh.pages[0].tags[256].value, imgh.pages[0].tags[257].value
            print(dims0)
    
    with open(args.roifile, 'r') as tempfile:
        info = json.load(tempfile)
    
    is_contour = False
    if 'location' in info['0'] and 'location' in info['1'] and 'size' in info['0'] and 'size' in info['1']:
        icoords = int(dims0[0] * info['0']['location']), int(dims0[1] * info['1']['location'])
        size = int(dims0[0] * info['0']['size']), int(dims0[1] * info['1']['size'])
    elif 'points' in info['0'] and 'points' in info['1']:
        is_contour = True
        icoords = int(dims0[0] * min(info['0']['points'])), int(dims0[1] * min(info['1']['points']))
        size = int(dims0[0] * (max(info['0']['points']) - min(info['0']['points']))), int(dims0[1] * (max(info['1']['points']) - min(info['1']['points'])))
    else:
        raise ValueError("Invalid ROI file format")
    print(dims0, '\t', icoords, '\t', size)
    
    if args.wholeside:
        sizegp = round(dims0[0] * dims0[1] / 10**6)
    else:
        sizegp = round(size[0] * size[1] / 10**6)

    with open(args.sizefile, 'w') as tempfile:
        tempfile.write(str(sizegp))

    def filterOutsideContour(img, info, dims0, icoords, q=0.95):
        mask = np.zeros(img.shape[:2], dtype=np.uint8)  # uint8, not bool
        points0 = (np.array(info['0']['points']) * dims0[0] - icoords[0]).astype(np.int32)
        points1 = (np.array(info['1']['points']) * dims0[1] - icoords[1]).astype(np.int32)
        points = np.stack([points0, points1], axis=-1).reshape((-1, 1, 2))  # (N, 1, 2)
        cv2.fillPoly(mask, [points], 255)  # list of contours, fill value 255
        img[mask == 0] = 255  # Set outside contour to white (255)
        t2 = img[::8, ::8, :].copy()
        t2[t2==255] = 0
        qv = np.quantile(t2, q, axis=(0, 1)).astype(np.uint8)
        print(f"Quantile value for outside contour: {qv}")
        img[(img==255).all(axis=-1)] = qv
        return img
        
    if args.extract=="True":
        print('Extracting ROI image')
        try:
            img = slide.read_region(location=icoords, level=0, size=size).convert('RGB')
            tifffile.imwrite(args.outfile, filterOutsideContour(np.array(img), info, dims0, icoords) if is_contour else np.array(img), bigtiff=True)
            img.close()
        except Exception as exception:
            print(exception)
            # If the slide is too large openslide may fail to read
            img = tifffile.imread(fileslide)[icoords[1]:icoords[1]+size[1],icoords[0]:icoords[0]+size[0],:]
            tifffile.imwrite(args.outfile, filterOutsideContour(img, info, dims0) if is_contour else img, bigtiff=True)
            del img

    exit(0)
    