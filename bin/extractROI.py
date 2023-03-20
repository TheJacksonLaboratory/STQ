import os
import openslide
import json
import tifffile
from tifffile import TiffFile
import numpy as np
import argparse

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--fileslide", type=str, required=True, help="")
    parser.add_argument("--roifile", type=str, required=True, help="")
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
        temptif = TiffFile(fileslide)
        dims0 = temptif.pages[0].shape
        temptif.close()
    
    with open(args.roifile, 'r') as tempfile:
        info = json.load(tempfile)
    
    icoords = int(dims0[0] * info['0']['location']), int(dims0[1] * info['1']['location'])
    size = int(dims0[0] * info['0']['size']), int(dims0[1] * info['1']['size'])
    print(dims0, '\t', icoords, '\t', size)
    
    with open(args.sizefile, 'w') as tempfile:
        sizegp = round(size[0] * size[1] / 10**6)
        tempfile.write(str(sizegp))
    
    if args.extract=="True":
        print('Extracting ROI image')
        try:
            img = slide.read_region(location=icoords, level=0, size=size).convert('RGB')
            tifffile.imwrite(args.outfile, np.array(img), bigtiff=True)
            img.close()
        except Exception as exception:
            print(exception)
            # If the slide is too large openslide may fail to read
            img = tifffile.imread(fileslide)[icoords[1]:icoords[1]+size[1],icoords[0]:icoords[0]+size[0],:]
            tifffile.imwrite(args.outfile, img, bigtiff=True)
            del img

    exit(0)
    