import openslide
import json
import tifffile
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

    slide = openslide.open_slide(args.fileslide)
    
    with open(args.roifile, 'r') as tempfile:
        info = json.load(tempfile)
    
    icoords = int(slide.dimensions[0] * info['0']['location']), int(slide.dimensions[1] * info['1']['location'])
    size = int(slide.dimensions[0] * info['0']['size']), int(slide.dimensions[1] * info['1']['size'])
    print(slide.dimensions, '\t', icoords, '\t', size)
    
    with open(args.sizefile, 'w') as tempfile:
        sizegp = round(size[0] * size[1] / 10**6)
        tempfile.write(str(sizegp))
    
    if args.extract=="True":
        print('Extracting ROI image')
        img = slide.read_region(location=icoords, level=0, size=size).convert('RGB')
        tifffile.imwrite(args.outfile, np.array(img), bigtiff=True)
        img.close()

    exit(0)
    