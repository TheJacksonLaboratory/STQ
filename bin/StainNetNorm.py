# Prepared by Domanskyi
# from https://github.com/khtao/StainNet repository
# Here I use their pretrained net to run image stain normalization

import os
import argparse
import imageio
import tifffile
import numpy as np
import torch
import torch.nn as nn
from tqdm import tqdm
from PIL import Image
from torch.utils.data import dataset, DataLoader
from glob import glob

import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

class StainNet(nn.Module):
    def __init__(self, input_nc=3, output_nc=3, n_layer=3, n_channel=32, kernel_size=1):
        super(StainNet, self).__init__()
        model_list = []
        model_list.append(nn.Conv2d(input_nc, n_channel, kernel_size=kernel_size, bias=True, padding=kernel_size // 2))
        model_list.append(nn.ReLU(True))
        for n in range(n_layer - 2):
            model_list.append(
                nn.Conv2d(n_channel, n_channel, kernel_size=kernel_size, bias=True, padding=kernel_size // 2))
            model_list.append(nn.ReLU(True))
        model_list.append(nn.Conv2d(n_channel, output_nc, kernel_size=kernel_size, bias=True, padding=kernel_size // 2))

        self.rgb_trans = nn.Sequential(*model_list)

    def forward(self, x):
        return self.rgb_trans(x)

def list_file_tree(path, file_type="tif"):
    if file_type.find("*") < 0:
        file_type = "*" + file_type
    image_list = glob(os.path.join(path, "*" + file_type), recursive=True)
    return image_list

class SingleImage(dataset.Dataset):
    def __init__(self, data_path, transform=None, augment=None):
        self.data_path = data_path
        self.transform = transform
        self.augment = augment
        self.image_list = list_file_tree(os.path.join(data_path), "png")
        self.image_list += list_file_tree(os.path.join(data_path), "jpg")
        self.image_list += list_file_tree(os.path.join(data_path), "tif")
        self.image_list += list_file_tree(os.path.join(data_path), "tiff")
        self.image_list.sort()

    def __len__(self):
        return len(self.image_list)

    def __getitem__(self, item):
        img = Image.open(self.image_list[item])
        img = (np.array(img, dtype=np.float32) / 255.0).transpose((2, 0, 1))
        return img

def process_images(opt, model, s = 4096):
    dataset = SingleImage(opt.source_dir)
    dataloader = DataLoader(dataset, batch_size=1, num_workers=1, drop_last=False)
    file_list = dataset.image_list
    num = 0
    for imgs in dataloader:
        print(imgs.shape)
        imgs_corrected = imgs.numpy()

        dims = imgs.shape[2], imgs.shape[3]
        r = [np.append(s*np.array(range(0, int(np.floor(dims[i]/s))+1)), [dims[i]]) for i in range(2)]
        coords = []
        for i in range(len(r[0])-1):
            for j in range(len(r[1])-1):
                coords.append((i,j))

        for i, j in tqdm(coords):
            imgs_temp = imgs[:, :, r[0][i]:r[0][i+1], r[1][j]:r[1][j+1]]
            print(imgs_temp.shape)
            if (imgs_temp.shape[2]!=0) and (imgs_temp.shape[3]!=0):
                with torch.no_grad():
                    imgs_temp = imgs_temp.cpu()
                    imgs_temp = (imgs_temp - 0.5) * 2
                    outputs = (model(imgs_temp) * 0.5 + 0.5).clamp(0, 1).detach().cpu().numpy()
        
                imgs_corrected[:, :, r[0][i]:r[0][i+1], r[1][j]:r[1][j+1]] = outputs

        for out in imgs_corrected:
            file_path = file_list[num]
            file_path = os.path.join(os.path.join(opt.save_dir), os.path.split(file_path)[1])
            os.makedirs(os.path.split(file_path)[0], exist_ok=True)
            print('\n', file_path)
            ext = os.path.splitext(file_path)[1]
            tifffile.imwrite(file_path[:-len(ext)] + ".tiff", np.array(np.array(Image.fromarray((out * 255).astype(np.uint8).transpose((1, 2, 0))))), bigtiff=True) # v2            
            #Image.fromarray((out * 255).astype(np.uint8).transpose((1, 2, 0))).save(file_path[:-len(ext)] + ".tiff", compression='raw') # v1
            #imageio.imwrite(file_path[:-len(ext)] + ".tiff", (out * 255).astype(np.uint8).transpose((1, 2, 0))) # v0
            num += 1

    return

def run_normalization(opt):
    model = StainNet()
    model = model.cpu()
    checkpoint = torch.load(opt.model_path, map_location=torch.device('cpu'))
    model.load_state_dict(checkpoint)
    model.eval()
    process_images(opt, model)
    return

if __name__ == '__main__':

    # python norm.py --source_dir "input_images/" --save_dir "output_images/" --model_path "StainNet-Public_layer3_ch32.pth"

    parser = argparse.ArgumentParser()
    parser.add_argument("--source_dir", type=str, required=True, help="path to source images")
    parser.add_argument("--save_dir", type=str, required=True, help="path to save images")
    parser.add_argument('--model_path', type=str, required=True, help='models path to load')
    args = parser.parse_args()

    run_normalization(args)
