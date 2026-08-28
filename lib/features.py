import os
import argparse
import json
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import openslide
import PIL.Image
import torch
from tqdm import tqdm

from forfeatures import load_positions, resolve_expansion_and_downsample, _try_read_tile, save_features
from model_loaders import MODEL_CHOICES, GPU_ONLY_MODELS, get_loader

PIL.Image.MAX_IMAGE_PIXELS = None

def parse_args():
    p = argparse.ArgumentParser(description='Compute features of each tile')
    p.add_argument('--wsi-file', dest='wsi_file', required=True,
                    help="Path to the WSI file (e.g., svs or ndpi), readable by openslide.")
    p.add_argument('--model-name', dest='model_name', required=True, choices=MODEL_CHOICES,
                    help="Name of the pre-trained model to use: " + ', '.join(MODEL_CHOICES))
    p.add_argument('--model-checkpoint-path', dest='modelPath', required=True,
                    help="Path to the model weights. For timm-based models "
                         "(uni, uni2, hoptimus0, hoptimus1, h0mini, "
                         "ctranspath, mocov3, conch, virchow, virchow2, "
                         "provgigapath, provgigapathflash, inception) this is "
                         "a single weight file (.bin/.pth/.safetensors). For "
                         "transformers-based models (phikon, phikon2, "
                         "kaikomidnight12k) this is either the directory "
                         "containing config.json / model.safetensors "
                         "(or pytorch_model.bin) / preprocessor_config.json, "
                         "or a path to one of those weight files inside it.")
    p.add_argument('--positions-list-file', dest='positions_list_file', required=True,
                    help="spaceranger positions_list.csv: one row per spot, tissue flag + x/y pixel coords.")
    p.add_argument('--scalefactors-json-file', dest='scalefactors_json_file', required=True,
                    help="spaceranger scalefactors_json.json defining the spot diameter in full resolution.")
    p.add_argument('--output-path', dest='output_path', required=True,
                    help="Name of output CSV file (rows=tiles, cols=features). Compressed if named *.gz")
    p.add_argument('--model-suffix', dest='modelSuffix', default=None, required=False,
                    help="Suffix used to name feature columns. Defaults to --model-name.")
    p.add_argument('--tile-mask', dest='tile_mask', default=None, required=False,
                    help="Path to tile mask file.")
    p.add_argument('--downsample-expanded', dest='downsample', default=True, required=False,
                    help="If expansion factor > 1, downsample tiles back to the input size.")
    p.add_argument('--expansion-factor', dest='expansion', required=True,
                    help="Expansion factor, 1 means no expansion.")
    p.add_argument('--cuda-visible-devices', dest='cuda_visible_devices', default="", required=False,
                    help="List of GPUs to use.")
    p.add_argument('--num-workers', dest='num_workers', type=int, default=8, required=False,
                    help="Number of threads used to read tiles from the slide concurrently.")
    p.add_argument('--batch-size', dest='batch_size_sf', type=int, default=6*10**7, required=False,
                    help="Number of tiles to process in a single batch")
    p.add_argument('--use-conch-normalizer', dest='useCONCHnormalizer', default=False, required=False,
                    help="For --model-name=conch, use CONCH's own normalizer instead of the default imagenet one.")
    return p.parse_args()

def load_model(model_name, model_path, use_conch_normalizer=False):
    """Build the requested model, its input transform, and its forward/
    postprocess functions.

    Returns (model, destination, transform, forward_fn, postprocess_fn).
    forward_fn(model, batch) returns the raw (not-yet-pooled-or-sliced)
    feature tensor for a batch of images; postprocess_fn(np.ndarray) turns
    that (already moved to CPU/numpy) into the final (batch, num_features)
    array. Models that already return a plain pooled vector use an identity
    postprocess_fn.

    Per-model details (architecture kwargs, custom nn.Module subclasses,
    extra package imports) live in model_loaders/<model_name>.py so that,
    e.g., importing `timm` for one model can't clash with the `vits` module
    mocov3 needs on sys.path, and models needing an extra package (conch,
    prov-gigapath-flash) only require it when actually requested.
    """
    if model_name in GPU_ONLY_MODELS and not torch.cuda.is_available():
        raise RuntimeError(f"Model '{model_name}' requires a GPU, but none is available.")

    destination = "cuda" if torch.cuda.is_available() else "cpu"
    loader_fn = get_loader(model_name)
    result = loader_fn(model_path, destination, use_conch_normalizer=use_conch_normalizer)

    return (result['model'], destination, result['transform'],
            result['forward_fn'], result.get('postprocess_fn'))


def extract_features(slide, pos, model, destination, w, h, expansion, downsample,
                      batch_size, num_batches, transform, forward_fn, postprocess_fn,
                      num_workers=8):
    if postprocess_fn is None:
        postprocess_fn = lambda x: x

    features = []
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        for ibatch in tqdm(range(num_batches)):
            # executor.map preserves input order, so results come back in the same
            # row_idx-ascending order as the original sequential loop.
            tasks = ((slide, pos, indx + ibatch * batch_size, w, h, expansion, downsample)
                      for indx in range(batch_size))
            images = [img for img in executor.map(_try_read_tile, tasks) if img is not None]
            print('Number of tiles:', len(images))

            if len(images) > 0:
                batch = torch.stack([transform(PIL.Image.fromarray(im)) for im in images], 0)
                if destination == 'cuda':
                    batch = batch.to(destination, non_blocking=True)
                with torch.inference_mode():
                    features_ = forward_fn(model, batch).cpu().numpy()
                    features_ = postprocess_fn(features_)
                    features.append(features_)
    return np.vstack(features)


def extract():
    args = parse_args()
    os.environ["CUDA_VISIBLE_DEVICES"] = args.cuda_visible_devices

    expansion = float(args.expansion)
    downsample = args.downsample == 'true'
    expansion, downsample = resolve_expansion_and_downsample(expansion, downsample)

    pos, _ = load_positions(args.positions_list_file, args.tile_mask)

    with open(args.scalefactors_json_file) as f:
        scalefactors_tbl = json.load(f)
    spot_diameter_fullres = scalefactors_tbl['spot_diameter_fullres']
    spot_diameter_wsi = round(spot_diameter_fullres * 1)  # scale_factor is always 1

    if downsample:
        num_rows = num_cols = round(spot_diameter_wsi)
    else:
        num_rows = num_cols = round(spot_diameter_wsi * expansion)

    use_conch_normalizer = str(args.useCONCHnormalizer) == 'true'
    model, destination, transform, forward_fn, postprocess_fn = load_model(
        args.model_name, args.modelPath, use_conch_normalizer=use_conch_normalizer)

    num_images = len(pos)
    batch_size = int(int(args.batch_size_sf) / (expansion * expansion * num_cols * num_rows))
    num_batches = int(np.ceil(num_images / batch_size))

    print('Reading and processing tiles:', num_images)
    print('Batch size:', batch_size)
    print('Number of batches:', num_batches)

    slide = openslide.open_slide(args.wsi_file)

    features = extract_features(slide, pos, model, destination,
                                 num_cols, num_rows, expansion, downsample,
                                 batch_size, num_batches, transform, forward_fn,
                                 postprocess_fn, num_workers=args.num_workers)

    model_suffix = args.modelSuffix if args.modelSuffix else args.model_name
    save_features(features, pos, model_suffix, args.output_path)


if __name__ == '__main__':

    extract()
