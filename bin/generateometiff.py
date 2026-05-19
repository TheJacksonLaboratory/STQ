# Author: Kevin Anderson and Sergii Domanskyi
# Organization: The Jackson Laboratory
# Date: 2026-05-19
# Description: This module generates a pyramidal OME-TIFF file from a multichannel TIFF image
# (MIF) or an H&E RGB TIFF image, incorporating physical pixel size and optional metadata.
# Pixel size (mpp, or microns per pixel) is inferred from the input TIFF if not provided.
# The script auto-detects whether the image is H&E (RGB/RGBA, axes YXS) or MIF (multi-channel
# grayscale, axes CYX) and applies the correct photometric interpretation and thumbnail method.
#
# H&E images are kept in interleaved YXS layout (channels last), which is required for JPEG
# compression. MIF images use CYX layout (channels first).
#
# Example (H&E):
# python generatetiff.py \
#   --inputImagePath input.tiff \
#   --outputImagePath output.ome.tif \
#   --compression jpeg \
#   --jpegQuality 90 \
#   --tileSize 512
#
# Example (MIF):
# python generatetiff.py \
#   --inputImagePath input.tiff \
#   --outputImagePath output.ome.tif \
#   --tileSize 512
#
# License: MIT License
# https://github.com/TheJacksonLaboratory/spatial-omics-tools

import argparse
import os
import numpy as np
import tifffile


def detect_image_type(img):
    """
    Detect whether the image is H&E (RGB/RGBA) or MIF (multi-channel grayscale).

    For H&E the input is expected to already be in YXS layout (channels last,
    S == 3 or 4). MIF images have channels first (CYX) or are single-plane (YX).

    Returns:
        'he'  — shape YXS, last dim is 3 or 4
        'mif' — shape CYX or YX
    """
    if img.ndim == 3 and img.shape[2] in (3, 4):
        return 'he'
    elif img.ndim == 3:
        return 'mif'
    elif img.ndim == 2:
        return 'mif'
    else:
        raise ValueError(
            f"Unexpected image shape {img.shape}. "
            "Expected YXS (H&E) or CYX / YX (MIF)."
        )


def make_thumbnail_he(img, factor, dtype):
    """
    Build an interleaved RGB uint8 thumbnail for H&E images (shape YXS).
    Alpha channel is dropped if present.
    """
    thumb = img[::factor, ::factor, :3]   # YXS, drop alpha
    if dtype == np.uint8:
        return thumb.astype(np.uint8)
    elif dtype == np.uint16:
        return (thumb >> 8).astype(np.uint8)
    else:
        raise TypeError(f"Unsupported H&E thumbnail dtype: {dtype}")


def make_thumbnail_mif(img, factor, dtype):
    """
    Build a single-channel uint8 thumbnail for MIF images (shape CYX).
    Uses channel 0 (typically DAPI or the first marker).
    """
    thumb = img[0, ::factor, ::factor]
    if dtype == np.uint8:
        return thumb.astype(np.uint8)
    elif dtype == np.uint16:
        # Assumes 10-bit data packed in uint16; shift by 2 to get uint8
        return (thumb >> 2).astype(np.uint8)
    else:
        raise TypeError(f"Unsupported MIF thumbnail dtype: {dtype}")


def compute_subresolutions(img, image_type, tile_size, pyramid_scale):
    """
    Determine how many sub-resolution pyramid levels to write.
    Extracts spatial dimensions correctly for both YXS and CYX layouts.
    """
    if image_type == 'he':
        h, w = img.shape[0], img.shape[1]   # YXS — spatial dims are 0 and 1
    else:
        h, w = img.shape[-2], img.shape[-1]  # CYX or YX — spatial dims are last two

    min_dim = min(h, w)
    subresolutions = 0
    while min_dim / (pyramid_scale ** subresolutions) >= tile_size:
        subresolutions += 1
    return max(subresolutions - 1, 0)


def downsample(img, mag, image_type):
    """Return a downsampled view of `img` by integer factor `mag`."""
    if image_type == 'he':
        return img[::mag, ::mag, :]       # YXS — downsample spatial dims only
    else:
        return img[..., ::mag, ::mag]     # CYX or YX


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Generate a pyramidal OME-TIFF from an H&E or MIF TIFF image.'
    )
    parser.add_argument('--inputImagePath',   type=str,   required=True,
                        help='Path to input TIFF image')
    parser.add_argument('--outputImagePath',  type=str,   required=True,
                        help='Path to output OME-TIFF image (must end with .ome.tif)')
    parser.add_argument('--tileSize',         type=int,   default=512,
                        help='Tile size for OME-TIFF (default: 512)')
    parser.add_argument('--mpp',              type=float, default=None,
                        help='Physical pixel size in micrometers (inferred from TIFF tags if omitted)')
    parser.add_argument('--compression',      type=str,   default='zlib',
                        help='Compression method: zlib, lzw, jpeg, etc. (default: zlib)')
    parser.add_argument('--compressionLevel', type=int,   default=8,
                        help='zlib/lzw compression level 0-9 (default: 8); ignored for JPEG')
    parser.add_argument('--jpegQuality',      type=int,   default=90,
                        help='JPEG quality 0-100 (default: 90); only used when --compression jpeg')
    parser.add_argument('--thumbnailFactor',  type=int,   default=8,
                        help='Thumbnail downsampling factor (default: 8)')
    parser.add_argument('--pyramidScale',     type=int,   default=2,
                        help='Pyramid downsampling factor between levels (default: 2)')
    args = parser.parse_args()

    # ── Validate arguments ────────────────────────────────────────────────────
    if not args.outputImagePath.endswith('.ome.tif') and not args.outputImagePath.endswith('.ome.tiff'):
        raise ValueError("Output path must end with '.ome.tiff' or '.ome.tif'")

    if not os.path.exists(args.inputImagePath):
        raise FileNotFoundError(f"File not found: {args.inputImagePath}")

    # ── Load image ────────────────────────────────────────────────────────────
    img = tifffile.imread(args.inputImagePath)
    print(f'Image shape: {img.shape}  dtype: {img.dtype}')

    image_type = detect_image_type(img)
    print(f'Detected image type: {"H&E (RGB, YXS)" if image_type == "he" else "MIF (multi-channel, CYX)"}')

    # ── Set photometric and axes based on image type ──────────────────────────
    # H&E stays in YXS (channels last, interleaved) — required for JPEG compression.
    # MIF stays in CYX (channels first).
    if image_type == 'he':
        photometric = 'rgb'
        axes        = 'YXS'
    else:
        photometric = 'minisblack'
        axes        = 'CYX' if img.ndim == 3 else 'YX'

    # ── Compute pyramid levels ────────────────────────────────────────────────
    subresolutions = compute_subresolutions(
        img, image_type, args.tileSize, args.pyramidScale
    )
    print(f'Number of pyramid levels: {subresolutions + 1}')

    # ── Resolve physical pixel size ───────────────────────────────────────────
    with tifffile.TiffFile(args.inputImagePath) as tif:
        page = list(tif.series[0])[0]

    if args.mpp is None:
        try:
            x_res = page.tags['XResolution'].value
            mpp   = 1e4 / (x_res[0] / x_res[1])
            print(f'Inferred physical pixel size: {mpp:.4f} µm')
        except (KeyError, ZeroDivisionError):
            mpp = 1.0
            print('Warning: XResolution tag not found; defaulting to mpp=1.0 µm')
    else:
        mpp = args.mpp
        print(f'Using provided physical pixel size: {mpp:.4f} µm')

    # ── Build compression args ────────────────────────────────────────────────
    # imagecodecs uses 'level' for both zlib and JPEG (where level == quality 0-100).
    if args.compression.lower() in ('jpeg', 'jpg'):
        if image_type == 'mif':
            raise ValueError(
                "JPEG compression is not supported for MIF images (uint16 multi-channel). "
                "Use --compression zlib or --compression lzw instead."
            )
        compression_args = {'level': args.jpegQuality}
    else:
        compression_args = {'level': args.compressionLevel}

    # ── Build write options ───────────────────────────────────────────────────
    options = dict(
        photometric=photometric,
        tile=(args.tileSize, args.tileSize),
        compression=args.compression,
        compressionargs=compression_args,
        resolutionunit='CENTIMETER',
    )

    # JPEG requires interleaved (contig) layout; data is already YXS so this
    # is consistent. Setting explicitly avoids tifffile guessing wrong.
    if image_type == 'he':
        options['planarconfig'] = 'contig'

    metadata = {
        'axes': axes,
        'PhysicalSizeX': mpp,
        'PhysicalSizeXUnit': 'µm',
        'PhysicalSizeY': mpp,
        'PhysicalSizeYUnit': 'µm',
    }

    # ── Write OME-TIFF ────────────────────────────────────────────────────────
    print('Generating OME-TIFF image…')

    with tifffile.TiffWriter(args.outputImagePath, bigtiff=True) as tif:

        # Full-resolution base level
        tif.write(
            data=img,
            subifds=subresolutions,
            resolution=(1e4 / mpp, 1e4 / mpp),
            metadata=metadata,
            **options,
        )

        # Sub-resolution pyramid levels
        for level in range(subresolutions):
            mag = args.pyramidScale ** (level + 1)
            tif.write(
                downsample(img, mag, image_type),
                subfiletype=1,
                resolution=(1e4 / mag / mpp, 1e4 / mag / mpp),
                **options,
            )

        # Thumbnail — interleaved RGB YXS for H&E, single-channel for MIF
        if image_type == 'he':
            thumbnail = make_thumbnail_he(img, args.thumbnailFactor, img.dtype.type)
            tif.write(thumbnail, photometric='rgb', metadata={'Name': 'thumbnail'})
        else:
            thumbnail = make_thumbnail_mif(img, args.thumbnailFactor, img.dtype.type)
            tif.write(thumbnail, metadata={'Name': 'thumbnail'})

    print(f'OME-TIFF image saved to: {args.outputImagePath}')