"""Prov-GigaPath tile-encoder loaders
(University of Washington / Microsoft Research / Providence).

prov-gigapath: registered timm architecture name is
'vit_giant_patch14_dinov2', but per the shipped config.json's model_args the
actual weights use patch_size=16 (the "patch14" label is historical/kept
from the base architecture family) and a SwiGLU MLP (mlp_ratio ~2.66667*2,
same convention as UNI2). No register tokens. global_pool="token" means
forward() already returns pooled CLS-token features -- no postprocessing
needed.

prov-gigapath-flash: NOT a built-in timm architecture. It's registered by
importing `gigapath.tile_encoder` from the `prov-gigapath` pip package
(pip install -e 'git+https://github.com/prov-gigapath/prov-gigapath'), which
calls timm's `@register_model` under the hood for the name
'gigapath_tile_enc_dinov2s'. That import is isolated inside load_flash()
below so it's never required (and never triggered) unless
prov-gigapath-flash is actually requested. Per its config.json,
global_pool="token" here too -- plain pooled output, no postprocessing.
"""
import torch.nn as nn

from .common import build_resize_crop_normalizer, identity_postprocess


def load(model_path, destination, **kwargs):
    import timm  # local import: avoid clashing with mocov3's `vits` module

    model = timm.create_model(
        'vit_giant_patch14_dinov2', pretrained=False, num_classes=0,
        img_size=224, patch_size=16, embed_dim=1536, depth=40, num_heads=24,
        init_values=1e-5, mlp_ratio=5.33334,
        mlp_layer=timm.layers.SwiGLUPacked, act_layer=nn.SiLU,
        checkpoint_path=model_path)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': build_resize_crop_normalizer('imagenet', resize=256, crop=224),
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': identity_postprocess,
    }


def load_flash(model_path, destination, **kwargs):
    import timm  # local import: avoid clashing with mocov3's `vits` module

    try:
        import gigapath.tile_encoder  # noqa: F401  registers 'gigapath_tile_enc_dinov2s'
    except ImportError as exc:
        raise ImportError(
            "prov-gigapath-flash requires the `gigapath` package: "
            "pip install -e 'git+https://github.com/prov-gigapath/prov-gigapath'"
        ) from exc

    model = timm.create_model(
        'gigapath_tile_enc_dinov2s', pretrained=False, num_classes=0,
        checkpoint_path=model_path)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': build_resize_crop_normalizer('imagenet', resize=256, crop=224),
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': identity_postprocess,
    }
