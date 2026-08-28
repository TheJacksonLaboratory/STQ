"""Bioptimus H-optimus / H0-mini family loaders.

H-optimus-0 and H-optimus-1 share the exact same architecture
(vit_giant_patch14_reg4_dinov2, confirmed against both models' shipped
config.json) and the exact same preprocessing stats -- only the checkpoint
weights differ, so one loader parameterized by model_path serves both
--model-name=hoptimus0 and --model-name=hoptimus1.

H0-mini is architecturally different (ViT-Base/14, distilled from
H-optimus-0, SwiGLU MLP, 4 register tokens) but shares the same
normalization stats and lineage, hence living in the same file.

Per config_H-optimus-1.json / config_H-optimus-mini.json:
  - H-optimus-0 / H-optimus-1: global_pool="token" -> forward() already
    returns pooled CLS-token features directly. No postprocessing needed.
  - H0-mini: global_pool="" in model_args -> forward() returns the raw
    [CLS, reg x4, patch...] token sequence. Bioptimus report using the CLS
    token alone for all their embedders (including H-optimus-1), so we
    slice out index 0 and discard the register + patch tokens.
"""
import torch

from .common import build_normalizer, cls_token_postprocess, identity_postprocess

_NORM = 'hoptimus'


def load_hoptimus(model_path, destination, **kwargs):
    """Shared loader for H-optimus-0 and H-optimus-1."""
    import timm  # local import: avoid clashing with mocov3's `vits` module

    model = timm.create_model(
        'vit_giant_patch14_reg4_dinov2', pretrained=False, num_classes=0,
        img_size=224, checkpoint_path=model_path)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': build_normalizer(_NORM),
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': identity_postprocess,
    }


def load_h0mini(model_path, destination, **kwargs):
    import timm  # local import: avoid clashing with mocov3's `vits` module

    model = timm.create_model(
        'vit_base_patch14_reg4_dinov2', pretrained=False, num_classes=0,
        img_size=224, init_values=1e-5, reg_tokens=4, mlp_ratio=5.33334,
        global_pool='', dynamic_img_size=True,
        mlp_layer=timm.layers.SwiGLUPacked, act_layer=torch.nn.SiLU,
        checkpoint_path=model_path)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': build_normalizer(_NORM),
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': cls_token_postprocess,
    }
