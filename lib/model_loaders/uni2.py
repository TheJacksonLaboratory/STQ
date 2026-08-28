"""MahmoodLab UNI2-h tile encoder loader (https://huggingface.co/MahmoodLab/UNI2-h)."""
import torch

from .common import build_normalizer, identity_postprocess


def load(model_path, destination, **kwargs):
    import timm  # local import: avoid clashing with mocov3's `vits` module

    model = timm.create_model(
        pretrained=False, model_name='vit_giant_patch14_224',
        img_size=224, patch_size=14, depth=24, num_heads=24,
        init_values=1e-5, embed_dim=1536, mlp_ratio=2.66667 * 2,
        num_classes=0, no_embed_class=True,
        mlp_layer=timm.layers.SwiGLUPacked, act_layer=torch.nn.SiLU,
        reg_tokens=8, dynamic_img_size=True)
    model.load_state_dict(torch.load(model_path, map_location="cpu"), strict=True)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': build_normalizer('imagenet'),
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': identity_postprocess,
    }
