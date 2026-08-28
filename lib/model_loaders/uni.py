"""MahmoodLab UNI tile encoder loader (https://huggingface.co/MahmoodLab/UNI)."""
import torch

from .common import build_normalizer, identity_postprocess


def load(model_path, destination, **kwargs):
    import timm  # local import: avoid clashing with mocov3's `vits` module

    model = timm.create_model(
        "vit_large_patch16_224", img_size=224, patch_size=16,
        init_values=1e-5, num_classes=0, dynamic_img_size=True)
    model.load_state_dict(torch.load(model_path, map_location="cpu"), strict=True)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': build_normalizer('imagenet'),
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': identity_postprocess,
    }
