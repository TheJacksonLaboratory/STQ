"""MahmoodLab CONCH loader.

Requires the `conch` package (pip install git+https://github.com/mahmoodlab/CONCH.git),
which no other model in this codebase needs -- isolated here.
"""
from .common import build_normalizer, identity_postprocess


def load(model_path, destination, use_conch_normalizer=False, **kwargs):
    from conch.open_clip_custom import create_model_from_pretrained

    model, conch_normalizer = create_model_from_pretrained(
        "conch_ViT-B-16", checkpoint_path=model_path)
    model.to(destination)
    model.eval()

    transform = conch_normalizer if use_conch_normalizer else build_normalizer('imagenet')

    return {
        'model': model,
        'transform': transform,
        'forward_fn': lambda m, batch: m.encode_image(batch, proj_contrast=False, normalize=False),
        'postprocess_fn': identity_postprocess,
    }
