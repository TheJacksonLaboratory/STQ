"""Kaiko-ai Midnight-12k loader (transformers Dinov2Model, MIT license).

Per config_kaikoMidnight12k.json: Dinov2Model, hidden_size=1536, patch14,
no register tokens (register-token support isn't present in this config,
unlike the timm-side H-optimus/H0-mini family).

Unlike Phikon/Phikon-v2, Midnight's own model-card usage snippet does not
actually call AutoImageProcessor -- it builds a plain torchvision v2
pipeline with (0.5, 0.5, 0.5) mean/std and Resize(224)/CenterCrop(224), so
we replicate that directly via build_normalizer('midnight') rather than
relying on a possibly-mismatched preprocessor_config.json.

Feature convention (from the model card): concat(CLS token, mean of patch
tokens) -> 2 * hidden_size = 3072-dim output.
"""
from .common import resolve_hf_checkpoint_dir, build_normalizer, concat_cls_mean_patch_postprocess


def load(model_path, destination, **kwargs):
    from transformers import AutoModel

    ckpt_dir = resolve_hf_checkpoint_dir(model_path)
    model = AutoModel.from_pretrained(ckpt_dir)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': build_normalizer('midnight'),
        'forward_fn': lambda m, batch: m(pixel_values=batch).last_hidden_state,
        'postprocess_fn': concat_cls_mean_patch_postprocess(num_prefix_tokens=1),
    }
