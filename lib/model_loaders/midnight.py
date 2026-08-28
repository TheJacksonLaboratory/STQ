"""Kaiko-ai Midnight-12k loader (transformers Dinov2Model, MIT license).

Per config_kaikoMidnight12k.json: Dinov2Model, hidden_size=1536, patch14.
Midnight is trained with DINOv2-style register tokens, so its
Dinov2Model/Dinov2WithRegistersModel config carries a `num_register_tokens`
attribute; last_hidden_state is laid out as [CLS, reg_0..reg_{n-1}, patch...].
We read that count off the loaded model's config (defaulting to 0 for
configs without it) so the CLS/patch split below correctly skips the
register tokens instead of averaging one of them into the patch mean.

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

    num_register_tokens = getattr(model.config, 'num_register_tokens', 0)

    return {
        'model': model,
        'transform': build_normalizer('midnight'),
        'forward_fn': lambda m, batch: m(pixel_values=batch).last_hidden_state,
        'postprocess_fn': concat_cls_mean_patch_postprocess(
            num_prefix_tokens=1 + num_register_tokens),
    }
