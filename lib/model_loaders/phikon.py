"""Owkin Phikon / Phikon-v2 loaders (transformers-based, not timm).

Phikon-v2 is shipped as a full HF-style directory (config.json + weights +
preprocessor_config.json), so `model_path` is resolved to that directory
via resolve_hf_checkpoint_dir. `AutoModel.from_pretrained` /
`AutoImageProcessor.from_pretrained` handle `model.safetensors`
transparently -- no manual safetensors handling needed.

Phikon (v1) is loaded differently: rather than requiring model_path to be
an HF-style directory (config.json + pytorch_model.bin colocated, which
transformers' from_pretrained requires by filename), we build the model
directly from its known architecture (ViTModel, patch16, hidden_size=768,
12 layers/12 heads, image_size=224 -- this matches the config.json shipped
alongside the v1 checkpoint) and load the raw state dict via torch.load.
This lets --model-checkpoint-path point straight at the single weight file
Owkin ships (pytorch_model.bin), with no repacking into a directory. Its
own usage snippet loads ViTModel explicitly with add_pooling_layer=False
and reads the CLS token from last_hidden_state; preprocessing uses the
standard ViTImageProcessor defaults (resize 224, normalize mean=std=0.5)
since Owkin doesn't ship a custom preprocessor_config.json for v1.

Phikon-v2: DINOv2-family ViT-Large, hidden_size=1024, no register tokens.
Loaded via AutoModel (its config.json's model_type resolves automatically);
also CLS-token-only per its own model card.
"""
import torch
from torchvision import transforms

from .common import resolve_hf_checkpoint_dir, cls_token_postprocess

# Owkin Phikon (v1) ViTModel architecture, per its shipped config.json.
_PHIKON_V1_CONFIG = dict(
    hidden_size=768,
    num_hidden_layers=12,
    num_attention_heads=12,
    intermediate_size=3072,
    image_size=224,
    patch_size=16,
    num_channels=3,
    qkv_bias=True,
    layer_norm_eps=1e-6,
    encoder_stride=16,
)


def _hf_processor_transform(processor):
    def _transform(image):
        return processor(images=image, return_tensors='pt')['pixel_values'][0]
    return _transform


def load_phikon(model_path, destination, **kwargs):
    from transformers import ViTConfig, ViTModel

    config = ViTConfig(**_PHIKON_V1_CONFIG)
    model = ViTModel(config, add_pooling_layer=False)
    state_dict = torch.load(model_path, map_location='cpu', weights_only=True)
    # add_pooling_layer=False means the model has no `pooler.*` params, but
    # the shipped checkpoint includes the original pooler weights -- ignore
    # those extras rather than failing on them.
    model.load_state_dict(state_dict, strict=False)
    model.to(destination)
    model.eval()

    # Standard ViTImageProcessor defaults (Owkin ships no custom
    # preprocessor_config.json for v1): resize to 224, normalize mean=std=0.5.
    transform = transforms.Compose([
        transforms.Resize((224, 224)),
        transforms.ToTensor(),
        transforms.Normalize(mean=(0.5, 0.5, 0.5), std=(0.5, 0.5, 0.5)),
    ])

    return {
        'model': model,
        'transform': transform,
        'forward_fn': lambda m, batch: m(pixel_values=batch).last_hidden_state,
        'postprocess_fn': cls_token_postprocess,
    }


def load_phikon_v2(model_path, destination, **kwargs):
    from transformers import AutoImageProcessor, AutoModel

    ckpt_dir = resolve_hf_checkpoint_dir(model_path)
    processor = AutoImageProcessor.from_pretrained(ckpt_dir)
    model = AutoModel.from_pretrained(ckpt_dir)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': _hf_processor_transform(processor),
        'forward_fn': lambda m, batch: m(pixel_values=batch).last_hidden_state,
        'postprocess_fn': cls_token_postprocess,
    }
