"""Shared utilities for per-model feature-extractor loaders.

Each model-specific loader module (uni.py, hoptimus.py, phikon.py, ...)
exposes one or more `load_xxx(model_path, destination, **kwargs)` functions
that return a dict with the keys expected by `load_model()` in
compute_features.py:

    {
        'model':          nn.Module, already .to(destination) and .eval()
        'transform':      callable(PIL.Image) -> torch.Tensor
        'forward_fn':     callable(model, batch: torch.Tensor) -> torch.Tensor
        'postprocess_fn': callable(np.ndarray) -> np.ndarray   (optional)
    }

`postprocess_fn` runs *after* forward_fn's output has been moved to CPU and
converted to a numpy array. It's how each model turns whatever its own
forward pass returns (a pooled vector, a full token sequence, CLS + patch
tokens, ...) into a flat (batch, num_features) array. Models that already
return a pooled vector can omit it (defaults to identity).

This module holds only the pieces that are genuinely shared across many
model families. Anything specific to one model/family (custom nn.Module
subclasses, one-off package imports, architecture kwargs) lives in that
model's own file so dependencies stay isolated.
"""
import os

import numpy as np
from torchvision import transforms

NORMALIZATION_STATS = {
    'imagenet': {'mean': (0.485, 0.456, 0.406), 'std': (0.229, 0.224, 0.225)},
    'hoptimus': {'mean': (0.707223, 0.578729, 0.703617),
                 'std': (0.211883, 0.230117, 0.177517)},
    'midnight': {'mean': (0.5, 0.5, 0.5), 'std': (0.5, 0.5, 0.5)},
    'inception': {'mean': (0.5, 0.5, 0.5), 'std': (0.5, 0.5, 0.5)},
}


def build_normalizer(normalization, size=224):
    """Resize -> ToTensor -> Normalize, using one of NORMALIZATION_STATS."""
    stats = NORMALIZATION_STATS[normalization]
    return transforms.Compose([
        transforms.Resize(size),
        transforms.ToTensor(),
        transforms.Normalize(mean=stats['mean'], std=stats['std']),
    ])


def build_resize_crop_normalizer(normalization, resize=256, crop=224):
    """Resize -> CenterCrop -> ToTensor -> Normalize.

    prov-gigapath and prov-gigapath-flash are documented with an explicit
    Resize(256)/CenterCrop(224) pair rather than a plain Resize(224). If your
    tiling pipeline already emits exactly `crop` x `crop` tiles, the
    CenterCrop below is a no-op, so this is safe to use unconditionally.
    """
    stats = NORMALIZATION_STATS[normalization]
    return transforms.Compose([
        transforms.Resize(resize, interpolation=transforms.InterpolationMode.BICUBIC),
        transforms.CenterCrop(crop),
        transforms.ToTensor(),
        transforms.Normalize(mean=stats['mean'], std=stats['std']),
    ])


def resolve_hf_checkpoint_dir(path):
    """Resolve a user-supplied --model-checkpoint-path into an HF-style
    model directory (config.json + weights + preprocessor_config.json).

    Transformers-family models (Phikon, Phikon-v2, Kaiko-Midnight) are
    shipped as a directory, not a single file. The rest of this codebase's
    CLI takes one --model-checkpoint-path, so we accept either:
      - a directory containing those files directly, or
      - a path to one specific weight file inside that directory
        (model.safetensors, pytorch_model.bin, *.bin, *.pth) -- in which
        case we use its parent directory.
    """
    if os.path.isdir(path):
        candidate = path
    else:
        candidate = os.path.dirname(path) or '.'

    if os.path.isfile(os.path.join(candidate, 'config.json')):
        return candidate

    raise FileNotFoundError(
        f"--model-checkpoint-path={path!r} does not point at a HuggingFace-"
        f"style model directory (no config.json found in {candidate!r}). "
        f"For transformers-based models, pass either the directory "
        f"containing config.json / model.safetensors (or pytorch_model.bin) "
        f"/ preprocessor_config.json, or a path to one of those weight "
        f"files inside it.")


def identity_postprocess(features):
    """No-op: the model's forward_fn already returns (batch, num_features)."""
    return features


def cls_token_postprocess(hidden_states):
    """hidden_states: (batch, seq_len, hidden) -> CLS token only.

    Used for models whose forward_fn returns the full token sequence (CLS,
    optionally followed by register tokens, then patch tokens) but whose own
    convention is "the CLS token is the embedding" -- e.g. Phikon,
    Phikon-v2, H0-mini.
    """
    return hidden_states[:, 0, :]


def concat_cls_mean_patch_postprocess(num_prefix_tokens=1):
    """Build a postprocess fn: concat(CLS token, mean of patch tokens).

    `num_prefix_tokens` is the number of non-patch tokens at the start of
    the sequence (1 for CLS-only; 1 + num_register_tokens for models that
    also use DINOv2-style register tokens, which should then be excluded
    from both the CLS slice and the patch-token mean).
    """
    def _postprocess(hidden_states):
        cls = hidden_states[:, 0, :]
        patches = hidden_states[:, num_prefix_tokens:, :]
        return np.concatenate([cls, patches.mean(axis=1)], axis=-1)
    return _postprocess
