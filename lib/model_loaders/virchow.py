"""Paige Virchow / Virchow2 loaders."""
import torch


def _build(model_path, destination, add_reg=False):
    import timm  # local import: avoid clashing with mocov3's `vits` module
    from timm.data import resolve_data_config
    from timm.data.transforms_factory import create_transform
    from timm.layers import SwiGLUPacked

    kwargs = dict(pretrained=False,
        checkpoint_path=model_path,
        mlp_layer=SwiGLUPacked,
        act_layer=torch.nn.SiLU,
        img_size=224,
        init_values=1e-5,
        num_classes=0,
        mlp_ratio=5.3375,      # non-default ratio used for the SwiGLU hidden dim
        global_pool="",
        dynamic_img_size=True)

    if add_reg:
        kwargs['reg_tokens'] = 4

    model = timm.create_model(
        "vit_huge_patch14_224",
        **kwargs)

    transform = create_transform(**resolve_data_config(model.pretrained_cfg, model=model))
    model.to(destination)
    model.eval()
    return model, transform


def load_virchow(model_path, destination, **kwargs):
    model, transform = _build(model_path, destination, add_reg=False)

    def postprocess(features):
        # Compute mean over tokens.
        return features.mean(axis=1)

    return {
        'model': model,
        'transform': transform,
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': postprocess,
    }


def load_virchow2(model_path, destination, **kwargs):
    model, transform = _build(model_path, destination, add_reg=True)

    def postprocess(features):
        # Token layout is [CLS, reg_0..reg_3, patch_0, ...]: CLS is index 0
        # and the 4 register tokens occupy indices 1-4, so patch tokens
        # start at index 5 (per Virchow2's official model-card recipe).
        # The previous features[:, 4:] slice still included the last
        # register token (index 4) in the mean -- off by one.
        return features[:, 5:].mean(axis=1)

    return {
        'model': model,
        'transform': transform,
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': postprocess,
    }
