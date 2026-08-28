"""XiyueWang CTransPath loader (Swin-Tiny with a custom convolutional stem)."""
import torch
import torch.nn as nn

from .common import build_normalizer, identity_postprocess


def load(model_path, destination, **kwargs):
    import timm  # local import: avoid clashing with mocov3's `vits` module
    from timm.models.layers.helpers import to_2tuple

    class ConvStem(nn.Module):
        def __init__(self, img_size=224, patch_size=4, in_chans=3, embed_dim=768,
                     norm_layer=None, flatten=True):
            super().__init__()

            assert patch_size == 4
            assert embed_dim % 8 == 0

            img_size = to_2tuple(img_size)
            patch_size = to_2tuple(patch_size)
            self.img_size = img_size
            self.patch_size = patch_size
            self.grid_size = (img_size[0] // patch_size[0], img_size[1] // patch_size[1])
            self.num_patches = self.grid_size[0] * self.grid_size[1]
            self.flatten = flatten

            stem = []
            input_dim, output_dim = 3, embed_dim // 8
            for _ in range(2):
                stem.append(nn.Conv2d(input_dim, output_dim, kernel_size=3, stride=2, padding=1, bias=False))
                stem.append(nn.BatchNorm2d(output_dim))
                stem.append(nn.ReLU(inplace=True))
                input_dim = output_dim
                output_dim *= 2
            stem.append(nn.Conv2d(input_dim, embed_dim, kernel_size=1))
            self.proj = nn.Sequential(*stem)
            self.norm = norm_layer(embed_dim) if norm_layer else nn.Identity()

        def forward(self, x):
            B, C, H, W = x.shape
            assert H == self.img_size[0] and W == self.img_size[1], \
                f"Input image size ({H}*{W}) doesn't match model ({self.img_size[0]}*{self.img_size[1]})."
            x = self.proj(x)
            if self.flatten:
                x = x.flatten(2).transpose(1, 2)  # BCHW -> BNC
            x = self.norm(x)
            return x

    model = timm.create_model('swin_tiny_patch4_window7_224', embed_layer=ConvStem, pretrained=False)
    model.head = nn.Identity()
    model.load_state_dict(torch.load(model_path)['model'], strict=True)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': build_normalizer('imagenet'),
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': identity_postprocess,
    }
