"""Google Inception V3 loader, timm/torch (`timm/inception_v3.tv_in1k`).

Per the shipped config.json (torchvision `tv_in1k` weights):
  - input_size=(3, 299, 299), mean=std=0.5 (not imagenet stats)
  - num_classes=0 -> timm's default average-pool head returns a plain
    (batch, 2048) pooled vector, so no postprocessing is needed.

This replaces the earlier TensorFlow/Keras InceptionV3 implementation
(imagenet weights, GlobalAveragePooling2D) that used to live here and ran
via a dedicated TF container/CLI (featuresinception.py,
bin/run-extract-inception.py). Loading real tv_in1k weights through timm
lets Inception go through the same torch pipeline (run-extract.py /
compute_features.py) as every other model in this package.
"""
from .common import build_normalizer, identity_postprocess


def load(model_path, destination, **kwargs):
    import timm  # local import: avoid clashing with mocov3's `vits` module

    # The checkpoint is the full tv_in1k classifier (1000-way fc head), so it
    # must be loaded with num_classes=1000 to match state_dict keys. Drop the
    # head afterward via reset_classifier(0) to get plain pooled features.
    model = timm.create_model(
        'inception_v3', pretrained=False, num_classes=1000,
        checkpoint_path=model_path)
    model.reset_classifier(0)
    model.to(destination)
    model.eval()

    return {
        'model': model,
        'transform': build_normalizer('inception', size=299),
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': identity_postprocess,
    }
