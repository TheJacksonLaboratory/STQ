"""Mahmood/TransPath MoCoV3 (ViT-Small) loader. GPU-only.

Needs `/TransPath/vits.py` on the path (shipped inside the TransPath
container image). This is exactly why the original codebase kept this
import local rather than at module top level: `vits` clashes with plain
`timm` model registration used elsewhere in the pipeline.
"""
import sys
from functools import partial

import torch
import torch.nn as nn

from .common import build_normalizer, identity_postprocess

REQUIRES_GPU = True


def load(model_path, destination, **kwargs):
    sys.path.append('/TransPath/')
    import vits  # noqa: local module shipped inside the TransPath container image

    class MoCo(nn.Module):
        """
        Build a MoCo model with a base encoder, a momentum encoder, and two MLPs
        https://arxiv.org/abs/1911.05722
        """
        def __init__(self, base_encoder, dim=256, mlp_dim=4096, T=1.0):
            super(MoCo, self).__init__()
            self.T = T

            self.base_encoder = base_encoder(num_classes=mlp_dim)
            self.momentum_encoder = base_encoder(num_classes=mlp_dim)

            self._build_projector_and_predictor_mlps(dim, mlp_dim)

            for param_b, param_m in zip(self.base_encoder.parameters(), self.momentum_encoder.parameters()):
                param_m.data.copy_(param_b.data)
                param_m.requires_grad = False

        def _build_mlp(self, num_layers, input_dim, mlp_dim, output_dim, last_bn=True):
            mlp = []
            for l in range(num_layers):
                dim1 = input_dim if l == 0 else mlp_dim
                dim2 = output_dim if l == num_layers - 1 else mlp_dim

                mlp.append(nn.Linear(dim1, dim2, bias=False))

                if l < num_layers - 1:
                    mlp.append(nn.BatchNorm1d(dim2))
                    mlp.append(nn.ReLU(inplace=True))
                elif last_bn:
                    mlp.append(nn.BatchNorm1d(dim2, affine=False))

            return nn.Sequential(*mlp)

        def _build_projector_and_predictor_mlps(self, dim, mlp_dim):
            pass

        @torch.no_grad()
        def _update_momentum_encoder(self, m):
            for param_b, param_m in zip(self.base_encoder.parameters(), self.momentum_encoder.parameters()):
                param_m.data = param_m.data * m + param_b.data * (1. - m)

        def contrastive_loss(self, q, k):
            q = nn.functional.normalize(q, dim=1)
            k = nn.functional.normalize(k, dim=1)
            k = concat_all_gather(k)
            logits = torch.einsum('nc,mc->nm', [q, k]) / self.T
            N = logits.shape[0]
            labels = (torch.arange(N, dtype=torch.long) + N * torch.distributed.get_rank()).cuda()
            return nn.CrossEntropyLoss()(logits, labels) * (2 * self.T)

        def forward(self, x1):
            return self.base_encoder(x1)

        def forward_train(self, x1, x2, m):
            q1 = self.predictor(self.base_encoder(x1))
            q2 = self.predictor(self.base_encoder(x2))

            with torch.no_grad():
                self._update_momentum_encoder(m)
                k1 = self.momentum_encoder(x1)
                k2 = self.momentum_encoder(x2)

            return self.contrastive_loss(q1, k2) + self.contrastive_loss(q2, k1)

    class MoCo_ViT(MoCo):
        def _build_projector_and_predictor_mlps(self, dim, mlp_dim):
            hidden_dim = self.base_encoder.head.weight.shape[1]
            del self.base_encoder.head, self.momentum_encoder.head

            self.base_encoder.head = self._build_mlp(3, hidden_dim, mlp_dim, dim)
            self.momentum_encoder.head = self._build_mlp(3, hidden_dim, mlp_dim, dim)
            self.predictor = self._build_mlp(2, dim, mlp_dim, dim)

    @torch.no_grad()
    def concat_all_gather(tensor):
        """
        Performs all_gather operation on the provided tensors.
        *** Warning ***: torch.distributed.all_gather has no gradient.
        """
        tensors_gather = [torch.ones_like(tensor)
                           for _ in range(torch.distributed.get_world_size())]
        torch.distributed.all_gather(tensors_gather, tensor, async_op=False)
        return torch.cat(tensors_gather, dim=0)

    model = MoCo_ViT(partial(vits.__dict__['vit_small'], stop_grad_conv1=True))
    model = nn.DataParallel(model).cuda()
    model.load_state_dict(torch.load(model_path)['state_dict'], strict=True)
    model.eval()

    return {
        'model': model,
        'transform': build_normalizer('imagenet'),
        'forward_fn': lambda m, batch: m(batch),
        'postprocess_fn': identity_postprocess,
    }
