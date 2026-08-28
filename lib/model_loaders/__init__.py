"""Registry of tile-feature-extractor loaders.

Each entry maps a --model-name CLI choice to (module path, function name,
gpu_only). Modules are imported lazily inside `get_loader()` so that:

  - importing `timm` for one model doesn't clash with the custom `vits`
    module mocov3 needs on sys.path,
  - a model requiring an extra package (prov-gigapath-flash needs
    `gigapath`; conch needs `conch`) only forces that import when it's
    actually requested, and
  - heavy imports (torch, transformers, timm) stay out of any module that
    doesn't need them.

To add a new model: write a `model_loaders/<name>.py` with a
`load(model_path, destination, **kwargs) -> dict` (or several such
functions, if the file covers a family of related models -- see
hoptimus.py and virchow.py), then add one line below.
"""
import importlib

_REGISTRY = {
    # model_name          (module,                     function,          gpu_only)
    'inception':           ('model_loaders.inception',   'load',           False),
    'ctranspath':          ('model_loaders.ctranspath',  'load',           False),
    'mocov3':              ('model_loaders.mocov3',      'load',           True),
    'conch':               ('model_loaders.conch',       'load',           False),
    'uni':                 ('model_loaders.uni',         'load',           False),
    'uni2':                ('model_loaders.uni2',        'load',           False),
    'hoptimus0':           ('model_loaders.hoptimus',    'load_hoptimus',  False),
    'hoptimus1':           ('model_loaders.hoptimus',    'load_hoptimus',  False),
    'h0mini':              ('model_loaders.hoptimus',    'load_h0mini',    False),
    'virchow':             ('model_loaders.virchow',     'load_virchow',   False),
    'virchow2':            ('model_loaders.virchow',     'load_virchow2',  False),
    'provgigapath':        ('model_loaders.gigapath',    'load',           False),
    'provgigapathflash':   ('model_loaders.gigapath',    'load_flash',     False),
    'phikon':              ('model_loaders.phikon',      'load_phikon',    False),
    'phikon2':             ('model_loaders.phikon',      'load_phikon_v2', False),
    'kaikomidnight12k':    ('model_loaders.midnight',    'load',           False),
}

MODEL_CHOICES = sorted(_REGISTRY.keys())
GPU_ONLY_MODELS = {name for name, (_, _, gpu_only) in _REGISTRY.items() if gpu_only}


def get_loader(model_name):
    """Return the `load(model_path, destination, **kwargs) -> dict` callable
    for the given --model-name, importing its module on first use only.
    """
    if model_name not in _REGISTRY:
        raise ValueError(f"Unknown model name: {model_name!r}. "
                          f"Choices: {MODEL_CHOICES}")
    module_path, func_name, _ = _REGISTRY[model_name]
    module = importlib.import_module(module_path)
    return getattr(module, func_name)
