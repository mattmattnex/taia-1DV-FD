import torch
import numpy as np

def bgc_var(name, scale=1.0, label=None):
    """Animator source for an Eco3M tracer, by name.

    Resolves the index lazily so panel order does not depend on when the
    coupler was created. scale converts to readable units, e.g. 1e6 for
    mol/L -> umol/L.
    """
    def _get(model):
        j = model.config.bgc_names.index(name)
        return model.state.bgc[:, :, j] * scale
    _get.__name__ = label or name
    return _get


def eparz(bgc, floor=1e-6):
    """Animator source for the PAR profile, (ne, nz)."""
    def _get(model):
        e = np.maximum(bgc.get_eparz(), floor)      # log axis needs > 0
        return torch.from_numpy(e).to(model.state.bgc.device)
    _get.__name__ = "E_PARZ"
    return _get