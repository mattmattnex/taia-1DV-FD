from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Sequence, Union

import numpy as np
import torch

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

Source = Union[str, Callable[[Any], torch.Tensor]]


# =========================
# --- Plot declarations ---
# =========================
@dataclass
class Var:
    source: Source
    label: str | None = None
    color: str | None = None

    def __post_init__(self):
        if self.label is None:
            self.label = self.source if isinstance(self.source, str) \
                else getattr(self.source, "__name__", "var")

    def get(self, model) -> torch.Tensor:
        if callable(self.source):
            return self.source(model)
        obj = model.state
        for attr in self.source.split("."):
            obj = getattr(obj, attr)
        return obj


VarLike = Union[Var, str, tuple, Callable]


@dataclass
class Panel:
    vars: Sequence[VarLike]
    title: str | None = None
    xlabel: str = ""
    xlim: tuple[float, float] | None = None   # None -> follows the data
    xscale: str = "linear"
    std: bool = True                          # draw the ensemble spread envelope
    legend: bool = True

    def __post_init__(self):
        if isinstance(self.vars, (str, Var)) or callable(self.vars):
            self.vars = [self.vars]
        self.vars = [_as_var(v) for v in self.vars]
        if self.title is None:
            self.title = ", ".join(v.label for v in self.vars)


def _as_var(v: VarLike) -> Var:
    if isinstance(v, Var):
        return v
    if isinstance(v, tuple):
        return Var(*v)
    return Var(v)


def _as_panel(p) -> Panel:
    if isinstance(p, Panel):
        return p
    return Panel(list(p) if isinstance(p, list) else [p])


# =================
# --- Animator  ---
# =================
class Animator:
    def __init__(self, model, panels: Sequence[Panel | VarLike], every: int = 1,
                 figsize: tuple[float, float] | None = None):
        grid = model.grid
        self.z_c = grid.z_c[0, :, 0].cpu().numpy()
        self.z_w = grid.z_w[0, :, 0].cpu().numpy()
        self.every = max(1, int(every))
        self.panels = [_as_panel(p) for p in panels]

        n = len(self.panels)
        plt.ion()
        self.fig, axes = plt.subplots(1, n, figsize=figsize or (2.8 * n + 1, 5),
                                      sharey=True, squeeze=False)
        self.axes = axes[0]

        # One list of (var, line, poly) per panel
        self.artists: list[list[tuple[Var, Any, Any]]] = []
        for ax, panel in zip(self.axes, self.panels):
            ax.set_title(panel.title)
            ax.set_xlabel(panel.xlabel)
            ax.grid(True)
            if panel.xlim is not None:
                ax.set_xlim(*panel.xlim)
            if panel.xscale != "linear":
                ax.set_xscale(panel.xscale)
            arts = []
            for v in panel.vars:
                line, poly = self._gen_line_poly(ax, v, panel.std)
                arts.append((v, line, poly))
            self.artists.append(arts)

            if panel.legend:
                ax.legend(loc="lower right", fontsize="small")

        self.axes[0].set_ylabel("Depth (m)")
        z_all = np.concatenate([self.z_c, self.z_w])
        self.axes[0].set_ylim(z_all.min(), z_all.max())
        self.fig.suptitle("t = 0 d")          # reserve room before laying out
        self.fig.tight_layout()
        plt.show(block=False)

    # ------------------------------------------------------------------
    def update(self, model, kt: int | None = None) -> None:
        """Redraw. If `kt` is given, only every `self.every` steps."""
        if kt is not None and kt % self.every:
            return

        t = float(model.state.time)
        self.fig.suptitle(f"t = {t / 86400:.2f} d  ({t:.0f} s)")

        for ax, panel, arts in zip(self.axes, self.panels, self.artists):
            lo, hi = np.inf, -np.inf
            for var, line, poly in arts:
                m, s = self._get_stats_np(var.get(model), var.label)
                z = self._z_for(m.size, var.label)
                line.set_data(m, z)
                if poly is not None:
                    poly_x = np.concatenate([m - s, (m + s)[::-1]])
                    poly_y = np.concatenate([z, z[::-1]])
                    poly.set_verts([np.column_stack([poly_x, poly_y])])
                    lo, hi = min(lo, np.nanmin(m - s)), max(hi, np.nanmax(m + s))
                else:
                    lo, hi = min(lo, np.nanmin(m)), max(hi, np.nanmax(m))

            if panel.xlim is None and np.isfinite(lo) and np.isfinite(hi):
                if panel.xscale == "log":
                    lo = max(lo, 1e-30)
                    ax.set_xlim(lo * 0.5, hi * 2.0)
                else:
                    pad = 0.05 * (hi - lo) if hi > lo else 0.05 * max(abs(hi), 1e-12)
                    ax.set_xlim(lo - pad, hi + pad)

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    update_plot = update  # backward compatibility

    def close(self, block: bool = True) -> None:
        """Call after the loop: keeps the window open (block=True) or closes it."""
        plt.ioff()
        if block:
            plt.show()
        else:
            plt.close(self.fig)

    # ==========================
    # --- Internal functions ---
    # ==========================
    def _gen_line_poly(self, ax, var: Var, with_std: bool):
        nan = np.full_like(self.z_c, np.nan)
        (line,) = ax.plot(nan, self.z_c, label=var.label, color=var.color)
        poly = None
        if with_std:
            poly = ax.fill_betweenx(self.z_c, nan, nan, alpha=0.3, color=line.get_color())
        return line, poly

    def _z_for(self, nz: int, label: str) -> np.ndarray:
        if nz == self.z_c.size:
            return self.z_c
        if nz == self.z_w.size:
            return self.z_w
        raise ValueError(f"'{label}' has {nz} vertical points; expected "
                         f"{self.z_c.size} (z_c) or {self.z_w.size} (z_w).")

    @staticmethod
    def _get_stats_np(var: torch.Tensor, label: str = ""):
        var = var.detach()
        if var.is_complex():
            raise TypeError(f"'{label}' is complex: plot '{label}.real' / '{label}.imag' "
                            f"or pass a callable, e.g. lambda m: m.state.{label}.abs()")
        var = var.reshape(var.shape[0], -1).float()          # (ne, nz)
        corr = 1 if var.shape[0] > 1 else 0                   # avoid NaN std when ne == 1
        std, mean = torch.std_mean(var, dim=0, correction=corr)
        return mean.cpu().numpy(), std.cpu().numpy()