"""Gráficas de bandas electrónicas, DOS y bandas con peso orbital (matplotlib). Energías respecto a E_F."""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.electronic_structure.core import Spin

_SPIN_COLOR = {Spin.up: "C0", Spin.down: "C3"}


def _ticks(bs):
    """Posiciones y etiquetas de los puntos de alta simetría (une 'A|B' en los saltos del camino)."""
    ticks, labels = [], []
    for br in bs.branches:
        s, e = br["start_index"], br["end_index"]
        a, b = br["name"].split("-")
        for pos, lab in ((bs.distance[s], a), (bs.distance[e], b)):
            if ticks and abs(ticks[-1] - pos) < 1e-6:
                if labels[-1] != lab:
                    labels[-1] = f"{labels[-1]}|{lab}"
            else:
                ticks.append(pos); labels.append(lab)
    fix = {"GAMMA": "Γ", "\\Gamma": "Γ", "G": "Γ"}
    return ticks, [fix.get(l, l) for l in labels]


def plot_bands(bs, ax=None, emin: float = -6.0, emax: float = 6.0, weights: np.ndarray | None = None,
               weight_label: str = "peso orbital"):
    """Bandas con E_F en 0. Con ``weights`` (nbandas, nk) se dibujan bandas con peso (color/tamaño)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4.5))
    d = np.asarray(bs.distance)
    sc = None
    for spin, e in bs.bands.items():
        e = np.asarray(e) - bs.efermi
        for ib in range(e.shape[0]):
            if weights is None or spin != Spin.up:
                ax.plot(d, e[ib], color=_SPIN_COLOR.get(spin, "C0"), lw=1.1,
                        ls="-" if spin == Spin.up else "--")
            else:
                ax.plot(d, e[ib], color="0.75", lw=0.6, zorder=1)
                sc = ax.scatter(d, e[ib], c=weights[ib], s=6 + 24 * weights[ib], cmap="viridis",
                                vmin=0, vmax=1, zorder=2)
    if sc is not None:
        ax.figure.colorbar(sc, ax=ax, label=weight_label, pad=0.02)
    ticks, labels = _ticks(bs)
    for t in ticks:
        ax.axvline(t, color="0.7", lw=0.6)
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xticks(ticks); ax.set_xticklabels(labels)
    ax.set_xlim(d[0], d[-1]); ax.set_ylim(emin, emax)
    ax.set_ylabel(r"$E - E_F$ (eV)")
    return ax


def plot_dos(dos, ax=None, emin: float = -6.0, emax: float = 6.0, horizontal: bool = False):
    """DOS total (por espín) con E_F en 0."""
    if ax is None:
        _, ax = plt.subplots(figsize=(3, 4.5))
    e = np.asarray(dos.energies) - dos.efermi
    for spin, rho in dos.densities.items():
        rho = np.asarray(rho) * (1 if spin == Spin.up else -1 if len(dos.densities) > 1 else 1)
        if horizontal:
            ax.fill_betweenx(e, 0, rho, color=_SPIN_COLOR.get(spin, "C1"), alpha=0.5)
        else:
            ax.fill_between(e, 0, rho, color=_SPIN_COLOR.get(spin, "C1"), alpha=0.5)
    if horizontal:
        ax.set_ylim(emin, emax); ax.set_xlabel("DOS (estados/eV)"); ax.axhline(0, color="k", lw=0.6, ls=":")
    else:
        ax.set_xlim(emin, emax); ax.set_xlabel(r"$E - E_F$ (eV)"); ax.set_ylabel("DOS (estados/eV)")
        ax.axvline(0, color="k", lw=0.6, ls=":")
    return ax


def plot_bands_dos(bs, dos, emin: float = -6.0, emax: float = 6.0, title: str | None = None):
    fig, (a1, a2) = plt.subplots(1, 2, sharey=True, figsize=(7, 4.5),
                                 gridspec_kw={"width_ratios": [3, 1], "wspace": 0.12})
    plot_bands(bs, a1, emin, emax); plot_dos(dos, a2, emin, emax, horizontal=True); a2.set_ylabel("")
    if title:
        fig.suptitle(title, fontsize=10)
    return fig
