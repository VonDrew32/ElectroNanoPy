"""Gráficas de bandas de fonones, DOS y combinadas (matplotlib)."""
from __future__ import annotations
import matplotlib.pyplot as plt
from .analysis import BandData, DosData, StabilityReport, THZ_TO_CM1


def _unit(unit: str) -> tuple[float, str]:
    if unit == "THz":
        return 1.0, "Frecuencia (THz)"
    if unit in ("cm-1", "cm^-1"):
        return THZ_TO_CM1, "Frecuencia (cm$^{-1}$)"
    raise ValueError("unit debe ser 'THz' o 'cm-1'")


def plot_bands(bands: BandData, ax=None, unit: str = "THz", color: str = "C0"):
    k, label = _unit(unit)
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    ax.plot(bands.distance, bands.frequencies * k, color=color, lw=1.2)
    for t in bands.ticks:
        ax.axvline(t, color="0.7", lw=0.6)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(bands.ticks); ax.set_xticklabels(bands.tick_labels)
    ax.set_xlim(bands.distance[0], bands.distance[-1])
    ax.set_ylabel(label)
    lo, _ = ax.get_ylim()
    if lo < 0:
        ax.axhspan(lo, 0, color="red", alpha=0.07)   # región de frecuencias imaginarias
    return ax


def plot_dos(dos: DosData, ax=None, unit: str = "THz", horizontal: bool = False, color: str = "C1"):
    k, label = _unit(unit)
    if ax is None:
        _, ax = plt.subplots(figsize=(3, 4))
    if horizontal:
        ax.fill_betweenx(dos.frequencies * k, 0, dos.dos / k, color=color, alpha=0.5)
        ax.set_xlabel("DOS"); ax.set_ylabel(label); ax.set_xlim(left=0)
    else:
        ax.fill_between(dos.frequencies * k, 0, dos.dos / k, color=color, alpha=0.5)
        ax.set_xlabel(label); ax.set_ylabel("DOS")
    return ax


def plot_bands_dos(bands: BandData, dos: DosData, stability: StabilityReport | None = None,
                   unit: str = "THz"):
    """Bandas a la izquierda y DOS a la derecha con el eje de frecuencias compartido."""
    fig, (a1, a2) = plt.subplots(1, 2, sharey=True, figsize=(7, 4),
                                 gridspec_kw={"width_ratios": [3, 1], "wspace": 0.12})
    plot_bands(bands, a1, unit); plot_dos(dos, a2, unit, horizontal=True)
    a2.set_ylabel("")
    if stability is not None:
        fig.suptitle(("Estable" if stability.stable else "Inestabilidad dinámica")
                     + f" (f_min = {(0.0 if abs(stability.min_frequency) < 5e-4 else stability.min_frequency):.3f} THz)", fontsize=10)
    return fig
