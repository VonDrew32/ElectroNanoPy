"""Visualización 3D de la celda y de superceldas con matplotlib."""
from __future__ import annotations
import itertools
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from .descriptors import bond_table, is_two_dimensional

# Colores por elemento (subconjunto; el resto usa un gris).
COLORS = {"H": "#FFFFFF", "C": "#555555", "N": "#3050F8", "O": "#FF0D0D", "S": "#FFFF30",
          "Mn": "#9C7AC7", "Sr": "#00FF00", "Sb": "#9E63B5", "Mo": "#54B5B5", "W": "#2194D6",
          "Se": "#FFA100", "Ga": "#C28F8F", "B": "#FFB5B5", "Fe": "#E06633"}
RADII = {"default": 60}  # tamaño base de los marcadores


def _cell_edges(lattice):
    m = lattice.matrix
    corners = {c: sum(ci * v for ci, v in zip(c, m)) for c in itertools.product((0, 1), repeat=3)}
    return [(corners[a], corners[b]) for a, b in itertools.combinations(corners, 2)
            if sum(x != y for x, y in zip(a, b)) == 1]


def plot_structure(structure: Structure, supercell: tuple[int, int, int] = (1, 1, 1),
                   show_bonds: bool = True, tolerance: float = 1.15, ax=None, view: str = "auto"):
    """Dibuja la estructura. Devuelve la figura (para incrustarla en la interfaz).

    ``view``: ``"auto"`` (superior si el material es 2D, isométrica si no), ``"top"`` (mirando a lo largo de c),
    ``"side"`` (mirando a lo largo de b) o ``"iso"``.
    """
    sc = structure.copy()
    sc.make_supercell(supercell)
    if ax is None:
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111, projection="3d")
    else:
        fig = ax.figure
    for p, q in _cell_edges(sc.lattice):
        ax.plot(*zip(p, q), color="0.5", lw=0.8)
    for el in sc.composition.elements:
        pts = np.array([s.coords for s in sc if s.specie == el])
        ax.scatter(*pts.T, s=RADII["default"], color=COLORS.get(el.symbol, "#BBBBBB"),
                   edgecolors="k", linewidths=0.4, label=el.symbol, depthshade=False)
    ghosts: dict[str, list] = {}                                    # vecinos fuera de la celda (imágenes periódicas)
    if show_bonds:
        for b in bond_table(sc, tolerance):
            p = sc[b["site"]].coords
            frac_q = sc.frac_coords[b["neighbor_index"]] + np.array(b["image"])
            q = sc.lattice.get_cartesian_coords(frac_q)
            ax.plot(*zip(p, q), color="0.3", lw=1.0)
            if ((frac_q < -1e-3) | (frac_q > 1 + 1e-3)).any():
                ghosts.setdefault(b["neighbor"], []).append(q)
    for sym, pts in ghosts.items():                                 # el enlace no termina en el vacío: se marca el átomo imagen
        pts = np.unique(np.round(np.array(pts), 4), axis=0)
        ax.scatter(*pts.T, s=0.55 * RADII["default"], color=COLORS.get(sym, "#BBBBBB"), alpha=0.35,
                   edgecolors="0.4", linewidths=0.3, depthshade=False)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    ax.set_title(sc.composition.reduced_formula)
    ax.legend(loc="upper right")
    allpts = np.array([p for e in _cell_edges(sc.lattice) for p in e] + [q for v in ghosts.values() for q in v])
    ext = np.ptp(allpts, axis=0)
    lo, hi = allpts.min(axis=0), allpts.max(axis=0)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    if view == "auto":
        view = "top" if is_two_dimensional(structure) else "iso"
    if view == "top":
        ax.view_init(elev=90, azim=-90); ax.set_zticks([]); ax.set_zlabel("")
        ext = np.array([ext[0], ext[1], 0.05 * max(ext[0], ext[1])])           # la altura no resta espacio al plano
    elif view == "side":
        ax.view_init(elev=0, azim=-90); ax.set_yticks([]); ax.set_ylabel("")
    ax.set_box_aspect(ext)
    return fig


_XLABEL = {"lattice": "Parámetro de red a (Å)", "kpoints": "Malla k (N×N×1)", "encut": "ENCUT (eV)",
           "vacuum": "Vacío (Å)", "other": "Parámetro"}


def _bm_of(series, rows):
    from .series import series_summary
    return series_summary(series, rows).get("bm")


def plot_series(series, rows, ax=None, table_ax=None):
    """Energía relativa (meV) frente al parámetro de la serie; en la serie de red añade el esfuerzo σ_xx.

    La referencia es la energía mínima (serie de red) o la del último cálculo, el más convergido (las demás).
    Con ``table_ax`` dibuja además la tabla de la serie.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    fig = ax.figure
    ok = [r for r in rows if r["energy"] is not None and r["value"] is not None]
    if not ok:
        ax.text(0.5, 0.5, "La serie no tiene energías legibles (falta OSZICAR).", ha="center", va="center",
                transform=ax.transAxes); ax.set_axis_off(); return fig
    x = np.array([r["value"] for r in ok]); e = np.array([r["energy"] for r in ok])
    ref = e.min() if series.kind == "lattice" else e[-1]
    de = (e - ref) * 1000
    ax.plot(x, de, "o-", color="tab:blue")
    d = series.default_member.name
    sel = [i for i, r in enumerate(ok) if r["name"] == d]
    if sel:
        ax.plot(x[sel], de[sel], "o", ms=11, mfc="none", mec="tab:red", label="cálculo por defecto"); ax.legend(fontsize=8, loc="upper center")
    bm = series.kind == "lattice" and _bm_of(series, rows)
    if bm:
        ax.plot(bm["curve_a"], (bm["curve_e"] - ref) * 1000, "-", color="tab:orange", lw=1.4, label="Birch–Murnaghan")
        ax.axvline(bm["a0"], color="tab:orange", ls=":", lw=1)
        ax.plot([bm["a0"]], [(bm["e0"] - ref) * 1000], "*", color="tab:orange", ms=11, label=f"a₀ = {bm['a0']:.4f} Å")
        ax.legend(fontsize=8, loc="upper center")
    ax.set_xlabel(_XLABEL.get(series.kind, "Parámetro")); ax.set_ylabel("E − E_ref (meV)")
    ax.set_title(f"{series.name}: energía relativa")
    if series.kind == "lattice":
        st = [(r["value"], r["stress"][0]) for r in ok if r["stress"] is not None]
        if st:
            ax2 = ax.twinx(); ax2.plot(*zip(*st), "s--", color="tab:green", ms=4)
            ax2.axhline(0, color="tab:green", lw=0.6); ax2.set_ylabel("σ_xx (kB)", color="tab:green")
    if table_ax is not None:
        table_ax.set_axis_off()
        head = ["valor", "a (Å)", "c (Å)", "d_enlace (Å)", "vacío (Å)", "E (eV)"]
        cells = [[r["name"], *(f"{r[k]:.4f}" if r.get(k) is not None else "—" for k in ("a", "c", "bond_min", "vacuum")),
                  f"{r['energy']:.5f}" if r["energy"] is not None else "—"] for r in rows]
        t = table_ax.table(cellText=cells, colLabels=head, loc="upper center", cellLoc="center")
        t.auto_set_font_size(False); t.set_fontsize(7); t.scale(1, 1.0)
    return fig
