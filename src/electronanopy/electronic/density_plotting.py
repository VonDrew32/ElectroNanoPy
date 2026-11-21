"""Gráficas de densidad de carga: corte 2D, promedio planar e isosuperficie."""
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from .density import slice_2d, planar_average, isosurface


def _plane_basis(lattice, axis):
    i, j = [k for k in range(3) if k != axis]
    a, b = lattice.matrix[i], lattice.matrix[j]
    e1 = a / np.linalg.norm(a)
    e2 = b - np.dot(b, e1) * e1
    return a, b, e1, e2 / np.linalg.norm(e2)


def plot_density_slice(chg, axis: int = 2, frac: float = 0.5, repeats: tuple[int, int] = (2, 2),
                       component: str = "total", ax=None, cmap: str = "magma", atoms_tol: float = 0.3):
    """Mapa de color del corte 2D en coordenadas cartesianas (celdas oblicuas incluidas).

    Se repite la celda ``repeats`` veces en el plano. Se marcan los átomos a menos de ``atoms_tol`` Å del plano.
    """
    s = chg.structure
    Z = slice_2d(chg, axis, frac, component)
    n1, n2 = Z.shape
    r1, r2 = repeats
    Zt = np.tile(Z, (r1, r2))
    f1 = np.arange(n1 * r1) / n1
    f2 = np.arange(n2 * r2) / n2
    F1, F2 = np.meshgrid(f1, f2, indexing="ij")
    a, b, e1, e2 = _plane_basis(s.lattice, axis)
    P = F1[..., None] * a + F2[..., None] * b
    X, Y = P @ e1, P @ e2
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4.5))
    m = ax.pcolormesh(X, Y, Zt, shading="gouraud", cmap=cmap)
    ax.figure.colorbar(m, ax=ax, label=r"$\rho$ (e/Å$^3$)")
    height = np.linalg.norm(s.lattice.matrix[axis])
    for site in s:
        d = (site.frac_coords[axis] - frac) % 1.0
        d = min(d, 1 - d) * height
        if d <= atoms_tol:
            i, j = [k for k in range(3) if k != axis]
            for u in range(r1):
                for v in range(r2):
                    p = (site.frac_coords[i] + u) * a + (site.frac_coords[j] + v) * b
                    ax.plot(p @ e1, p @ e2, "o", mfc="none", mec="w", ms=7)
    ax.set_aspect("equal"); ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)")
    return ax


def plot_planar_average(chg, axis: int = 2, component: str = "total", ax=None):
    z, prof = planar_average(chg, axis, component)
    if ax is None:
        _, ax = plt.subplots(figsize=(4.5, 3.5))
    ax.plot(z, prof, color="C0")
    ax.set_xlabel("posición (Å)"); ax.set_ylabel(r"$\bar\rho$ (e/Å$^3$)")
    return ax


def plot_isosurface(chg, level: float, component: str = "total", ax=None, color: str = "C1", alpha: float = 0.6):
    verts, faces = isosurface(chg, level, component)
    if ax is None:
        fig = plt.figure(figsize=(5.5, 5)); ax = fig.add_subplot(111, projection="3d")
    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], color=color, alpha=alpha, linewidth=0)
    pts = np.array([s.coords for s in chg.structure])
    ax.scatter(*pts.T, color="k", s=25)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    ax.set_box_aspect(np.ptp(verts, axis=0) + 1e-9)
    return ax
