"""Grafeno SINTÉTICO (modelo tight-binding de primeros vecinos) para pruebas y demostraciones.

E(k) = ± t |1 + e^{-i k·a1} + e^{i k·a2}| con t = 2.7 eV. No es un cálculo DFT: sirve para verificar el
análisis (brecha, punto de Dirac, velocidad de Fermi) contra resultados analíticos conocidos:
ħ v_F = (√3/2) a t, con cono de Dirac en K y brecha nula.
"""
from __future__ import annotations
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import Dos

A, T = 2.46, 2.7        # Å, eV


def graphene_structure(vacuum: float = 15.0) -> Structure:
    return Structure(Lattice.hexagonal(A, vacuum), ["C", "C"], [[1/3, 2/3, .5], [2/3, 1/3, .5]])


def _energies(kfrac: np.ndarray) -> np.ndarray:
    # vecinos B de un átomo A: desplazamientos de celda (0,0), (-1,0), (0,1) con la base hexagonal de 120°
    f = 1 + np.exp(-2j * np.pi * kfrac[:, 0]) + np.exp(2j * np.pi * kfrac[:, 1])
    e = T * np.abs(f)
    return np.vstack([-e, e])                       # (2 bandas, nk)


def graphene_bands(npts: int = 60, gap: float = 0.0, with_projections: bool = True) -> BandStructureSymmLine:
    """Bandas Γ-M-K-Γ; ``gap`` (eV) abre una brecha simétrica (simula una subred rota, tipo h-BN)."""
    s = graphene_structure()
    G, M, K = np.zeros(3), np.array([.5, 0, 0]), np.array([1/3, 1/3, 0])
    seg = [(G, M, "\\Gamma", "M"), (M, K, "M", "K"), (K, G, "K", "\\Gamma")]
    kpts, labels = [], {}
    for a, b, la, lb in seg:
        pts = np.linspace(a, b, npts)
        labels.setdefault(la, a.copy()); labels[lb] = b.copy()
        kpts.append(pts)
    k = np.vstack(kpts)
    e = _energies(k)
    e = np.sign(e) * np.sqrt(e ** 2 + (gap / 2) ** 2)
    proj = None
    if with_projections:     # los dos estados π son orbitales pz (índice 2) mitad en cada átomo
        p = np.zeros((2, len(k), 9, 2)); p[:, :, 2, :] = 0.5
        proj = {Spin.up: p}
    return BandStructureSymmLine(k, {Spin.up: e}, s.lattice.reciprocal_lattice, 0.0, labels,
                                 structure=s, projections=proj)


def graphene_dos(n: int = 90, sigma: float = 0.08) -> Dos:
    """DOS por muestreo de la zona de Brillouin y suavizado gaussiano (2 estados por celda)."""
    u = (np.arange(n) + 0.5) / n
    k = np.array([(x, y, 0) for x in u for y in u])
    e = _energies(k).ravel()
    grid = np.linspace(-9.5, 9.5, 800)
    rho = np.exp(-((grid[:, None] - e[None, :]) / sigma) ** 2 / 2).sum(axis=1) / (sigma * np.sqrt(2 * np.pi) * len(k))
    return Dos(0.0, grid, {Spin.up: rho})


def write_outcar_magnetization(path, moments: list[float]):
    """Escribe un OUTCAR mínimo con el bloque ``magnetization (x)`` (formato de VASP)."""
    rows = "\n".join(f"{i + 1:>5d}   0.000   0.000   {m:5.3f}   {m:5.3f}" for i, m in enumerate(moments))
    tot = sum(moments)
    text = f""" vasp.6.x
 ...
 magnetization (x)

# of ion       s       p       d       tot
------------------------------------------
{rows}
------------------------------------------
tot            0.000   0.000   {tot:5.3f}   {tot:5.3f}

 other text
"""
    from pathlib import Path
    Path(path).write_text(text)
    return path


def graphene_chgcar(ngrid=(24, 24, 100), electrons: float = 8.0, sigma_xy: float = 0.6, sigma_z: float = 0.5,
                    spin_moments: tuple[float, float] | None = None):
    """CHGCAR sintético: gaussianas en cada átomo, normalizado a ``electrons`` (8 = 2 C x 4 e de valencia).

    ``spin_moments`` (µB por átomo) añade el componente 'diff' (magnetización) para pruebas de espín.
    No es una densidad DFT.
    """
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.io.vasp.outputs import Chgcar
    s = graphene_structure()
    n1, n2, n3 = ngrid
    f = np.stack(np.meshgrid(np.arange(n1) / n1, np.arange(n2) / n2, np.arange(n3) / n3, indexing="ij"), axis=-1)

    def blob(site):
        rho = np.zeros(ngrid)
        for u in (-1, 0, 1):
            for v in (-1, 0, 1):
                d = (f - (site.frac_coords + np.array([u, v, 0]))) @ s.lattice.matrix          # (n1,n2,n3,3) en Å
                r2 = (d[..., 0] ** 2 + d[..., 1] ** 2) / sigma_xy ** 2 + d[..., 2] ** 2 / sigma_z ** 2
                rho += np.exp(-r2 / 2)
        return rho / (rho.sum() * (s.volume / rho.size))          # integra 1 en la celda (e o µB por unidad)

    blobs = [blob(site) for site in s]
    data = {"total": sum(b * electrons / len(s) for b in blobs) * s.volume}      # VASP guarda rho * V_celda
    if spin_moments is not None:
        data["diff"] = sum(m * b for m, b in zip(spin_moments, blobs)) * s.volume
    return Chgcar(Poscar(s), data)
