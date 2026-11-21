"""Cálculos de fonones SINTÉTICOS (modelos de resortes) para pruebas y demostraciones.

No representan ningún material real: sirven para verificar el flujo de lectura, el cálculo de bandas,
la DOS y el criterio de estabilidad sin necesitar salidas de VASP.
"""
from __future__ import annotations
import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms


def _build(cell, supercell, mass_symbol, pair_fc, onsite=0.0):
    ph = Phonopy(PhonopyAtoms(symbols=[mass_symbol], cell=cell, scaled_positions=[[0, 0, 0]]),
                 supercell_matrix=supercell, primitive_matrix=np.eye(3))
    sc = ph.supercell
    n = len(sc)
    frac = sc.scaled_positions
    fc = np.zeros((n, n, 3, 3))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            df = frac[j] - frac[i]
            df -= np.round(df)
            m = pair_fc(df @ sc.cell)
            if m is not None:
                fc[i, j] = m
    for i in range(n):
        fc[i, i] = -fc[i].sum(axis=0) - onsite * np.eye(3)
    ph.force_constants = fc
    ph.synthetic_fc = fc          # copia completa (n, n, 3, 3) para generar fuerzas sintéticas
    return ph


def cubic(k1: float = 1.0, k2: float = 0.5, a: float = 3.0, onsite: float = 0.0) -> Phonopy:
    """Red cúbica simple (1 átomo de C), resortes centrales a 1.º (k1) y 2.º vecinos (k2) en eV/Å^2.
    ``onsite`` > 0 ablanda todos los modos por igual (rompe la regla de suma acústica)."""
    def pair(d):
        r = np.linalg.norm(d)
        k = k1 if abs(r - a) < 1e-3 else k2 if abs(r - a * 2 ** 0.5) < 1e-3 else 0.0
        return None if k == 0.0 else -k * np.outer(d / r, d / r)
    return _build(a * np.eye(3), [3, 3, 3], "C", pair, onsite)


def triangular_2d(k: float = 1.0, kz: float = 0.3, a: float = 2.5, vacuum: float = 15.0) -> Phonopy:
    """Red triangular 2D (1 átomo de C, vacío en c): resortes centrales en el plano y rigidez fuera del plano."""
    cell = np.array([[a, 0, 0], [-a / 2, a * 3 ** 0.5 / 2, 0], [0, 0, vacuum]])
    def pair(d):
        r = np.linalg.norm(d)
        if abs(r - a) > 1e-3:
            return None
        m = -k * np.outer(d / r, d / r)
        m[2, 2] = -kz
        return m
    return _build(cell, [4, 4, 1], "C", pair)


def save_as_files(ph: Phonopy, folder, layout: str = "params"):
    """Escribe el cálculo en disco como lo haría Phonopy con VASP ('params' o 'disp+forcesets')."""
    from pathlib import Path
    from phonopy.file_IO import write_FORCE_SETS
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    ph.generate_displacements(distance=0.01)
    fc = ph.synthetic_fc
    forces = []
    for d in ph.dataset["first_atoms"]:
        u = np.zeros((len(ph.supercell), 3)); u[d["number"]] = d["displacement"]
        forces.append(-np.einsum("ijab,jb->ia", fc, u))
    ph.forces = np.array(forces)
    if layout == "params":
        ph.save(str(folder / "phonopy_params.yaml"))
    else:
        ph.save(str(folder / "phonopy_disp.yaml"), settings={"force_sets": False})
        write_FORCE_SETS(ph.dataset, str(folder / "FORCE_SETS"))
    return folder
