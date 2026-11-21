"""Densidad de carga (CHGCAR): electrones totales, promedio planar, cortes 2D e isosuperficies.

VASP guarda en el CHGCAR el producto ρ(r)·V_celda; aquí todo se devuelve ya en e/Å^3. Las mallas son
periódicas, por lo que la interpolación y las isosuperficies tratan los bordes como periódicos.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
from pymatgen.io.vasp.outputs import Chgcar


class DensityReadError(ValueError):
    """No se pudo leer el CHGCAR."""


def load_chgcar(path: str | Path) -> Chgcar:
    """Lee un CHGCAR (o CHG). Lanza ``DensityReadError`` con un mensaje claro si falla."""
    path = Path(path)
    if path.is_dir():
        path = path / "CHGCAR"
    if not path.is_file():
        raise DensityReadError(f"No existe el archivo: {path}")
    if path.stat().st_size == 0:
        raise DensityReadError(f"El archivo está vacío: {path}")
    try:
        chg = Chgcar.from_file(str(path))
    except Exception as exc:
        raise DensityReadError(f"No se pudo interpretar {path.name} como CHGCAR: {exc}") from exc
    if "total" not in chg.data:
        raise DensityReadError(f"{path.name} no contiene la densidad total.")
    return chg


def _raw(chg: Chgcar, component: str) -> np.ndarray:
    if component not in chg.data:
        raise ValueError(f"El CHGCAR no tiene el componente '{component}' (disponibles: {list(chg.data)}). "
                         "La magnetización ('diff') solo existe en cálculos polarizados en espín.")
    return np.asarray(chg.data[component])


def density(chg: Chgcar, component: str = "total") -> np.ndarray:
    """Densidad en e/Å^3 sobre la malla (nx, ny, nz). ``component``: 'total' o 'diff' (magnetización)."""
    return _raw(chg, component) / chg.structure.volume


def total_electrons(chg: Chgcar) -> float:
    """Número de electrones integrados en la celda (debe igualar los electrones de valencia del cálculo)."""
    raw = _raw(chg, "total")
    return float(raw.sum() / raw.size)


def net_magnetization(chg: Chgcar) -> float | None:
    """Magnetización neta (µB) integrada, o None si el cálculo no es polarizado en espín."""
    if "diff" not in chg.data:
        return None
    raw = np.asarray(chg.data["diff"])
    return float(raw.sum() / raw.size)


def planar_average(chg: Chgcar, axis: int = 2, component: str = "total") -> tuple[np.ndarray, np.ndarray]:
    """Promedio de la densidad en cada plano perpendicular a ``axis``: (posición en Å, e/Å^3).

    La posición es ``fracción * |a_axis|`` (supone que ese eje es perpendicular al plano, como en un
    material 2D con vacío a lo largo de c). La integral ``∑ ρ̄ dz * Área`` da el número de electrones.
    """
    rho = density(chg, component)
    other = tuple(i for i in range(3) if i != axis)
    prof = rho.mean(axis=other)
    n = rho.shape[axis]
    z = np.arange(n) / n * np.linalg.norm(chg.structure.lattice.matrix[axis])
    return z, prof


def density_at(chg: Chgcar, frac_coords, component: str = "total") -> np.ndarray:
    """Densidad (e/Å^3) en coordenadas fraccionarias, con interpolación lineal periódica."""
    from scipy.ndimage import map_coordinates
    rho = density(chg, component)
    f = np.atleast_2d(np.asarray(frac_coords, dtype=float))
    idx = (f * np.array(rho.shape)).T                      # coordenadas de índice (con borde periódico)
    return map_coordinates(rho, idx, order=1, mode="grid-wrap")


def slice_2d(chg: Chgcar, axis: int = 2, frac: float = 0.5, component: str = "total") -> np.ndarray:
    """Corte 2D de la densidad en el plano ``frac`` (fraccionario) perpendicular a ``axis``."""
    rho = density(chg, component)
    n = rho.shape[axis]
    pos = (frac % 1.0) * n
    i0 = int(np.floor(pos)) % n
    i1 = (i0 + 1) % n
    w = pos - np.floor(pos)
    return (1 - w) * np.take(rho, i0, axis=axis) + w * np.take(rho, i1, axis=axis)


def isosurface(chg: Chgcar, level: float, component: str = "total", max_faces: int = 1_000_000):
    """Isosuperficie ρ = ``level`` (e/Å^3) como malla de triángulos en coordenadas cartesianas (Å).

    Devuelve (vértices (n,3), caras (m,3)). Requiere scikit-image. La malla es periódica: se calcula
    con un borde envuelto y se descartan los vértices fuera de la celda original solo al graficar.

    Si la malla resultante supera ``max_faces`` triángulos (un nivel muy bajo en una malla densa: gráficas de varios
    GB y decenas de segundos), se recalcula tomando un punto de cada dos hasta quedar por debajo.
    """
    try:
        from skimage.measure import marching_cubes
    except ImportError as exc:
        raise ImportError("Las isosuperficies requieren scikit-image: pip install scikit-image") from exc
    rho = density(chg, component)
    if not (rho.min() < level < rho.max()):
        raise ValueError(f"level={level} fuera del rango de la densidad [{rho.min():.4g}, {rho.max():.4g}] e/Å^3")
    while True:
        padded = np.pad(rho, 1, mode="wrap")
        verts, faces, _, _ = marching_cubes(padded, level)
        if len(faces) <= max_faces or min(rho.shape) < 24:
            break
        rho = rho[::2, ::2, ::2]                                  # la malla periódica sigue siendo válida
    frac = (verts - 1.0) / np.array(rho.shape)
    return frac @ chg.structure.lattice.matrix, faces
