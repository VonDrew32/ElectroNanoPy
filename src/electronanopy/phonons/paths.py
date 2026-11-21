"""Caminos de alta simetría para bandas de fonones (2D con camino propio; 3D con seekpath)."""
from __future__ import annotations
import warnings
import numpy as np

# Un camino es una lista de segmentos: (etiqueta_inicio, etiqueta_fin, q_inicio, q_fin) en coordenadas
# fraccionarias de la celda primitiva. ``connected[i]`` indica si el segmento i continúa en el i+1.
Segment = tuple[str, str, np.ndarray, np.ndarray]


def reciprocal_lattice(cell: np.ndarray) -> np.ndarray:
    """Filas = vectores recíprocos (1/Å, con el factor 2π)."""
    return 2 * np.pi * np.linalg.inv(np.asarray(cell, dtype=float)).T


def in_plane_class(cell: np.ndarray, tol: float = 0.02) -> str:
    """Clasifica la red 2D (vectores a, b de la celda): 'hex', 'square', 'rect' u 'oblique'."""
    a, b = np.asarray(cell, dtype=float)[:2]
    la, lb = np.linalg.norm(a), np.linalg.norm(b)
    gamma = np.degrees(np.arccos(np.dot(a, b) / (la * lb)))
    equal = abs(la - lb) / max(la, lb) < tol
    if equal and (abs(gamma - 120) < 1 or abs(gamma - 60) < 1):
        return "hex"
    if abs(gamma - 90) < 1:
        return "square" if equal else "rect"
    return "oblique"


def path_2d(cell: np.ndarray) -> tuple[list[Segment], list[bool]]:
    """Camino estándar en el plano (kz = 0) según el tipo de red."""
    cell = np.asarray(cell, dtype=float)
    kind = in_plane_class(cell)
    G = np.zeros(3)
    if kind == "hex":
        rec = reciprocal_lattice(cell)
        M = np.array([0.5, 0, 0])
        # K es el vértice de la zona: |K| = |b1| / sqrt(3); se prueba (1/3,1/3) y (1/3,-1/3)
        target = np.linalg.norm(rec[0]) / np.sqrt(3)
        K = min((np.array([1/3, 1/3, 0]), np.array([1/3, -1/3, 0])),
                key=lambda q: abs(np.linalg.norm(q @ rec) - target))
        return [("Γ", "M", G, M), ("M", "K", M, K), ("K", "Γ", K, G)], [True, True, False]
    X, Y = np.array([0.5, 0, 0]), np.array([0, 0.5, 0])
    S = np.array([0.5, 0.5, 0])
    if kind == "square":
        return [("Γ", "X", G, X), ("X", "M", X, S), ("M", "Γ", S, G)], [True, True, False]
    return [("Γ", "X", G, X), ("X", "S", X, S), ("S", "Y", S, Y), ("Y", "Γ", Y, G)], [True, True, True, False]


def _to_user_basis(L_std: np.ndarray, rot: np.ndarray, L_user: np.ndarray, tol: float = 5e-3):
    """Matriz que pasa de coordenadas fraccionarias recíprocas de la celda estándar de seekpath a las de la celda
    del cálculo, o None si no describen la misma red.

    ``rot`` es la rotación cartesiana (columnas) de la celda original a la estándar (convención de spglib), así que
    un vector estándar ``v`` (fila) vuelve al marco del usuario como ``v @ rot``. Las celdas describen la misma red si,
    tras esa rotación, los vectores del usuario son combinaciones enteras de los estándar (matriz unimodular).
    También se admite la identidad por si la celda ya estaba en el marco estándar.
    """
    for R in (rot, np.eye(3)):
        Ls = L_std @ R                                       # vectores estándar (filas) en el marco del usuario
        M = L_user @ np.linalg.inv(Ls)                        # L_user = M · Ls
        if np.abs(M - np.round(M)).max() < tol and abs(abs(np.linalg.det(np.round(M))) - 1) < tol:
            rec_std, rec_user = reciprocal_lattice(L_std), reciprocal_lattice(L_user)
            # fracción estándar -> cartesiano (marco estándar) -> marco del usuario -> fracción del usuario
            return rec_std @ R @ np.linalg.inv(rec_user)
    return None


def path_3d(primitive_cell: np.ndarray, numbers: np.ndarray, scaled_positions: np.ndarray
            ) -> tuple[list[Segment], list[bool]]:
    """Camino 3D con seekpath, expresado en la base de la celda primitiva del cálculo.

    seekpath estandariza la celda (otra base de la misma red, a veces rotada); los puntos q se transforman a la
    base del usuario. Si las redes no coinciden se lanza ``ValueError`` para que se indique el camino a mano.
    """
    import seekpath
    L_p = np.asarray(primitive_cell, dtype=float)
    sp = seekpath.get_path((L_p, np.asarray(scaled_positions), np.asarray(numbers)))
    T = _to_user_basis(np.asarray(sp["primitive_lattice"], dtype=float),
                       np.asarray(sp.get("rotation_matrix", np.eye(3)), dtype=float), L_p)
    if T is None:
        raise ValueError("La celda primitiva de seekpath no corresponde a la del cálculo de fonones; "
                         "indique el camino de alta simetría manualmente (parámetro path=).")
    pts = {k: np.asarray(v, dtype=float) @ T for k, v in sp["point_coords"].items()}
    segs, conn = [], []
    prev_end = None
    for a, b in sp["path"]:
        segs.append((a.replace("GAMMA", "Γ"), b.replace("GAMMA", "Γ"), pts[a], pts[b]))
        if prev_end is not None:
            conn.append(prev_end == a)
        prev_end = b
    conn.append(False)
    return segs, conn


def interpolate(segments: list[Segment], npoints: int) -> list[np.ndarray]:
    """Puntos q (npoints por segmento, extremos incluidos)."""
    return [np.linspace(qs, qe, npoints) for _, _, qs, qe in segments]
