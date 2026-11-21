"""Descriptores estructurales: red, simetría, enlaces, coordinación y dimensionalidad."""
from __future__ import annotations
import itertools
import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def lattice_parameters(structure: Structure) -> dict:
    """Parámetros de red (Å, grados) y volumen (Å^3)."""
    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles
    return {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma,
            "volume": structure.volume}


def space_group(structure: Structure, symprec: float = 1e-3) -> tuple[str, int]:
    """(símbolo, número) del grupo espacial. ``symprec`` (Å) es la tolerancia de simetría."""
    sga = SpacegroupAnalyzer(structure, symprec=symprec)
    return sga.get_space_group_symbol(), sga.get_space_group_number()


def vacuum_thickness(structure: Structure, axis: int = 2) -> float:
    """Mayor región vacía (Å) a lo largo de ``axis`` medida perpendicular a los otros dos vectores.

    Para un cristal 3D da el mayor hueco entre capas de átomos; para un material
    2D con región de vacío, da el espesor del vacío.
    """
    lat = structure.lattice.matrix
    others = [i for i in range(3) if i != axis]
    area = np.linalg.norm(np.cross(lat[others[0]], lat[others[1]]))
    height = structure.volume / area                      # altura perpendicular de la celda
    z = np.sort(structure.frac_coords[:, axis] % 1.0)
    gaps = np.diff(np.concatenate([z, [z[0] + 1.0]]))     # incluye el hueco periódico
    return float(gaps.max() * height)


def is_two_dimensional(structure: Structure, min_vacuum: float = 8.0, axis: int = 2) -> bool:
    """Heurística: hay vacío mayor que ``min_vacuum`` Å a lo largo de ``axis``."""
    return vacuum_thickness(structure, axis) >= min_vacuum


def _radius(el) -> float:
    return float(el.atomic_radius or 1.5)


def _cutoff(ei, ej, tolerance: float, mixed: bool = False) -> float:
    """Distancia máxima de enlace entre dos elementos: tolerance·(r_i + r_j); 0 si no hay enlace.

    Entre dos metales no se cuentan enlaces cuando el compuesto contiene algún no metal (LiF, SrMnSb2: la distancia
    catión–catión no es un enlace y inflaría la coordinación). En una estructura solo de metales la tolerancia entre
    metales se limita a 1.0.
    """
    if ei.is_metal and ej.is_metal:
        if mixed:
            return 0.0
        tolerance = min(tolerance, 1.0)
    return tolerance * (_radius(ei) + _radius(ej))


def _has_nonmetal(structure: Structure) -> bool:
    return any(not e.is_metal for e in structure.composition.elements)


def bond_table(structure: Structure, tolerance: float = 1.15) -> list[dict]:
    """Enlaces únicos por sitio: dos átomos están enlazados si d <= tolerance*(r_i + r_j).

    Los radios son los atómicos de pymatgen; ``tolerance`` es ajustable desde la interfaz (entre dos metales se
    usa como máximo 1.0, ver ``_cutoff``).
    Devuelve una lista de dicts con índices del sitio, elemento vecino y distancia (Å).
    """
    max_r = tolerance * 2 * max(_radius(sp) for sp in structure.composition.elements)
    out = []
    mixed = _has_nonmetal(structure)
    for i, nbrs in enumerate(structure.get_all_neighbors(max_r)):
        for n in nbrs:
            if n.nn_distance <= _cutoff(structure[i].specie, n.specie, tolerance, mixed):
                out.append({"site": i, "element": structure[i].specie.symbol,
                            "neighbor": n.specie.symbol, "distance": float(n.nn_distance),
                            "neighbor_index": int(n.index), "image": tuple(int(x) for x in n.image)})
    return out


def coordination_numbers(structure: Structure, tolerance: float = 1.15) -> list[int]:
    """Número de coordinación de cada sitio."""
    counts = [0] * len(structure)
    for b in bond_table(structure, tolerance):
        counts[b["site"]] += 1
    return counts


def bond_angles(structure: Structure, tolerance: float = 1.15) -> list[dict]:
    """Ángulos de enlace (grados) para cada par de vecinos de cada sitio."""
    max_r = tolerance * 2 * max(_radius(sp) for sp in structure.composition.elements)
    angles = []
    mixed = _has_nonmetal(structure)
    for i, nbrs in enumerate(structure.get_all_neighbors(max_r)):
        bonded = [n for n in nbrs if n.nn_distance <= _cutoff(structure[i].specie, n.specie, tolerance, mixed)]
        center = structure[i].coords
        for n1, n2 in itertools.combinations(bonded, 2):
            v1, v2 = n1.coords - center, n2.coords - center
            cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            angles.append({"center": i, "element": structure[i].specie.symbol,
                           "pair": (n1.specie.symbol, n2.specie.symbol),
                           "angle": float(np.degrees(np.arccos(np.clip(cosang, -1, 1))))})
    return angles


def summary(structure: Structure, tolerance: float = 1.15) -> dict:
    """Resumen estructural listo para mostrar en la interfaz."""
    symbol, number = space_group(structure)
    bonds = bond_table(structure, tolerance)
    d = [b["distance"] for b in bonds]
    return {"formula": structure.composition.reduced_formula, "n_atoms": len(structure),
            "lattice": lattice_parameters(structure), "space_group": symbol,
            "space_group_number": number, "vacuum": vacuum_thickness(structure),
            "two_dimensional": is_two_dimensional(structure),
            "coordination": coordination_numbers(structure, tolerance),
            "bond_min": min(d) if d else None, "bond_max": max(d) if d else None}
