"""Lectura de estructuras desde archivos de VASP (POSCAR/CONTCAR)."""
from __future__ import annotations
from pathlib import Path
from pymatgen.core import Structure


class StructureReadError(ValueError):
    """El archivo no existe, está vacío o no tiene formato POSCAR/CONTCAR válido."""


def load_structure(path: str | Path) -> Structure:
    """Lee un POSCAR/CONTCAR y devuelve un ``pymatgen.core.Structure``.

    Lanza ``StructureReadError`` con un mensaje claro si el archivo es inválido,
    para que la interfaz gráfica pueda mostrarlo al usuario.
    """
    path = Path(path)
    if not path.is_file():
        raise StructureReadError(f"No existe el archivo: {path}")
    if path.stat().st_size == 0:
        raise StructureReadError(f"El archivo está vacío: {path}")
    try:
        structure = Structure.from_file(path)
    except Exception as exc:  # pymatgen lanza distintos tipos según el fallo
        raise StructureReadError(
            f"No se pudo interpretar {path.name} como POSCAR/CONTCAR: {exc}") from exc
    if len(structure) == 0:
        raise StructureReadError(f"La estructura de {path.name} no contiene átomos.")
    return structure
