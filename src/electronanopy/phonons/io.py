"""Lectura de cálculos de fonones hechos con Phonopy a partir de salidas de VASP.

Se aceptan dos disposiciones de archivos (las dos habituales de Phonopy):

1. ``phonopy_params.yaml`` (contiene estructura, desplazamientos y fuerzas).
2. ``phonopy_disp.yaml`` + ``FORCE_SETS`` (se crea con ``phonopy-init -f disp-*/vasprun.xml`` en Phonopy ≥ 4, o ``phonopy -f ...`` en 2 y 3).

Si hay un archivo ``BORN`` en la misma carpeta se usa la corrección no analítica (NAC).
"""
from __future__ import annotations
from pathlib import Path
import phonopy
from phonopy import Phonopy


class PhononReadError(ValueError):
    """No se pudo construir el cálculo de fonones con los archivos dados."""


def _find_inputs(path: Path) -> tuple[Path, Path | None]:
    """Devuelve (yaml, FORCE_SETS o None) a partir de una carpeta o de un archivo yaml."""
    if path.is_dir():
        for name in ("phonopy_params.yaml", "phonopy_params.yaml.xz", "phonopy_params.yaml.gz"):
            if (path / name).is_file():
                return path / name, None
        disp = path / "phonopy_disp.yaml"
        if disp.is_file():
            fs = path / "FORCE_SETS"
            if not fs.is_file():
                raise PhononReadError(
                    f"Se encontró phonopy_disp.yaml pero falta FORCE_SETS en {path}. "
                    "Cree FORCE_SETS con: phonopy-init -f disp-*/vasprun.xml (Phonopy >= 4) o phonopy -f disp-*/vasprun.xml (2 y 3)")
            return disp, fs
        raise PhononReadError(
            f"No hay phonopy_params.yaml ni phonopy_disp.yaml + FORCE_SETS en {path}")
    if not path.is_file():
        raise PhononReadError(f"No existe: {path}")
    fs = path.parent / "FORCE_SETS"
    return path, (fs if fs.is_file() and "disp" in path.name else None)


def load_phonopy(path: str | Path, supercell_matrix=None) -> Phonopy:
    """Carga un cálculo de fonones y construye las constantes de fuerza.

    ``path`` puede ser una carpeta o el archivo yaml. ``supercell_matrix`` solo es necesario si
    el yaml no lo incluye. Lanza ``PhononReadError`` con un mensaje accionable ante datos faltantes.
    """
    path = Path(path)
    yaml_file, force_sets = _find_inputs(path)
    born = yaml_file.parent / "BORN"
    try:
        ph = phonopy.load(
            phonopy_yaml=str(yaml_file),
            supercell_matrix=supercell_matrix,
            force_sets_filename=str(force_sets) if force_sets else None,
            born_filename=str(born) if born.is_file() else None,
            is_nac=born.is_file(),
            is_compact_fc=False,
            log_level=0,
        )
    except Exception as exc:
        raise PhononReadError(f"No se pudo leer {yaml_file.name}: {exc}") from exc
    if ph.force_constants is None:
        raise PhononReadError(
            "El cálculo no tiene fuerzas ni constantes de fuerza. Verifique que FORCE_SETS "
            "(o el yaml) contenga las fuerzas de todas las supercedas desplazadas.")
    return ph
