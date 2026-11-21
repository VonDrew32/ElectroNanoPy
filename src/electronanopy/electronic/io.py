"""Lectura de resultados electrónicos de VASP (vasprun.xml, KPOINTS, OUTCAR)."""
from __future__ import annotations
from pathlib import Path
from pymatgen.io.vasp.outputs import Vasprun


class ElectronicReadError(ValueError):
    """No se pudo leer o interpretar el archivo de VASP."""


def load_vasprun(path: str | Path, parse_projected: bool = False) -> Vasprun:
    """Lee un vasprun.xml. ``parse_projected=True`` carga las proyecciones (más lento; necesita LORBIT)."""
    path = Path(path)
    if path.is_dir():
        path = path / "vasprun.xml"
    if not path.is_file():
        raise ElectronicReadError(f"No existe el archivo: {path}")
    try:
        return Vasprun(str(path), parse_projected_eigen=parse_projected, parse_potcar_file=False)
    except Exception as exc:
        raise ElectronicReadError(
            f"No se pudo leer {path.name} (¿cálculo incompleto o XML dañado?): {exc}") from exc


def _fermi_energy(efermi, default: float) -> float:
    """``efermi``: None (usar el del vasprun), un número (eV) o la ruta a otro vasprun (p. ej. el SCF)."""
    if efermi is None:
        return default
    if isinstance(efermi, (int, float)):
        return float(efermi)
    return float(load_vasprun(efermi).efermi)


def band_structure_from_vasprun(vr: Vasprun, kpoints: str | Path, efermi=None):
    """Estructura de bandas (modo línea) a partir de un ``Vasprun`` ya leído y su KPOINTS."""
    kp = Path(kpoints)
    if not kp.is_file():
        raise ElectronicReadError(
            f"Falta el KPOINTS en modo línea ({kp}). Se necesita para las etiquetas del camino k.")
    try:
        return vr.get_band_structure(kpoints_filename=str(kp), efermi=_fermi_energy(efermi, vr.efermi),
                                     line_mode=True)
    except Exception as exc:
        raise ElectronicReadError(f"No se pudo construir la estructura de bandas: {exc}") from exc


def load_band_structure(vasprun: str | Path, kpoints: str | Path | None = None, efermi=None,
                        parse_projected: bool = False):
    """Estructura de bandas de un cálculo en modo línea (camino de alta simetría).

    ``kpoints`` es el KPOINTS en modo línea (con las etiquetas); si no se da, se busca junto al vasprun.
    ``efermi`` permite fijar la energía de Fermi: un número, o la ruta al vasprun del cálculo SCF
    (recomendado en cálculos no autoconsistentes de bandas).
    """
    vr = load_vasprun(vasprun, parse_projected)
    vrpath = Path(vasprun)
    kp = Path(kpoints) if kpoints else (vrpath if vrpath.is_dir() else vrpath.parent) / "KPOINTS"
    return band_structure_from_vasprun(vr, kp, efermi)


def read_scf_info(vasprun: str | Path) -> tuple[float, float | None]:
    """(energía de Fermi, NELECT) de un vasprun.xml, sin cargar eigenvalores ni DOS (rápido)."""
    path = Path(vasprun)
    if path.is_dir():
        path = path / "vasprun.xml"
    if not path.is_file():
        raise ElectronicReadError(f"No existe el archivo: {path}")
    try:
        vr = Vasprun(str(path), parse_dos=False, parse_eigen=False, parse_potcar_file=False)
    except Exception as exc:
        raise ElectronicReadError(f"No se pudo leer {path.name}: {exc}") from exc
    nelect = vr.parameters.get("NELECT")
    efermi = vr.efermi
    if efermi is None:                    # sin <dos> pymatgen no la lee: se pide solo la DOS (sin eigenvalores)
        try:
            efermi = Vasprun(str(path), parse_dos=True, parse_eigen=False, parse_potcar_file=False).efermi
        except Exception:
            efermi = None
    if efermi is None:
        efermi = _efermi_from_outcar(path.parent / "OUTCAR")
    if efermi is None:
        raise ElectronicReadError(f"{path.name} no contiene la energía de Fermi (ni hay OUTCAR con «E-fermi»).")
    return float(efermi), (float(nelect) if nelect is not None else None)


def _efermi_from_outcar(path: Path) -> float | None:
    """Último valor «E-fermi :» del OUTCAR (eV), o None."""
    import re
    try:
        hits = re.findall(r"E-fermi\s*:\s*(-?\d+\.\d+)", path.read_text(errors="replace"))
    except OSError:
        return None
    return float(hits[-1]) if hits else None


def mesh_band_structure(vr: Vasprun, efermi=None):
    """Bandas sobre la malla k uniforme de un vasprun (SCF o DOS), sin camino de alta simetría.

    Sirve para estimar la brecha cuando no hay un cálculo de bandas en modo línea; no tiene etiquetas k.
    """
    try:
        return vr.get_band_structure(efermi=_fermi_energy(efermi, vr.efermi), line_mode=False)
    except Exception as exc:
        raise ElectronicReadError(f"No se pudieron construir las bandas de la malla k: {exc}") from exc


def load_dos(vasprun: str | Path, parse_projected: bool = True):
    """Devuelve (DOS total, DOS completa o None). La completa existe si el cálculo usó LORBIT."""
    vr = load_vasprun(vasprun, parse_projected=parse_projected)
    return vr.tdos, getattr(vr, "complete_dos", None)
