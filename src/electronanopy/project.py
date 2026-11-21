"""Capa de integración (OE4): detección de archivos de VASP/Phonopy en una carpeta y carga de una sesión.

No depende de la interfaz gráfica: la ventana solo usa estas funciones. La carga es tolerante a fallos:
si un módulo falla (archivo dañado, faltante), los demás se cargan igual y el error queda en ``Session.errors``.
"""
from __future__ import annotations
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from typing import Callable

MODULES = ("estructura", "fonones", "electronico", "densidad")


@dataclass
class Series:
    """Serie de cálculos hermanos que varían un parámetro (malla k, ENCUT, parámetro de red, vacío…).

    Cada miembro es una carpeta de cálculo. ``default`` es el índice del miembro que se usa por defecto:
    el de menor energía en una serie de parámetro de red, y el de mayor valor (el más convergido) en las demás.
    """
    name: str                      # nombre de la carpeta contenedora, p. ej. «3-Lattice_Parameter»
    kind: str                      # lattice | kpoints | encut | vacuum | other
    folder: Path
    members: list[Path]
    values: list[float | None]
    default: int = 0

    @property
    def default_member(self) -> Path:
        return self.members[self.default]


# Palabras (tokens del nombre de la carpeta, separados por no alfanuméricos) que indican el parámetro que varía.
_KIND_TOKENS = (("lattice", ("lattice", "latt", "lat", "red", "parametro", "parameter", "a0")),
                ("kpoints", ("kpoint", "kpoints", "kpts", "kp", "kmesh", "mesh", "kgrid")),
                ("encut", ("encut", "cutoff", "ecut")), ("vacuum", ("vac", "vacio", "vacuum")))


_KIND_COMPACT = (("lattice", ("lattice",)), ("kpoints", ("kpoint", "kpts", "kmesh", "kgrid")),
                 ("encut", ("encut", "cutoff")), ("vacuum", ("vacuum", "vacio")))


def _series_kind(name: str) -> str:
    compact = re.sub(r"[^a-z0-9]", "", name.lower())            # «k-points» -> «kpoints»
    for kind, words in _KIND_COMPACT:
        if any(w in compact for w in words):
            return kind
    tokens = [t for t in re.split(r"[^a-z0-9áéíóú]+", name.lower()) if t]
    for kind, words in _KIND_TOKENS:
        if any(t in words or (kind in ("encut", "vacuum") and t.startswith(words[0])) for t in tokens):
            return kind
    return "other"


_NUM = re.compile(r"^[+-]?(\d+([.,]\d*)?|[.,]\d+)([eE][+-]?\d+)?$")
_MESH = re.compile(r"^(\d+)\s*[x×]\s*\d+(\s*[x×]\s*\d+)?$", re.I)


def _to_float(text: str) -> float | None:
    """Número de un texto (admite coma decimal); rechaza 'inf', 'nan' y '1_000'."""
    text = text.strip()
    if not _NUM.match(text):
        return None
    try:
        return float(text.replace(",", "."))
    except ValueError:
        return None


def _name_value(name: str) -> float | None:
    """Valor de una carpeta de serie: el número del nombre, o el primer entero de una malla '6x6x1'."""
    v = _to_float(name)
    if v is None:
        m = _MESH.match(name.strip())
        v = float(m.group(1)) if m else None
    return v


def _natural_key(text: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", text)]


def last_energy(folder: Path) -> float | None:
    """Energía E0 (σ→0, eV) del último paso iónico según OSZICAR; ``None`` si no se puede leer."""
    try:
        text = (folder / "OSZICAR").read_text(errors="replace")
    except OSError:
        return None
    hits = re.findall(r"E0=\s*(\S+)", text)
    return _to_float(hits[-1]) if hits else None


_HINTS = ("relax", "scf", "band", "dos", "phonon", "fonon")


def _find_series(vruns: list[Path]) -> list[Series]:
    """Agrupa los vasprun.xml por carpeta contenedora; ≥3 cálculos sin pistas de otro tipo forman una serie."""
    groups: dict[Path, list[Path]] = {}
    for v in vruns:
        low = v.parent.name.lower()
        kp = v.parent / "KPOINTS"
        if any(h in low for h in _HINTS) or (kp.is_file() and is_line_mode_kpoints(kp)):
            continue
        groups.setdefault(v.parent.parent, []).append(v.parent)
    out = []
    for folder, members in groups.items():
        if len(members) < 3:
            continue
        vals = {m: _name_value(m.name) for m in members}
        members.sort(key=lambda m: (vals[m] is None, vals[m] if vals[m] is not None else 0.0, _natural_key(m.name)))
        kind = _series_kind(folder.name)
        values = [vals[m] for m in members]
        default = len(members) - 1
        if kind == "lattice":
            energies = [last_energy(m) for m in members]
            known = [(e, i) for i, e in enumerate(energies) if e is not None]
            if known:
                default = min(known)[1]
        out.append(Series(folder.name, kind, folder, members, values, default))
    order = {"lattice": 0, "kpoints": 1, "encut": 2, "vacuum": 3, "other": 4}
    out.sort(key=lambda s_: (order[s_.kind], s_.name))
    return out


@dataclass
class ProjectFiles:
    root: Path
    structure: Path | None = None        # CONTCAR (preferido) o POSCAR
    phonons: Path | None = None          # carpeta con phonopy_params.yaml o phonopy_disp.yaml + FORCE_SETS
    bands_vasprun: Path | None = None
    bands_kpoints: Path | None = None
    dos_vasprun: Path | None = None      # DOS (o el SCF si no hay una carpeta de DOS)
    scf_vasprun: Path | None = None
    outcar: Path | None = None
    chgcar: Path | None = None
    notes: list[str] = field(default_factory=list)
    series: list[Series] = field(default_factory=list)

    def found(self) -> dict[str, Path | None]:
        return {"Estructura": self.structure, "Fonones": self.phonons, "Bandas": self.bands_vasprun,
                "KPOINTS (bandas)": self.bands_kpoints, "DOS": self.dos_vasprun, "SCF": self.scf_vasprun,
                "OUTCAR": self.outcar, "CHGCAR": self.chgcar}


def is_line_mode_kpoints(path: Path) -> bool:
    """¿El KPOINTS está en modo línea (camino de alta simetría)? La 3.ª línea empieza por 'L'."""
    try:
        lines = path.read_text(errors="replace").splitlines()
    except OSError:
        return False
    return len(lines) >= 3 and lines[2].strip().lower().startswith("l")


def _walk(root: Path, max_depth: int):
    root = root.resolve()
    for dirpath, dirnames, filenames in os.walk(root):
        depth = len(Path(dirpath).relative_to(root).parts)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith(".") and not d.startswith("disp-"))
        if depth >= max_depth:
            dirnames[:] = []
        yield Path(dirpath), filenames


def scan_folder(root: str | Path, max_depth: int = 4) -> ProjectFiles:
    """Busca los archivos de interés en ``root`` (hasta ``max_depth`` niveles); ignora las carpetas ``disp-*``.

    Los cálculos hermanos que varían un parámetro (malla k, ENCUT, parámetro de red, vacío) se agrupan en
    ``ProjectFiles.series`` en lugar de competir entre sí.
    """
    root = Path(root)
    if not root.is_dir():
        raise NotADirectoryError(f"No es una carpeta: {root}")
    return _scan(root, max_depth)


def _supercell_factor(structure, phonon) -> int:
    """n > 1 si ``structure`` es una supercelda (n celdas unidad) de la celda de los fonones; 0 si no."""
    try:
        nu, ns = len(phonon.unitcell.symbols), len(structure)
        vu = abs(np.linalg.det(np.asarray(phonon.unitcell.cell)))
        if nu and ns > nu and ns % nu == 0:
            n = ns // nu
            if abs(structure.volume / vu - n) < 0.02 * n:
                return n
    except Exception:
        pass
    return 0


def _scan(root: Path, max_depth: int) -> ProjectFiles:
    pf = ProjectFiles(root=root)
    rroot = root.resolve()
    contcars, poscars, params, disps, vruns, outcars, chgcars = [], [], [], [], [], [], []
    ph_struct, ph_outcars = [], []          # de carpetas de Phonopy: son de la supercelda desplazada
    for d, files in _walk(root, max_depth):
        names = set(files)
        in_phonopy = bool(names & {"phonopy_disp.yaml", "phonopy_params.yaml", "phonopy_params.yaml.xz",
                                   "phonopy_params.yaml.gz", "FORCE_SETS"})
        for fname, dest in (("CONTCAR", contcars), ("POSCAR", poscars), ("POSCAR-unitcell", poscars)):
            if fname in names:
                (ph_struct if in_phonopy else dest).append(d / fname)
        if names & {"phonopy_params.yaml", "phonopy_params.yaml.xz", "phonopy_params.yaml.gz"}: params.append(d)
        if "phonopy_disp.yaml" in names and "FORCE_SETS" in names: disps.append(d)
        if "vasprun.xml" in names and not (names & {"phonopy_disp.yaml", "phonopy_params.yaml"}):
            vruns.append(d / "vasprun.xml")
        if "OUTCAR" in names: (ph_outcars if in_phonopy else outcars).append(d / "OUTCAR")
        if "CHGCAR" in names: chgcars.append(d / "CHGCAR")

    pf.series = _find_series(vruns)
    # Ni los cálculos de una serie ni los archivos plantilla de su carpeta contenedora compiten con el resto;
    # la serie solo se usa (como un único cálculo coherente) si no hay nada más.
    in_series = {m for s in pf.series for m in s.members} | {s.folder for s in pf.series}
    outside = lambda paths: [p for p in paths if p.parent not in in_series]
    vruns, contcars, poscars = outside(vruns), outside(contcars), outside(poscars)
    outcars, chgcars = outside(outcars), outside(chgcars)

    def rel(p: Path) -> str:
        """Ruta del archivo relativa a la raíz del proyecto, en minúsculas: las pistas (scf, dos, band…) solo se
        buscan dentro del proyecto y no en los directorios que lo contienen."""
        try:
            parts = p.resolve().relative_to(rroot).parts[:-1]
            return "/".join(parts).lower() if parts else rroot.name.lower()
        except ValueError:
            return p.parent.name.lower()

    def prefer(paths, *tokens):
        paths = sorted(paths, key=lambda p: _natural_key(str(p)))
        for t in tokens:
            for p in paths:
                if t in rel(p):
                    return p
        return paths[0] if paths else None

    pf.structure = prefer(contcars, "relax") or prefer(poscars, "relax") or prefer(ph_struct)
    ph = (params or disps)
    pf.phonons = ph[0] if ph else None
    bands, scf, dos, nohint, relax = [], [], [], [], []
    for v in vruns:
        kp = v.parent / "KPOINTS"
        low = rel(v)
        if "band" in low or (kp.is_file() and is_line_mode_kpoints(kp)):
            bands.append(v)
        elif "dos" in low:
            dos.append(v)
        elif "scf" in low:
            scf.append(v)
        elif "relax" in low:
            relax.append(v)                         # una relajación solo sirve de SCF si no hay otro cálculo
        else:
            nohint.append(v); pf.notes.append(f"{v.parent.name}: sin pista de si es SCF o DOS; se usará como SCF y DOS.")
    scf = scf + nohint + relax
    if bands:
        pf.bands_vasprun = bands[0]
        kp = bands[0].parent / "KPOINTS"
        pf.bands_kpoints = kp if kp.is_file() else None
    pf.scf_vasprun = scf[0] if scf else None
    pf.dos_vasprun = dos[0] if dos else pf.scf_vasprun
    pf.outcar = prefer(outcars, "scf", "relax") or prefer(ph_outcars)
    pf.chgcar = prefer(chgcars, "scf")

    if pf.series:
        for s in pf.series:
            vals = [v for v in s.values if v is not None]
            rango = f"{min(vals):g}–{max(vals):g}" if vals else "sin valores numéricos"
            pf.notes.append(f"Serie «{s.name}»: {len(s.members)} cálculos ({rango}); por defecto «{s.default_member.name}».")
        # Un CONTCAR de una carpeta de Phonopy es de la supercelda desplazada: no cuenta como estructura propia.
        struct_ph = pf.structure is not None and pf.structure in ph_struct
        if (pf.structure is None or struct_ph) and pf.scf_vasprun is None:
            # Un único cálculo coherente: todos los archivos salen de la misma carpeta.
            m = pf.series[0].default_member
            for attr, name in (("structure", "CONTCAR"), ("scf_vasprun", "vasprun.xml"), ("outcar", "OUTCAR"),
                               ("chgcar", "CHGCAR")):
                if (m / name).is_file() and not (attr == "structure" and struct_ph):
                    setattr(pf, attr, m / name)
            if pf.structure is None and (m / "POSCAR").is_file():
                pf.structure = m / "POSCAR"
            pf.dos_vasprun = pf.scf_vasprun
            if pf.series[0].kind == "lattice":
                pf.notes.append("El mínimo de la serie de red no es el óptimo exacto: ajuste la curva E(a) para obtenerlo.")
        if pf.chgcar is None:
            # El CHGCAR puede estar solo en otra serie (p. ej. la de ENCUT): se usa el del cálculo por defecto de esa serie.
            for ser in pf.series:
                if (ser.default_member / "CHGCAR").is_file():
                    pf.chgcar = ser.default_member / "CHGCAR"
                    pf.notes.append(f"CHGCAR tomado de la serie «{ser.name}» («{ser.default_member.name}»).")
                    break
    return pf


@dataclass
class Session:
    """Datos cargados de un proyecto. Los módulos que no se pudieron cargar quedan en None con su error."""
    files: ProjectFiles | None = None
    structure: object = None             # pymatgen Structure
    phonon: object = None                # phonopy.Phonopy
    bands: object = None                 # pymatgen BandStructureSymmLine
    mesh_bands: object = None            # bandas en la malla k uniforme del SCF/DOS (sin camino); solo para la brecha
    tdos: object = None
    cdos: object = None
    chg: object = None                   # pymatgen Chgcar
    magnetization: dict | None = None
    series: list = field(default_factory=list)   # [(Series, filas de structural.series_table)]
    scf_efermi: float | None = None
    bands_efermi_native: float | None = None   # E_F del propio vasprun de bandas (antes de reemplazarla por la del SCF)
    nelect: float | None = None
    discarded_structure: str | None = None      # CONTCAR/POSCAR que era la supercelda de los fonones
    structure_source: str | None = None         # si no hubo CONTCAR/POSCAR: de dónde se tomó la estructura
    lsorbit: bool | None = None                 # espín-órbita (LSORBIT) del cálculo de bandas o del OUTCAR
    errors: dict[str, str] = field(default_factory=dict)

    def available(self) -> dict[str, bool]:
        return {"estructura": self.structure is not None, "fonones": self.phonon is not None,
                "electronico": self.bands is not None or self.tdos is not None,   # lo que se puede dibujar
                "densidad": self.chg is not None}


def load_session(files: ProjectFiles, progress: Callable[[str], None] | None = None) -> Session:
    """Carga todos los módulos posibles; cada fallo se registra en ``errors`` sin detener los demás."""
    from .structural import load_structure
    from .phonons import load_phonopy
    from .electronic import (load_dos, load_vasprun, band_structure_from_vasprun, read_scf_info, load_chgcar,
                             read_magnetization, mesh_band_structure)
    s = Session(files=files)
    say = progress or (lambda m: None)

    def guard(key, fn):
        try:
            return fn()
        except Exception as exc:                  # un módulo dañado no debe bloquear a los demás
            s.errors[key] = f"{type(exc).__name__}: {exc}"
            return None

    if files.structure:
        say("Leyendo la estructura…")
        s.structure = guard("estructura", lambda: load_structure(files.structure))
    if files.series:
        say("Leyendo las series de cálculos…")
        from .structural import series_table
        for ser in files.series:
            rows = guard("series", lambda ser=ser: series_table(ser))
            if rows is not None:
                s.series.append((ser, rows))
    if files.phonons:
        say("Leyendo los fonones…")
        s.phonon = guard("fonones", lambda: load_phonopy(files.phonons))
    if s.structure is not None and s.phonon is not None:
        n_super = _supercell_factor(s.structure, s.phonon)
        if n_super:
            s.discarded_structure = f"{files.structure.name} de {files.structure.parent.name}: supercelda de {len(s.structure)} átomos"
            s.structure = None
    if files.scf_vasprun:
        say("Leyendo el SCF…")
        info = guard("electronico", lambda: read_scf_info(files.scf_vasprun))
        if info:
            s.scf_efermi, s.nelect = info
    if files.bands_vasprun and files.bands_kpoints:
        say("Leyendo las bandas…")
        def _bands():
            vr = load_vasprun(files.bands_vasprun, parse_projected=True)
            s.bands_efermi_native = float(vr.efermi)
            soc = vr.parameters.get("LSORBIT")
            s.lsorbit = bool(soc) if soc is not None else s.lsorbit
            return band_structure_from_vasprun(vr, files.bands_kpoints, s.scf_efermi)
        s.bands = guard("electronico", _bands)
    elif files.bands_vasprun:
        s.errors["electronico"] = "Se encontraron bandas pero falta el KPOINTS en modo línea."
    if files.dos_vasprun:
        say("Leyendo la DOS…")
        def _dos():
            vr = load_vasprun(files.dos_vasprun, parse_projected=True)
            s.tdos, s.cdos = vr.tdos, getattr(vr, "complete_dos", None)
            soc = vr.parameters.get("LSORBIT")
            if s.lsorbit is None and soc is not None:
                s.lsorbit = bool(soc)
            if s.bands is None:                       # sin camino k: brecha aproximada en la malla del cálculo
                try:
                    s.mesh_bands = mesh_band_structure(vr, s.scf_efermi)
                except Exception:
                    s.mesh_bands = None
        guard("electronico", _dos)
    if s.bands is None and s.mesh_bands is None and files.scf_vasprun:
        try:
            s.mesh_bands = mesh_band_structure(load_vasprun(files.scf_vasprun), s.scf_efermi)
        except Exception:
            s.mesh_bands = None
    if files.chgcar:
        say("Leyendo la densidad de carga…")
        s.chg = guard("densidad", lambda: load_chgcar(files.chgcar))
    if files.outcar:
        s.magnetization = guard("electronico", lambda: read_magnetization(files.outcar, s.structure))
        if s.lsorbit is None and s.magnetization:
            s.lsorbit = s.magnetization.get("lsorbit")
    if s.structure is None:                        # estructura alternativa desde otro módulo
        for name, src in (("la densidad de carga", s.chg), ("las bandas", s.bands), ("la malla k del SCF/DOS", s.mesh_bands)):
            st = getattr(src, "structure", None)
            if st is not None:
                s.structure, s.structure_source = st, name
                break
    if s.structure is None and s.phonon is not None:
        try:
            from pymatgen.core import Structure
            u = s.phonon.unitcell
            s.structure = Structure(u.cell, list(u.symbols), u.scaled_positions)
            s.structure_source = "la celda del cálculo de fonones" + (
                f" (se descartó {s.discarded_structure}; la celda unidad es de {len(s.structure)} átomos)"
                if getattr(s, "discarded_structure", None) else "")
        except Exception:
            pass
    return s
