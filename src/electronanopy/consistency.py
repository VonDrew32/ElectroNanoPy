"""Consistencia entre cálculos (OE4): ¿los cálculos de un mismo proyecto son comparables entre sí?

Las bandas, la DOS y los fonones solo tienen sentido si se calcularon con el mismo corte de energía, el mismo
funcional y los mismos pseudopotenciales que el SCF y la relajación. Este módulo lee esos ajustes del
``vasprun.xml`` de cada cálculo (sin cargarlo completo) y avisa de las diferencias. No usa ML: son comparaciones directas.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re
import xml.etree.ElementTree as ET


@dataclass
class CalcSettings:
    source: str
    encut: float | None = None            # eV: corte efectivo (ENMAX del vasprun)
    ispin: int | None = None
    nelect: float | None = None
    functional: str | None = None         # derivado de los pseudopotenciales y de GGA/METAGGA
    potcars: tuple = ()
    lhfcalc: bool | None = None
    ldau: bool | None = None
    lsorbit: bool | None = None
    kmesh: tuple | None = None            # divisiones de la malla; None en modo línea o con lista explícita
    natoms: int | None = None             # átomos de la celda (la supercelda de fonones tiene más que la celda unidad)


def _bool(text):
    return None if text is None else text.strip().upper().startswith("T")


def _functional(potcars, gga, metagga, lhf):
    kinds = sorted({re.sub(r"^PAW_?", "", p.split()[0]) if p else "?" for p in potcars})
    base = "/".join(kinds) if kinds else "desconocido"
    extra = []
    if gga and gga.strip() not in ("--", "", "PE"):
        extra.append(f"GGA={gga.strip()}")
    if metagga and metagga.strip() not in ("--", "None", ""):
        extra.append(f"METAGGA={metagga.strip()}")
    if lhf:
        extra.append("híbrido")
    return base + (" + " + ", ".join(extra) if extra else "")


def read_settings(vasprun: str | Path, label: str | None = None) -> CalcSettings:
    """Lee los ajustes de un ``vasprun.xml`` deteniéndose al terminar ``atominfo`` (no carga los pasos iónicos)."""
    p = Path(vasprun)
    out = CalcSettings(source=label or str(p))
    params: dict[str, str] = {}
    incar: dict[str, str] = {}
    for _, el in ET.iterparse(str(p), events=("end",)):
        if el.tag == "incar":
            incar = {i.get("name"): (i.text or "").strip() for i in el.iter("i")}
        elif el.tag == "kpoints" and out.kmesh is None:
            gen = el.find("generation")
            if gen is not None:
                v = gen.find("v[@name='divisions']")
                if v is not None and v.text:
                    out.kmesh = tuple(int(float(x)) for x in v.text.split())
        elif el.tag == "parameters":
            params = {i.get("name"): (i.text or "").strip() for i in el.iter("i")}
        elif el.tag == "atominfo":
            pots = []
            for arr in el.iter("array"):
                if arr.get("name") == "atomtypes":
                    pots = [rc.findall("c")[-1].text.strip() for rc in arr.iter("rc")]
            out.potcars = tuple(pots)
            na = el.find("atoms")
            out.natoms = int(na.text) if na is not None and na.text and na.text.strip().isdigit() else None
            break
    def num(key, cast=float):
        for d in (params, incar):
            if key in d:
                try:
                    return cast(d[key])
                except ValueError:
                    return None
        return None
    out.encut = num("ENMAX") or num("ENCUT")
    out.ispin = num("ISPIN", int)
    out.nelect = num("NELECT")
    out.lhfcalc = _bool(params.get("LHFCALC"))
    out.ldau = _bool(params.get("LDAU"))
    out.lsorbit = _bool(params.get("LSORBIT"))
    out.functional = _functional(out.potcars, params.get("GGA") or incar.get("GGA"),
                                 params.get("METAGGA") or incar.get("METAGGA"), out.lhfcalc)
    return out


def compare_settings(named: dict[str, CalcSettings], encut_tol: float = 1.0) -> list[str]:
    """Compara los ajustes entre cálculos y devuelve las advertencias (lista vacía si todo es coherente)."""
    names = list(named)
    if len(names) < 2:
        return []
    warns: list[str] = []

    def values(attr):
        return {n: getattr(named[n], attr) for n in names if getattr(named[n], attr) not in (None, ())}

    def differ(attr, fmt=str):
        v = values(attr)
        if len({fmt(x) for x in v.values()}) > 1:
            return ", ".join(f"{n}: {fmt(x)}" for n, x in v.items())
        return None

    enc = values("encut")
    if enc and max(enc.values()) - min(enc.values()) > encut_tol:
        warns.append("Corte de energía distinto entre cálculos (" +
                     ", ".join(f"{n}: {x:.0f} eV" for n, x in enc.items()) +
                     "): las energías y las bandas no son directamente comparables.")
    for attr, text in (("functional", "Funcional distinto entre cálculos"),
                       ("potcars", "Pseudopotenciales distintos entre cálculos")):
        d = differ(attr, (lambda x: " ".join(x)) if attr == "potcars" else str)
        if d:
            warns.append(f"{text} ({d}).")
    for attr, text in (("ispin", "ISPIN distinto"), ("nelect", "Número de electrones (NELECT) distinto"),
                       ("lsorbit", "Acoplamiento espín-órbita distinto"), ("ldau", "LDAU (corrección U) distinto")):
        if attr == "nelect":
            nat = {n: named[n].natoms for n in names if named[n].nelect is not None}
            if nat and all(nat.values()) and len(set(nat.values())) > 1:
                # celdas de distinto tamaño (p. ej. supercelda de fonones): se compara por átomo
                per = {n: named[n].nelect / nat[n] for n in nat}
                if max(per.values()) - min(per.values()) > 1e-3:
                    warns.append("Número de electrones por átomo (NELECT/N) distinto entre cálculos (" +
                                 ", ".join(f"{n}: {x:.2f}" for n, x in per.items()) + ").")
                continue
        d = differ(attr)
        if d:
            warns.append(f"{text} entre cálculos ({d}).")
    return warns


def phonon_displacement_vasprun(folder: str | Path) -> Path | None:
    """Primer ``vasprun.xml`` de desplazamiento (``disp-*``) bajo la carpeta de fonones, si existe."""
    f = Path(folder)
    for pattern in ("disp-*/vasprun.xml", "*/disp-*/vasprun.xml", "disp-*/*/vasprun.xml", "vasprun.xml"):
        hits = sorted(f.glob(pattern))
        if hits:
            return hits[0]
    return None


def gather_settings(files) -> dict[str, CalcSettings]:
    """Ajustes de cada cálculo detectado (``ProjectFiles``); omite los archivos ilegibles y los repetidos."""
    candidates = []
    ph_dir = Path(files.phonons).resolve() if files.phonons is not None else None
    if files.structure is not None and Path(files.structure).resolve().parent != ph_dir:   # no la supercelda de fonones
        candidates.append(("Relajación", Path(files.structure).parent / "vasprun.xml"))
    candidates += [("SCF", files.scf_vasprun), ("Bandas", files.bands_vasprun), ("DOS", files.dos_vasprun)]
    if files.phonons is not None:
        candidates.append(("Fonones", phonon_displacement_vasprun(files.phonons)))
    named, seen = {}, set()
    for label, path in candidates:
        if path is None or not Path(path).is_file():
            continue
        key = Path(path).resolve()
        if key in seen:
            continue
        seen.add(key)
        try:
            named[label] = read_settings(path, label)
        except Exception:
            continue
    return named
