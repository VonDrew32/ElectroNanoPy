"""Series de cálculos (barridos de malla k, ENCUT, parámetro de red, vacío): tabla estructural y resumen.

Describe los datos que ya están en las carpetas (estructura, energía, esfuerzo).
El ajuste de Birch–Murnaghan está en ``eos.py`` y se aplica a la serie de red.
"""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np
from .io import load_structure
from . import descriptors as D
from .eos import fit_birch_murnaghan

_KB = re.compile(r"in kB\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")


def read_stress(folder: Path) -> tuple[float, float, float] | None:
    """Esfuerzo (kB) xx, yy, zz del último paso iónico, leído de la línea «in kB» del OUTCAR."""
    try:
        text = (folder / "OUTCAR").read_text(errors="replace")
    except OSError:
        return None
    hits = _KB.findall(text)
    if not hits:
        return None
    try:
        return tuple(float(x) for x in hits[-1][:3])
    except ValueError:
        return None


def series_table(series, tolerance: float = 1.15) -> list[dict]:
    """Una fila por cálculo de la serie: parámetro, estructura, a, c, d_enlace, vacío, energía y esfuerzo.

    Un cálculo que no se puede leer produce una fila con ``error`` en lugar de detener la tabla.
    """
    from ..project import last_energy
    rows = []
    for member, value in zip(series.members, series.values):
        row = {"name": member.name, "value": value, "path": member, "structure": None, "error": None,
               "energy": last_energy(member), "stress": read_stress(member)}
        f = next((member / n for n in ("CONTCAR", "POSCAR") if (member / n).is_file()), None)
        try:
            if f is None:
                raise FileNotFoundError("sin CONTCAR/POSCAR")
            s = load_structure(f)
            d = [b["distance"] for b in D.bond_table(s, tolerance)]
            lat = D.lattice_parameters(s)
            row.update(structure=s, a=lat["a"], b=lat["b"], c=lat["c"], gamma=lat["gamma"], volume=lat["volume"],
                       bond_min=min(d) if d else None, vacuum=D.vacuum_thickness(s))
        except Exception as exc:
            row["error"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)
    return rows


def series_summary(series, rows: list[dict]) -> dict:
    """Cifras clave de la serie para el reporte (energías relativas en meV).

    * ``kpoints``/``encut``/``vacuum``: diferencia entre los dos últimos puntos y respecto al más convergido.
    * ``lattice``: mínimo de la serie, valor de ``a`` donde el esfuerzo in-plane σ_xx se anula (interpolación lineal
      entre vecinos) y ajuste de Birch–Murnaghan (``bm``, ver ``structural.eos``; ``None`` con < 5 puntos).
    """
    ok = [r for r in rows if r["energy"] is not None and r["value"] is not None]
    out = {"n": len(rows), "n_energy": len(ok), "kind": series.kind}
    if not ok:
        return out
    e = np.array([r["energy"] for r in ok]); v = np.array([r["value"] for r in ok])
    if series.kind == "lattice":
        i = int(e.argmin())
        out.update(a_min=float(v[i]), e_min=float(e[i]), e_span_mev=float((e.max() - e.min()) * 1000))
        st = [(r["value"], r["stress"][0]) for r in ok if r["stress"] is not None]
        for (a1, s1), (a2, s2) in zip(st, st[1:]):
            if s1 == 0:                                        # el esfuerzo se anula justo en un punto de la malla
                out["a_zero_stress"] = float(a1); break
            if s1 * s2 < 0:
                out["a_zero_stress"] = float(a1 + (a2 - a1) * s1 / (s1 - s2)); break
        else:
            if st and st[-1][1] == 0:
                out["a_zero_stress"] = float(st[-1][0])
        first = next((r["structure"] for r in rows if r.get("structure") is not None), None)
        out["is_2d"] = bool(first is not None and D.is_two_dimensional(first))
        out["bm"] = fit_birch_murnaghan(rows, is_2d=out["is_2d"])
    else:
        out["last_step_mev"] = float((e[-1] - e[-2]) * 1000) if len(e) > 1 else None
        out["spread_last3_mev"] = float((e[-3:].max() - e[-3:].min()) * 1000)
        out["monotonic"] = bool(np.all(np.diff(e) <= 1e-6) or np.all(np.diff(e) >= -1e-6))
    return out
