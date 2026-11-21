"""Análisis electrónico: brecha, puntos de Dirac, velocidad de Fermi, DOS en E_F, pesos orbitales y magnetismo."""
from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin

HBAR_EV_S = 6.582119569e-16     # eV·s


@dataclass
class GapReport:
    gap: float                  # eV (0 si es metálico o semimetálico)
    direct: bool | None
    is_metal: bool
    vbm: float | None           # eV (respecto a E_F de la estructura de bandas)
    cbm: float | None
    transition: str | None      # p. ej. 'Γ-X'
    message: str = ""


def _k_name(bs, idx: int) -> str:
    lab = _k_label(bs, idx)
    return lab if lab else "(" + ",".join(f"{x:.3f}" for x in bs.kpoints[idx].frac_coords) + ")"


def _gap_from_occupation(bs, nelect: float, gap_tol: float, magnetization: float | None = None,
                         soc: bool = False) -> GapReport | None:
    """Brecha por conteo de bandas ocupadas (no depende de dónde caiga E_F). ``None`` si no se puede aplicar.

    - Sin polarización de espín: NELECT/2 bandas ocupadas (un número impar de electrones indica metal).
    - Con espín-órbita (``soc``): un solo canal con NELECT bandas ocupadas.
    - Con polarización de espín: N↑ = (N + M)/2 y N↓ = (N − M)/2, lo que exige el momento total M (µB); sin M se
      devuelve ``None`` (se usa el método de E_F). Si N↑ o N↓ no son enteros el sistema es metálico.
    La brecha es min(banda de conducción) − max(banda de valencia) sobre todos los k y canales.
    """
    spins = list(bs.bands)
    chans = []                                               # (energías ordenadas (nbandas, nk), nocupadas)
    if len(spins) == 2:
        if magnetization is None:
            return None
        counts = [(nelect + magnetization) / 2.0, (nelect - magnetization) / 2.0]
        if any(abs(c - round(c)) > 0.05 for c in counts):
            return GapReport(0.0, None, True, None, None, None,
                             "Sin brecha: " + ("con NELECT impar y sin momento magnético las bandas no pueden estar llenas (metal)."
                                              if abs(magnetization) < 0.05 else
                                              "el momento magnético no corresponde a bandas llenas (metal)."))
        for sp, c in zip((Spin.up, Spin.down), counts):
            chans.append((np.sort(np.asarray(bs.bands[sp]), axis=0), int(round(c))))
    else:
        e = np.sort(np.asarray(bs.bands[spins[0]]), axis=0)
        nocc_f = nelect if soc else nelect / 2.0
        if abs(nocc_f - round(nocc_f)) > 1e-6:
            return GapReport(0.0, None, True, None, None, None,
                             "Sin brecha: número impar de electrones (metal).")
        chans.append((e, int(round(nocc_f))))
    vb = [(e[n - 1].max(), int(e[n - 1].argmax())) for e, n in chans if 1 <= n <= e.shape[0]]
    cb = [(e[n].min(), int(e[n].argmin())) for e, n in chans if 0 <= n < e.shape[0]]
    if not vb or not cb:
        return None                                          # no hay banda de valencia o de conducción calculada
    (vbm, iv), (cbm, ic) = max(vb), min(cb)
    gap = float(cbm - vbm)
    if gap <= gap_tol:
        what = "se tocan" if gap > -gap_tol else "se solapan"
        return GapReport(0.0, None, True, None, None, None,
                         f"Brecha nula (< {gap_tol} eV): las bandas {what} (semimetal o metal). "
                         "Revise los puntos de Dirac.")
    direct = bool(np.allclose(bs.kpoints[iv].frac_coords, bs.kpoints[ic].frac_coords, atol=1e-4))
    trans = f"{_k_name(bs, iv)}-{_k_name(bs, ic)}"
    return GapReport(gap, direct, False, float(vbm) - bs.efermi, float(cbm) - bs.efermi, trans,
                     f"Brecha {'directa' if direct else 'indirecta'} de {gap:.3f} eV ({trans}).")


def band_gap(bs, gap_tol: float = 0.01, nelect: float | None = None, magnetization: float | None = None,
             soc: bool = False) -> GapReport:
    """Brecha de energía (la del cálculo de origen: no se aplican correcciones).

    Con ``nelect`` (electrones de valencia, de NELECT) y espín no polarizado la brecha se calcula contando bandas
    ocupadas, lo que no depende de dónde caiga E_F: es el método recomendado (en el grafeno E_F coincide con el
    punto de Dirac con una diferencia de ~10⁻⁵ eV y el método de E_F es inestable). Sin ``nelect`` se usa E_F:
    si hay bandas cruzando E_F o la brecha es menor que ``gap_tol`` eV (p. ej. bandas que solo se tocan), se
    informa brecha nula y ``is_metal=True`` (metal o semimetal).

    ``magnetization`` (momento total, µB) permite el conteo por canal de espín en cálculos polarizados (ISPIN = 2) y
    ``soc=True`` indica un cálculo con espín-órbita (una sola banda por estado, NELECT bandas ocupadas).
    """
    if nelect is not None:
        g = _gap_from_occupation(bs, nelect, gap_tol, magnetization, soc)
        if g is not None:
            return g
    if bs.is_metal():
        return GapReport(0.0, None, True, None, None, None,
                         "Sin brecha: hay bandas cruzando el nivel de Fermi (metal o semimetal).")
    g = bs.get_band_gap()
    if g["energy"] < gap_tol:
        return GapReport(0.0, None, True, None, None, None,
                         f"Brecha nula (< {gap_tol} eV): las bandas se tocan en E_F (semimetal o metal). "
                         "Revise los puntos de Dirac.")
    vbm, cbm = bs.get_vbm()["energy"], bs.get_cbm()["energy"]
    kind = "directa" if g["direct"] else "indirecta"
    return GapReport(float(g["energy"]), bool(g["direct"]), False, float(vbm), float(cbm),
                     g.get("transition"), f"Brecha {kind} de {g['energy']:.3f} eV ({g.get('transition')}).")


def _k_label(bs, idx: int) -> str | None:
    for lab, kp in bs.labels_dict.items():
        if np.allclose(kp.frac_coords, bs.kpoints[idx].frac_coords, atol=1e-4):
            return lab
    return None


def dirac_points(bs, tol: float = 0.05, ef_window: float = 0.1, max_degenerate: float = 0.3) -> list[dict]:
    """Puntos del camino donde dos bandas consecutivas se tocan (separación < tol) cerca de E_F.

    Devuelve una lista de dicts con el índice, la etiqueta (si es un punto etiquetado), las coordenadas
    fraccionarias, la energía del cruce (respecto a E_F) y la separación. Sirve para validar el cono de Dirac
    del grafeno en K. ``tol`` y ``ef_window`` en eV.
    """
    out = []
    for spin, e in bs.bands.items():
        e = np.sort(np.asarray(e), axis=0)                           # (nbandas, nk)
        gaps = np.diff(e, axis=0)                                    # (nbandas-1, nk)
        mid = (e[:-1] + e[1:]) / 2 - bs.efermi
        hit = (gaps < tol) & (np.abs(mid) < ef_window)               # (nbandas-1, nk)
        hit &= ((gaps < tol).mean(axis=1) <= max_degenerate)[:, None]   # descarta pares degenerados en todo el camino
        ks = np.nonzero(hit.any(axis=0))[0]
        if ks.size == 0:
            continue
        groups = np.split(ks, np.nonzero(np.diff(ks) > 1)[0] + 1)    # índices k contiguos = un mismo punto
        for grp in groups:
            band = int(np.argmin(np.where(hit[:, grp], gaps[:, grp], np.inf).min(axis=1)))
            sub = gaps[band, grp]
            i = int(grp[np.argmin(sub)])
            out.append({"spin": spin.name, "k_index": i, "label": _k_label(bs, i) or _k_label(bs, int(grp[0])),
                        "frac_coords": tuple(float(x) for x in bs.kpoints[i].frac_coords),
                        "energy": float(mid[band, i]), "separation": float(gaps[band, i]), "band": band})
    return out


def fermi_velocity(bs, dirac: dict, window: float = 0.1) -> dict:
    """Pendiente lineal de las bandas alrededor de un punto de Dirac.

    Ajusta |E - E_D| = ħ·v_F·|Δk| (recta por el origen) con los puntos a menos de ``window`` (1/Å) del
    punto de Dirac, en las dos bandas y a ambos lados. Devuelve ħv_F en eV·Å y v_F en m/s.
    Promedia las direcciones del camino que llegan al punto (no distingue anisotropía).
    """
    spin = Spin[dirac["spin"]]
    e = np.sort(np.asarray(bs.bands[spin]), axis=0)
    b, i0 = dirac["band"], dirac["k_index"]
    d = np.asarray(bs.distance)
    ed = (e[b, i0] + e[b + 1, i0]) / 2
    x = np.abs(d - d[i0])
    sel = (x > 1e-9) & (x <= window)
    if sel.sum() < 4:
        raise ValueError("Muy pocos puntos k cerca del punto de Dirac; aumente window o la densidad del camino.")
    xs = np.concatenate([x[sel], x[sel]])
    ys = np.concatenate([e[b + 1, sel] - ed, ed - e[b, sel]])
    slope = float(np.sum(xs * ys) / np.sum(xs * xs))
    return {"hbar_vF_eV_A": slope, "vF_m_s": slope * 1e-10 / HBAR_EV_S, "n_points": int(sel.sum()) * 2}


def dos_at_fermi(dos, window: float = 0.1) -> float:
    """Densidad de estados media (estados/eV por celda) en E_F ± ``window`` eV (suma de espines)."""
    e = np.asarray(dos.energies) - dos.efermi
    tot = sum(np.asarray(v) for v in dos.densities.values())
    m = np.abs(e) <= window
    return float(tot[m].mean()) if m.any() else float("nan")


def dos_energy_step(dos) -> float:
    """Paso medio de la rejilla de energía de la DOS (eV). Con NEDOS = 301 en un rango de ~36 eV es ~0.12 eV."""
    e = np.asarray(dos.energies)
    return float(np.diff(e).mean()) if e.size > 1 else float("nan")


def classify(bs, dos=None, gap_tol: float = 0.01, dos_tol: float = 0.1, nelect: float | None = None,
             max_dos_step: float = 0.05, magnetization: float | None = None, soc: bool = False) -> str:
    """Clasificación orientativa: semiconductor/aislante, semimetal de Dirac o metal.

    Regla heurística: brecha > gap_tol -> semiconductor/aislante; si no, con puntos de Dirac y (si se da la
    DOS) densidad en E_F < dos_tol estados/eV -> semimetal de Dirac; en otro caso metal. La DOS solo se usa si su
    rejilla de energía es más fina que ``max_dos_step`` eV (con una rejilla gruesa el valor en E_F no es fiable).
    """
    g = band_gap(bs, gap_tol, nelect, magnetization, soc)
    if not g.is_metal and g.gap > gap_tol:
        return "semiconductor/aislante"
    dos_ok = dos is not None and dos_energy_step(dos) <= max_dos_step
    if dirac_points(bs) and (not dos_ok or dos_at_fermi(dos) < dos_tol):
        return "semimetal de Dirac"
    return "metal"


_ORB = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2"]     # orden de VASP/pymatgen


def orbital_weights(bs, orbitals=("pz",), element: str | None = None, structure: Structure | None = None,
                    spin: Spin = Spin.up, normalize: bool = True) -> np.ndarray:
    """Peso (0-1) de ciertos orbitales (y elemento) en cada banda y punto k: matriz (nbandas, nk).

    Requiere un cálculo con LORBIT y proyecciones cargadas (``parse_projected=True``).

    VASP proyecta solo dentro de las esferas PAW, así que la suma de todas las proyecciones de una banda es < 1
    (≈ 0.5 en el grafeno). Con ``normalize=True`` (por defecto) el peso es relativo al total proyectado de esa banda
    y punto k (carácter orbital); con ``False`` se devuelve la proyección cruda de VASP.
    """
    if not getattr(bs, "projections", None):
        raise ValueError("La estructura de bandas no tiene proyecciones (use LORBIT y parse_projected=True).")
    proj = np.asarray(bs.projections[spin])                         # (nbandas, nk, norb, nion)
    norb = proj.shape[2]
    names = _ORB if norb >= len(_ORB) else ["s", "p", "d"][:norb] if norb <= 3 else _ORB[:norb]
    idx = []
    for o in orbitals:
        if o not in _ORB and o not in ("p", "d"):
            raise ValueError(f"Orbital desconocido: {o}. Opciones: {_ORB}")
        if o not in names:
            raise ValueError(f"Las proyecciones solo contienen {names} (¿LORBIT = 10?); para «{o}» calcule con LORBIT = 11.")
        idx.append(names.index(o))
    s = structure or bs.structure
    ions = [i for i, site in enumerate(s) if element is None or site.specie.symbol == element]
    if not ions:
        raise ValueError(f"No hay átomos de {element} en la estructura.")
    w = proj[:, :, idx, :][:, :, :, ions].sum(axis=(2, 3))
    if normalize:
        total = proj.sum(axis=(2, 3))
        w = np.divide(w, total, out=np.zeros_like(w), where=total > 1e-9)
    return w


def read_magnetization(outcar: str | Path, structure: Structure | None = None) -> dict:
    """Momentos magnéticos del último bloque ``magnetization (x)`` de un OUTCAR (µB).

    Devuelve ``{'spin_polarized', 'total', 'per_atom'}``. Si el cálculo no es polarizado en espín
    (no hay bloque), ``spin_polarized`` es False y los valores son None.
    """
    path = Path(outcar)
    if not path.is_file():
        raise FileNotFoundError(f"No existe el archivo: {path}")
    lines = path.read_text(errors="replace").splitlines()
    text = "\n".join(lines)
    m_ispin = re.search(r"\bISPIN\s*=\s*(\d)", text)
    m_soc = re.search(r"\bLSORBIT\s*=\s*([TF])", text)
    ispin = int(m_ispin.group(1)) if m_ispin else None
    lsorbit = (m_soc.group(1) == "T") if m_soc else None
    # «number of electron  N magnetization  M» (ISPIN = 2, aparece sin LORBIT); sin polarizar el valor va vacío
    tot_lines = re.findall(r"number of electron\s+\S+\s+magnetization\s+(-?\d+\.\d+)", text)
    starts = [i for i, ln in enumerate(lines) if re.match(r"\s*magnetization \(x\)", ln)]
    if not starts:
        if tot_lines or ispin == 2:
            return {"spin_polarized": True, "total": float(tot_lines[-1]) if tot_lines else None,
                    "per_atom": None, "ispin": ispin, "lsorbit": lsorbit}
        return {"spin_polarized": False, "total": None, "per_atom": None, "ispin": ispin, "lsorbit": lsorbit}
    def block(start: int):
        """Filas por átomo y total de un bloque «magnetization (x|y|z)» que empieza en la línea ``start``."""
        i = start + 1
        while i < len(lines) and not lines[i].lstrip().startswith("# of ion"):
            i += 1
        header = lines[i].replace("# of ion", "").split()            # ['s','p','d','tot']
        i += 2                                                       # salta la línea de guiones
        rows = []
        while i < len(lines) and not lines[i].startswith("---"):
            parts = lines[i].split()
            if parts and not parts[0].isdigit():
                break                                # otra sección (con un solo átomo VASP no escribe la suma)
            if parts:
                rows.append({"index": int(parts[0]) - 1, **{h: float(v) for h, v in zip(header, parts[1:])}})
            i += 1
        tot = lines[i + 1].split() if i < len(lines) and lines[i].startswith("---") and i + 1 < len(lines) else []
        total = float(tot[-1]) if tot and tot[0] == "tot" else float(sum(r["tot"] for r in rows))
        return rows, total

    rows, total = block(starts[-1])
    noncollinear = False
    comps = {}
    for axis in ("y", "z"):                                          # no colineal: bloques (y) y (z) tras el (x)
        nxt = [i for i, ln in enumerate(lines) if re.match(rf"\s*magnetization \({axis}\)", ln) and i > starts[-1]]
        if nxt and nxt[0] - starts[-1] < 4 * (len(rows) + 12):
            comps[axis] = block(nxt[0])
    if comps:
        noncollinear = True
        vec = [total, comps.get("y", ([], 0.0))[1], comps.get("z", ([], 0.0))[1]]
        for k, r in enumerate(rows):
            m = [r["tot"]] + [comps[a][0][k]["tot"] if a in comps and k < len(comps[a][0]) else 0.0 for a in ("y", "z")]
            r["mx"], r["my"], r["mz"] = m
            r["tot"] = float(np.linalg.norm(m))              # módulo del momento local (el signo lo dan mx, my, mz)
            r["tot_signed"] = m[int(np.argmax(np.abs(m)))]   # componente dominante, con signo
        total = float(np.linalg.norm(vec))
    if structure is not None and len(structure) >= len(rows):       # si no, el OUTCAR es de otra celda
        for r in rows:
            r["element"] = structure[r["index"]].specie.symbol
    return {"spin_polarized": True, "total": total, "per_atom": rows, "ispin": ispin, "lsorbit": lsorbit,
            "noncollinear": noncollinear}
