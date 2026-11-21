"""Ajuste de la ecuación de estado de Birch–Murnaghan (3.er orden) a una serie de parámetro de red.

E(V) = E0 + 9·V0·B0/16 · { [(V0/V)^(2/3) − 1]^3 · B0' + [(V0/V)^(2/3) − 1]^2 · [6 − 4·(V0/V)^(2/3)] }

Se ajusta E frente al volumen de la celda leído de cada CONTCAR. El parámetro de red óptimo a0 se obtiene con la
relación entre a y V de la propia serie (exponente n = d ln V / d ln a: 3 si la celda escala por igual en las tres
direcciones; 2 si solo cambia el plano, como en una monocapa con c fijo).
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import curve_fit

EV_A3_TO_GPA = 160.21766208


def birch_murnaghan(v, e0, v0, b0, b0p):
    eta = (v0 / np.asarray(v, dtype=float)) ** (2.0 / 3.0)
    return e0 + 9.0 * v0 * b0 / 16.0 * ((eta - 1) ** 3 * b0p + (eta - 1) ** 2 * (6 - 4 * eta))


def fit_birch_murnaghan(rows: list[dict], min_points: int = 5, is_2d: bool = False) -> dict | None:
    """Ajusta E(V) a las filas de ``series_table``; ``None`` si no hay suficientes puntos para 4 parámetros.

    Devuelve ``a0`` (valor del parámetro de la serie en el mínimo; es el que se escribe en el POSCAR), ``lattice0``
    (a, b, c de los CONTCAR interpolados en V0), ``v0``, ``e0`` (eV), ``b0_gpa`` (``None`` si ``is_2d`` —con vacío fijo no es un módulo físico— o si la celda no escala por igual en las tres
    direcciones, |n − 3| > 0.15: entonces dV/da no es 3V/a y B₀ no sería el módulo volumétrico),
    ``b0p``, ``rmse_mev``, ``n_points``, ``exponent``, ``extrapolated`` (V0 fuera del rango de la serie) y la curva
    ajustada (``curve_a``, ``curve_e``) para dibujarla frente a ``a``.
    """
    ok = [r for r in rows if r.get("energy") is not None and r.get("volume") and r.get("a")]
    if len(ok) < min_points:
        return None
    ok.sort(key=lambda r: r["volume"])
    # Eje del parámetro: el valor de la carpeta (lo que se escribe en el POSCAR); si falta, el a del CONTCAR.
    # Es el mismo eje en el que se dibujan los puntos, así la curva no se desplaza cuando la celda cambia de
    # forma al relajarse (a del CONTCAR ≠ parámetro de la carpeta).
    a = np.array([r["value"] if r.get("value") is not None else r["a"] for r in ok], dtype=float)
    v = np.array([r["volume"] for r in ok]); e = np.array([r["energy"] for r in ok])
    n = float(np.polyfit(np.log(a), np.log(v), 1)[0])
    c = np.polyfit(v, e, 2)                                    # semilla: parábola
    v0 = -c[1] / (2 * c[0]) if c[0] > 0 else v[e.argmin()]
    p0 = [e.min(), v0, max(2 * c[0] * v0, 1e-3), 4.0]
    try:
        p, _ = curve_fit(birch_murnaghan, v, e, p0=p0, maxfev=50000)
    except (RuntimeError, ValueError):
        return None
    e0, v0, b0, b0p = (float(x) for x in p)
    if not (v0 > 0 and b0 > 0):
        return None
    rmse = float(np.sqrt(np.mean((birch_murnaghan(v, *p) - e) ** 2)) * 1000)
    k = np.polyfit(np.log(v), np.log(a), 1)                    # ln a = k0·ln V + k1 (interpolación entre vecinos)
    a_of_v = lambda vv: np.exp(np.polyval(k, np.log(vv)))
    vg = np.linspace(v.min(), v.max(), 200)
    lattice0 = {}                                             # a, b, c del CONTCAR interpolados en V0
    for key in ("a", "b", "c"):
        vals = np.array([r.get(key) or np.nan for r in ok], dtype=float)
        if np.all(np.isfinite(vals)) and np.all(vals > 0):
            lattice0[key] = float(np.exp(np.polyval(np.polyfit(np.log(v), np.log(vals), 1), np.log(v0))))
    return {"a0": float(a_of_v(v0)), "lattice0": lattice0, "v0": v0, "e0": e0, "b0_gpa": None if (is_2d or abs(n - 3.0) > 0.15) else b0 * EV_A3_TO_GPA, "b0p": b0p,
            "rmse_mev": rmse, "n_points": len(ok), "exponent": n, "extrapolated": bool(not v.min() <= v0 <= v.max()),
            "curve_a": a_of_v(vg), "curve_e": birch_murnaghan(vg, *p)}
