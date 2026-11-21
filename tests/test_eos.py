"""Ajuste de Birch–Murnaghan con parámetros conocidos."""
import numpy as np
import pytest
from electronanopy.project import scan_folder, load_session
from electronanopy.report import analyze
from electronanopy.structural import birch_murnaghan, fit_birch_murnaghan
from electronanopy.structural.eos import EV_A3_TO_GPA
from test_series import _calc

E0, V0, B0, B0P = -10.0, 100.0, 0.6, 4.2          # B0 = 0.6 eV/Å³ ≈ 96 GPa


def _filas(vs, expo=3, c=1.0):
    """Filas tipo series_table para una celda que escala con a (cúbica: V = a³; monocapa de c fijo: V = a²·c)."""
    out = []
    for v in vs:
        a = v ** (1 / expo) if expo == 3 else (v / c) ** 0.5
        out.append({"value": a, "a": a, "volume": v, "energy": float(birch_murnaghan(v, E0, V0, B0, B0P)), "stress": None})
    return out


def test_recupera_los_parametros_en_3d():
    bm = fit_birch_murnaghan(_filas(np.linspace(88, 114, 9)))
    assert bm["v0"] == pytest.approx(V0, rel=1e-6) and bm["e0"] == pytest.approx(E0, abs=1e-9)
    assert bm["b0_gpa"] == pytest.approx(B0 * EV_A3_TO_GPA, rel=1e-6) and bm["b0p"] == pytest.approx(B0P, rel=1e-5)
    assert bm["a0"] == pytest.approx(V0 ** (1 / 3), rel=1e-6) and bm["exponent"] == pytest.approx(3.0, abs=1e-6)
    assert bm["rmse_mev"] < 1e-3 and not bm["extrapolated"]


def test_monocapa_con_c_fijo_usa_exponente_2_y_no_da_b0():
    bm = fit_birch_murnaghan(_filas(np.linspace(88, 114, 9), expo=2, c=20.0), is_2d=True)
    assert bm["exponent"] == pytest.approx(2.0, abs=1e-6) and bm["b0_gpa"] is None
    assert bm["a0"] == pytest.approx((V0 / 20.0) ** 0.5, rel=1e-6)


def test_pocos_puntos_o_sin_energia_devuelve_none():
    assert fit_birch_murnaghan(_filas(np.linspace(90, 110, 4))) is None
    filas = _filas(np.linspace(88, 114, 8)); [f.update(energy=None) for f in filas]
    assert fit_birch_murnaghan(filas) is None


def test_marca_la_extrapolacion_si_el_minimo_esta_fuera_del_rango():
    bm = fit_birch_murnaghan(_filas(np.linspace(104, 130, 9)))      # solo la rama de expansión
    assert bm is None or bm["extrapolated"]


def _serie_red(tmp, a_cero_esfuerzo):
    for a in np.arange(2.40, 2.5401, 0.02):
        v = np.sqrt(3) / 2 * a * a * 20.0
        _calc(tmp / "Graphene" / "3-Lattice_Parameter" / f"{a:.3f}", a=a,
              energy=float(birch_murnaghan(v, -18.45, np.sqrt(3) / 2 * 2.468 ** 2 * 20.0, 0.35, 4.3)),
              sxx=(a_cero_esfuerzo - a) * 900)
    return tmp / "Graphene"


def test_reporte_con_a0_coherente_con_el_esfuerzo(tmp_path):
    m = analyze(load_session(scan_folder(_serie_red(tmp_path, 2.468)))).modules["estructura"]
    assert any("a₀ = 2.4680 Å" in h for h in m.highlights)
    assert not any("difieren" in w or "fuera del rango" in w for w in m.warnings)
    assert m.params["a0"] == pytest.approx(2.468, abs=1e-4)
    assert any("Repita el cálculo" in s for s in m.suggestions)


def test_reporte_avisa_si_a0_y_esfuerzo_no_coinciden(tmp_path):
    m = analyze(load_session(scan_folder(_serie_red(tmp_path, 2.50)))).modules["estructura"]
    assert any("difieren" in w for w in m.warnings)


def test_la_curva_usa_el_parametro_de_la_carpeta_aunque_el_contcar_cambie_de_forma():
    """Con ISIF=4 o c relajado, a del CONTCAR ≠ valor de la carpeta: la curva no debe desplazarse (a0 = 4.4 aquí)."""
    v0, filas = 300.0, []
    for x in np.linspace(4.2, 4.6, 9):
        v = v0 * (x / 4.4) ** 3
        filas.append({"value": x, "a": x * 1.03, "b": x * 1.03, "c": x * 0.9, "volume": v, "stress": None,
                      "energy": float(birch_murnaghan(v, -50, v0, 0.4, 4.0))})
    bm = fit_birch_murnaghan(filas)
    assert bm["a0"] == pytest.approx(4.4, abs=1e-6)
    assert bm["curve_a"][np.argmin(bm["curve_e"])] == pytest.approx(4.4, abs=0.01)
    assert bm["lattice0"]["a"] == pytest.approx(4.4 * 1.03, abs=1e-6) and bm["lattice0"]["c"] == pytest.approx(4.4 * 0.9, abs=1e-6)
