import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")
from phonopy.physical_units import get_physical_units
VaspToTHz = get_physical_units().DefaultToTHz
from electronanopy.phonons import (load_phonopy, PhononReadError, compute_bands, compute_dos,
                                   check_stability, plot_bands, plot_bands_dos)
from electronanopy.phonons import _synthetic as S
from electronanopy.phonons.paths import path_2d, in_plane_class


def test_frecuencias_analiticas_en_X():
    """En la red cúbica con resortes k1 (1.º vecinos) y k2 (2.º): D_xx = 4k1+8k2 y D_yy = D_zz = 4k2 en X."""
    k1, k2 = 1.0, 0.5
    ph = S.cubic(k1=k1, k2=k2)
    m = ph.primitive.masses[0]
    ph.run_qpoints([[0.5, 0, 0]])
    f = np.sort(ph.qpoints.frequencies[0])
    esperado = np.sort(VaspToTHz * np.sqrt(np.array([4 * k2, 4 * k2, 4 * k1 + 8 * k2]) / m))
    assert f == pytest.approx(esperado, abs=1e-3)


def test_bandas_3d_estable():
    ph = S.cubic()
    b = compute_bands(ph)
    assert not b.two_d and b.tick_labels[0] == "Γ" and "X" in b.tick_labels
    assert b.frequencies.shape[1] == 3
    assert b.frequencies[0] == pytest.approx(0, abs=1e-3)     # acústicos nulos en Γ
    r = check_stability(ph, bands=b)
    assert r.stable and r.imaginary == []


def test_dos_normalizada_a_3N():
    for ph in (S.cubic(), S.triangular_2d()):
        d = compute_dos(ph)
        area = d.dos.sum() * (d.frequencies[1] - d.frequencies[0])
        assert area == pytest.approx(3.0, rel=0.02)            # 3 ramas por átomo


def test_material_2d_camino_gamma_m_k_gamma():
    ph = S.triangular_2d()
    b = compute_bands(ph)
    assert b.two_d and b.tick_labels == ["Γ", "M", "K", "Γ"]
    seg = path_2d(np.asarray(ph.primitive.cell))[0]
    rec = 2 * np.pi * np.linalg.inv(np.asarray(ph.primitive.cell)).T
    K = seg[1][3]
    assert np.linalg.norm(K @ rec) == pytest.approx(4 * np.pi / (3 * 2.5), rel=1e-3)   # |K| = 4π/(3a)
    assert check_stability(ph, bands=b).stable


def test_clasificacion_de_redes_2d():
    assert in_plane_class(np.array([[2.5, 0, 0], [-1.25, 2.165, 0], [0, 0, 15]])) == "hex"
    assert in_plane_class(np.diag([3.0, 3.0, 15])) == "square"
    assert in_plane_class(np.diag([3.0, 4.0, 15])) == "rect"


def test_inestabilidad_detectada():
    r = check_stability(S.cubic(k2=-0.2))
    assert not r.stable and r.imaginary and r.min_frequency < -1.0
    assert "inestabilidad" in r.message.lower()


def test_modos_acusticos_cerca_de_gamma_se_ignoran():
    r = check_stability(S.cubic(onsite=0.02))                  # f(Γ) ≈ -0.64 THz solo en el entorno de Γ
    assert r.stable and r.n_ignored_near_gamma > 0 and r.min_frequency < -0.05
    r2 = check_stability(S.cubic(onsite=0.02), gamma_radius=0.0)   # sin exclusión se marca inestable
    assert not r2.stable


@pytest.mark.parametrize("layout", ["params", "disp+forcesets"])
def test_lectura_desde_disco(tmp_path, layout):
    S.save_as_files(S.cubic(), tmp_path, layout)
    leido = compute_bands(load_phonopy(tmp_path))
    original = compute_bands(S.cubic())
    assert leido.frequencies == pytest.approx(original.frequencies, abs=1e-3)


def test_errores_de_lectura(tmp_path):
    with pytest.raises(PhononReadError):
        load_phonopy(tmp_path / "no_existe")
    with pytest.raises(PhononReadError):
        load_phonopy(tmp_path)                                  # carpeta vacía
    S.save_as_files(S.cubic(), tmp_path, "disp+forcesets")
    (tmp_path / "FORCE_SETS").unlink()
    with pytest.raises(PhononReadError, match="FORCE_SETS"):
        load_phonopy(tmp_path)


def test_graficas_y_unidades():
    ph = S.triangular_2d()
    b, d = compute_bands(ph), compute_dos(ph)
    fig = plot_bands_dos(b, d, check_stability(ph))
    assert len(fig.axes) == 2
    ax = plot_bands(b, unit="cm-1")
    assert ax.get_ylabel().startswith("Frecuencia")
    with pytest.raises(ValueError):
        plot_bands(b, unit="meV")


def test_nacl_real_de_phonopy_con_vasp():
    """Regresión con un cálculo real (VASP + Phonopy, NaCl): se lee con NAC y da fonones estables."""
    from pathlib import Path
    ph = load_phonopy(Path(__file__).parent / "data" / "nacl")
    assert ph.nac_params is not None                           # se aplicó la corrección no analítica (BORN)
    b = compute_bands(ph)
    assert b.frequencies.shape[1] == 6 and not b.two_d
    assert b.frequencies.max() == pytest.approx(7.4, abs=0.2)
    near_gamma = np.sort(b.frequencies[1])
    assert near_gamma[-1] - near_gamma[-2] > 2.0               # separación LO-TO (efecto de la NAC)
    assert check_stability(ph, bands=b).stable
    d = compute_dos(ph)
    assert d.dos.sum() * (d.frequencies[1] - d.frequencies[0]) == pytest.approx(6.0, rel=0.02)


# --- diagnóstico de calidad (regla de suma acústica y tamaño de supercelda) ---
from electronanopy.phonons import phonon_quality, apply_asr


def test_calidad_sin_ruido_cumple_la_regla_de_suma_acustica():
    q = phonon_quality(S.cubic())
    assert q.asr_max < 1e-3 and q.periodic_dims == 3
    assert q.min_cell_length == pytest.approx(9.0)                  # 3 × 3 Å
    q2 = phonon_quality(S.triangular_2d())
    assert q2.periodic_dims == 2 and q2.min_cell_length == pytest.approx(10.0)   # c (15 Å) no cuenta


def test_ruido_en_las_fuerzas_viola_la_asr_y_apply_asr_la_corrige_sin_tocar_el_original():
    ph = S.cubic(onsite=0.02)
    assert phonon_quality(ph).asr_max > 0.5                        # ≈ 0.64 THz en Γ
    fixed = apply_asr(ph)
    assert phonon_quality(fixed).asr_max < 0.05
    assert fixed is not ph and phonon_quality(ph).asr_max > 0.5    # el original sigue igual
