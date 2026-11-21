import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")
from electronanopy.electronic import (load_chgcar, DensityReadError, density, total_electrons, net_magnetization,
    planar_average, density_at, slice_2d, isosurface, plot_density_slice, plot_planar_average, plot_isosurface)
from electronanopy.electronic import _synthetic as S


@pytest.fixture
def chg_file(tmp_path):
    f = tmp_path / "CHGCAR"
    S.graphene_chgcar().write_file(str(f))
    return f


def test_lectura_y_electrones_totales(chg_file):
    chg = load_chgcar(chg_file)
    assert total_electrons(chg) == pytest.approx(8.0, abs=1e-6)        # 2 C x 4 e de valencia
    assert chg.structure.composition.reduced_formula == "C"
    assert density(chg).shape == (24, 24, 100)


def test_carpeta_con_chgcar(chg_file):
    assert total_electrons(load_chgcar(chg_file.parent)) == pytest.approx(8.0, abs=1e-6)


def test_promedio_planar(chg_file):
    chg = load_chgcar(chg_file)
    z, prof = planar_average(chg)
    assert z[np.argmax(prof)] == pytest.approx(7.5, abs=0.2)           # plano del grafeno a c/2
    assert prof[0] < 1e-6                                              # vacío
    area = np.linalg.norm(np.cross(*chg.structure.lattice.matrix[:2]))
    assert prof.sum() * (z[1] - z[0]) * area == pytest.approx(8.0, abs=1e-6)


def test_densidad_en_atomo_mayor_que_en_el_hueco(chg_file):
    chg = load_chgcar(chg_file)
    en_atomo, en_hueco = density_at(chg, [[1 / 3, 2 / 3, 0.5], [0.0, 0.0, 0.5]])
    assert en_atomo > 2 * en_hueco > 0
    # en un nodo de la malla coincide con el valor de la malla
    assert density_at(chg, [[0, 0, 0.5]])[0] == pytest.approx(density(chg)[0, 0, 50])


def test_corte_2d(chg_file):
    chg = load_chgcar(chg_file)
    plano, vacio = slice_2d(chg, 2, 0.5), slice_2d(chg, 2, 0.0)
    assert plano.shape == (24, 24) and plano.max() > 1e3 * vacio.max()
    medio = slice_2d(chg, 2, 0.5 + 0.5 / 100)                          # entre dos planos de la malla
    assert min(slice_2d(chg, 2, 0.5).max(), slice_2d(chg, 2, 0.51).max()) <= medio.max() <= max(
        slice_2d(chg, 2, 0.5).max(), slice_2d(chg, 2, 0.51).max())


def test_isosuperficie_envuelve_el_plano(chg_file):
    pytest.importorskip("skimage")
    chg = load_chgcar(chg_file)
    v, f = isosurface(chg, 0.5 * density(chg).max())
    assert v.shape[1] == 3 and f.shape[1] == 3 and len(v) > 100
    assert abs(v[:, 2].mean() - 7.5) < 0.1 and v[:, 2].min() > 5.5 and v[:, 2].max() < 9.5
    with pytest.raises(ValueError, match="fuera del rango"):
        isosurface(chg, 1e6)


def test_magnetizacion(tmp_path):
    assert net_magnetization(S.graphene_chgcar()) is None
    af = S.graphene_chgcar(spin_moments=(3.0, -3.0))
    assert net_magnetization(af) == pytest.approx(0.0, abs=1e-9)
    assert net_magnetization(S.graphene_chgcar(spin_moments=(2.0, 1.0))) == pytest.approx(3.0, abs=1e-9)
    assert density(af, "diff").max() > 0 > density(af, "diff").min()
    with pytest.raises(ValueError, match="solo existe en cálculos polarizados"):
        density(S.graphene_chgcar(), "diff")


def test_errores_de_lectura(tmp_path):
    with pytest.raises(DensityReadError):
        load_chgcar(tmp_path / "no_existe")
    vacio = tmp_path / "CHGCAR"; vacio.write_text("")
    with pytest.raises(DensityReadError, match="vacío"):
        load_chgcar(vacio)
    malo = tmp_path / "malo"; malo.write_text("esto no es un CHGCAR\n1.0\n")
    with pytest.raises(DensityReadError):
        load_chgcar(malo)


def test_graficas(chg_file):
    pytest.importorskip("skimage")
    chg = load_chgcar(chg_file)
    ax = plot_density_slice(chg, 2, 0.5)
    assert ax.get_xlabel() == "x (Å)" and len(ax.lines) >= 8          # átomos marcados en las celdas repetidas
    assert plot_planar_average(chg).get_ylabel().startswith(r"$\bar")
    assert len(plot_isosurface(chg, 0.5 * density(chg).max()).collections) >= 1


def test_isosuperficie_limita_el_numero_de_caras():
    from electronanopy.electronic import _synthetic as SE
    from electronanopy.electronic.density import isosurface, density
    chg = SE.graphene_chgcar(ngrid=(48, 48, 100))
    lvl = 0.05 * float(density(chg).max())
    _, caras_todas = isosurface(chg, lvl, max_faces=10 ** 9)
    verts, caras = isosurface(chg, lvl, max_faces=len(caras_todas) // 2)
    assert len(caras) <= len(caras_todas) // 2 < len(caras_todas)
    assert np.isfinite(verts).all()
