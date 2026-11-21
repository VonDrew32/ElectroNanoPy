import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")
from electronanopy.electronic import (band_gap, dirac_points, fermi_velocity, dos_at_fermi, classify,
    orbital_weights, read_magnetization, plot_bands, plot_dos, plot_bands_dos, load_vasprun, ElectronicReadError)
from electronanopy.electronic import _synthetic as S

HBAR_VF = 3 ** 0.5 / 2 * S.A * S.T      # ħ v_F analítico del grafeno tight-binding (eV·Å)


def test_grafeno_brecha_nula_y_punto_de_dirac_en_K():
    bs = S.graphene_bands()
    g = band_gap(bs)
    assert g.is_metal and g.gap == 0.0 and "nula" in g.message
    dp = dirac_points(bs)
    assert len(dp) == 1 and dp[0]["label"] == "K"
    assert dp[0]["frac_coords"][:2] == pytest.approx((1 / 3, 1 / 3), abs=1e-6)
    assert dp[0]["energy"] == pytest.approx(0.0, abs=1e-6)


def test_velocidad_de_fermi_coincide_con_el_valor_analitico():
    bs = S.graphene_bands(npts=80)
    v = fermi_velocity(bs, dirac_points(bs)[0])
    assert v["hbar_vF_eV_A"] == pytest.approx(HBAR_VF, rel=0.02)
    assert v["vF_m_s"] == pytest.approx(HBAR_VF * 1e-10 / 6.582119569e-16, rel=0.02)


def test_velocidad_con_pocos_puntos_falla_con_mensaje():
    bs = S.graphene_bands(npts=60)
    with pytest.raises(ValueError, match="Muy pocos"):
        fermi_velocity(bs, dirac_points(bs)[0], window=0.001)


def test_dos_en_fermi_y_clasificacion():
    bs, dos = S.graphene_bands(), S.graphene_dos()
    area = dos.densities[next(iter(dos.densities))].sum() * (dos.energies[1] - dos.energies[0])
    assert area == pytest.approx(2.0, rel=0.02)                  # 2 estados π por celda
    assert dos_at_fermi(dos) < 0.1                               # DOS -> 0 en el punto de Dirac
    assert classify(bs, dos) == "semimetal de Dirac"


def test_brecha_abierta_es_semiconductor_directo_en_K():
    bs = S.graphene_bands(gap=0.5)
    g = band_gap(bs)
    assert not g.is_metal and g.gap == pytest.approx(0.5, abs=1e-3) and g.direct
    assert dirac_points(bs) == []
    assert classify(bs) == "semiconductor/aislante"


def test_pesos_orbitales():
    bs = S.graphene_bands()
    w = orbital_weights(bs, ("pz",), element="C")
    assert w.shape == (2, len(bs.kpoints)) and np.allclose(w, 1.0)
    assert np.allclose(orbital_weights(bs, ("s", "px")), 0.0)
    with pytest.raises(ValueError, match="Orbital desconocido"):
        orbital_weights(bs, ("zz",))
    with pytest.raises(ValueError, match="No hay átomos"):
        orbital_weights(bs, ("pz",), element="Fe")
    with pytest.raises(ValueError, match="proyecciones"):
        orbital_weights(S.graphene_bands(with_projections=False))


def test_magnetizacion_desde_outcar(tmp_path):
    f = S.write_outcar_magnetization(tmp_path / "OUTCAR", [3.5, -3.5])
    m = read_magnetization(f, S.graphene_structure())
    assert m["spin_polarized"] and m["total"] == pytest.approx(0.0, abs=1e-9)
    assert [a["tot"] for a in m["per_atom"]] == [3.5, -3.5]
    assert m["per_atom"][0]["element"] == "C" and m["per_atom"][0]["d"] == 3.5


def test_magnetizacion_usa_el_ultimo_bloque(tmp_path):
    a = S.write_outcar_magnetization(tmp_path / "a", [1.0]).read_text()
    b = S.write_outcar_magnetization(tmp_path / "b", [2.0]).read_text()
    (tmp_path / "OUTCAR").write_text(a + b)
    assert read_magnetization(tmp_path / "OUTCAR")["total"] == pytest.approx(2.0)


def test_sin_polarizacion_de_espin(tmp_path):
    (tmp_path / "OUTCAR").write_text(" vasp.6.x\n sin bloque de magnetizacion\n")
    m = read_magnetization(tmp_path / "OUTCAR")
    assert m["spin_polarized"] is False and m["total"] is None and m["per_atom"] is None
    with pytest.raises(FileNotFoundError):
        read_magnetization(tmp_path / "no_existe")


def test_lectura_invalida(tmp_path):
    with pytest.raises(ElectronicReadError):
        load_vasprun(tmp_path / "no_existe.xml")
    malo = tmp_path / "vasprun.xml"; malo.write_text("<modeling><esto no es un vasprun")
    with pytest.raises(ElectronicReadError):
        load_vasprun(malo)


def test_graficas():
    bs, dos = S.graphene_bands(), S.graphene_dos()
    fig = plot_bands_dos(bs, dos, title="grafeno")
    assert len(fig.axes) == 2
    ax = plot_bands(bs, weights=orbital_weights(bs, ("pz",)))
    labels = [t.get_text() for t in ax.get_xticklabels()]
    assert labels == ["Γ", "M", "K", "Γ"]
    assert plot_dos(dos).get_xlabel().startswith("$E")


# --- casos hallados con un vasprun.xml real de grafeno (E_F coincide con la energía del punto de Dirac) ---
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import Dos
from electronanopy.electronic import dos_energy_step


def _grafeno_con_fermi(delta):
    """Bandas del grafeno cuya E_F queda ``delta`` eV por debajo del cruce de Dirac (en el real fue 4·10⁻⁵ eV)."""
    bs = S.graphene_bands()
    e = {sp: np.asarray(v) + delta for sp, v in bs.bands.items()}        # bajar E_F = subir las bandas
    return BandStructureSymmLine([k.frac_coords for k in bs.kpoints], e, bs.structure.lattice.reciprocal_lattice, 0.0,
                                 {l: k.frac_coords for l, k in bs.labels_dict.items()}, structure=bs.structure)


@pytest.mark.parametrize("delta", [-1e-3, -4e-5, 0.0, 4e-5, 1e-3])
def test_brecha_del_grafeno_no_depende_de_donde_caiga_E_F(delta):
    bs = _grafeno_con_fermi(delta)
    g = band_gap(bs, nelect=2.0)                                       # 2 bandas π: 1 ocupada (en el cálculo real, 8 e⁻ y 8 bandas)
    assert g.gap == 0.0 and g.is_metal and "nula" in g.message
    assert classify(bs, nelect=2.0) == "semimetal de Dirac"


def test_brecha_por_ocupacion_en_un_semiconductor_y_numero_impar_de_electrones():
    bs = S.graphene_bands(gap=1.0)                                      # 2 bandas π: 1 ocupada → NELECT = 2
    g = band_gap(bs, nelect=2.0)
    assert g.gap == pytest.approx(1.0, abs=0.05) and g.direct and not g.is_metal
    assert band_gap(bs, nelect=3.0).is_metal                            # número impar de electrones → metal


def test_dos_con_rejilla_gruesa_no_se_usa_para_clasificar():
    bs = S.graphene_bands()
    grid = np.linspace(-20, 15, 301)                                    # NEDOS = 301 → paso 0.117 eV
    gruesa = Dos(0.0, grid, {Spin.up: np.full(301, 0.5)})               # «DOS alta» pero con rejilla gruesa
    assert dos_energy_step(gruesa) == pytest.approx(35 / 300)
    assert classify(bs, gruesa) == "semimetal de Dirac"                 # se ignora: no es fiable
    fina = Dos(0.0, np.linspace(-20, 15, 3501), {Spin.up: np.full(3501, 0.5)})
    assert classify(bs, fina) == "metal"                                # con rejilla fina sí cuenta


def test_pesos_orbitales_relativos_y_crudos():
    bs = S.graphene_bands(with_projections=True)
    k = bs.projections[Spin.up]
    k[:] = k * 0.525                                                    # VASP solo proyecta dentro de las esferas PAW
    assert orbital_weights(bs, ("pz",)).max() == pytest.approx(1.0)
    assert orbital_weights(bs, ("pz",), normalize=False).max() == pytest.approx(0.525)
    assert orbital_weights(bs, ("s",)).max() == 0.0
