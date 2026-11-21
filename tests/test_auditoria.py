"""Pruebas de los errores hallados en la revisión independiente (v0.11)."""
import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin
from electronanopy.electronic import band_gap, classify, orbital_weights, read_magnetization
from electronanopy.project import scan_folder, load_session, _series_kind, _to_float, _name_value, Session
from electronanopy.report import analyze
from electronanopy.structural import fit_birch_murnaghan, birch_murnaghan
from test_series import _calc
from test_eos import _filas, E0, V0, B0, B0P


def _bandas(up, down=None, efermi=0.0, proj=None):
    """Bandas con k a lo largo de x: ``up``/``down`` son arreglos (nbandas, nk)."""
    s = Structure(Lattice.cubic(4.0), ["Si"], [[0, 0, 0]])
    nk = np.asarray(up).shape[1]
    k = np.array([[x, 0, 0] for x in np.linspace(0, .5, nk)])
    bands = {Spin.up: np.asarray(up, float)}
    if down is not None:
        bands[Spin.down] = np.asarray(down, float)
    return BandStructureSymmLine(k, bands, s.lattice.reciprocal_lattice, efermi,
                                 {"\\Gamma": k[0], "X": k[-1]}, structure=s, projections=proj)


# Semiconductor de 1 eV: 4 bandas de valencia (≤ 0) y 4 de conducción (≥ 1)
VAL = np.array([[-1.0 - i * 0.5 - 0.3 * np.cos(t) for t in np.linspace(0, np.pi, 5)] for i in range(4)])
VAL = VAL - VAL.max()                                         # VBM = 0
CON = np.array([[1.0 + i * 0.5 + 0.3 * (1 - np.cos(t)) for t in np.linspace(0, np.pi, 5)] for i in range(4)])
SEMI = np.vstack([VAL, CON])


@pytest.mark.parametrize("efermi", [0.0, -0.0005, -0.05, 0.4])
def test_brecha_con_espin_no_depende_de_donde_caiga_ef(efermi):
    bs = _bandas(SEMI, SEMI, efermi=efermi)                   # ISPIN = 2 sin momento neto
    g = band_gap(bs, nelect=8, magnetization=0.0)
    assert not g.is_metal and g.gap == pytest.approx(1.0, abs=1e-6)


def test_brecha_con_espin_y_momento_distintos_por_canal():
    up = SEMI                                                  # 4 bandas ocupadas
    down = np.vstack([VAL[:3], CON])                           # 3 ocupadas (N = 7, M = 1 µB)
    g = band_gap(_bandas(up, down), nelect=7, magnetization=1.0)
    assert not g.is_metal and g.gap == pytest.approx(1.0, abs=1e-6)


def test_momento_no_entero_es_metal_y_sin_momento_se_usa_ef():
    bs = _bandas(SEMI, SEMI)
    assert band_gap(bs, nelect=8, magnetization=0.37).is_metal
    assert band_gap(bs, nelect=8).is_metal in (True, False)      # sin momento: recurre a E_F, no falla


def test_brecha_con_espin_orbita_cuenta_nelect_bandas():
    # con SOC hay un solo canal y cada banda ya contiene los dos espines: 8 e⁻ = 8 bandas ocupadas
    val = np.vstack([VAL, VAL - 5.0]); con = np.vstack([CON, CON + 5.0])
    bs = _bandas(np.vstack([val, con]))
    g = band_gap(bs, nelect=8, soc=True)
    assert not g.is_metal and g.gap == pytest.approx(1.0, abs=1e-6)
    assert band_gap(bs, nelect=8).gap != pytest.approx(1.0, abs=1e-6)        # sin soc cuenta mal: 4 bandas


def test_clasifica_semiconductor_con_espin():
    assert "semiconductor" in classify(_bandas(SEMI, SEMI), nelect=8, magnetization=0.0).lower()


def test_magnetizacion_sin_bloque_por_atomo(tmp_path):
    f = tmp_path / "OUTCAR"
    f.write_text("   ISPIN  =      2    spin polarized calculation?\n   LSORBIT =      F\n"
                 " number of electron      16.0000000 magnetization       3.2000000\n"
                 " number of electron      16.0000000 magnetization       3.9000000\n")
    m = read_magnetization(f)
    assert m["spin_polarized"] and m["total"] == pytest.approx(3.9) and m["per_atom"] is None and m["ispin"] == 2
    f.write_text("   ISPIN  =      1\n number of electron       8.0000001 magnetization \n")
    assert read_magnetization(f)["spin_polarized"] is False


def test_lorbit_10_avisa_en_vez_de_leer_otra_columna():
    proj = np.zeros((8, 5, 3, 1)); proj[:, :, 1, :] = 1.0          # s, p, d: todo el peso en p
    bs = _bandas(SEMI, proj={Spin.up: proj})
    assert orbital_weights(bs, ("p",), normalize=False).max() == pytest.approx(1.0)
    with pytest.raises(ValueError, match="LORBIT"):
        orbital_weights(bs, ("pz",))


def test_pistas_solo_dentro_del_proyecto(tmp_path):
    r = tmp_path / "bandgap_study" / "graphene"
    for sub in ("1-relax", "2-scf", "3-dos"):
        (r / sub).mkdir(parents=True); (r / sub / "vasprun.xml").write_text("<modeling/>")
    pf = scan_folder(r)
    assert pf.scf_vasprun.parent.name == "2-scf" and pf.dos_vasprun.parent.name == "3-dos"
    assert pf.bands_vasprun is None


def test_resultado_determinista_y_natural(tmp_path):
    for n in ("12x12x1", "6x6x1", "18x18x1", "3x3x1"):
        _calc(tmp_path / "kpoints" / n, energy=-18.0 - 1 / int(n.split("x")[0]))
    s = scan_folder(tmp_path).series[0]
    assert [m.name for m in s.members] == ["3x3x1", "6x6x1", "12x12x1", "18x18x1"]
    assert s.values == [3.0, 6.0, 12.0, 18.0] and s.default_member.name == "18x18x1"


@pytest.mark.parametrize("name, kind", [("3-Lattice_Parameter", "lattice"), ("Reduced_cell_tests", "other"),
                                        ("supercell", "other"), ("Parameters", "other"), ("ordered_runs", "other"),
                                        ("2-Encuts", "encut"), ("1-KPOINTS", "kpoints"), ("4-Vacuum", "vacuum"),
                                        ("k-points", "kpoints"), ("a0_scan", "lattice")])
def test_tipo_de_serie_por_palabras(name, kind):
    assert _series_kind(name) == kind


def test_valores_numericos_estrictos():
    assert _to_float("2,46") == 2.46 and _to_float("1e-3") == 1e-3
    assert _to_float("inf") is None and _to_float("nan") is None and _to_float("1_000") is None
    assert _name_value("6x6x1") == 6.0 and _name_value("abc") is None


def test_b0_solo_si_la_celda_escala_en_las_tres_direcciones():
    # 3D con solo a variando y c fijo: exponente 2, dV/da ≠ 3V/a -> sin B0
    bm = fit_birch_murnaghan(_filas(np.linspace(88, 114, 9), expo=2, c=20.0), is_2d=False)
    assert bm["exponent"] == pytest.approx(2.0, abs=1e-6) and bm["b0_gpa"] is None
    assert fit_birch_murnaghan(_filas(np.linspace(88, 114, 9)))["b0_gpa"] is not None


def test_corte_de_densidad_con_lamina_en_el_borde_de_la_celda():
    from electronanopy.electronic import _synthetic as SE
    from electronanopy.report import _density
    s = Session(chg=SE.graphene_chgcar())
    from electronanopy.electronic import density
    from electronanopy.report import Thresholds
    r = _density(s, Thresholds())
    assert r.params["slice_frac"] == pytest.approx(0.5, abs=0.01)
    assert r.params["iso_level"] == pytest.approx(0.1 * density(s.chg).max())
    # lámina que cruza el borde de la celda (z = 0.02 y 0.98): el plano de los átomos es z ≈ 0, no 0.5
    st = s.chg.structure.copy()
    st.translate_sites(range(len(st)), [0, 0, 0.5 - st.frac_coords[:, 2].mean() - 0.5 + 0.0], frac_coords=True)
    s.chg.structure = Structure(st.lattice, ["C", "C"], [[1/3, 2/3, 0.02], [2/3, 1/3, 0.98]])
    assert _density(s, Thresholds()).params["slice_frac"] in (pytest.approx(0.0, abs=0.01), pytest.approx(1.0, abs=0.01))
