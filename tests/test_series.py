"""Series de cálculos (barridos de malla k, ENCUT, parámetro de red, vacío) con una disposición tipo Graphene."""
import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from electronanopy.project import scan_folder, load_session, last_energy, Session
from electronanopy.report import analyze
from electronanopy.structural import series_table, series_summary, read_stress, plot_series
from electronanopy.electronic.io import _efermi_from_outcar


def _calc(folder, a=2.46, c=20.0, energy=-18.45, sxx=0.0, vasprun=True, contcar=True):
    """Carpeta de un cálculo con CONTCAR, OSZICAR y OUTCAR mínimos (sin POTCAR)."""
    folder.mkdir(parents=True, exist_ok=True)
    if contcar:
        Structure(Lattice.hexagonal(a, c), ["C", "C"], [[1/3, 2/3, 0.5], [2/3, 1/3, 0.5]]).to(
            filename=str(folder / "CONTCAR"), fmt="poscar")
    if vasprun:
        (folder / "vasprun.xml").write_text("<modeling/>")
    (folder / "OSZICAR").write_text(f"   1 F= {energy - 0.0002:.8E} E0= {energy:.8E}  d E =-.3E-03\n")
    (folder / "OUTCAR").write_text(f"  in kB   {sxx:.5f}  {sxx:.5f}  -0.50000  0.00000  0.00000  0.00000\n"
                                   "  E-fermi :  -2.7087     XC(G=0):\n")


def _grafeno(tmp_path):
    """Réplica sin datos con licencia de la carpeta de pruebas: 4 series + plantillas en el contenedor."""
    r = tmp_path / "Graphene"
    energias = {2.40: -18.30, 2.44: -18.43, 2.46: -18.4495, 2.48: -18.4480, 2.52: -18.39}
    for a, e in energias.items():                                     # σ_xx cambia de signo cerca de 2.47
        _calc(r / "3-Lattice_Parameter" / f"{a:.3f}", a=a, energy=e, sxx=(2.47 - a) * 900)
    for n, e in {2: -16.3, 4: -18.48, 10: -18.4500, 8: -18.4490, 18: -18.4458}.items():
        _calc(r / "1-KPOINTS" / str(n), energy=e)
    for n, e in {400: -18.4636, 500: -18.4456, 600: -18.4500, 700: -18.4594}.items():
        _calc(r / "2-Encuts" / str(n), energy=e, c=10.0)
    for n, e in {10: -18.451567, 15: -18.451767, 20: -18.451832, 25: -18.45189}.items():
        _calc(r / "4-Vacuum" / f"{n:.2f}", energy=e, c=float(n))
    for sub in ("1-KPOINTS", "2-Encuts"):                             # plantillas sueltas en el contenedor
        (r / sub / "POSCAR").write_text("plantilla")
        (r / sub / "INCAR").write_text('incar = """\nENCUT = 500\n')
    return r


def test_las_cuatro_series_se_detectan_y_se_ordenan_numericamente(tmp_path):
    pf = scan_folder(_grafeno(tmp_path))
    assert [s.kind for s in pf.series] == ["lattice", "kpoints", "encut", "vacuum"]
    k = pf.series[1]
    assert k.values == [2, 4, 8, 10, 18]                              # 10 va después de 8 (no orden alfabético)
    assert [m.name for m in k.members] == ["2", "4", "8", "10", "18"]


def test_serie_de_red_usa_el_minimo_de_energia_y_las_demas_la_ultima(tmp_path):
    pf = scan_folder(_grafeno(tmp_path))
    assert pf.series[0].default_member.name == "2.460"
    assert [s.default_member.name for s in pf.series[1:]] == ["18", "700", "25.00"]


def test_un_solo_calculo_coherente_y_notas_resumidas(tmp_path):
    pf = scan_folder(_grafeno(tmp_path))
    assert pf.structure.parent == pf.scf_vasprun.parent == pf.outcar.parent
    assert pf.structure.parent.name == "2.460"
    assert pf.chgcar is None                                          # no se toma de otro cálculo
    assert len(pf.notes) == 5                                         # 4 series + aviso del mínimo de red
    assert not any("sin pista" in n for n in pf.notes)


def test_plantillas_del_contenedor_no_se_toman_como_estructura(tmp_path):
    pf = scan_folder(_grafeno(tmp_path))
    assert pf.structure.name == "CONTCAR"


def test_profundidad_suficiente_para_la_carpeta_madre(tmp_path):
    _grafeno(tmp_path)
    assert scan_folder(tmp_path).series                              # Graphene/<serie>/<calculo> a 3 niveles


def test_menos_de_tres_calculos_o_con_pistas_de_otro_tipo_no_son_serie(tmp_path):
    for n in ("a", "b"):
        _calc(tmp_path / "x" / n)
    for n in ("01_relax", "02_scf", "04_dos"):
        _calc(tmp_path / "y" / n)
    pf = scan_folder(tmp_path)
    assert pf.series == []


def test_lectura_de_energia_y_esfuerzo(tmp_path):
    _calc(tmp_path / "c", energy=-18.4496, sxx=9.569)
    assert last_energy(tmp_path / "c") == pytest.approx(-18.4496)
    assert read_stress(tmp_path / "c") == pytest.approx((9.569, 9.569, -0.5))
    assert last_energy(tmp_path / "no_existe") is None and read_stress(tmp_path / "no_existe") is None


def test_tabla_de_la_serie_y_tolerancia_a_cálculos_danados(tmp_path):
    r = _grafeno(tmp_path)
    (r / "3-Lattice_Parameter" / "2.400" / "CONTCAR").write_text("basura")
    ser = scan_folder(r).series[0]
    rows = series_table(ser)
    assert rows[0]["error"] and rows[0]["structure"] is None            # la fila dañada no detiene la tabla
    assert all(x["error"] is None for x in rows[1:])
    assert rows[2]["a"] == pytest.approx(2.46) and rows[2]["bond_min"] == pytest.approx(2.46 / np.sqrt(3), abs=1e-3)
    assert rows[2]["vacuum"] == pytest.approx(20.0, abs=1e-2)


def test_resumen_de_la_serie_de_red_interpola_el_esfuerzo(tmp_path):
    ser = scan_folder(_grafeno(tmp_path)).series[0]
    sm = series_summary(ser, series_table(ser))
    assert sm["a_min"] == pytest.approx(2.46) and sm["a_zero_stress"] == pytest.approx(2.47, abs=1e-6)


def test_resumen_de_convergencia(tmp_path):
    pf = scan_folder(_grafeno(tmp_path))
    enc = series_summary(pf.series[2], series_table(pf.series[2]))
    assert enc["spread_last3_mev"] == pytest.approx(13.8, abs=0.1) and not enc["monotonic"]
    vac = series_summary(pf.series[3], series_table(pf.series[3]))
    assert vac["spread_last3_mev"] < 1 and vac["monotonic"]


def test_serie_sin_energias_no_rompe(tmp_path):
    for n in (1, 2, 3):
        d = tmp_path / "k" / str(n); _calc(d); (d / "OSZICAR").unlink()
    ser = scan_folder(tmp_path).series[0]
    rows = series_table(ser)
    assert series_summary(ser, rows)["n_energy"] == 0
    fig = plot_series(ser, rows, ax=None)
    assert fig is not None


def test_reporte_marca_lo_no_convergido_y_la_no_monotonia(tmp_path):
    s = load_session(scan_folder(_grafeno(tmp_path)))
    assert len(s.series) == 4 and "series" not in s.errors
    m = analyze(s).modules["estructura"]
    assert any("2-Encuts" in w and "sin converger" in w for w in m.warnings)
    assert any("1-KPOINTS" in w and "monótona" in w for w in m.warnings)
    assert not any("4-Vacuum" in w for w in m.warnings)
    assert any("3-Lattice_Parameter" in h and "2.470" in h for h in m.highlights)


def test_efermi_desde_el_outcar(tmp_path):
    _calc(tmp_path / "c")
    assert _efermi_from_outcar(tmp_path / "c" / "OUTCAR") == pytest.approx(-2.7087)
    assert _efermi_from_outcar(tmp_path / "no_existe") is None
