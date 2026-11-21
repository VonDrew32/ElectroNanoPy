"""Regresiones halladas al probar con salidas reales y públicas de VASP/Phonopy (v0.13).

Los casos reales (pymatgen-test-files, repositorio de Phonopy) no se distribuyen con el proyecto: aquí se reproducen
con texto y estructuras mínimas que contienen el mismo rasgo.
"""
from pathlib import Path
import shutil
import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from electronanopy.electronic import read_magnetization
from electronanopy.electronic import _synthetic as SE
from electronanopy.phonons import _synthetic as SP
from electronanopy.project import Session, scan_folder, load_session
from electronanopy.report import analyze
from electronanopy.structural import coordination_numbers, bond_table
from test_auditoria import _bandas, SEMI
from test_project import _estructura_de_la_celda_de_fonones

DATA = Path(__file__).parent / "data"


def _outcar(tmp_path, cuerpo, cabecera="   ISPIN  =      2    spin polarized calculation?\n"):
    p = tmp_path / "OUTCAR"
    p.write_text(cabecera + cuerpo)
    return p


# --- OUTCAR: un solo átomo (VASP no escribe la línea de suma) ---
def test_magnetizacion_un_solo_atomo_sin_linea_de_suma(tmp_path):
    cuerpo = (" magnetization (x)\n\n# of ion     s       p       d       tot\n"
              "----------------------------------------\n    1       0.000   0.000   0.000   0.000\n"
              "\n BZINTS: Fermi energy:  5.0\n")
    m = read_magnetization(_outcar(tmp_path, cuerpo))
    assert m["total"] == pytest.approx(0.0) and len(m["per_atom"]) == 1


# --- OUTCAR no colineal (LSORBIT): momentos en (x), (y), (z); el AFM tiene total 0 pero momentos locales ±1.67 µB ---
def _bloque(eje, v):
    filas = "".join(f"  {i + 1}   0.000   0.000   {x:.3f}   {x:.3f}\n" for i, x in enumerate(v))
    return (f" magnetization ({eje})\n \n# of ion     s       p       d       tot\n"
            "----------------------------------------\n" + filas +
            "------------------------------------------------\n" f"tot   0.000   0.000   {sum(v):.3f}   {sum(v):.3f}\n \n\n")


def test_magnetizacion_no_colineal_lee_xyz(tmp_path):
    cuerpo = _bloque("x", [0, 0]) + _bloque("y", [0, 0]) + _bloque("z", [1.67, -1.67])
    m = read_magnetization(_outcar(tmp_path, cuerpo, "   ISPIN  =      1\n   LNONCOLLINEAR =      T\n   LSORBIT =      T\n"))
    assert m["noncollinear"] and m["lsorbit"]
    assert m["total"] == pytest.approx(0.0, abs=1e-9)
    assert [r["tot"] for r in m["per_atom"]] == pytest.approx([1.67, 1.67])         # módulo
    assert [r["tot_signed"] for r in m["per_atom"]] == pytest.approx([1.67, -1.67])


def test_informe_antiferromagnetico_no_colineal(tmp_path):
    cuerpo = _bloque("x", [0, 0]) + _bloque("y", [0, 0]) + _bloque("z", [1.67, -1.67])
    m = read_magnetization(_outcar(tmp_path, cuerpo, "   ISPIN  =      1\n   LNONCOLLINEAR =      T\n   LSORBIT =      T\n"))
    r = analyze(Session(structure=_estructura_de_la_celda_de_fonones(), magnetization=m)).modules["electronico"]
    txt = " ".join(r.highlights)
    assert r.available and "antiferromagn" in txt and "+1.67" in txt and "-1.67" in txt


def test_momento_cero_no_se_muestra_como_negativo():
    m = {"spin_polarized": True, "total": -0.0001, "per_atom": None, "ispin": 2, "lsorbit": False}
    r = analyze(Session(structure=_estructura_de_la_celda_de_fonones(), magnetization=m)).modules["electronico"]
    assert any("0.000 µB" in h and "-0.000" not in h for h in r.highlights)


# --- LiF: la distancia Li–Li no es un enlace ---
def test_halogenuro_alcalino_no_cuenta_cation_cation():
    lif = Structure.from_spacegroup("Fm-3m", Lattice.cubic(4.03), ["Li", "F"], [[0, 0, 0], [.5, .5, .5]])
    assert set(coordination_numbers(lif)) == {6}
    assert {(b["element"], b["neighbor"]) for b in bond_table(lif)} == {("Li", "F"), ("F", "Li")}


# --- Brecha en la malla k cuando no hay bandas en línea ---
def test_brecha_aproximada_en_malla_del_scf():
    s = Session(structure=_estructura_de_la_celda_de_fonones(), mesh_bands=_bandas(SEMI, SEMI), nelect=8.0,
                magnetization={"spin_polarized": True, "total": 0.0, "per_atom": None, "ispin": 2, "lsorbit": False})
    r = analyze(s).modules["electronico"]
    assert r.available and r.summary["gap"] == pytest.approx(1.0, abs=1e-6)
    assert r.summary["gap_source"] == "malla" and "aproximada" in " ".join(r.highlights)
    assert not s.available()["electronico"]                      # no hay nada que dibujar en la pestaña


def test_con_bandas_en_linea_no_se_usa_la_malla():
    r = analyze(Session(structure=_estructura_de_la_celda_de_fonones(), bands=SE.graphene_bands(),
                        mesh_bands=_bandas(SEMI), scf_efermi=0.0, nelect=2.0)).modules["electronico"]
    assert "gap_source" not in r.summary


def test_metal_con_nelect_impar_explica_el_motivo():
    from electronanopy.electronic import band_gap
    bs = _bandas(np.array([[-1.0, -0.5, -0.2], [0.2, 0.6, 1.0]]), np.array([[-1.0, -0.5, -0.2], [0.2, 0.6, 1.0]]))
    g = band_gap(bs, nelect=3.0, magnetization=0.0)
    assert g.is_metal and "NELECT impar" in g.message


# --- Fonones: seekpath puede cambiar la celda; la estabilidad no debe perderse ---
def test_fonones_sin_camino_evalua_en_la_malla(monkeypatch):
    import electronanopy.phonons as P

    def _falla(*a, **k):
        raise ValueError("La celda primitiva de seekpath difiere de la del cálculo de fonones")
    monkeypatch.setattr(P, "compute_bands", _falla)
    r = analyze(Session(phonon=SP.triangular_2d())).modules["fonones"]
    assert r.summary["stable"] is True and "f_max" in r.summary
    assert any("solo en la malla q" in w for w in r.warnings)


# --- Sin CONTCAR/POSCAR la estructura sale de otro módulo y se avisa ---
def test_estructura_desde_la_celda_de_fonones(tmp_path):
    for f in ("phonopy_disp.yaml", "FORCE_SETS", "BORN"):
        shutil.copy(DATA / "nacl" / f, tmp_path / f)
    s = load_session(scan_folder(tmp_path))
    assert s.structure is not None and "fonones" in s.structure_source
    m = analyze(s).modules["estructura"]
    assert any("se tomó de la celda del cálculo de fonones" in w for w in m.warnings)


# --- Carpeta de Phonopy con el CONTCAR de la supercelda desplazada (así sale de un cálculo real de grafeno 5×5×1) ---
def _proyecto_con_supercelda(tmp_path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    """05_fonones/: phonopy_disp.yaml + FORCE_SETS de NaCl y un CONTCAR que es la supercelda 2×1×1 de esa celda."""
    from electronanopy.phonons import load_phonopy
    carpeta = tmp_path / "05_fonones"
    carpeta.mkdir()
    for f in ("phonopy_disp.yaml", "FORCE_SETS", "BORN"):
        shutil.copy(DATA / "nacl" / f, carpeta / f)
    ph = load_phonopy(carpeta)
    u = ph.unitcell
    cell = Structure(u.cell, list(u.symbols), u.scaled_positions)
    cell.make_supercell([2, 1, 1])
    cell.to(filename=str(carpeta / "CONTCAR"), fmt="poscar")
    return carpeta, len(u.symbols)


def test_contcar_de_la_supercelda_de_fonones_no_es_la_estructura(tmp_path):
    _, n = _proyecto_con_supercelda(tmp_path)
    s = load_session(scan_folder(tmp_path))
    assert len(s.structure) == n and "supercelda" in s.discarded_structure and "fonones" in s.structure_source
    r = analyze(s).modules["consistencia"]
    assert not any("no coincide con la estructura" in w for w in r.warnings)


def test_estructura_propia_gana_a_la_supercelda_de_fonones(tmp_path):
    carpeta, n = _proyecto_con_supercelda(tmp_path)
    from electronanopy.phonons import load_phonopy
    u = load_phonopy(carpeta).unitcell
    relax = tmp_path / "01_relax"
    relax.mkdir()
    Structure(u.cell, list(u.symbols), u.scaled_positions).to(filename=str(relax / "CONTCAR"), fmt="poscar")
    s = load_session(scan_folder(tmp_path))
    assert s.files.structure == relax / "CONTCAR" and s.structure_source is None and s.discarded_structure is None


def test_carpeta_de_fonones_no_oculta_la_serie_de_cálculos(tmp_path):
    """Con una serie de red y una carpeta de Phonopy, el SCF sale de la serie y la estructura de los fonones."""
    from test_project import _touch
    for a in ("2.40", "2.45", "2.50", "2.55"):
        _touch(tmp_path / "3-Lattice_Parameter" / a / "vasprun.xml")
        _touch(tmp_path / "3-Lattice_Parameter" / a / "CONTCAR")
    (tmp_path / "5-Fonones").mkdir()
    _proyecto_con_supercelda(tmp_path / "5-Fonones")
    f = scan_folder(tmp_path)
    assert f.scf_vasprun is not None and "3-Lattice_Parameter" in str(f.scf_vasprun)


def test_nelect_se_compara_por_atomo_entre_celda_unidad_y_supercelda():
    from electronanopy.consistency import CalcSettings, compare_settings
    uni = CalcSettings(source="SCF", encut=500, nelect=8.0, natoms=2, functional="PBE", potcars=("PAW_PBE C",))
    sup = CalcSettings(source="Fonones", encut=500, nelect=200.0, natoms=50, functional="PBE", potcars=("PAW_PBE C",))
    assert compare_settings({"SCF": uni, "Fonones": sup}) == []
    sup.nelect = 196.0
    assert any("por átomo" in w for w in compare_settings({"SCF": uni, "Fonones": sup}))


def test_el_vasprun_junto_a_phonopy_disp_se_reconoce_como_desplazamiento(tmp_path):
    from electronanopy.consistency import phonon_displacement_vasprun
    (tmp_path / "vasprun.xml").write_text("<modeling/>")
    assert phonon_displacement_vasprun(tmp_path) == tmp_path / "vasprun.xml"


# --- Validación con el cálculo real de grafeno 5×5×1 del usuario (PBE, ENCUT 500, un desplazamiento) ---
LIT_GGA = {  # cm^-1, Maultzsch et al., PRL 92, 075501 (2004), Tabla 1, columna GGA; orden ascendente
    "Γ": [0, 0, 0, 884, 1569, 1569], "M": [476, 634, 640, 1338, 1346, 1396], "K": [539, 539, 1004, 1221, 1221, 1289]}
THZ = 33.35641


def test_fonones_de_grafeno_reales_concuerdan_con_la_literatura():
    from electronanopy.phonons import load_phonopy
    ph = load_phonopy(DATA / "grafeno_5x5")
    qs = {"Γ": [0, 0, 0], "M": [.5, 0, 0], "K": [1 / 3, 1 / 3, 0]}
    ph.run_qpoints(list(qs.values()))
    f = ph.get_qpoints_dict()["frequencies"] * THZ
    for row, (q, lit) in zip(f, LIT_GGA.items()):
        lit = np.array(lit, float)
        if q == "Γ":
            assert np.abs(row[:3]).max() < 1.0                                  # acústicos ≈ 0
        # TO en K (+57 cm⁻¹ con 5×5×1): K = 1/3 no es conmensurable con la supercelda y el modo tiene anomalía de Kohn
        tol = 70.0 if q == "K" else 25.0
        assert np.abs(row[3:] - lit[3:]).max() < tol, (q, row, lit)
        if q != "K":
            assert np.abs(row[:3] - lit[:3]).max() < 25.0 or q == "Γ"
    assert abs(f[1][3:] - np.array(LIT_GGA["M"][3:])).max() < 15       # M converge mejor (≤ 15 cm⁻¹)


def test_informe_de_los_fonones_de_grafeno_reales():
    s = Session(phonon=__import__("electronanopy.phonons", fromlist=["load_phonopy"]).load_phonopy(DATA / "grafeno_5x5"))
    r = analyze(s).modules["fonones"]
    assert r.summary["stable"] and r.summary["two_d"] and 47.0 < r.summary["f_max"] < 48.5       # THz


def test_frecuencia_minima_no_se_imprime_como_negativa_cero():
    from electronanopy.phonons.analysis import _z
    assert f"{_z(-0.0002):.3f}" == "0.000" and _z(-0.3) == -0.3


# --- Ruta 3D de seekpath en la base y orientación de la celda del usuario ---
def _rotacion(eje, angulo):
    eje = np.asarray(eje, float) / np.linalg.norm(eje)
    K = np.array([[0, -eje[2], eje[1]], [eje[2], 0, -eje[0]], [-eje[1], eje[0], 0]])
    return np.eye(3) + np.sin(angulo) * K + (1 - np.cos(angulo)) * K @ K


def test_ruta_3d_en_celda_rotada_y_otra_base():
    """NaCl (Oh) en otra base y girada: X debe quedar sobre un eje cúbico y L sobre una diagonal (tras deshacer el giro)."""
    from electronanopy.phonons.paths import path_3d, reciprocal_lattice
    a = 5.64
    L0 = a / 2 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)
    Q = _rotacion([1, 2, 3], 0.7)
    M = np.array([[1, 0, 0], [1, 1, 0], [0, 1, 1]], float)           # otra base de la misma red
    L = (M @ L0) @ Q.T
    frac = np.array([[0, 0, 0], [.5, .5, .5]])                        # Na, Cl en la base original (no afecta a la red)
    segs, _ = path_3d(L, np.array([11, 17]), frac @ np.linalg.inv(M) @ np.eye(3))
    rec = reciprocal_lattice(L)
    puntos = {}
    for la, lb, qa, qb in segs:
        puntos[la], puntos[lb] = qa @ rec @ Q, qb @ rec @ Q           # al marco cúbico original
    X, Lp = puntos["X"], puntos["L"]
    assert np.sort(np.abs(X))[:2].max() < 1e-6 and abs(np.linalg.norm(X) - 2 * np.pi / a) < 1e-6          # X = (2π/a)(1,0,0)
    assert np.allclose(np.abs(Lp), np.abs(Lp)[0], atol=1e-6) and abs(np.linalg.norm(Lp) - np.pi / a * np.sqrt(3)) < 1e-6     # L = (π/a)(1,1,1)


def test_bandas_iguales_con_otra_base_de_la_celda_primitiva():
    import phonopy
    from electronanopy.phonons import compute_bands
    from electronanopy.phonons.paths import reciprocal_lattice
    d = str(DATA / "nacl") + "/"
    pm = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    M = np.array([[1, 0, 0], [1, 1, 0], [0, 1, 1]])
    cargar = lambda p: phonopy.load(phonopy_yaml=d + "phonopy_disp.yaml", force_sets_filename=d + "FORCE_SETS", is_nac=False,
                                    primitive_matrix=np.asarray(p).tolist(), produce_fc=True, log_level=0)
    A, B = cargar(pm), cargar(pm @ M)
    bA, bB = compute_bands(A, npoints=7), compute_bands(B, npoints=7)
    cA = bA.qpoints @ reciprocal_lattice(np.asarray(A.primitive.cell))
    cB = bB.qpoints @ reciprocal_lattice(np.asarray(B.primitive.cell))
    assert not np.allclose(bA.qpoints, bB.qpoints)                    # la base cambió, así que hubo que transformar
    assert np.allclose(cA, cB, atol=1e-6) and np.allclose(bA.frequencies, bB.frequencies, atol=1e-6)


def test_acusticos_en_gamma_no_se_imprimen_como_negativo_cero():
    from electronanopy.phonons import load_phonopy
    r = analyze(Session(phonon=load_phonopy(DATA / "grafeno_5x5"))).modules["fonones"]
    linea = next(h for h in r.highlights if h.startswith("Modos acústicos en Γ"))
    assert "-0.00" not in linea


def test_chgcar_en_otra_serie_se_usa_si_falta_en_la_principal(tmp_path):
    from test_project import _touch
    for a in ("2.40", "2.45", "2.50"):
        _touch(tmp_path / "3-Lattice_Parameter" / a / "vasprun.xml")
    for e in ("400", "500", "600"):
        _touch(tmp_path / "2-Encuts" / e / "vasprun.xml")
        _touch(tmp_path / "2-Encuts" / e / "CHGCAR")
    f = scan_folder(tmp_path)
    assert f.chgcar is not None and "2-Encuts" in str(f.chgcar)
    assert any("CHGCAR tomado de la serie" in n for n in f.notes)
