import numpy as np
import pytest
from electronanopy.project import scan_folder, load_session, Session, ProjectFiles, is_line_mode_kpoints
from electronanopy.report import analyze, Thresholds
from electronanopy.phonons import _synthetic as SP
from electronanopy.electronic import _synthetic as SE
from electronanopy.structural.io import StructureReadError

KPOINTS_LINEA = "Camino\n60\nLine-mode\nReciprocal\n0 0 0 ! G\n0.5 0 0 ! M\n"
KPOINTS_MALLA = "Malla\n0\nGamma\n24 24 1\n"


def _touch(p, text="x"):
    p.parent.mkdir(parents=True, exist_ok=True); p.write_text(text); return p


def _arbol(tmp_path):
    for n in ("01_relax", "02_scf", "03_bandas", "04_dos"):
        (tmp_path / n).mkdir()
    _touch(tmp_path / "01_relax/CONTCAR"); _touch(tmp_path / "01_relax/POSCAR")
    _touch(tmp_path / "02_scf/vasprun.xml"); _touch(tmp_path / "02_scf/OUTCAR"); _touch(tmp_path / "02_scf/CHGCAR")
    _touch(tmp_path / "02_scf/KPOINTS", KPOINTS_MALLA)
    _touch(tmp_path / "03_bandas/vasprun.xml"); _touch(tmp_path / "03_bandas/KPOINTS", KPOINTS_LINEA)
    _touch(tmp_path / "04_dos/vasprun.xml")
    SP.save_as_files(SP.triangular_2d(), tmp_path / "05_fonones", "params")
    _touch(tmp_path / "05_fonones/disp-001/vasprun.xml")            # no debe confundirse con bandas/SCF
    return tmp_path


def test_kpoints_en_modo_linea(tmp_path):
    assert is_line_mode_kpoints(_touch(tmp_path / "a", KPOINTS_LINEA))
    assert not is_line_mode_kpoints(_touch(tmp_path / "b", KPOINTS_MALLA))
    assert not is_line_mode_kpoints(tmp_path / "no_existe")


def test_deteccion_de_archivos_en_la_disposicion_recomendada(tmp_path):
    pf = scan_folder(_arbol(tmp_path))
    assert pf.structure == (tmp_path / "01_relax/CONTCAR").resolve()
    assert pf.scf_vasprun.parent.name == "02_scf" and pf.dos_vasprun.parent.name == "04_dos"
    assert pf.bands_vasprun.parent.name == "03_bandas" and pf.bands_kpoints.name == "KPOINTS"
    assert pf.outcar.parent.name == "02_scf" and pf.chgcar.parent.name == "02_scf"
    assert pf.phonons.name == "05_fonones"


def test_deteccion_por_kpoints_sin_pistas_en_el_nombre(tmp_path):
    _touch(tmp_path / "a/vasprun.xml"); _touch(tmp_path / "a/KPOINTS", KPOINTS_LINEA)
    _touch(tmp_path / "b/vasprun.xml"); _touch(tmp_path / "b/KPOINTS", KPOINTS_MALLA)
    pf = scan_folder(tmp_path)
    assert pf.bands_vasprun.parent.name == "a"
    assert pf.scf_vasprun.parent.name == "b" and pf.dos_vasprun == pf.scf_vasprun and pf.notes


def test_carpeta_inexistente_y_vacia(tmp_path):
    with pytest.raises(NotADirectoryError):
        scan_folder(tmp_path / "no_existe")
    pf = scan_folder(tmp_path)
    assert all(v is None for k, v in pf.found().items())


def test_sesion_tolerante_a_fallos(tmp_path):
    _arbol(tmp_path)
    s = load_session(scan_folder(tmp_path))
    assert s.phonon is not None                                    # los fonones sí se leen
    assert "estructura" in s.errors                                # CONTCAR de relleno: error aislado
    assert "electronico" in s.errors and "densidad" in s.errors    # vasprun/CHGCAR falsos
    assert s.available()["fonones"] and not s.available()["electronico"]


def test_sesion_con_estructura_real(tmp_path):
    from pymatgen.core import Structure
    SE.graphene_structure().to(filename=str(tmp_path / "CONTCAR"), fmt="poscar")
    s = load_session(scan_folder(tmp_path))
    assert s.structure is not None and s.errors == {}


def _estructura_de_la_celda_de_fonones():
    """Grafeno con la celda del modelo de fonones sintético (a = 2.5 Å), para que ambas sean coherentes."""
    from pymatgen.core import Lattice, Structure
    return Structure(Lattice.hexagonal(2.5, 15.0), ["C", "C"], [[1/3, 2/3, .5], [2/3, 1/3, .5]])


def _sesion_sintetica(**kw):
    d = dict(structure=_estructura_de_la_celda_de_fonones(), bands=SE.graphene_bands(), tdos=SE.graphene_dos(),
             chg=SE.graphene_chgcar(), scf_efermi=0.0, bands_efermi_native=0.0, nelect=8.0,
             phonon=SP.triangular_2d(), magnetization={"spin_polarized": False, "total": None, "per_atom": None})
    d.update(kw)
    return Session(**d)


def test_reporte_resalta_lo_relevante():
    r = analyze(_sesion_sintetica())
    txt = " ".join(h for m in r.modules.values() for h in m.highlights)
    assert "P6/mmm" in txt and "Material bidimensional" in txt
    assert "Punto(s) de Dirac en K" in txt and "ħv_F" in txt
    assert "Electrones integrados en la celda: 8.000" in txt
    assert r.modules["estructura"].params["supercell"] == (3, 3, 1)
    assert r.modules["densidad"].params["slice_frac"] == pytest.approx(0.5)
    assert r.all_warnings() == []                                # sesión coherente: ninguna advertencia


def test_reporte_detecta_incoherencias():
    r = analyze(_sesion_sintetica(scf_efermi=0.3))
    assert any("difieren en 0.300 eV" in w for w in r.modules["electronico"].warnings)
    r = analyze(_sesion_sintetica(nelect=10.0))
    assert any("NELECT" in w for w in r.modules["densidad"].warnings)


def test_reporte_con_inestabilidad_y_supercelda_pequena():
    r = analyze(_sesion_sintetica(phonon=SP.cubic(k2=-0.2)))
    w = r.modules["fonones"].warnings
    assert any("frecuencias imaginarias" in x for x in w)
    assert not any("Supercelda" in x for x in w)                 # 3×3×3 cumple el mínimo de 3 celdas
    r2 = analyze(_sesion_sintetica(phonon=SP.triangular_2d()), Thresholds(min_supercell=5))
    assert any("Supercelda 4×4×1" in x for x in r2.modules["fonones"].warnings)


def test_modulos_ausentes_dan_mensajes_claros():
    r = analyze(Session(files=ProjectFiles(root=None)))
    assert not any(m.available for m in r.modules.values())
    assert all(m.warnings for k, m in r.modules.items() if k not in ("consistencia", "integracion"))


def test_reporte_avisa_de_asr_y_de_supercelda_corta():
    w = analyze(_sesion_sintetica(phonon=SP.cubic(onsite=0.02))).modules["fonones"].warnings
    assert any("regla de suma acústica" in x for x in w)
    assert any("9.0 Å" in x for x in w)                              # 3×3×3 de 3 Å: cumple el nº de celdas, no los 10 Å
    ok = analyze(_sesion_sintetica(phonon=SP.triangular_2d())).modules["fonones"]
    assert not any("suma acústica" in x or "Å (<" in x for x in ok.warnings)   # 4×4 de 2.5 Å = 10 Å
    assert any("Modos acústicos en Γ" in h for h in ok.highlights)
