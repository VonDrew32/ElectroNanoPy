import numpy as np
from electronanopy.consistency import read_settings, compare_settings, gather_settings, phonon_displacement_vasprun
from electronanopy.project import ProjectFiles, Session
from electronanopy.report import analyze
from electronanopy.phonons import _synthetic as SP
from electronanopy.electronic import _synthetic as SE

XML = """<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
 <incar><i type="string" name="PREC">normal</i><i name="ENCUT">{encut}</i></incar>
 <kpoints><generation param="Gamma"><v type="int" name="divisions">{k} {k} 1</v></generation></kpoints>
 <parameters><separator name="electronic"><i name="ENMAX">{encut}</i><i type="int" name="ISPIN">{ispin}</i>
 <i name="NELECT">8.0</i></separator><separator name="electronic exchange-correlation"><i name="GGA">{gga}</i>
 <i type="logical" name="LHFCALC">F</i></separator></parameters>
 <atominfo><array name="atomtypes"><set><rc><c>2</c><c>C </c><c>12.0</c><c>4.0</c><c>  {pot}  </c></rc></set></array></atominfo>
</modeling>"""


def _vr(path, encut=600, k=18, ispin=1, gga="--", pot="PAW_PBE C 08Apr2002"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(XML.format(encut=encut, k=k, ispin=ispin, gga=gga, pot=pot))
    return path


def test_lee_los_ajustes(tmp_path):
    s = read_settings(_vr(tmp_path / "vasprun.xml"))
    assert s.encut == 600 and s.ispin == 1 and s.nelect == 8.0
    assert s.functional == "PBE" and s.kmesh == (18, 18, 1) and s.lhfcalc is False


def test_calculos_iguales_no_avisan(tmp_path):
    a = read_settings(_vr(tmp_path / "a/vasprun.xml"), "SCF")
    b = read_settings(_vr(tmp_path / "b/vasprun.xml", k=30), "Bandas")      # la malla k puede diferir
    assert compare_settings({"SCF": a, "Bandas": b}) == []


def test_detecta_encut_funcional_pseudopotencial_y_espin(tmp_path):
    a = read_settings(_vr(tmp_path / "a/vasprun.xml"), "SCF")
    b = read_settings(_vr(tmp_path / "b/vasprun.xml", encut=400, ispin=2, pot="PAW_LDA C 08Apr2002"), "Bandas")
    w = " ".join(compare_settings({"SCF": a, "Bandas": b}))
    assert "Corte de energía" in w and "400" in w
    assert "Funcional distinto" in w and "Pseudopotenciales" in w and "ISPIN" in w


def test_tolerancia_del_corte(tmp_path):
    a = read_settings(_vr(tmp_path / "a/vasprun.xml", encut=600), "A")
    b = read_settings(_vr(tmp_path / "b/vasprun.xml", encut=600.5), "B")
    assert compare_settings({"A": a, "B": b}, encut_tol=1.0) == []


def test_un_solo_calculo_no_se_compara(tmp_path):
    a = read_settings(_vr(tmp_path / "a/vasprun.xml"), "A")
    assert compare_settings({"A": a}) == []


def test_reune_los_calculos_del_proyecto_sin_repetir(tmp_path):
    scf = _vr(tmp_path / "02_scf/vasprun.xml")
    bands = _vr(tmp_path / "03_bandas/vasprun.xml", encut=400)
    disp = _vr(tmp_path / "05_fonones/disp-001/vasprun.xml")
    files = ProjectFiles(root=tmp_path, scf_vasprun=scf, dos_vasprun=scf, bands_vasprun=bands,
                         phonons=tmp_path / "05_fonones")
    named = gather_settings(files)
    assert set(named) == {"SCF", "Bandas", "Fonones"}               # la DOS repite el SCF
    assert phonon_displacement_vasprun(tmp_path / "05_fonones") == disp


def test_reporte_incluye_consistencia(tmp_path):
    files = ProjectFiles(root=tmp_path, scf_vasprun=_vr(tmp_path / "a/vasprun.xml"),
                         bands_vasprun=_vr(tmp_path / "b/vasprun.xml", encut=400))
    rep = analyze(Session(files=files, structure=SE.graphene_structure()))
    m = rep.modules["consistencia"]
    assert m.available and any("Corte de energía" in w for w in m.warnings)
    assert any("Consistencia" in w for w in rep.all_warnings())


def test_reporte_coherente_resalta(tmp_path):
    files = ProjectFiles(root=tmp_path, scf_vasprun=_vr(tmp_path / "a/vasprun.xml"),
                         bands_vasprun=_vr(tmp_path / "b/vasprun.xml"))
    m = analyze(Session(files=files)).modules["consistencia"]
    assert not m.warnings and m.highlights


def test_celda_de_fonones_distinta_de_la_estructura():
    ph, st = SP.triangular_2d(), SE.graphene_structure()
    rep = analyze(Session(structure=st, phonon=ph))
    ls = np.array(st.lattice.abc); lp = np.linalg.norm(np.asarray(ph.unitcell.cell), axis=1)
    hay = any("celda de los fonones" in w for w in rep.modules["consistencia"].warnings)
    assert hay == bool((np.abs(ls - lp) > 0.005).any())


def test_sin_archivos_no_falla():
    m = analyze(Session()).modules["consistencia"]
    assert not m.available and not m.warnings
