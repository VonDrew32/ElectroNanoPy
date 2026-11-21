"""Reglas entre módulos (OE4), Dirac en 3D y enlaces en metales (v0.12)."""
import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from electronanopy.electronic import dirac_points
from electronanopy.electronic import _synthetic as SE
from electronanopy.phonons import _synthetic as SP
from electronanopy.project import Session
from electronanopy.report import analyze, Thresholds
from electronanopy.structural import bond_table, coordination_numbers
from test_auditoria import _bandas, SEMI
from test_project import _estructura_de_la_celda_de_fonones, _sesion_sintetica


def _graphene_session(**kw):
    d = dict(structure=_estructura_de_la_celda_de_fonones(), bands=SE.graphene_bands(), tdos=SE.graphene_dos(),
             scf_efermi=0.0, bands_efermi_native=0.0, nelect=2.0)
    d.update(kw)
    return Session(**d)


def test_grafeno_semimetal_2d_con_dirac():
    m = analyze(_graphene_session()).modules["integracion"]
    assert any("semimetálico con cono de Dirac en" in h for h in m.highlights)
    assert m.params["tab"] == "electronico"


def test_estable_y_clasificacion():
    m = analyze(_graphene_session(phonon=SP.triangular_2d())).modules["integracion"]
    assert any(h.startswith("Dinámicamente estable") for h in m.highlights)


def test_inestable_sugiere_fonones_y_causas():
    ph = SP.triangular_2d(kz=-0.3)                           # constante de fuerza fuera del plano negativa
    rep = analyze(_graphene_session(phonon=ph))
    m = rep.modules["integracion"]
    assert rep.modules["fonones"].summary["stable"] is False
    assert m.params["tab"] == "fonones" and any("frecuencias imaginarias" in w for w in m.warnings)


def test_sugiere_usar_a0_si_la_estructura_difiere():
    rep = analyze(_graphene_session())
    r = rep.modules["integracion"]
    rep.modules["estructura"].params["a0"] = 2.60            # simulado: el mínimo de la serie
    from electronanopy.report import _integration
    r2 = _integration(_graphene_session(), Thresholds(), rep.modules)
    assert any("use a₀" in t for t in r2.suggestions)


def test_sin_datos_no_sugiere_nada():
    m = analyze(Session()).modules["integracion"]
    assert not m.available and "tab" not in m.params


def test_la_ventana_abre_la_pestana_sugerida(qtbot=None):
    from PySide6.QtWidgets import QApplication
    from electronanopy.app.main_window import MainWindow
    app = QApplication.instance() or QApplication([])
    w = MainWindow()
    w.set_session(_graphene_session())
    assert w.tabs.currentIndex() == 2 and "sugerida" in w.statusBar().currentMessage()


# --- Dirac en 3D: una degeneración en todo el camino no es un punto de Dirac ---
def test_degeneracion_general_no_es_dirac():
    # bandas dobles (Kramers): pares degenerados en todos los k, justo en E_F
    deg = np.vstack([SEMI[:3] - 3.0, np.tile([[-0.001], [0.001]], (1, 5)), SEMI[4:] + 3.0])
    assert dirac_points(_bandas(deg, efermi=0.0)) == []


def test_cruce_aislado_si_es_dirac():
    e = np.array([[-1, -0.5, 0.0, -0.5, -1.0], [1, 0.5, 0.0, 0.5, 1.0]], float)
    pts = dirac_points(_bandas(e, efermi=0.0))
    assert len(pts) == 1 and pts[0]["k_index"] == 2


# --- Enlaces: dos metales no se enlazan a 1.15·(r_i + r_j) ---
def test_metales_vecinos_lejanos_no_son_enlace():
    st = Structure(Lattice.tetragonal(4.5, 8.0), ["Sr", "Sr"], [[0, 0, 0], [0.5, 0.5, 0.5]])   # Sr–Sr a 4.5 Å
    assert bond_table(st) == []
    cu = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.61), ["Cu"], [[0, 0, 0]])
    assert set(coordination_numbers(cu)) == {12}                       # Cu–Cu a 2.55 Å sí es enlace


def test_enlace_covalente_conserva_tolerancia():
    gr = _estructura_de_la_celda_de_fonones()
    assert set(coordination_numbers(gr)) == {3}
