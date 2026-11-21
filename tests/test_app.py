import time
import pytest
pytest.importorskip("PySide6")
import matplotlib
from PySide6.QtWidgets import QApplication
from electronanopy.app.main_window import MainWindow, TABS
from electronanopy.project import Session
from electronanopy.phonons import _synthetic as SP
from electronanopy.electronic import _synthetic as SE


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _sesion(**kw):
    d = dict(structure=SE.graphene_structure(), bands=SE.graphene_bands(), tdos=SE.graphene_dos(),
             chg=SE.graphene_chgcar(), scf_efermi=0.0, bands_efermi_native=0.0, nelect=8.0, phonon=SP.triangular_2d(),
             magnetization={"spin_polarized": False, "total": None, "per_atom": None})
    d.update(kw)
    return Session(**d)


@pytest.fixture
def ventana(app):
    w = MainWindow(); w.show()
    yield w
    w.close()


def _axes(w, i):
    return [p.figure for p in (w.p_struct, w.p_ph, w.p_el, w.p_den)][i].axes


def test_sin_sesion_las_pestanas_estan_deshabilitadas(ventana):
    assert not any(ventana.tabs.isTabEnabled(i) for i in range(4))
    assert "Abra la carpeta" in ventana.summary.toPlainText()


def test_todas_las_pestanas_dibujan(ventana):
    ventana.set_session(_sesion())
    assert all(ventana.tabs.isTabEnabled(i) for i in range(4))
    for i in range(4):
        ventana.tabs.setCurrentIndex(i)
        assert _axes(ventana, i), f"la pestaña {TABS[i][1]} no dibujó nada"
    txt = ventana.summary.toPlainText()
    assert "Punto(s) de Dirac en K" in txt and "Electrones integrados" in txt and "Sugerencia" in txt


def test_valores_por_defecto_vienen_de_las_sugerencias(ventana):
    ventana.set_session(_sesion())
    assert [w.value() for w in ventana.sc] == [3, 3, 1]                  # supercelda 2D sugerida
    assert ventana.d_frac.value() == pytest.approx(0.5)
    assert ventana.d_level.value() > 0


def test_opciones_de_visualizacion(ventana):
    ventana.set_session(_sesion())
    ventana.tabs.setCurrentIndex(0)
    ventana.sc[0].setValue(2); ventana.bonds.setChecked(False)           # redibuja sin error
    assert _axes(ventana, 0)
    ventana.tabs.setCurrentIndex(1)
    ventana.ph_unit.setCurrentText("cm-1")
    assert "cm" in _axes(ventana, 1)[0].get_ylabel()
    ventana.ph_dos.setChecked(False)
    assert len(_axes(ventana, 1)) == 1
    ventana.tabs.setCurrentIndex(2)
    assert ventana.e_orb.isEnabled()
    ventana.e_orb.setCurrentText("pz"); ventana.e_dos.setChecked(False)
    assert len(_axes(ventana, 2)) >= 1
    ventana.tabs.setCurrentIndex(3)
    for v in (1, 2, 0):
        ventana.d_view.setCurrentIndex(v)
        assert _axes(ventana, 3)


def test_sin_proyecciones_se_deshabilita_el_peso_orbital(ventana):
    ventana.set_session(_sesion(bands=SE.graphene_bands(with_projections=False)))
    assert not ventana.e_orb.isEnabled()


def test_modulos_faltantes_quedan_deshabilitados_con_motivo(ventana):
    ventana.set_session(_sesion(bands=None, tdos=None, chg=None))
    assert ventana.tabs.isTabEnabled(0) and ventana.tabs.isTabEnabled(1)
    assert not ventana.tabs.isTabEnabled(2) and not ventana.tabs.isTabEnabled(3)
    assert ventana.tabs.tabToolTip(3) != ""


def test_error_al_dibujar_no_rompe_la_ventana(ventana):
    ventana.set_session(_sesion())
    ventana.session.phonon = object()                                    # objeto inválido
    ventana._dirty.add(1); ventana.tabs.setCurrentIndex(1); ventana._redraw_current()
    assert "Error al dibujar" in ventana.statusBar().currentMessage()


def test_cargar_carpeta_real(ventana, tmp_path):
    SE.graphene_structure().to(filename=str(tmp_path / "CONTCAR"), fmt="poscar")
    SP.save_as_files(SP.triangular_2d(), tmp_path / "fonones", "params")
    ventana.load_folder(tmp_path, background=False)
    items = [ventana.file_list.item(i).text() for i in range(ventana.file_list.count())]
    assert any(t.startswith("✓ Estructura") for t in items) and any(t.startswith("✓ Fonones") for t in items)
    assert any(t.startswith("✗ CHGCAR") for t in items)
    assert ventana.tabs.isTabEnabled(0) and ventana.tabs.isTabEnabled(1) and not ventana.tabs.isTabEnabled(3)


def test_carga_en_segundo_plano(ventana, tmp_path):
    SE.graphene_structure().to(filename=str(tmp_path / "CONTCAR"), fmt="poscar")
    ventana.load_folder(tmp_path)
    t0 = time.time()
    while (ventana._loading or ventana.session is None) and time.time() - t0 < 30:
        QApplication.processEvents(); time.sleep(0.01)
    assert ventana.session is not None and ventana.tabs.isTabEnabled(0)


def test_carpeta_invalida_no_falla(ventana, tmp_path):
    ventana.load_folder(tmp_path / "no_existe")
    assert "No es una carpeta" in ventana.statusBar().currentMessage()


def test_series_de_calculos_en_la_pestana_estructura(ventana, tmp_path):
    from test_series import _grafeno
    from electronanopy.project import scan_folder, load_session
    ventana.set_session(load_session(scan_folder(_grafeno(tmp_path))))
    w = ventana
    assert w.s_view.isEnabled() and w.s_series.count() == 4
    assert w.s_member.count() == 1 + 5 and "por defecto" in w.s_member.itemText(1 + 2)      # serie de red: 2.460 es el 3.º
    w.s_view.setCurrentIndex(1)                                                              # vista de serie
    w.tabs.setCurrentIndex(0); w._redraw_current()
    assert len(_axes(w, 0)) >= 3                                                             # E(a), σ y tabla
    w.s_view.setCurrentIndex(0); w.s_series.setCurrentIndex(1); w.s_member.setCurrentIndex(2)
    w._redraw_current()
    assert len(_axes(w, 0)) == 1 and "No se pudo" not in w.statusBar().currentMessage()      # estructura 3D del cálculo


def test_sin_series_los_selectores_quedan_deshabilitados(ventana):
    ventana.set_session(_sesion())
    assert not ventana.s_view.isEnabled() and not ventana.s_series.isEnabled()


def test_menu_ver_reabre_los_paneles_cerrados(ventana):
    ventana.dock_files.close(); ventana.dock_summary.close()
    assert not ventana.dock_files.isVisible() and not ventana.dock_summary.isVisible()
    ver = next(a for a in ventana.menuBar().actions() if a.text().replace('&', '') == 'Ver')
    acciones = {a.text(): a for a in ver.menu().actions()}
    acciones["Archivos detectados"].trigger(); acciones["Resumen"].trigger()
    assert ventana.dock_files.isVisible() and ventana.dock_summary.isVisible()


def test_selector_de_camino_y_asr_en_la_pestana_fonones(ventana):
    ventana.set_session(_sesion(phonon=SP.cubic(onsite=0.02)))
    ventana.tabs.setCurrentIndex(1); ventana._redraw_current()
    assert "inestab" not in ventana.p_ph.figure._suptitle.get_text().lower()
    ventana.ph_asr.setChecked(True)
    ventana.ph_path.setCurrentIndex(2)                                # 3D (seekpath)
    ventana._redraw_current()
    assert "ASR impuesta" in ventana.p_ph.figure._suptitle.get_text()
    assert "No se pudo dibujar" not in ventana.statusBar().currentMessage()
    ventana.ph_path.setCurrentIndex(1)                                # forzar 2D sobre un material 3D no debe romper
    ventana._redraw_current()
    assert "No se pudo dibujar" not in ventana.statusBar().currentMessage()


def test_exportar_figura_y_reportes(ventana, tmp_path):
    ventana.set_session(_sesion())
    for ext in ("png", "svg", "pdf"):
        out = ventana.export_figure(str(tmp_path / f"fig.{ext}"))
        assert out is not None and out.stat().st_size > 100
    assert ventana.export_report(str(tmp_path / "r.json")).exists()
    informe = ventana.export_full_report(str(tmp_path / "informe.html"))
    html = informe.read_text(encoding="utf-8")
    assert html.count("data:image/png;base64,") == 4          # estructura, fonones, electrónico, densidad
    assert "Consistencia entre cálculos" in html


def test_exportar_sin_sesion_no_falla(ventana, tmp_path):
    assert ventana.export_figure(str(tmp_path / "f.png")) is None
    assert ventana.export_report(str(tmp_path / "r.json")) is None
    assert ventana.export_full_report(str(tmp_path / "i.html")) is None
