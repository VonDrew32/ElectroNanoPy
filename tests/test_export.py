import json
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from electronanopy.export import (report_to_dict, report_to_json, report_to_html, save_report, save_figure,
                                  figure_png_bytes)
from electronanopy.project import Session
from electronanopy.report import analyze
from electronanopy.electronic import _synthetic as SE


@pytest.fixture
def reporte():
    return analyze(Session(structure=SE.graphene_structure(), bands=SE.graphene_bands(), tdos=SE.graphene_dos(),
                           scf_efermi=0.0, bands_efermi_native=0.0, nelect=8.0))


def test_json_es_valido_y_completo(reporte):
    d = json.loads(report_to_json(reporte, "proyecto"))
    assert d["proyecto"] == "proyecto" and set(d["modulos"]) >= {"estructura", "electronico", "consistencia"}
    assert d["modulos"]["estructura"]["resaltados"]


def test_html_escapa_e_incrusta_figuras(reporte):
    reporte.modules["estructura"].warnings.append("a < b & <script>")
    fig, ax = plt.subplots(); ax.plot([0, 1])
    html = report_to_html(reporte, "p", {"Estructura": figure_png_bytes(fig, 50)})
    assert "&lt;script&gt;" in html and "<script>" not in html
    assert "data:image/png;base64," in html


def test_guardar_segun_extension(reporte, tmp_path):
    assert json.loads(save_report(reporte, tmp_path / "r.json").read_text(encoding="utf-8"))
    assert "<html" in save_report(reporte, tmp_path / "r.html").read_text(encoding="utf-8")
    with pytest.raises(ValueError):
        save_report(reporte, tmp_path / "r.txt")


@pytest.mark.parametrize("ext", ["png", "svg", "pdf"])
def test_guardar_figura(tmp_path, ext):
    fig, ax = plt.subplots(); ax.plot([0, 1])
    assert save_figure(fig, tmp_path / f"f.{ext}", dpi=80).stat().st_size > 100


def test_figura_formato_no_admitido(tmp_path):
    with pytest.raises(ValueError):
        save_figure(plt.figure(), tmp_path / "f.bmp")


def test_valores_numpy_y_no_finitos_se_serializan(reporte):
    import numpy as np
    reporte.modules["estructura"].summary["x"] = {"a": np.float64(1.5), "b": np.array([1, 2]), "c": float("nan")}
    d = json.loads(report_to_json(reporte))
    assert d["modulos"]["estructura"]["resumen"]["x"] == {"a": 1.5, "b": [1, 2], "c": None}


def test_resumen_en_texto_de_una_carpeta(tmp_path, capsys, monkeypatch):
    from test_series import _grafeno
    from electronanopy.export import summarize_folder
    txt = summarize_folder(_grafeno(tmp_path))
    assert "== Estructura ==" in txt and "Birch–Murnaghan" in txt and "== Consistencia entre cálculos ==" in txt
    import sys
    from electronanopy.app.main_window import run
    monkeypatch.setattr(sys, "argv", ["electronanopy", "--resumen", str(tmp_path / "Graphene")])
    run()
    assert "Birch–Murnaghan" in capsys.readouterr().out


def test_estilos_cambian_el_tamano_y_restauran_la_figura(tmp_path):
    from PIL import Image
    from electronanopy.export import STYLES, styled
    fig, ax = plt.subplots(figsize=(6, 4)); ax.plot([0, 1], [0, 1], lw=1.5); ax.set_xlabel("x", fontsize=10)
    dims = {}
    for st in STYLES:
        p = save_figure(fig, tmp_path / f"{st}.png", style=st)
        dims[st] = Image.open(p).size
    assert dims["publicacion_1col"][0] < dims["publicacion_2col"][0]
    assert 6.0 < dims["publicacion_1col"][0] / 600 * 2.54 <= 8.6                  # ≤ 8.5 cm a 600 dpi (bbox tight recorta márgenes)
    assert 28.0 < dims["poster"][0] / 300 * 2.54 <= 35.5
    assert tuple(fig.get_size_inches()) == (6, 4)
    assert ax.xaxis.label.get_fontsize() == 10 and ax.lines[0].get_linewidth() == 1.5
    with pytest.raises(ValueError):
        save_figure(fig, tmp_path / "x.png", style="inexistente")


def test_la_ventana_exporta_con_el_estilo_elegido(tmp_path):
    from PIL import Image
    from PySide6.QtWidgets import QApplication
    from electronanopy.app.main_window import MainWindow
    from test_integracion import _graphene_session
    app = QApplication.instance() or QApplication([])
    w = MainWindow(); w.set_session(_graphene_session())
    a = Image.open(w.export_figure(str(tmp_path / "a.png"))).size
    w.export_style = "poster"
    b = Image.open(w.export_figure(str(tmp_path / "b.png"))).size
    assert b != a
    stylem = [m for m in w.menuBar().findChildren(type(w.menuBar().addMenu("x"))) if m.title().startswith("Estilo")]
    assert stylem
