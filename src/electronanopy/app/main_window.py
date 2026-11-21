"""Ventana principal: carpeta del proyecto -> archivos detectados -> pestañas de visualización + resumen inteligente."""
from __future__ import annotations
import os
import sys
from pathlib import Path

os.environ.setdefault("QT_API", "pyside6")           # matplotlib debe usar PySide6

from PySide6.QtCore import Qt, QObject, QThread, Signal
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDockWidget, QDoubleSpinBox, QFileDialog,
                               QFormLayout, QHBoxLayout, QLabel, QListWidget, QMainWindow, QSpinBox, QTabWidget,
                               QTextBrowser, QVBoxLayout, QWidget)
from matplotlib.ticker import MaxNLocator
from PySide6.QtGui import QAction, QActionGroup, QKeySequence
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

from ..project import scan_folder, load_session, Session
from ..report import analyze, Report

TABS = [("estructura", "Estructura"), ("fonones", "Fonones"), ("electronico", "Electrónico"), ("densidad", "Densidad")]


class _Worker(QObject):
    finished = Signal(object)
    progress = Signal(str)

    def __init__(self, files):
        super().__init__(); self.files = files

    def run(self):
        self.finished.emit(load_session(self.files, self.progress.emit))


class _Panel(QWidget):
    """Controles a la izquierda y lienzo de matplotlib (con barra de herramientas) a la derecha."""

    def __init__(self, on_change):
        super().__init__()
        self.on_change = on_change
        self.form = QFormLayout()
        side = QWidget(); side.setFixedWidth(270)
        sl = QVBoxLayout(side); sl.addLayout(self.form); sl.addStretch(1)
        self.figure = Figure(figsize=(6, 4.5), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        right = QVBoxLayout(); right.addWidget(NavigationToolbar2QT(self.canvas, self)); right.addWidget(self.canvas)
        lay = QHBoxLayout(self); lay.addWidget(side); lay.addLayout(right, 1)
        self.message = QLabel(""); self.message.setWordWrap(True); sl.addWidget(self.message)

    def add(self, label, widget, signal=None):
        self.form.addRow(label, widget)
        (signal or getattr(widget, "valueChanged", None) or getattr(widget, "currentIndexChanged", None)
         or getattr(widget, "toggled")).connect(lambda *_: self.on_change())
        return widget

    def show_error(self, text):
        self.figure.clear()
        self.figure.text(0.5, 0.5, text, ha="center", va="center", wrap=True, color="firebrick")
        self.canvas.draw_idle()


def _spin(lo, hi, val, step=1, decimals=None):
    w = QDoubleSpinBox() if decimals is not None else QSpinBox()
    if decimals is not None:
        w.setDecimals(decimals)
    w.setRange(lo, hi); w.setSingleStep(step); w.setValue(val)
    return w


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ElectroNanoPy")
        self.resize(1250, 780)
        self.session: Session | None = None
        self.report: Report | None = None
        self._loading = False
        self._dirty = set()

        # --- acciones ---
        open_act = QAction("Abrir carpeta…", self); open_act.setShortcut(QKeySequence.Open)
        open_act.triggered.connect(self.open_dialog)
        reload_act = QAction("Recargar", self); reload_act.setShortcut(QKeySequence.Refresh)
        reload_act.triggered.connect(lambda: self.root and self.load_folder(self.root))
        tb = self.addToolBar("Proyecto"); tb.addAction(open_act); tb.addAction(reload_act)
        self.root: Path | None = None

        # --- panel de archivos y resumen ---
        self.file_list = QListWidget()
        dock = QDockWidget("Archivos detectados", self); dock.setWidget(self.file_list)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.summary = QTextBrowser(); self.summary.setOpenExternalLinks(False)
        sdock = QDockWidget("Resumen", self); sdock.setWidget(self.summary)
        self.addDockWidget(Qt.BottomDockWidgetArea, sdock)
        filem = self.menuBar().addMenu("&Archivo")
        filem.addAction(open_act)
        from ..export import STYLES
        self.export_style = "pantalla"
        stylem = filem.addMenu("Estilo de exportación de figuras")
        group = QActionGroup(self)
        for key, st in STYLES.items():
            act = QAction(st["label"], self, checkable=True, checked=(key == self.export_style))
            act.triggered.connect(lambda _=False, k=key: setattr(self, "export_style", k))
            group.addAction(act); stylem.addAction(act)
        for text, slot, key in (("Exportar figura actual…", self.export_figure, "Ctrl+E"),
                                ("Exportar reporte (JSON / HTML)…", self.export_report, None),
                                ("Exportar informe con figuras (HTML)…", self.export_full_report, None)):
            act = QAction(text, self); act.triggered.connect(slot)
            if key:
                act.setShortcut(key)
            filem.addAction(act)
        view = self.menuBar().addMenu("&Ver")                       # permite reabrir los paneles cerrados
        view.addAction(dock.toggleViewAction()); view.addAction(sdock.toggleViewAction())
        self.dock_files, self.dock_summary = dock, sdock
        self.summary.setHtml("<p>Abra la carpeta de un proyecto (Ctrl+O) con los resultados de VASP y Phonopy.</p>")

        # --- pestañas ---
        self.tabs = QTabWidget(); self.setCentralWidget(self.tabs)
        self._build_structure_tab(); self._build_phonon_tab(); self._build_electronic_tab(); self._build_density_tab()
        self.tabs.currentChanged.connect(self._redraw_current)
        for i in range(self.tabs.count()):
            self.tabs.setTabEnabled(i, False)

    # ------------------------------------------------------------------ pestañas
    def _build_structure_tab(self):
        p = self.p_struct = _Panel(lambda: self._mark(0))
        self.sc = [p.add(f"Supercelda {ax}", _spin(1, 8, 1)) for ax in "abc"]
        self.bonds = p.add("Mostrar enlaces", QCheckBox(), None); self.bonds.setChecked(True)
        self.tol = p.add("Tolerancia de enlaces", _spin(0.8, 1.8, 1.15, 0.05, 2))
        self.s_angle = p.add("Ángulo", QComboBox()); self.s_angle.addItems(["Automático", "Superior", "Lateral", "Isométrico"])
        self.s_view = p.add("Vista", QComboBox()); self.s_view.addItems(["Estructura", "Serie"])
        self.s_series = p.add("Serie", QComboBox()); self.s_member = p.add("Cálculo", QComboBox())
        self.s_series.currentIndexChanged.connect(lambda *_: self._fill_members())
        for w in (self.s_view, self.s_series, self.s_member):
            w.setEnabled(False)
        self.tabs.addTab(p, TABS[0][1])

    def _build_phonon_tab(self):
        p = self.p_ph = _Panel(lambda: self._mark(1))
        self.ph_unit = p.add("Unidad", QComboBox()); self.ph_unit.addItems(["THz", "cm-1"])
        self.ph_npts = p.add("Puntos por segmento", _spin(21, 401, 101, 10))
        self.ph_mesh = p.add("Malla de la DOS", _spin(8, 60, 24, 2))
        self.ph_dos = p.add("Mostrar DOS", QCheckBox()); self.ph_dos.setChecked(True)
        self.ph_path = p.add("Camino", QComboBox()); self.ph_path.addItems(["Automático", "2D (kz = 0)", "3D (seekpath)"])
        self.ph_asr = p.add("Imponer regla de suma acústica", QCheckBox())
        self.ph_asr.setToolTip("Aplica la regla de suma acústica a las constantes de fuerza (copia; el cálculo original no cambia).")
        self.tabs.addTab(p, TABS[1][1])

    def _build_electronic_tab(self):
        p = self.p_el = _Panel(lambda: self._mark(2))
        self.e_min = p.add("E mínima (eV)", _spin(-30, 0, -6, 1, 1)); self.e_max = p.add("E máxima (eV)", _spin(0, 30, 6, 1, 1))
        self.e_dos = p.add("Mostrar DOS", QCheckBox()); self.e_dos.setChecked(True)
        self.e_orb = p.add("Peso orbital", QComboBox()); self.e_orb.addItems(["Ninguno", "pz", "s", "px + py"])
        self.tabs.addTab(p, TABS[2][1])

    def _build_density_tab(self):
        p = self.p_den = _Panel(lambda: self._mark(3))
        self.d_view = p.add("Vista", QComboBox()); self.d_view.addItems(["Corte 2D", "Promedio planar", "Isosuperficie"])
        self.d_axis = p.add("Eje perpendicular", QComboBox()); self.d_axis.addItems(["a", "b", "c"]); self.d_axis.setCurrentIndex(2)
        self.d_frac = p.add("Posición (fracción)", _spin(0.0, 1.0, 0.5, 0.01, 3))
        self.d_level = p.add("Nivel isosuperficie (e/Å³)", _spin(0.0, 1000.0, 0.5, 0.05, 3))
        self.d_comp = p.add("Componente", QComboBox()); self.d_comp.addItems(["total"])
        self.d_cmap = p.add("Mapa de color", QComboBox()); self.d_cmap.addItems(["magma", "viridis", "inferno", "coolwarm"])
        self.tabs.addTab(p, TABS[3][1])

    # ------------------------------------------------------------------ carga
    def open_dialog(self):
        d = QFileDialog.getExistingDirectory(self, "Carpeta del proyecto", str(self.root or Path.home()))
        if d:
            self.load_folder(d)

    def load_folder(self, path, background: bool = True):
        """Detecta los archivos y carga la sesión (en un hilo aparte salvo ``background=False``)."""
        try:
            files = scan_folder(path)
        except NotADirectoryError as exc:
            self.statusBar().showMessage(str(exc)); return
        self.root = Path(path)
        self._show_files(files)
        if not background:
            self.set_session(load_session(files)); return
        if self._loading:
            return
        self._loading = True
        self.statusBar().showMessage("Cargando…")
        self._thread = QThread(self); self._worker = _Worker(files)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.progress.connect(self.statusBar().showMessage)
        self._worker.finished.connect(self._on_loaded)
        self._thread.start()

    def _on_loaded(self, session):
        self._thread.quit(); self._thread.wait(); self._loading = False
        self.set_session(session)

    def _show_files(self, files):
        self.file_list.clear()
        for name, p in files.found().items():
            shown = "no encontrado" if p is None else str(Path(p).relative_to(files.root.resolve())
                                                          if files.root.resolve() in Path(p).parents else p)
            self.file_list.addItem(("✓ " if p else "✗ ") + f"{name}: {shown}")
        for n in files.notes:
            self.file_list.addItem("• " + n)

    def set_session(self, session: Session):
        """Muestra una sesión ya cargada (también sirve con datos inyectados, p. ej. en pruebas)."""
        self.session = session
        self.report = analyze(session)
        avail = session.available()
        for i, (key, _) in enumerate(TABS):
            self.tabs.setTabEnabled(i, avail[key])
            m = self.report.modules[key]
            self.tabs.setTabToolTip(i, "" if avail[key] else (m.warnings[0] if m.warnings else "No disponible"))
        self._apply_suggestions()
        self.summary.setHtml(self._summary_html())
        first = next((i for i, (k, _) in enumerate(TABS) if avail[k]), None)
        sug = self.report.modules["integracion"].params               # pestaña sugerida por las reglas de integración
        reason = None
        if sug.get("tab") and avail.get(sug["tab"]):
            first = [k for k, _ in TABS].index(sug["tab"])
            reason = sug.get("reason")
        self._dirty = set(range(4))
        if first is not None:
            self.tabs.setCurrentIndex(first)
            self._redraw_current()
        self.statusBar().showMessage((f"Listo. Pestaña sugerida: {TABS[first][1]} ({reason})." if reason else "Listo.")
                                     if first is not None else "No se pudo cargar ningún módulo.")

    def _apply_suggestions(self):
        """Valores por defecto de los controles según las sugerencias del reporte."""
        r = self.report.modules
        for w, v in zip(self.sc, r["estructura"].params.get("supercell", (1, 1, 1))):
            w.blockSignals(True); w.setValue(v); w.blockSignals(False)
        self._fill_series()
        d = r["densidad"].params
        for w, key in ((self.d_frac, "slice_frac"), (self.d_level, "iso_level")):
            if key in d:
                w.blockSignals(True); w.setValue(d[key]); w.blockSignals(False)
        self.d_comp.blockSignals(True); self.d_comp.clear(); self.d_comp.addItem("total")
        if self.session.chg is not None and "diff" in self.session.chg.data:
            self.d_comp.addItem("diff")
        self.d_comp.blockSignals(False)
        has_proj = bool(r["electronico"].params.get("has_projections"))
        self.e_orb.setEnabled(has_proj)
        if not has_proj:
            self.e_orb.setCurrentIndex(0)

    def _fill_series(self):
        """Selectores de serie y de cálculo (solo si el proyecto tiene series de cálculos)."""
        has = bool(self.session.series)
        for w in (self.s_view, self.s_series, self.s_member):
            w.setEnabled(has)
        self.s_series.blockSignals(True); self.s_series.clear()
        self.s_series.addItems([ser.name for ser, _ in self.session.series]); self.s_series.blockSignals(False)
        self.s_view.blockSignals(True); self.s_view.setCurrentIndex(0); self.s_view.blockSignals(False)
        self._fill_members()

    def _fill_members(self):
        """Cálculos de la serie elegida; el primero es la estructura principal de la sesión."""
        self.s_member.blockSignals(True); self.s_member.clear()
        self.s_member.addItem("Principal (sesión)")
        i = self.s_series.currentIndex()
        if self.session is not None and 0 <= i < len(self.session.series):
            ser, rows = self.session.series[i]
            for k, r in enumerate(rows):
                self.s_member.addItem(r["name"] + (" (por defecto)" if k == ser.default else ""))
        self.s_member.blockSignals(False)
        self._mark(0)

    def _summary_html(self) -> str:
        parts = []
        for m in self.report.modules.values():
            parts.append(f"<h4 style='margin:6px 0 2px'>{m.name}</h4><ul style='margin:0'>")
            parts += [f"<li>{h}</li>" for h in m.highlights]
            parts += [f"<li style='color:#b45309'><b>Advertencia:</b> {w}</li>" for w in m.warnings]
            parts += [f"<li style='color:#1d4ed8'>Sugerencia: {s}</li>" for s in m.suggestions]
            parts.append("</ul>")
        return "".join(parts)

    # ------------------------------------------------------------------ dibujo
    def _mark(self, i):
        self._dirty.add(i)
        if self.tabs.currentIndex() == i:
            self._redraw_current()

    def _redraw_current(self, *_):
        i = self.tabs.currentIndex()
        if self.session is None or i < 0 or i not in self._dirty or not self.tabs.isTabEnabled(i):
            return
        self._draw_tab(i)

    def _draw_tab(self, i):
        panel = [self.p_struct, self.p_ph, self.p_el, self.p_den][i]
        try:
            [self._draw_structure, self._draw_phonons, self._draw_electronic, self._draw_density][i](panel)
            panel.canvas.draw_idle()
            self._dirty.discard(i)
        except Exception as exc:
            panel.show_error(f"No se pudo dibujar: {exc}")
            self.statusBar().showMessage(f"Error al dibujar {TABS[i][1]}: {exc}")

    # ------------------------------------------------------------------ exportar
    def _project_name(self):
        return str(self.root) if self.root else None

    def _ask_path(self, title, filters, default):
        path, chosen = QFileDialog.getSaveFileName(self, title, str(Path(self.root or ".") / default), filters)
        return path, chosen

    def export_figure(self, path: str | None = None):
        """Guarda la figura de la pestaña activa (PNG, SVG o PDF)."""
        from ..export import save_figure
        i = self.tabs.currentIndex()
        if self.session is None or i < 0 or not self.tabs.isTabEnabled(i):
            self.statusBar().showMessage("No hay ninguna figura para exportar."); return None
        if i in self._dirty:
            self._draw_tab(i)
        panel = [self.p_struct, self.p_ph, self.p_el, self.p_den][i]
        if path is None:
            path, _ = self._ask_path("Exportar figura", "PNG (*.png);;SVG (*.svg);;PDF (*.pdf)", f"{TABS[i][0]}.png")
        if not path:
            return None
        try:
            out = save_figure(panel.figure, path, style=self.export_style)
        except Exception as exc:
            self.statusBar().showMessage(f"No se pudo exportar: {exc}"); return None
        self.statusBar().showMessage(f"Figura guardada en {out}")
        return out

    def export_report(self, path: str | None = None):
        """Guarda el reporte de texto (resaltados, advertencias y sugerencias) en JSON o HTML."""
        from ..export import save_report
        if self.report is None:
            self.statusBar().showMessage("No hay reporte para exportar."); return None
        if path is None:
            path, _ = self._ask_path("Exportar reporte", "JSON (*.json);;HTML (*.html)", "reporte.json")
        if not path:
            return None
        try:
            out = save_report(self.report, path, self._project_name())
        except Exception as exc:
            self.statusBar().showMessage(f"No se pudo exportar: {exc}"); return None
        self.statusBar().showMessage(f"Reporte guardado en {out}")
        return out

    def export_full_report(self, path: str | None = None):
        """HTML con el reporte y las figuras de todas las pestañas disponibles."""
        from ..export import save_report, figure_png_bytes
        if self.report is None:
            self.statusBar().showMessage("No hay reporte para exportar."); return None
        if path is None:
            path, _ = self._ask_path("Exportar informe", "HTML (*.html)", "informe.html")
        if not path:
            return None
        figures = {}
        panels = [self.p_struct, self.p_ph, self.p_el, self.p_den]
        for i, (_, title) in enumerate(TABS):
            if not self.tabs.isTabEnabled(i):
                continue
            try:
                if i in self._dirty:
                    self._draw_tab(i)
                figures[title] = figure_png_bytes(panels[i].figure)
            except Exception:
                continue
        try:
            out = save_report(self.report, path, self._project_name(), figures)
        except Exception as exc:
            self.statusBar().showMessage(f"No se pudo exportar: {exc}"); return None
        self.statusBar().showMessage(f"Informe guardado en {out}")
        return out

    def _draw_structure(self, panel):
        from ..structural import plot_structure, plot_series
        panel.figure.clear()
        k = self.s_series.currentIndex()
        if self.s_view.currentIndex() == 1 and 0 <= k < len(self.session.series):
            ser, rows = self.session.series[k]
            top, bottom = panel.figure.add_gridspec(2, 1, height_ratios=[2.2, 1.6])
            plot_series(ser, rows, ax=panel.figure.add_subplot(top), table_ax=panel.figure.add_subplot(bottom))
            return
        structure = self.session.structure
        m = self.s_member.currentIndex() - 1                         # 0 = principal
        if m >= 0 and 0 <= k < len(self.session.series):
            row = self.session.series[k][1][m]
            if row["structure"] is None:
                raise ValueError(f"No se pudo leer la estructura de «{row['name']}»: {row['error']}")
            structure = row["structure"]
        ax = panel.figure.add_subplot(111, projection="3d")
        view = ("auto", "top", "side", "iso")[self.s_angle.currentIndex()]
        plot_structure(structure, tuple(w.value() for w in self.sc), self.bonds.isChecked(), self.tol.value(), ax=ax,
                       view=view)

    def _draw_phonons(self, panel):
        from ..phonons import compute_bands, compute_dos, check_stability, plot_bands, plot_dos
        from ..phonons import apply_asr
        ph, unit = self.session.phonon, self.ph_unit.currentText()
        two_d = (None, True, False)[self.ph_path.currentIndex()]
        if self.ph_asr.isChecked():
            if getattr(self, "_ph_asr", (None, None))[0] is not ph:          # una copia por cálculo cargado
                self._ph_asr = (ph, apply_asr(ph))
            ph = self._ph_asr[1]
        bands = compute_bands(ph, self.ph_npts.value(), two_d=two_d)
        rep = check_stability(ph, bands=bands, two_d=two_d)
        panel.figure.clear()
        if self.ph_dos.isChecked():
            a1, a2 = panel.figure.subplots(1, 2, sharey=True, gridspec_kw={"width_ratios": [3, 1]})
            plot_bands(bands, a1, unit); plot_dos(compute_dos(ph, self.ph_mesh.value(), two_d=two_d), a2, unit, horizontal=True)
            a2.set_ylabel(""); a2.xaxis.set_major_locator(MaxNLocator(nbins=3))   # eje estrecho: pocas marcas
        else:
            plot_bands(bands, panel.figure.add_subplot(111), unit)
        panel.figure.suptitle(("Estable" if rep.stable else "Inestabilidad dinámica") + f" (f_min = {(0.0 if abs(rep.min_frequency) < 5e-4 else rep.min_frequency):.3f} THz)"
                              + (" · ASR impuesta" if self.ph_asr.isChecked() else ""),
                              fontsize=10)

    def _draw_electronic(self, panel):
        from ..electronic import plot_bands, plot_dos, orbital_weights
        s = self.session
        emin, emax = self.e_min.value(), self.e_max.value()
        panel.figure.clear()
        show_dos = self.e_dos.isChecked() and s.tdos is not None
        if s.bands is not None and show_dos:
            a1, a2 = panel.figure.subplots(1, 2, sharey=True, gridspec_kw={"width_ratios": [3, 1]})
        elif s.bands is not None:
            a1, a2 = panel.figure.add_subplot(111), None
        else:
            a1, a2 = None, panel.figure.add_subplot(111)
        if a1 is not None:
            orb = self.e_orb.currentText()
            w = None
            if orb != "Ninguno" and s.bands is not None:
                w = orbital_weights(s.bands, tuple(o.strip() for o in orb.split("+")))
            plot_bands(s.bands, a1, emin, emax, weights=w, weight_label=f"peso {orb}")
        if a2 is not None:
            plot_dos(s.tdos, a2, emin, emax, horizontal=a1 is not None)
            if a1 is not None:
                a2.set_ylabel(""); a2.xaxis.set_major_locator(MaxNLocator(nbins=3))   # eje estrecho: pocas marcas

    def _draw_density(self, panel):
        from ..electronic import plot_density_slice, plot_planar_average, plot_isosurface
        chg, view, axis = self.session.chg, self.d_view.currentIndex(), self.d_axis.currentIndex()
        comp = self.d_comp.currentText()
        panel.figure.clear()
        if view == 0:
            plot_density_slice(chg, axis, self.d_frac.value(), (2, 2), comp, panel.figure.add_subplot(111),
                               self.d_cmap.currentText())
        elif view == 1:
            plot_planar_average_ax = panel.figure.add_subplot(111)
            plot_planar_average(chg, axis, comp, plot_planar_average_ax)
        else:
            plot_isosurface(chg, self.d_level.value(), comp, panel.figure.add_subplot(111, projection="3d"))


def run(folder: str | None = None):
    """Abre la aplicación. Uso: ``electronanopy [CARPETA]``; ``electronanopy --resumen CARPETA`` imprime el reporte
    en la terminal sin abrir ventanas."""
    args = sys.argv[1:]
    if args and args[0] == "--resumen":
        if len(args) < 2:
            sys.exit("Uso: electronanopy --resumen CARPETA")
        from ..export import summarize_folder
        print(summarize_folder(args[1]))
        return
    if folder is None and args:
        folder = args[0]
    app = QApplication.instance() or QApplication(sys.argv)
    w = MainWindow(); w.show()
    if folder:
        w.load_folder(folder)
    sys.exit(app.exec())
