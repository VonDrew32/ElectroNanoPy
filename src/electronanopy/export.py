"""Exportar resultados (OE4): figuras con calidad de publicación y el reporte en JSON o HTML. No depende de Qt."""
from __future__ import annotations
import base64
from contextlib import contextmanager
import html
import io
import json
from datetime import datetime
from pathlib import Path
import numpy as np

FIGURE_FORMATS = ("png", "svg", "pdf")

# Estilos de exportación: ancho final (cm), tamaño de letra base (pt; el texto de la figura se escala en proporción:
# un rótulo de 10 pt pasa a ``font_pt``), factor de grosor de línea y resolución de los PNG.
STYLES = {
    "pantalla": {"label": "Pantalla (tamaño actual, 300 dpi)", "width_cm": None, "font_pt": None, "lw": 1.0, "dpi": 300},
    "publicacion_1col": {"label": "Publicación, 1 columna (8.5 cm, letra de 8 pt, 600 dpi)", "width_cm": 8.5,
                         "font_pt": 8.0, "lw": 1.0, "dpi": 600},
    "publicacion_2col": {"label": "Publicación, 2 columnas (17 cm, letra de 9 pt, 600 dpi)", "width_cm": 17.0,
                         "font_pt": 9.0, "lw": 1.0, "dpi": 600},
    "poster": {"label": "Póster (35 cm de ancho, letra de 24 pt, líneas gruesas, 300 dpi)", "width_cm": 35.0,
               "font_pt": 24.0, "lw": 2.2, "dpi": 300},
}


def _plain(obj):
    """Convierte lo que no es JSON (numpy, Path, objetos) en tipos simples."""
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _plain(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    return str(obj)


def report_to_dict(report, project: str | None = None) -> dict:
    from . import __version__
    return {"proyecto": project, "version": __version__, "fecha": datetime.now().isoformat(timespec="seconds"),
            "modulos": {key: {"nombre": m.name, "disponible": m.available, "resumen": _plain(m.summary),
                              "resaltados": list(m.highlights), "advertencias": list(m.warnings),
                              "sugerencias": list(m.suggestions)} for key, m in report.modules.items()}}


def report_to_json(report, project: str | None = None) -> str:
    return json.dumps(report_to_dict(report, project), ensure_ascii=False, indent=2)


def report_to_html(report, project: str | None = None, figures: dict | None = None) -> str:
    """HTML autocontenido. ``figures``: {título: bytes PNG}, incrustadas como imágenes."""
    d = report_to_dict(report, project)
    e = html.escape
    out = ["<!doctype html><html lang='es'><head><meta charset='utf-8'><title>Reporte ElectroNanoPy</title>",
           "<style>body{font-family:system-ui,sans-serif;max-width:900px;margin:2em auto;padding:0 1em;color:#1f2937}"
           "h1{font-size:1.5em}h2{border-bottom:1px solid #d1d5db;padding-bottom:.2em;margin-top:1.6em}"
           ".w{color:#b45309}.s{color:#1d4ed8}.m{color:#6b7280;font-size:.85em}img{max-width:100%;margin:.6em 0}</style>"
           "</head><body>", "<h1>Reporte de ElectroNanoPy</h1>",
           f"<p class='m'>Proyecto: {e(str(d['proyecto'] or '—'))} · versión {e(d['version'])} · {e(d['fecha'])}</p>"]
    for m in d["modulos"].values():
        out.append(f"<h2>{e(m['nombre'])}</h2>")
        if not m["disponible"] and not m["advertencias"]:
            out.append("<p class='m'>Sin datos.</p>")
        out.append("<ul>")
        out += [f"<li>{e(t)}</li>" for t in m["resaltados"]]
        out += [f"<li class='w'><b>Advertencia:</b> {e(t)}</li>" for t in m["advertencias"]]
        out += [f"<li class='s'>Sugerencia: {e(t)}</li>" for t in m["sugerencias"]]
        out.append("</ul>")
    for title, png in (figures or {}).items():
        b64 = base64.b64encode(png).decode()
        out.append(f"<h2>{e(title)}</h2><img alt='{e(title)}' src='data:image/png;base64,{b64}'>")
    out.append("</body></html>")
    return "".join(out)


def report_to_text(report, project: str | None = None) -> str:
    """Reporte en texto plano (para la terminal)."""
    d = report_to_dict(report, project)
    out = [f"ElectroNanoPy {d['version']} · proyecto: {d['proyecto'] or '—'}"]
    for m in d["modulos"].values():
        out.append(f"\n== {m['nombre']} ==")
        if not m["disponible"] and not m["advertencias"]:
            out.append("  (sin datos)")
        out += [f"  · {t}" for t in m["resaltados"]]
        out += [f"  ! {t}" for t in m["advertencias"]]
        out += [f"  → {t}" for t in m["sugerencias"]]
    return "\n".join(out)


def summarize_folder(folder: str | Path) -> str:
    """Examina una carpeta de proyecto y devuelve el reporte en texto, sin interfaz gráfica."""
    from .project import scan_folder, load_session
    from .report import analyze
    session = load_session(scan_folder(folder))
    return report_to_text(analyze(session), str(folder))


def save_report(report, path: str | Path, project: str | None = None, figures: dict | None = None) -> Path:
    """Guarda el reporte según la extensión: .json o .html."""
    p = Path(path)
    ext = p.suffix.lower()
    if ext == ".json":
        p.write_text(report_to_json(report, project), encoding="utf-8")
    elif ext in (".html", ".htm"):
        p.write_text(report_to_html(report, project, figures), encoding="utf-8")
    else:
        raise ValueError(f"Formato no admitido «{ext}»; use .json o .html.")
    return p


def figure_png_bytes(fig, dpi: int = 150) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    return buf.getvalue()


@contextmanager
def styled(fig, style: str = "pantalla"):
    """Aplica temporalmente un estilo de ``STYLES`` a la figura (tamaño, letra y grosor de línea) y la restaura al salir."""
    if style not in STYLES:
        raise ValueError(f"Estilo desconocido «{style}». Opciones: {', '.join(STYLES)}.")
    st = STYLES[style]
    from matplotlib.lines import Line2D
    from matplotlib.text import Text
    size0 = fig.get_size_inches().copy()
    texts = [(t, t.get_fontsize()) for t in fig.findobj(Text)]
    lines = [(ln, ln.get_linewidth()) for ln in fig.findobj(Line2D)]
    try:
        if st["width_cm"]:
            w = st["width_cm"] / 2.54
            fig.set_size_inches(w, size0[1] * w / size0[0])
        if st["font_pt"]:
            for t, fs in texts:
                t.set_fontsize(fs * st["font_pt"] / 10.0)
        if st["lw"] != 1.0:
            for ln, lw in lines:
                ln.set_linewidth(lw * st["lw"])
        yield st
    finally:
        fig.set_size_inches(size0)
        for t, fs in texts:
            t.set_fontsize(fs)
        for ln, lw in lines:
            ln.set_linewidth(lw)


def save_figure(fig, path: str | Path, dpi: int | None = None, style: str = "pantalla") -> Path:
    """Guarda una figura de matplotlib en PNG, SVG o PDF (según la extensión) con el estilo indicado.

    ``dpi`` solo afecta al PNG; si no se da, se usa el del estilo. La figura queda como estaba al terminar.
    """
    p = Path(path)
    ext = p.suffix.lower().lstrip(".")
    if ext not in FIGURE_FORMATS:
        raise ValueError(f"Formato no admitido «.{ext}»; use {', '.join('.' + f for f in FIGURE_FORMATS)}.")
    with styled(fig, style) as st:
        fig.savefig(p, format=ext, dpi=dpi or st["dpi"], bbox_inches="tight")
    return p
