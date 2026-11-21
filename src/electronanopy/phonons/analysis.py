"""Bandas de fonones, densidad de estados vibracionales y criterio de estabilidad dinámica."""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from phonopy import Phonopy
from .paths import path_2d, path_3d, interpolate, reciprocal_lattice, Segment

THZ_TO_CM1 = 33.35641   # 1 THz = 33.35641 cm^-1


def _vacuum_two_d(ph: Phonopy, min_vacuum: float = 8.0) -> bool:
    """¿La celda tiene un vacío >= min_vacuum Å a lo largo de c? (criterio de material 2D)."""
    cell = np.asarray(ph.primitive.cell)
    z = np.sort(np.asarray(ph.primitive.scaled_positions)[:, 2] % 1.0)
    gaps = np.diff(np.concatenate([z, [z[0] + 1.0]]))
    height = abs(np.linalg.det(cell)) / np.linalg.norm(np.cross(cell[0], cell[1]))
    return float(gaps.max() * height) >= min_vacuum


def _z(x: float) -> float:
    """0.0 en lugar de -0.000 al imprimir (ruido numérico de los acústicos en Γ)."""
    return 0.0 if abs(x) < 5e-4 else float(x)


@dataclass
class BandData:
    """Dispersión de fonones. Frecuencias en THz; las imaginarias se representan como negativas."""
    distance: np.ndarray          # (nq_total,) distancia acumulada en el espacio recíproco (1/Å)
    frequencies: np.ndarray       # (nq_total, nbands)
    qpoints: np.ndarray           # (nq_total, 3) fraccionarias
    ticks: list[float]            # posiciones de las etiquetas
    tick_labels: list[str]
    two_d: bool


@dataclass
class DosData:
    frequencies: np.ndarray       # THz
    dos: np.ndarray               # estados / THz


@dataclass
class StabilityReport:
    stable: bool
    min_frequency: float                       # THz, sobre todos los puntos y ramas
    min_frequency_considered: float            # THz, sin los modos acústicos cerca de Γ
    n_ignored_near_gamma: int                  # puntos acústicos cerca de Γ con f < -tol
    imaginary: list[dict] = field(default_factory=list)   # puntos que violan el criterio
    tol: float = 0.05
    gamma_radius: float = 0.1
    message: str = ""


def compute_bands(ph: Phonopy, npoints: int = 101, path: list[Segment] | None = None,
                  connected: list[bool] | None = None, two_d: bool | None = None) -> BandData:
    """Calcula la dispersión a lo largo de un camino de alta simetría.

    Si ``path`` es None se elige automáticamente: camino 2D (kz = 0) para materiales con vacío,
    o seekpath para 3D. ``two_d=None`` detecta el vacío en la celda.
    """
    cell = np.asarray(ph.primitive.cell)
    if two_d is None:
        two_d = _vacuum_two_d(ph)
    if path is None:
        if two_d:
            path, connected = path_2d(cell)
        else:
            path, connected = path_3d(cell, np.asarray(ph.primitive.numbers),
                                      np.asarray(ph.primitive.scaled_positions))
    if connected is None:
        connected = [False] * len(path)
    qsets = interpolate(path, npoints)
    ph.run_band_structure(qsets, path_connections=connected)
    freqs = ph.band_structure.frequencies
    rec = reciprocal_lattice(cell)
    xs, qs, fs, ticks, labels = [], [], [], [], []
    x0 = 0.0
    for i, ((la, lb, _, _), q, f) in enumerate(zip(path, qsets, freqs)):
        d = x0 + np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(q @ rec, axis=0), axis=1))])
        xs.append(d); qs.append(q); fs.append(np.asarray(f))
        if i == 0:
            ticks.append(d[0]); labels.append(la)
        # etiqueta del extremo final: si el siguiente segmento continúa, se comparte; si no, 'A|B'
        nxt = path[i + 1][0] if i + 1 < len(path) else None
        if connected[i] or nxt is None:
            ticks.append(d[-1]); labels.append(lb)
        else:
            ticks.append(d[-1]); labels.append(f"{lb}|{nxt}")
            # el salto entre segmentos no suma distancia
        x0 = d[-1]
    # si el segmento no está conectado, la etiqueta combinada ya cubre el salto; se evita duplicar el inicio
    return BandData(np.concatenate(xs), np.vstack(fs), np.vstack(qs), ticks, labels, bool(two_d))


def compute_dos(ph: Phonopy, mesh: int | tuple[int, int, int] = 24, sigma: float | None = None,
                two_d: bool | None = None, npoints: int = 501) -> DosData:
    """Densidad de estados vibracionales en una malla centrada en Γ.

    En 2D la malla es (n, n, 1) y se usa suavizado gaussiano (``sigma`` THz, 0.05 por defecto);
    en 3D se usa el método del tetraedro salvo que se indique ``sigma``.
    """
    if two_d is None:
        two_d = _vacuum_two_d(ph)
    if isinstance(mesh, int):
        mesh = (mesh, mesh, 1) if two_d else (mesh, mesh, mesh)
    ph.run_mesh(list(mesh), is_gamma_center=True)
    if two_d and sigma is None:
        sigma = 0.05
    freqs = ph.mesh.frequencies
    ph.run_total_dos(sigma=sigma, freq_min=float(freqs.min()) - 0.5, freq_max=float(freqs.max()) + 0.5,
                     freq_pitch=(float(freqs.max() - freqs.min()) + 1.0) / (npoints - 1),
                     use_tetrahedron_method=sigma is None)
    return DosData(np.asarray(ph.total_dos.frequency_points), np.asarray(ph.total_dos.dos))


@dataclass
class QualityReport:
    """Diagnóstico de la calidad de las constantes de fuerza (no decide la estabilidad)."""
    gamma_acoustic: tuple[float, float, float]   # THz, los 3 modos más bajos en Γ (deben ser ≈ 0)
    asr_max: float                               # THz, máximo valor absoluto de gamma_acoustic
    min_cell_length: float                       # Å, vector más corto de la supercelda en las direcciones periódicas
    periodic_dims: int                           # 2 si el material es 2D (se ignora c), 3 en otro caso
    n_displacements: int | None = None
    has_nac: bool = False


def phonon_quality(ph: Phonopy, two_d: bool | None = None) -> QualityReport:
    """Mide cuánto se aleja el cálculo de las condiciones ideales.

    * Regla de suma acústica (ASR): en Γ, los 3 modos acústicos deben tener f ≈ 0. Si no, las fuerzas tienen
      ruido (EDIFF, malla k, ENCUT, ISMEAR) y los modos acústicos cerca de Γ no son fiables.
    * Tamaño de la supercelda: longitud del vector más corto en las direcciones periódicas.
    """
    if two_d is None:
        two_d = _vacuum_two_d(ph)
    ph.run_qpoints([[0.0, 0.0, 0.0]], with_eigenvectors=False)
    f = np.sort(np.asarray(ph.get_qpoints_dict()["frequencies"])[0])
    gam = tuple(float(x) for x in f[:3])
    cell = np.asarray(ph.supercell.cell)
    lengths = np.linalg.norm(cell, axis=1)
    dims = 2 if two_d else 3
    ds = getattr(ph, "dataset", None) or {}
    nd = len(ds["first_atoms"]) if "first_atoms" in ds else None
    return QualityReport(gam, float(np.max(np.abs(gam))), float(lengths[:dims].min()), dims, nd,
                         getattr(ph, "nac_params", None) is not None)


def apply_asr(ph: Phonopy) -> Phonopy:
    """Copia de ``ph`` con la regla de suma acústica (y la simetría de permutación) impuesta a las constantes de
    fuerza. No modifica el original. Corrige el ruido numérico; no arregla una supercelda pequeña ni fuerzas mal
    convergidas."""
    import copy
    out = copy.deepcopy(ph)
    out.symmetrize_force_constants()
    return out


def check_stability(ph: Phonopy, mesh: int | tuple[int, int, int] = 24, tol: float = 0.05,
                    gamma_radius: float = 0.1, two_d: bool | None = None,
                    bands: BandData | None = None) -> StabilityReport:
    """Criterio de estabilidad dinámica: ¿hay frecuencias imaginarias?

    Se evalúan todas las ramas en la malla (y en el camino, si se da ``bands``). Una frecuencia es
    imaginaria si f < -``tol`` (THz). Los 3 modos acústicos más bajos se ignoran a menos de
    ``gamma_radius`` (1/Å; 0 = sin exclusión) de Γ, donde el ruido numérico suele dar valores pequeños negativos.
    """
    if two_d is None:
        two_d = _vacuum_two_d(ph)
    if isinstance(mesh, int):
        mesh = (mesh, mesh, 1) if two_d else (mesh, mesh, mesh)
    ph.run_mesh(list(mesh), is_gamma_center=True)
    q = np.asarray(ph.mesh.qpoints); f = np.asarray(ph.mesh.frequencies)
    if bands is not None:
        q = np.vstack([q, bands.qpoints]); f = np.vstack([f, bands.frequencies])
    rec = reciprocal_lattice(np.asarray(ph.primitive.cell))
    qw = q - np.round(q)                                   # imagen más cercana a Γ (coordenadas fraccionarias)
    dist = np.linalg.norm(qw @ rec, axis=1)
    near = dist < gamma_radius
    acoustic = np.zeros_like(f, dtype=bool); acoustic[:, :3] = True
    ignored = acoustic & near[:, None]
    bad = (f < -tol) & ~ignored
    n_ign = int(((f < -tol) & ignored).sum())
    considered = np.where(ignored, np.inf, f)
    imag = [{"q": tuple(float(x) for x in np.round(q[i], 4)), "band": int(j), "frequency": float(f[i, j])}
            for i, j in zip(*np.nonzero(bad))]
    imag.sort(key=lambda d: d["frequency"])
    stable = not imag
    fmin = round(float(f.min()), 4) + 0.0            # evita mostrar -0.000
    if stable:
        msg = ("Dinámicamente estable: no hay frecuencias imaginarias (umbral -" f"{tol} THz)."
               + (f" Se ignoraron {n_ign} valores acústicos negativos cerca de Γ (ruido numérico; "
                  f"mínimo {_z(fmin):.3f} THz)." if n_ign else f" Frecuencia mínima: {_z(fmin):.3f} THz."))
    else:
        msg = (f"Posible inestabilidad dinámica: {len(imag)} modos con frecuencia imaginaria "
               f"(mínimo {imag[0]['frequency']:.3f} THz en q = {imag[0]['q']}).")
    return StabilityReport(stable, fmin, float(considered.min()), n_ign, imag[:50], tol, gamma_radius, msg)
