"""Módulo de fonones (OE2): lectura de Phonopy, bandas, DOS y estabilidad dinámica."""
from .io import load_phonopy, PhononReadError
from .analysis import (compute_bands, compute_dos, check_stability, phonon_quality, apply_asr, BandData, DosData,
                       StabilityReport, QualityReport)
from .plotting import plot_bands, plot_dos, plot_bands_dos

__all__ = ["load_phonopy", "PhononReadError", "compute_bands", "compute_dos", "check_stability", "phonon_quality", "apply_asr", "QualityReport",
           "BandData", "DosData", "StabilityReport", "plot_bands", "plot_dos", "plot_bands_dos"]
