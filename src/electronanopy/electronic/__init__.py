"""Módulo electrónico (OE3): bandas, DOS, brecha, punto de Dirac, magnetismo y densidad de carga."""
from .io import (load_vasprun, load_band_structure, band_structure_from_vasprun, read_scf_info, load_dos,
                 mesh_band_structure, ElectronicReadError)
from .analysis import (band_gap, dirac_points, fermi_velocity, dos_at_fermi, dos_energy_step, classify, orbital_weights,
                       read_magnetization, GapReport)
from .plotting import plot_bands, plot_dos, plot_bands_dos
from .density import (load_chgcar, DensityReadError, density, total_electrons, net_magnetization, planar_average,
                      density_at, slice_2d, isosurface)
from .density_plotting import plot_density_slice, plot_planar_average, plot_isosurface

__all__ = ["load_vasprun", "load_band_structure", "band_structure_from_vasprun", "read_scf_info", "load_dos", "mesh_band_structure", "ElectronicReadError", "band_gap",
           "dirac_points", "fermi_velocity", "dos_at_fermi", "dos_energy_step", "classify", "orbital_weights",
           "read_magnetization", "GapReport", "plot_bands", "plot_dos", "plot_bands_dos",
           "load_chgcar", "DensityReadError", "density", "total_electrons", "net_magnetization", "planar_average",
           "density_at", "slice_2d", "isosurface", "plot_density_slice", "plot_planar_average", "plot_isosurface"]
