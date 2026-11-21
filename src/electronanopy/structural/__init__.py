"""Módulo estructural (OE1): lectura, descriptores y visualización de la estructura."""
from .io import load_structure
from .descriptors import (
    lattice_parameters, space_group, vacuum_thickness, is_two_dimensional,
    bond_table, coordination_numbers, bond_angles, summary,
)
from .plotting import plot_structure, plot_series
from .eos import fit_birch_murnaghan, birch_murnaghan
from .series import series_table, series_summary, read_stress

__all__ = ["load_structure", "lattice_parameters", "space_group", "vacuum_thickness",
           "is_two_dimensional", "bond_table", "coordination_numbers", "bond_angles",
           "summary", "plot_structure", "plot_series", "series_table", "series_summary", "read_stress", "fit_birch_murnaghan", "birch_murnaghan"]
