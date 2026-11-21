import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")      # pruebas de la interfaz sin pantalla
import pytest
from pymatgen.core import Lattice, Structure


@pytest.fixture
def graphene(tmp_path):
    """Grafeno de una celda (2 átomos de C) con 15 Å de vacío en c; escrito como POSCAR."""
    a = 2.46
    s = Structure(Lattice.hexagonal(a, 15.0), ["C", "C"], [[1/3, 2/3, 0.5], [2/3, 1/3, 0.5]])
    path = tmp_path / "POSCAR"
    s.to(filename=str(path), fmt="poscar")
    return path
