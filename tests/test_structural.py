import math
import pytest
from electronanopy.structural import (load_structure, lattice_parameters, space_group,
    vacuum_thickness, is_two_dimensional, bond_table, coordination_numbers, bond_angles, summary)
from electronanopy.structural.io import StructureReadError


def test_lectura_y_red(graphene):
    p = lattice_parameters(load_structure(graphene))
    assert p["a"] == pytest.approx(2.46, abs=1e-3) and p["b"] == pytest.approx(2.46, abs=1e-3)
    assert p["gamma"] == pytest.approx(120.0, abs=1e-3)
    assert p["volume"] == pytest.approx(math.sqrt(3) / 2 * 2.46**2 * 15.0, rel=1e-4)


def test_grupo_espacial(graphene):
    assert space_group(load_structure(graphene))[1] == 191  # P6/mmm


def test_vacio_y_dimensionalidad(graphene):
    s = load_structure(graphene)
    assert vacuum_thickness(s) == pytest.approx(15.0, abs=1e-3)
    assert is_two_dimensional(s)


def test_enlaces_coordinacion_angulos(graphene):
    s = load_structure(graphene)
    assert coordination_numbers(s) == [3, 3]
    distancias = [b["distance"] for b in bond_table(s)]
    assert len(distancias) == 6                                     # 3 enlaces por átomo x 2 átomos
    assert all(d == pytest.approx(2.46 / math.sqrt(3), abs=1e-3) for d in distancias)  # 1.420 Å
    assert all(a["angle"] == pytest.approx(120.0, abs=1e-2) for a in bond_angles(s))


def test_resumen(graphene):
    r = summary(load_structure(graphene))
    assert r["formula"] == "C" and r["n_atoms"] == 2 and r["two_dimensional"]


@pytest.mark.parametrize("contenido", ["", "esto no es un POSCAR\n1\n"])
def test_archivo_invalido(tmp_path, contenido):
    f = tmp_path / "POSCAR"; f.write_text(contenido)
    with pytest.raises(StructureReadError):
        load_structure(f)


def test_archivo_inexistente(tmp_path):
    with pytest.raises(StructureReadError):
        load_structure(tmp_path / "no_existe")


def test_vista_automatica_superior_en_2d_e_isometrica_en_3d(graphene):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pymatgen.core import Lattice, Structure
    from electronanopy.structural import load_structure, plot_structure
    fig = plot_structure(load_structure(graphene))
    assert fig.axes[0].elev == 90
    fig = plot_structure(load_structure(graphene), view="side")
    assert fig.axes[0].elev == 0
    nacl = Structure(Lattice.cubic(5.6), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    assert plot_structure(nacl).axes[0].elev != 90
    plt.close("all")


def test_vista_3d_marca_los_vecinos_fuera_de_la_celda():
    import matplotlib; matplotlib.use("Agg")
    from pymatgen.core import Lattice, Structure
    from electronanopy.structural import plot_structure
    st = Structure(Lattice.cubic(1.45), ["C"], [[0, 0, 0]])                  # red cúbica simple, enlaces C–C periódicos
    fig = plot_structure(st)
    ax = fig.axes[0]
    n_scatter = len(ax.collections)
    assert n_scatter >= 2                                   # átomos + imágenes fantasma de los vecinos periódicos
