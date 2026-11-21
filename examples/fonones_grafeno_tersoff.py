"""Fonones de grafeno con el Tersoff optimizado de Lindsay y Broido (PRB 81, 205441) vía LAMMPS + Phonopy.
Sirve para verificar la cadena de OE2 (fuerzas -> constantes de fuerza -> bandas -> estabilidad) sin VASP."""
import sys, os, numpy as np
sys.path.insert(0, "/tmp/claude-0/-home-claude/177a799e-8a3f-5121-9f3b-2f024a39804f/scratchpad/venv/lib/python3.13/site-packages")
import lammps
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
POT = os.path.join(os.path.dirname(lammps.__file__), "share/lammps/potentials/BNC.tersoff")


def lmp_run(cell, frac, mode="energy"):
    """cell: (3,3) filas = vectores de red (triclinica con a1 || x). Devuelve energia (eV) y fuerzas (eV/A)."""
    cell = np.asarray(cell, float); frac = np.asarray(frac, float)
    # rotar para que a1 || x y a2 en el plano xy (condicion de LAMMPS)
    a1, a2, a3 = cell
    ex = a1 / np.linalg.norm(a1); ez = np.cross(a1, a2); ez /= np.linalg.norm(ez); ey = np.cross(ez, ex)
    R = np.array([ex, ey, ez])
    c = cell @ R.T
    pos = frac @ cell @ R.T
    L = lammps.lammps(cmdargs=["-log", "none", "-screen", "none"])
    L.commands_list([
        "units metal", "boundary p p p", "atom_style atomic", "atom_modify map array",
        f"region box prism 0 {c[0,0]:.12f} 0 {c[1,1]:.12f} 0 {c[2,2]:.12f} {c[1,0]:.12f} {c[2,0]:.12f} {c[2,1]:.12f} units box",
        "create_box 1 box", "mass 1 12.011"])
    for p in pos:
        L.command(f"create_atoms 1 single {p[0]:.12f} {p[1]:.12f} {p[2]:.12f} units box")
    L.commands_list(["pair_style tersoff", f"pair_coeff * * {POT} C", "thermo_style custom pe", "run 0"])
    e = L.get_thermo("pe")
    f = np.array(L.gather_atoms("f", 1, 3)).reshape(-1, 3) @ R          # de vuelta al marco original
    L.close()
    return e, f


def cell_for(a, c=15.0):
    return np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])

FRAC = [[1/3, 2/3, .5], [2/3, 1/3, .5]]

if __name__ == "__main__":
    out = sys.argv[1]
    from scipy.optimize import minimize_scalar
    a0 = minimize_scalar(lambda a: lmp_run(cell_for(a), FRAC)[0], bounds=(2.3, 2.6), method="bounded", options={"xatol": 1e-7}).x
    print("a0 empirico =", a0)
    ph = Phonopy(PhonopyAtoms(["C", "C"], cell=cell_for(a0), scaled_positions=FRAC, masses=[12.011] * 2),
                 supercell_matrix=[5, 5, 1], primitive_matrix="auto")
    ph.generate_displacements(distance=0.01)
    forces = []
    for sc in ph.supercells_with_displacements:
        _, f = lmp_run(sc.cell, sc.scaled_positions)
        forces.append(f)
    ph.forces = forces
    ph.produce_force_constants()
    os.makedirs(out, exist_ok=True)
    ph.save(os.path.join(out, "phonopy_params.yaml"), settings={"force_constants": False})
    print("guardado en", out)
