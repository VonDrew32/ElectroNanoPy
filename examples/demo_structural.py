"""Demo: crea grafeno, lo guarda como POSCAR, lo lee y muestra resumen y figura."""
import sys, json
import matplotlib
from pymatgen.core import Lattice, Structure
from electronanopy.structural import load_structure, summary, plot_structure

out = sys.argv[1] if len(sys.argv) > 1 else "grafeno_demo.png"
s = Structure(Lattice.hexagonal(2.46, 15.0), ["C", "C"], [[1/3, 2/3, .5], [2/3, 1/3, .5]])
s.to(filename="POSCAR_grafeno", fmt="poscar")
s = load_structure("POSCAR_grafeno")
print(json.dumps(summary(s), indent=2))
matplotlib.use("Agg")
fig = plot_structure(s, supercell=(3, 3, 1)); fig.savefig(out, dpi=130)
print("Figura:", out)
