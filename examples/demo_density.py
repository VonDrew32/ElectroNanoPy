"""Demo con densidad de carga SINTÉTICA de grafeno (gaussianas en los átomos): corte 2D, perfil y isosuperficie.

Con datos reales: chg = load_chgcar("02_scf/CHGCAR")  (calcule el SCF con LCHARG = .TRUE.)
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from electronanopy.electronic import (density, total_electrons, plot_density_slice, plot_planar_average,
                                      plot_isosurface)
from electronanopy.electronic import _synthetic as S

out = sys.argv[1] if len(sys.argv) > 1 else "densidad_demo.png"
chg = S.graphene_chgcar()                                          # sustituir por load_chgcar(...)
print(f"Electrones en la celda: {total_electrons(chg):.3f}")
fig = plt.figure(figsize=(12, 4.2))
plot_density_slice(chg, axis=2, frac=0.5, ax=fig.add_subplot(131)); fig.axes[0].set_title("Corte en el plano (z = c/2)")
plot_planar_average(chg, ax=fig.add_subplot(132)); fig.axes[-1].set_title("Promedio planar")
plot_isosurface(chg, 0.4 * density(chg).max(), ax=fig.add_subplot(133, projection="3d")); fig.axes[-1].set_title("Isosuperficie")
fig.tight_layout(); fig.savefig(out, dpi=120)
print("Figura:", out)
