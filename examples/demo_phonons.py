"""Demo con un cálculo SINTÉTICO (red triangular 2D de resortes): bandas, DOS y estabilidad.

Con datos reales: reemplace `ph = ...` por `ph = load_phonopy("carpeta_con_phonopy_params")`.
"""
import sys
import matplotlib
matplotlib.use("Agg")
from electronanopy.phonons import compute_bands, compute_dos, check_stability, plot_bands_dos
from electronanopy.phonons import _synthetic as S

out = sys.argv[1] if len(sys.argv) > 1 else "fonones_demo.png"
ph = S.triangular_2d()                       # sustituir por load_phonopy(...) con datos reales
bands, dos = compute_bands(ph), compute_dos(ph)
rep = check_stability(ph, bands=bands)
print(rep.message)
plot_bands_dos(bands, dos, rep).savefig(out, dpi=130, bbox_inches="tight")
print("Figura:", out)
