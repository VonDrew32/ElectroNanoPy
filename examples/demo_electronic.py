"""Demo con grafeno SINTÉTICO (tight-binding): bandas, DOS, brecha, punto de Dirac y velocidad de Fermi.

Con datos reales (cálculo de bandas en modo línea + KPOINTS):
    bs = load_band_structure("bandas/vasprun.xml", efermi="scf/vasprun.xml")
    tdos, cdos = load_dos("dos/vasprun.xml")
"""
import sys
import matplotlib
matplotlib.use("Agg")
from electronanopy.electronic import (band_gap, dirac_points, fermi_velocity, classify, plot_bands_dos)
from electronanopy.electronic import _synthetic as S

out = sys.argv[1] if len(sys.argv) > 1 else "electronico_demo.png"
bs, dos = S.graphene_bands(npts=80), S.graphene_dos()          # sustituir por load_band_structure / load_dos
print(band_gap(bs).message)
dp = dirac_points(bs)
v = fermi_velocity(bs, dp[0])
print(f"Punto de Dirac en {dp[0]['label']}; ħv_F = {v['hbar_vF_eV_A']:.2f} eV·Å ({v['vF_m_s']:.2e} m/s)")
print("Clasificación:", classify(bs, dos))
plot_bands_dos(bs, dos, title="Grafeno (tight-binding sintético)").savefig(out, dpi=130, bbox_inches="tight")
print("Figura:", out)
