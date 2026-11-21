"""Dispersión de fonones del grafeno (carpeta con phonopy_disp.yaml + FORCE_SETS o phonopy_params.yaml) frente a la literatura.

Uso:  python examples/comparar_fonones_grafeno.py CARPETA [salida.png]

Referencia: J. Maultzsch, S. Reich, C. Thomsen, H. Requardt, P. Ordejón, Phys. Rev. Lett. 92, 075501 (2004), Tabla 1
(DFT-GGA y experimento; los cálculos son de grafito, cuyas ramas en el plano son las del grafeno).
Imprime una tabla Γ, M, K y guarda la figura (bandas + DOS, en cm⁻¹).
"""
import sys
import warnings
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from electronanopy.phonons import load_phonopy, compute_bands, compute_dos, check_stability

THZ = 33.35641
# cm^-1, orden ascendente en cada punto
GGA = {"Γ": [0, 0, 0, 884, 1569, 1569], "M": [476, 634, 640, 1338, 1346, 1396], "K": [539, 539, 1004, 1221, 1221, 1289]}
EXP = {"Γ": {3: 868, 4: 1582, 5: 1587}, "M": {4: 1323, 5: 1390}, "K": {4: 1194, 5: 1265}}   # IXS (M, K); Raman/IR (Γ)

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    carpeta = sys.argv[1]
    salida = sys.argv[2] if len(sys.argv) > 2 else "fonones_grafeno_vs_literatura.png"
    ph = load_phonopy(carpeta)
    ph.run_qpoints([[0, 0, 0], [.5, 0, 0], [1 / 3, 1 / 3, 0]])
    f = ph.get_qpoints_dict()["frequencies"] * THZ
    print(f"{'':3s}{'calculado':>10s}{'lit. GGA':>10s}{'dif.':>8s}{'dif. %':>8s}")
    for row, (q, lit) in zip(f, GGA.items()):
        for a, b in zip(row, lit):
            print(f"{q:3s}{a:10.1f}{b:10.0f}{a - b:+8.1f}" + (f"{100 * (a - b) / b:+8.1f}" if b else "       –"))
    bands, dos = compute_bands(ph), compute_dos(ph, mesh=48)
    rep = check_stability(ph, bands=bands)
    fig, (a1, a2) = plt.subplots(1, 2, sharey=True, figsize=(8, 4.6), constrained_layout=True,
                                 gridspec_kw={"width_ratios": [3, 1]})
    a1.plot(bands.distance, bands.frequencies * THZ, color="C0", lw=1.3)
    xs = dict(zip("ΓMK", [bands.ticks[0], bands.ticks[1], bands.ticks[2]]))
    for q, lit in GGA.items():
        a1.plot([xs[q]] * 6, lit, "s", mfc="none", mec="C3", ms=6, label="DFT-GGA (Maultzsch 2004)" if q == "Γ" else None)
    for q, d in EXP.items():
        a1.plot([xs[q]] * len(d), list(d.values()), "o", color="k", ms=4, label="Experimento" if q == "Γ" else None)
    for t in bands.ticks:
        a1.axvline(t, color="0.8", lw=0.6)
    a1.set_xticks(bands.ticks); a1.set_xticklabels(bands.tick_labels)
    a1.set_xlim(bands.distance[0], bands.distance[-1]); a1.set_ylabel("Frecuencia (cm$^{-1}$)")
    a1.legend(loc="upper right", fontsize=8, frameon=False)
    a2.fill_betweenx(dos.frequencies * THZ, 0, dos.dos / THZ, color="C1", alpha=0.5); a2.set_xlabel("DOS")
    fig.suptitle(f"Grafeno: {'estable' if rep.stable else 'INESTABLE'} (f_min = {(0.0 if abs(rep.min_frequency) < 5e-4 else rep.min_frequency):.3f} THz)", fontsize=10)
    fig.savefig(salida, dpi=180)
    print("figura:", salida)
