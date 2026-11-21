"""Valida el análisis electrónico (OE3) con un vasprun.xml real de grafeno, SIN cálculo de bandas en línea.

La malla uniforme 18x18x1 Γ-centrada contiene Γ, M y K, pero el vasprun solo guarda los k irreducibles. Se reconstruye
un camino Γ–M–K–Γ de baja densidad buscando, para cada punto del camino, su equivalente por simetría entre los k
calculados, y se pasa por las mismas funciones que usa la aplicación. Requiere un cálculo con LORBIT = 11.

Uso:  python examples/validar_oe3_malla_uniforme.py CARPETA_DEL_CALCULO     (con vasprun.xml)
"""
import sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from electronanopy.electronic import analysis as A
from electronanopy.electronic import _synthetic as S
from pymatgen.io.vasp import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin


def build(path, efermi=None):
    vr = Vasprun(path + "/vasprun.xml", parse_projected_eigen=True, parse_potcar_file=False)
    s = vr.final_structure
    k = np.array(vr.actual_kpoints)
    ev = vr.eigenvalues[Spin.up][:, :, 0]
    proj = np.array(vr.projected_eigenvalues[Spin.up])
    ops = {tuple(np.round(o.rotation_matrix, 6).ravel()) for o in SpacegroupAnalyzer(s, symprec=1e-4).get_symmetry_operations()}
    R = [np.array(o).reshape(3, 3) for o in ops]
    rec_ops = [np.linalg.inv(r).T for r in R] + [-np.linalg.inv(r).T for r in R]

    def find(target):
        for i, kk in enumerate(k):
            for m in rec_ops:
                d = m @ target - kk
                if np.linalg.norm(d - np.round(d)) < 1e-4:
                    return i
        raise LookupError(f"{target} no está entre los k calculados")

    G_, M_, K_ = np.zeros(3), np.array([.5, 0, 0]), np.array([1/3, 1/3, 0])
    seg = [[G_ + (M_ - G_) * t / 9 for t in range(10)], [M_ + (K_ - M_) * t / 3 for t in range(4)],
           [K_ + (G_ - K_) * t / 6 for t in range(7)]]
    pts = [p for ps in seg for p in ps]
    idx = np.array([find(p) for p in pts])
    ef = vr.efermi if efermi is None else efermi
    bs = BandStructureSymmLine(pts, {Spin.up: ev[idx].T}, s.lattice.reciprocal_lattice, ef,
                               {"Γ": G_, "M": M_, "K": K_}, structure=s,
                               projections={Spin.up: proj[idx].transpose(1, 0, 3, 2)})
    return vr, bs, pts


if __name__ == "__main__":
    vr, bs, pts = build(sys.argv[1])
    nelect = float(vr.parameters["NELECT"])
    print(f"camino reconstruido: {len(pts)} puntos k | E_F = {vr.efermi:.4f} eV | NELECT = {nelect:g}")
    g = A.band_gap(bs, nelect=nelect); print("\n[brecha]", g.message)
    dp = A.dirac_points(bs); print("[Dirac]", [(d["label"], round(d["energy"], 4), round(d["separation"], 5)) for d in dp])
    print("[clasificación]", A.classify(bs, vr.tdos, nelect=nelect), f"(DOS(E_F) = {A.dos_at_fermi(vr.tdos):.3f}, paso de la DOS = {A.dos_energy_step(vr.tdos):.3f} eV)")
    print("\n[v_F] método sobre la malla real vs. sobre el modelo analítico (t = 2.7 eV, ħv_F exacto = %.3f eV·Å) en los MISMOS puntos k:" % (3 ** 0.5 / 2 * S.A * S.T))
    tb = BandStructureSymmLine(pts, {Spin.up: S._energies(np.array(pts)).T if S._energies(np.array(pts)).ndim == 2 and S._energies(np.array(pts)).shape[0] == len(pts) else S._energies(np.array(pts))},
                               S.graphene_structure().lattice.reciprocal_lattice, 0.0, {"Γ": pts[0], "M": pts[9], "K": pts[13]}, structure=S.graphene_structure())
    for w in (0.6, 0.9):
        v = A.fermi_velocity(bs, dp[0], window=w)
        try:
            dpt = A.dirac_points(tb)[0]; vt = A.fermi_velocity(tb, dpt, window=w)["hbar_vF_eV_A"]
        except Exception as e: vt = float("nan")
        print(f"   ventana {w}: DFT ħv_F = {v['hbar_vF_eV_A']:.3f} eV·Å ({v['vF_m_s']:.2e} m/s) | modelo en los mismos k = {vt:.3f} → sesgo de muestreo = {vt / (3 ** 0.5 / 2 * S.A * S.T) - 1:+.1%}")
    iK = next(i for i, p in enumerate(pts) if np.allclose(p, [1/3, 1/3, 0]))
    near = np.argsort(np.abs(bs.bands[Spin.up][:, iK] - vr.efermi))[:2]
    pz, raw = A.orbital_weights(bs, ("pz",)), A.orbital_weights(bs, ("pz",), normalize=False)
    print("\n[pesos orbitales en K, bandas junto a E_F] pz relativo:", [round(float(pz[b, iK]), 3) for b in near], "| pz crudo de VASP:", [round(float(raw[b, iK]), 3) for b in near])
    e = bs.bands[Spin.up] - vr.efermi
    print("[energías en Γ]", np.round(e[:, 0], 2).tolist()); print("[energías en M]", np.round(e[:, 9], 2).tolist()); print("[energías en K]", np.round(e[:, iK], 3).tolist())
