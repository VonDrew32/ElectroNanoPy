"""Compara lo que calcula ElectroNanoPy con valores de la literatura para el grafeno (OE1 y OE3), con datos reales de VASP.

Uso:  python examples/comparar_literatura_grafeno.py CARPETA_GRAPHENE [salida.png]
      CARPETA_GRAPHENE contiene 3-Lattice_Parameter/ (serie de red) y 2-Encuts/600/ (vasprun con malla uniforme, LORBIT = 11).

Valores de referencia (con su fuente) en REF. Lo que NO se compara aquí por falta de cálculos: fonones (OE2).
"""
import sys, warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np
from electronanopy.project import scan_folder
from electronanopy.structural import series_table, series_summary
from electronanopy.electronic.analysis import HBAR_EV_S
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin

EV_A2_TO_N_M = 16.02176634
REF = {   # valor, fuente
    "a_exp": (2.459, "Å, experimento (Lindsay & Broido, PRB 81, 205441, 2010; arXiv:1003.2236)"),
    "nu": (0.149, "razón de Poisson DFT (revisión de mecánica de 2D, arXiv:1611.01555, ref. [46])"),
    "Y_dft": (345.0, "N/m, DFT (arXiv:1611.01555, Tabla I)"),
    "Y_exp": (340.0, "N/m, experimento (Lee et al., Science 321, 385, 2008; ±50 N/m)"),
    "vF_lda": (0.89e6, "m/s, DFT-LDA (arXiv:1006.3443)"),
    "vF_gga": (0.88e6, "m/s, DFT-GGA/PW91 (arXiv:1006.3443)"),
    "vF_exp": (1.0e6, "m/s, ARPES 1.0±0.05e6 (arXiv:1006.3443, Tabla I); otros experimentos hasta 1.1e6"),
}


def oe1(root):
    pf = scan_folder(root)
    ser = next(s for s in pf.series if s.kind == "lattice")
    rows = series_table(ser)
    sm = series_summary(ser, rows)
    bm = sm["bm"]
    a = np.array([r["value"] for r in rows]); e = np.array([r["energy"] for r in rows])
    out = {"a0_bm": bm["a0"], "a_sigma0": sm.get("a_zero_stress"), "bond": bm["a0"] / np.sqrt(3), "a": a, "e": e, "bm": bm}
    # C11 + C12 (biaxial) a partir de d²E/dε² = 2 (C11 + C12) A0 ;  A0 = √3/2 a0²  ->  C11+C12 = E''(a)/√3
    sel = np.abs(a - bm["a0"]) <= 0.08
    p = np.polyfit(a[sel], e[sel], 4)
    epp = np.polyval(np.polyder(p, 2), bm["a0"])
    out["c11c12_energy"] = epp / np.sqrt(3) * EV_A2_TO_N_M
    # idem con el esfuerzo: σ_xx (kB) · c (Å) · 0.1 = N/m ; (C11 + C12) = dσ/dε con ε = Δa/a0
    c = rows[0]["c"]
    s = np.array([r["stress"][0] for r in rows]); ss = np.abs(a - bm["a0"]) <= 0.05
    k = -np.polyfit(a[ss], s[ss], 1)[0] * bm["a0"]                    # VASP da la presión (σ = −esfuerzo): kB por unidad de deformación
    out["c11c12_stress"] = k * 0.1 * c * 0.1                          # kB -> GPa (×0.1) -> N/m·(c en nm: ×0.1·Å)
    for key in ("c11c12_energy", "c11c12_stress"):
        out["Y_" + key.split("_")[1]] = out[key] * (1 - REF["nu"][0])  # Y = (C11 + C12)(1 − ν)
    return out


def oe3(root, t_max=6, nk=18):
    vr = Vasprun(str(Path(root) / "vasprun.xml"), parse_eigen=True, parse_potcar_file=False)
    s = vr.final_structure
    k = np.array(vr.actual_kpoints)
    ev = vr.eigenvalues[Spin.up][:, :, 0]
    rec = s.lattice.reciprocal_lattice.matrix
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    ops = {tuple(np.round(o.rotation_matrix, 6).ravel()) for o in SpacegroupAnalyzer(s, symprec=1e-4).get_symmetry_operations()}
    R = [np.array(o).reshape(3, 3) for o in ops]
    rec_ops = [np.linalg.inv(r).T for r in R] + [-np.linalg.inv(r).T for r in R]

    def find(t):
        for i, kk in enumerate(k):
            for m in rec_ops:
                d = m @ t - kk
                if np.linalg.norm(d - np.round(d)) < 1e-4:
                    return i
        raise LookupError(t)

    K = np.array([1 / 3, 1 / 3, 0.0])
    nocc = int(round(float(vr.parameters["NELECT"]) / 2))
    ev_ = np.sort(ev, axis=1)
    iK = find(K)
    ed = (ev_[iK, nocc - 1] + ev_[iK, nocc]) / 2
    shells = {}
    for name, step in (("K→Γ", np.array([-1.0, -1.0, 0.0]) / nk), ("K→M", np.array([1.0, -2.0, 0.0]) / nk)):
        shells[name] = []
        for n in (1, 2, 3):
            kk = K + n * step
            q = float(np.linalg.norm((n * step) @ rec))
            i = find(kk)
            shells[name].append((q, ev_[i, nocc] - ed, ed - ev_[i, nocc - 1]))      # (q, π*−E_D, E_D−π)
    ts = np.arange(t_max + 1) / nk                                      # k = (t, t, 0): de Γ (0) a K (6/18)
    idx = [find(np.array([t, t, 0.0])) for t in ts]
    kc = np.array([np.array([t, t, 0.0]) @ rec for t in ts])
    q = np.linalg.norm(kc - kc[-1], axis=1)                              # distancia a K en 1/Å
    pi, pis = ev_[idx, nocc - 1], ev_[idx, nocc]
    return {"q": q, "pi": pi, "pis": pis, "ED": ed, "shells": shells, "efermi": vr.efermi, "gap_K": float(pis[-1] - pi[-1])}


def estimar_vf(shells):
    """ħv_F de las 'capas' {dirección: [(q, π*−E_D, E_D−π), ...]}: promedio de direcciones (cancela la deformación trigonal
    de primer orden) y Richardson (4v₁ − v₂)/3 (cancela el término q²)."""
    per = {n: [np.mean([up, dn]) / q for q, up, dn in sh] for n, sh in shells.items()}
    v1 = np.mean([p[0] for p in per.values()]); v2 = np.mean([p[1] for p in per.values()])
    return per, v1, (4 * v1 - v2) / 3


def capas_modelo_tb(nk=18):
    """Las mismas capas con el tight-binding de primeros vecinos (t = 2.7 eV): ħv_F exacto = (√3/2)·a·t."""
    from electronanopy.electronic import _synthetic as S
    rec = S.graphene_structure().lattice.reciprocal_lattice.matrix
    K = np.array([1 / 3, 1 / 3, 0.0])
    out = {}
    for name, step in (("K→Γ", np.array([-1.0, -1.0, 0.0]) / nk), ("K→M", np.array([1.0, -2.0, 0.0]) / nk)):
        out[name] = []
        for n in (1, 2, 3):
            e = S._energies((K + n * step)[None, :])[:, 0]
            out[name].append((float(np.linalg.norm((n * step) @ rec)), float(e[1]), float(-e[0])))
    return out, 3 ** 0.5 / 2 * S.A * S.T


def fit_vf(q, de, qmax):
    """|E − E_D| = ħv q + c q³ con los puntos q ≤ qmax (la corrección cúbica absorbe la curvatura de la banda)."""
    m = (q > 0) & (q <= qmax)
    A = np.vstack([q[m], q[m] ** 3]).T
    (v, c), *_ = np.linalg.lstsq(A, de[m], rcond=None)
    return float(v), float(c), int(m.sum())


if __name__ == "__main__":
    root = Path(sys.argv[1])
    r1 = oe1(root)
    r3 = oe3(root / "2-Encuts" / "600")
    print("== OE1: estructura ==")
    print(f"a0 (Birch–Murnaghan)      = {r1['a0_bm']:.4f} Å   | σ_xx = 0 en {r1['a_sigma0']:.4f} Å | exp. {REF['a_exp'][0]} Å "
          f"(+{(r1['a0_bm'] / REF['a_exp'][0] - 1) * 100:.2f} %; PBE sobreestima la red ~0.3 %)")
    print(f"enlace C–C                = {r1['bond']:.4f} Å   (exp. 1.42 Å)")
    for tag, lab in (("energy", "desde E(a)"), ("stress", "desde σ_xx(a)")):
        print(f"C11 + C12 {lab:14s}= {r1['c11c12_' + tag]:.0f} N/m  ->  Y2D = {r1['Y_' + tag]:.0f} N/m "
              f"(DFT {REF['Y_dft'][0]:.0f}, exp. {REF['Y_exp'][0]:.0f} ± 50 N/m; usa ν = {REF['nu'][0]})")
    print("\n== OE3: electrónico (malla uniforme 18×18) ==")
    print(f"E_F = {r3['efermi']:.3f} eV | E_Dirac en K = {r3['ED']:.3f} eV | separación π–π* en K = {r3['gap_K'] * 1000:.2f} meV")
    # ħv_F: en cada dirección E(q) = ħv q (1 ± α q + β q² …); el término α q cambia de signo entre K→Γ y K→M
    # (deformación trigonal), así que el promedio de las dos direcciones lo cancela; Richardson quita el q².
    per, v1, vr_ = estimar_vf(r3["shells"])
    for name, sh in r3["shells"].items():
        print(f"  {name}: ħv_n (n = 1, 2, 3 pasos de malla, q = {', '.join(f'{x[0]:.3f}' for x in sh)} 1/Å) = "
              + ", ".join(f"{v:.3f}" for v in per[name]) + " eV·Å")
    tb, exacto = capas_modelo_tb()
    _, tb1, tbr = estimar_vf(tb)
    print(f"  mismo estimador sobre el modelo tight-binding: {tbr:.3f} eV·Å (exacto {exacto:.3f}) -> sesgo del estimador {(tbr / exacto - 1) * 100:+.1f} %")
    for lab, v in (("promedio de direcciones, 1.ª capa", v1), ("con Richardson", vr_), ("con Richardson y corrección del sesgo", vr_ * exacto / tbr)):
        print(f"ħv_F ({lab}) = {v:.3f} eV·Å = {v * 1e-10 / HBAR_EV_S:.3e} m/s")
    print(f"literatura: LDA {REF['vF_lda'][0]:.2e}, GGA {REF['vF_gga'][0]:.2e} m/s (= {REF['vF_gga'][0] * HBAR_EV_S * 1e10:.2f} eV·Å); "
          f"experimento ≈ {REF['vF_exp'][0]:.1e} m/s")
    if len(sys.argv) > 2:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        bm = r1["bm"]
        ax[0].plot(r1["a"], (r1["e"] - bm["e0"]) * 1000, "o", label="VASP (serie de red)")
        ax[0].plot(bm["curve_a"], (bm["curve_e"] - bm["e0"]) * 1000, "-", label="Birch–Murnaghan")
        ax[0].axvline(REF["a_exp"][0], color="k", ls="--", label=f"experimento {REF['a_exp'][0]} Å")
        ax[0].set_xlabel("a (Å)"); ax[0].set_ylabel("E − E₀ (meV/celda)"); ax[0].set_ylim(-5, 120); ax[0].legend(fontsize=8)
        ax[0].set_title("OE1: parámetro de red")
        for (name, sh), mk in zip(r3["shells"].items(), ("^", "v")):          # cada dirección por separado (deformación trigonal)
            ax[1].plot([x[0] for x in sh], [x[1] for x in sh], mk, color="0.65", ms=5, label=f"π* {name}")
        sh = list(r3["shells"].values())
        qs = np.array([x[0] for x in sh[0]])
        prom = np.mean([[(x[1] + x[2]) / 2 for x in d] for d in sh], axis=0)
        ax[1].plot(qs, prom, "o", color="C0", label="promedio de direcciones y de π/π* (VASP)")
        qq = np.linspace(0, 0.9, 50)
        for lab, vf, ls in (("DFT-GGA lit. 0.88e6 m/s", REF["vF_gga"][0], "--"), ("experimento 1.0e6 m/s", REF["vF_exp"][0], ":")):
            ax[1].plot(qq, vf * HBAR_EV_S * 1e10 * qq, "k" + ls, label=lab)
        ax[1].set_xlabel("distancia a K (1/Å)"); ax[1].set_ylabel("|E − E_D| (eV)"); ax[1].legend(fontsize=8); ax[1].set_title("OE3: cono de Dirac (malla 18×18)")
        fig.tight_layout(); fig.savefig(sys.argv[2], dpi=160)
        print("figura:", sys.argv[2])
