"""Reporte inteligente (OE4): resaltados, advertencias y sugerencias de visualización por módulo.

«Inteligente» significa reglas automáticas y transparentes (no machine learning): detectar lo relevante en
los resultados (brecha, punto de Dirac, frecuencias imaginarias, momentos magnéticos), comprobar la
consistencia entre archivos (E_F del SCF vs. bandas, electrones integrados vs. NELECT, tamaño de supercelda) y
proponer qué graficar con parámetros sensatos (supercelda 2D, plano del corte de densidad, nivel de la isosuperficie).
Todos los umbrales son parámetros de ``Thresholds``.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from .project import Session


@dataclass
class Thresholds:
    min_vacuum_2d: float = 10.0          # Å: vacío mínimo recomendado en un material 2D
    asr_tol: float = 0.3                 # THz: |f| máxima admitida de los modos acústicos en Γ (regla de suma acústica)
    max_dos_step: float = 0.05           # eV: paso máximo de la rejilla de energía de la DOS para fiarse del valor en E_F
    min_phonon_cell: float = 10.0        # Å: vector más corto admitido de la supercelda de fonones
    min_supercell: int = 3               # celdas mínimas por dirección en la supercelda de fonones (en el plano)
    fermi_mismatch: float = 0.05         # eV: diferencia máxima entre E_F del SCF y de las bandas
    dirac_offset: float = 0.05           # eV: desplazamiento máximo del punto de Dirac respecto a E_F
    nelect_tol: float = 0.01             # electrones: diferencia máxima entre la carga integrada y NELECT
    stability_tol: float = 0.05          # THz
    gamma_radius: float = 0.1            # 1/Å
    lattice_agreement: float = 0.005     # Å: máxima diferencia entre a0 (Birch–Murnaghan) y a con σ_xx = 0
    encut_tol: float = 1.0               # eV: diferencia admitida del corte de energía entre cálculos
    conv_tol_mev: float = 5.0            # meV/celda: variación admitida entre los últimos puntos de una serie de convergencia


@dataclass
class ModuleReport:
    name: str
    available: bool = False
    summary: dict = field(default_factory=dict)
    highlights: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    suggestions: list[str] = field(default_factory=list)
    params: dict = field(default_factory=dict)       # parámetros sugeridos para la interfaz


@dataclass
class Report:
    modules: dict[str, ModuleReport]

    def all_warnings(self) -> list[str]:
        return [f"[{m.name}] {w}" for m in self.modules.values() for w in m.warnings]


def _structure(s: Session, th: Thresholds) -> ModuleReport:
    from .structural import summary
    r = ModuleReport("Estructura", available=s.structure is not None)
    if not r.available:
        r.warnings.append(s.errors.get("estructura", "No se encontró CONTCAR/POSCAR."))
        return r
    d = summary(s.structure)
    r.summary = d
    if getattr(s, "structure_source", None):
        r.warnings.append(f"No hay un CONTCAR/POSCAR de la celda unidad: la estructura se tomó de {s.structure_source}.")
    L = d["lattice"]
    r.highlights.append(f"{d['formula']}: {d['n_atoms']} átomos, grupo espacial {d['space_group']} (#{d['space_group_number']}).")
    r.highlights.append(f"Red: a = {L['a']:.3f} Å, b = {L['b']:.3f} Å, c = {L['c']:.3f} Å; "
                        f"α = {L['alpha']:.1f}°, β = {L['beta']:.1f}°, γ = {L['gamma']:.1f}°.")
    if d["bond_min"] is not None:
        r.highlights.append(f"Enlaces entre {d['bond_min']:.3f} y {d['bond_max']:.3f} Å; "
                            f"coordinación {sorted(set(d['coordination']))}.")
    if d["two_dimensional"]:
        r.highlights.append(f"Material bidimensional (vacío de {d['vacuum']:.1f} Å a lo largo de c).")
        r.params["supercell"] = (3, 3, 1)
        r.suggestions.append("Mostrar una supercelda 3×3×1 para ver la red completa.")
        if d["vacuum"] < th.min_vacuum_2d:
            r.warnings.append(f"El vacío ({d['vacuum']:.1f} Å) es menor que {th.min_vacuum_2d:.0f} Å: "
                              "posible interacción entre capas periódicas.")
    else:
        r.params["supercell"] = (2, 2, 2)
        r.suggestions.append("Mostrar una supercelda 2×2×2.")
    _series(s, th, r)
    return r


_SERIES_LABEL = {"lattice": "parámetro de red", "kpoints": "malla k", "encut": "ENCUT", "vacuum": "vacío", "other": "parámetro"}


def _series(s: Session, th: Thresholds, r: ModuleReport):
    """Resaltados y advertencias de las series de cálculos (convergencia y parámetro de red)."""
    from .structural import series_summary
    for ser, rows in s.series:
        sm = series_summary(ser, rows)
        label = f"«{ser.name}» ({_SERIES_LABEL[ser.kind]}, {sm['n']} cálculos)"
        r.summary.setdefault("series", {})[ser.name] = sm
        bad = [x["name"] for x in rows if x["error"]]
        if bad:
            r.warnings.append(f"Serie {label}: no se pudo leer la estructura de {', '.join(bad)}.")
        if sm["n_energy"] < 2:
            r.warnings.append(f"Serie {label}: faltan energías (OSZICAR) para evaluarla.")
            continue
        if ser.kind == "lattice":
            msg = f"Serie {label}: mínimo de energía en a = {sm['a_min']:.3f} Å"
            if "a_zero_stress" in sm:
                msg += f"; el esfuerzo σ_xx se anula en a ≈ {sm['a_zero_stress']:.3f} Å"
            r.highlights.append(msg + ".")
            bm = sm.get("bm")
            if bm is None:
                r.warnings.append(f"Serie {label}: se necesitan al menos 5 energías para ajustar Birch–Murnaghan.")
            else:
                txt = f"Birch–Murnaghan: a₀ = {bm['a0']:.4f} Å, E₀ = {bm['e0']:.5f} eV, V₀ = {bm['v0']:.2f} Å³"
                if bm["b0_gpa"] is not None:
                    txt += f", B₀ = {bm['b0_gpa']:.1f} GPa, B₀' = {bm['b0p']:.2f}"
                elif abs(bm["exponent"] - 3.0) > 0.15 and not sm.get("is_2d"):
                    txt += (f" (sin B₀: la celda no escala por igual en las tres direcciones, n = {bm['exponent']:.2f}; "
                            "para B₀ varíe a, b y c juntos, p. ej. con ISIF = 3)")
                lat = bm.get("lattice0") or {}
                if lat and any(abs(lat[k] - bm["a0"]) > 1e-3 for k in ("a",) if k in lat):
                    txt += ("; en V₀ la celda del CONTCAR tiene " + ", ".join(f"{k} = {v:.4f} Å" for k, v in lat.items()) +
                            " (no coincide con el parámetro de la serie: la celda cambió de forma al relajarse)")
                r.highlights.append(f"Serie {label}: {txt} (error de ajuste {bm['rmse_mev']:.2f} meV, {bm['n_points']} puntos).")
                r.params["a0"] = bm["a0"]
                if bm["extrapolated"]:
                    r.warnings.append(f"Serie {label}: el mínimo del ajuste está fuera del rango calculado; amplíe la serie.")
                za = sm.get("a_zero_stress")
                if za is not None and abs(bm["a0"] - za) > th.lattice_agreement:
                    r.warnings.append(f"Serie {label}: a₀ de Birch–Murnaghan ({bm['a0']:.4f} Å) y a con σ_xx = 0 ({za:.4f} Å) "
                                      f"difieren más de {th.lattice_agreement:g} Å; revise ENCUT/Pulay o el rango de la serie.")
                r.suggestions.append(f"Repita el cálculo (relajación, SCF, bandas) en a₀ = {bm['a0']:.4f} Å con los mismos parámetros.")
        else:
            span = sm["spread_last3_mev"]
            ok = span <= th.conv_tol_mev
            r.highlights.append(f"Serie {label}: variación de {span:.1f} meV entre los tres últimos puntos "
                                f"({'dentro' if ok else 'fuera'} de {th.conv_tol_mev:g} meV).")
            if not ok:
                r.warnings.append(f"Serie {label} sin converger: {span:.1f} meV entre los tres últimos puntos "
                                  f"(umbral {th.conv_tol_mev:g} meV). Revise parámetros comunes a todos los cálculos.")
            elif not sm["monotonic"]:
                r.warnings.append(f"Serie {label}: la energía no es monótona; con el punto de Dirac y suavizado "
                                  "gaussiano es habitual, pero conviene justificar la elección.")


def _phonons(s: Session, th: Thresholds) -> ModuleReport:
    from .phonons import check_stability, compute_bands, phonon_quality
    from .phonons.analysis import THZ_TO_CM1
    r = ModuleReport("Fonones", available=s.phonon is not None)
    if not r.available:
        if s.files is not None and s.files.phonons is None:
            r.warnings.append("No se encontraron archivos de Phonopy (phonopy_params.yaml o phonopy_disp.yaml + FORCE_SETS).")
        else:
            r.warnings.append(s.errors.get("fonones", "No disponible."))
        return r
    ph = s.phonon
    from .phonons.analysis import _vacuum_two_d
    try:
        bands = compute_bands(ph)
    except Exception as exc:                       # p. ej. seekpath cambia la celda: se evalúa solo en la malla q
        bands = None
        r.warnings.append(f"No se pudo construir el camino de alta simetría ({exc}). La estabilidad se evaluó solo en la "
                          "malla q; indique el camino manualmente en la pestaña Fonones.")
    two_d = bands.two_d if bands is not None else bool(_vacuum_two_d(ph))
    try:
        rep = check_stability(ph, tol=th.stability_tol, gamma_radius=th.gamma_radius, bands=bands, two_d=two_d)
    except Exception as exc:
        r.warnings.append(f"No se pudo analizar la estabilidad: {exc}")
        return r
    fmax = float((bands.frequencies if bands is not None else ph.mesh.frequencies).max())
    r.summary = {"stable": rep.stable, "f_min": rep.min_frequency, "f_max": fmax, "two_d": two_d}
    r.highlights.append(rep.message)
    r.highlights.append(f"Frecuencia máxima: {fmax:.2f} THz ({fmax * THZ_TO_CM1:.0f} cm⁻¹).")
    if not rep.stable:
        r.warnings.append("Hay frecuencias imaginarias: revise la convergencia, el tamaño de la supercelda y que la "
                          "estructura esté bien relajada antes de concluir que es inestable.")
    dims = np.diag(np.asarray(ph.supercell_matrix)).astype(int)
    check = dims[:2] if two_d else dims
    if check.min() < th.min_supercell:
        r.warnings.append(f"Supercelda {'×'.join(map(str, dims))}: menor que {th.min_supercell} celdas en alguna "
                          "dirección; las constantes de fuerza pueden no estar convergidas.")
    if check.min() >= th.min_supercell:                     # si ya avisó por el número de celdas, no se repite
        q = phonon_quality(ph, two_d=two_d)
        if q.min_cell_length < th.min_phonon_cell - 1e-6:
            r.warnings.append(f"La supercelda tiene un vector de {q.min_cell_length:.1f} Å (< {th.min_phonon_cell:g} Å): "
                              "las constantes de fuerza pueden no estar convergidas con el tamaño de la supercelda.")
    q = phonon_quality(ph, two_d=two_d)
    r.summary.update(asr_max=q.asr_max, min_cell_length=q.min_cell_length)
    r.highlights.append("Modos acústicos en Γ: " + ", ".join(f"{(0.0 if abs(x) < 5e-3 else x):.2f}" for x in q.gamma_acoustic) + " THz (deben ser ≈ 0).")
    if q.asr_max > th.asr_tol:
        r.warnings.append(f"No se cumple la regla de suma acústica: |f| = {q.asr_max:.2f} THz en Γ (umbral {th.asr_tol:g} THz). "
                          "Suele indicar ruido en las fuerzas (EDIFF, malla k, ENCUT, suavizado); las frecuencias cerca de Γ "
                          "no son fiables. Puede imponer la regla en la pestaña Fonones, pero conviene mejorar el cálculo.")
    if bands is not None:
        r.suggestions.append("Ver bandas y DOS con el camino " + ("Γ–M–K–Γ (2D)." if two_d else "de alta simetría (3D)."))
    else:
        r.suggestions.append("Ver la DOS de fonones; para las bandas indique el camino q a mano.")
    return r


def _electronic(s: Session, th: Thresholds) -> ModuleReport:
    from .electronic import band_gap, dirac_points, fermi_velocity, classify, dos_at_fermi, dos_energy_step
    mesh = getattr(s, "mesh_bands", None)
    r = ModuleReport("Electrónico", available=(s.bands is not None or mesh is not None or s.tdos is not None
                                               or s.magnetization is not None))
    if not r.available:
        r.warnings.append(s.errors.get("electronico", "No se encontraron resultados electrónicos (vasprun.xml)."))
        return r
    if s.bands is not None:
        mag, soc = (s.magnetization or {}).get("total"), bool(s.lsorbit)
        g = band_gap(s.bands, nelect=s.nelect, magnetization=mag, soc=soc)
        r.summary.update(gap=g.gap, direct=g.direct, is_metal=g.is_metal)
        r.highlights.append(g.message)
        dp = dirac_points(s.bands)
        if dp:
            labels = ", ".join(sorted({str(d["label"]) for d in dp}))
            r.highlights.append(f"Punto(s) de Dirac en {labels}.")
            try:
                v = fermi_velocity(s.bands, dp[0])
                r.summary["vF_m_s"] = v["vF_m_s"]
                r.highlights.append(f"ħv_F = {v['hbar_vF_eV_A']:.2f} eV·Å (v_F = {v['vF_m_s']:.2e} m/s).")
            except ValueError as exc:
                r.warnings.append(f"No se pudo estimar v_F: {exc}")
            if abs(dp[0]["energy"]) > th.dirac_offset:
                r.warnings.append(f"El punto de Dirac está a {dp[0]['energy']:+.3f} eV de E_F; "
                                  "revise la energía de Fermi (conviene tomarla del SCF).")
        r.summary["classification"] = classify(s.bands, s.tdos, nelect=s.nelect, magnetization=mag, soc=soc)
        r.summary["dirac"] = sorted({str(d["label"]) for d in dp})
        r.highlights.append("Clasificación orientativa: " + r.summary["classification"] + ".")
        if s.scf_efermi is not None and s.bands_efermi_native is not None:
            diff = abs(s.scf_efermi - s.bands_efermi_native)
            if diff > th.fermi_mismatch:
                r.warnings.append(f"E_F del SCF ({s.scf_efermi:.3f} eV) y del cálculo de bandas "
                                  f"({s.bands_efermi_native:.3f} eV) difieren en {diff:.3f} eV; se usó la del SCF.")
        elif s.scf_efermi is None:
            r.suggestions.append("Sin SCF: se usa la E_F del cálculo de bandas; con un SCF aparte sería más fiable.")
        if getattr(s.bands, "projections", None):
            r.suggestions.append("Hay proyecciones orbitales: puede ver bandas con peso orbital (p. ej. pz).")
            r.params["has_projections"] = True
        else:
            r.suggestions.append("Para bandas con peso orbital, calcule con LORBIT = 11.")
    elif mesh is not None:                         # sin bandas en línea: brecha aproximada en la malla k del SCF/DOS
        mag, soc = (s.magnetization or {}).get("total"), bool(s.lsorbit)
        try:
            g = band_gap(mesh, nelect=s.nelect, magnetization=mag, soc=soc)
            r.summary.update(gap=g.gap, direct=g.direct, is_metal=g.is_metal, gap_source="malla")
            r.highlights.append(g.message + " (aproximada: solo los k de la malla del SCF/DOS; el valor correcto "
                                "requiere bandas en el camino de alta simetría.)")
            r.summary["classification"] = classify(mesh, s.tdos, nelect=s.nelect, magnetization=mag, soc=soc)
            r.highlights.append("Clasificación orientativa (malla k): " + r.summary["classification"] + ".")
        except Exception as exc:
            r.warnings.append(f"No se pudo estimar la brecha en la malla k: {exc}")
        r.suggestions.append("Para la brecha real y las bandas, calcule en modo línea (KPOINTS con camino de alta simetría).")
    if s.tdos is not None:
        v = dos_at_fermi(s.tdos)
        r.summary["dos_ef"] = v
        r.highlights.append(f"DOS en E_F: {v:.3f} estados/eV por celda.")
        step = dos_energy_step(s.tdos)
        r.summary["dos_step"] = step
        if step > th.max_dos_step:
            r.warnings.append(f"La DOS tiene un paso de energía de {step:.3f} eV (NEDOS demasiado bajo): el valor en E_F "
                              f"no es fiable. Recalcule la DOS con NEDOS ≥ 2000 (paso < {th.max_dos_step:g} eV).")
    m = s.magnetization
    if m is not None:
        if m["spin_polarized"]:
            if m["total"] is not None:
                tot = 0.0 if abs(m["total"]) < 5e-4 else m["total"]                # evita «-0.000»
                r.highlights.append(f"Momento magnético total: {tot:.3f} µB"
                                    + (" (módulo del vector; cálculo no colineal)." if m.get("noncollinear") else "."))
                r.summary["magnetization"] = m["total"]
                local = [x for x in (m.get("per_atom") or []) if abs(x["tot"]) >= 0.1]
                if local:
                    local.sort(key=lambda x: -abs(x["tot"]))
                    txt = ", ".join(f"{x.get('element', '#' + str(x['index'] + 1))}{x['index'] + 1}: "
                                    f"{x.get('tot_signed', x['tot']):+.2f}" for x in local[:6])
                    r.highlights.append("Momentos locales (µB): " + txt + (" …" if len(local) > 6 else "") + ".")
                    r.summary["local_moments"] = [(x["index"], x.get("tot_signed", x["tot"])) for x in local]
                    if abs(sum(x.get("tot_signed", x["tot"]) for x in local)) < 0.05 * max(abs(x["tot"]) for x in local):
                        r.highlights.append("Los momentos locales se compensan: ordenamiento antiferromagnético "
                                            "(o ferrimagnético compensado).")
            else:
                r.highlights.append("Cálculo con polarización de espín (sin momento total en el OUTCAR).")
        else:
            r.highlights.append("Cálculo sin polarización de espín (sin momentos magnéticos).")
    return r


def _density(s: Session, th: Thresholds) -> ModuleReport:
    from .electronic import total_electrons, density
    r = ModuleReport("Densidad", available=s.chg is not None)
    if not r.available:
        r.warnings.append(s.errors.get("densidad", "No se encontró CHGCAR (calcule el SCF con LCHARG = .TRUE.)."))
        return r
    ne = total_electrons(s.chg)
    r.summary["electrons"] = ne
    r.highlights.append(f"Electrones integrados en la celda: {ne:.3f}.")
    if s.nelect is not None and abs(ne - s.nelect) > th.nelect_tol:
        r.warnings.append(f"Los electrones integrados ({ne:.3f}) difieren de NELECT ({s.nelect:.0f}): "
                          "el CHGCAR podría no corresponder al mismo cálculo.")
    rho = density(s.chg)
    r.highlights.append(f"Densidad máxima: {rho.max():.3f} e/Å³.")
    st = s.chg.structure
    r.params["slice_axis"] = 2
    ang = 2 * np.pi * st.frac_coords[:, 2]                      # media circular: la celda es periódica en c
    zc = np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()) / (2 * np.pi) % 1.0
    r.params["slice_frac"] = float(np.round(zc, 3))             # plano que contiene a los átomos (2D)
    r.params["iso_level"] = float(0.1 * rho.max())              # lejos de los núcleos: muestra el enlace y no esferas
    r.suggestions.append(f"Corte a lo largo de c en fracción {r.params['slice_frac']:.3f} (plano de los átomos) "
                         f"e isosuperficie en {r.params['iso_level']:.3f} e/Å³.")
    if s.magnetization and s.magnetization.get("spin_polarized"):
        r.suggestions.append("Cálculo con espín: también puede ver la magnetización ('diff').")
    return r


def _consistency(s: Session, th: Thresholds, a0: float | None = None) -> ModuleReport:
    """Coherencia entre cálculos: ENCUT, funcional, pseudopotenciales, ISPIN; y entre estructuras (celda de fonones)."""
    from .consistency import gather_settings, compare_settings
    r = ModuleReport("Consistencia entre cálculos")
    named = gather_settings(s.files) if s.files is not None else {}
    r.summary = {k: {"encut": v.encut, "funcional": v.functional, "ispin": v.ispin, "kmesh": v.kmesh}
                 for k, v in named.items()}
    if len(named) >= 2:
        r.available = True
        r.warnings += compare_settings(named, th.encut_tol)
        if not r.warnings:
            r.highlights.append("Cálculos comparables: " + ", ".join(named) +
                                " usan el mismo corte, funcional y pseudopotenciales.")
    ph, st = s.phonon, s.structure
    if ph is not None and st is not None and a0 is not None and getattr(s, "structure_source", None):
        # La estructura no es un CONTCAR propio sino la de otro cálculo (p. ej. un punto de la serie de red): la
        # referencia correcta para la celda de los fonones es el a₀ de la serie.
        try:
            a = float(np.linalg.norm(np.asarray(ph.unitcell.cell, dtype=float)[0]))
            r.available = True
            if abs(a - a0) > th.lattice_agreement:
                r.warnings.append(f"La celda de los fonones (a = {a:.4f} Å) difiere del a₀ de la serie de red "
                                  f"({a0:.4f} Å): los fonones deben calcularse con la celda relajada.")
            else:
                r.highlights.append(f"La celda de los fonones (a = {a:.4f} Å) coincide con el a₀ de la serie de red "
                                    f"({a0:.4f} Å).")
        except Exception:
            pass
    elif ph is not None and st is not None:
        try:
            cell = np.asarray(ph.unitcell.cell, dtype=float)
            lp = np.linalg.norm(cell, axis=1)
            ls = np.array(st.lattice.abc)
            bad = np.abs(lp - ls) > th.lattice_agreement
            if bad.any():
                r.available = True
                r.warnings.append(
                    "La celda de los fonones (" + ", ".join(f"{x:.3f}" for x in lp) + " Å) no coincide con la estructura "
                    "cargada (" + ", ".join(f"{x:.3f}" for x in ls) + " Å): los fonones deben calcularse con la celda relajada.")
        except Exception:
            pass
    return r


def _integration(s: Session, th: Thresholds, m: dict) -> ModuleReport:
    """Reglas que relacionan varios módulos (OE4). Cada regla es una condición explícita sobre los resultados de otros
    módulos; no hay aprendizaje automático. También propone qué pestaña mostrar primero (``params['tab']``)."""
    r = ModuleReport("Integración entre módulos")
    st, ph, el = m["estructura"], m["fonones"], m["electronico"]
    two_d = bool(st.summary.get("two_dimensional"))
    unstable = ph.available and ph.summary.get("stable") is False
    stable = ph.available and ph.summary.get("stable") is True

    # R1: carácter electrónico + estabilidad dinámica
    cls = el.summary.get("classification")
    if cls and el.summary.get("is_metal") is not None:
        dirac = el.summary.get("dirac")
        if el.summary["is_metal"] and dirac and two_d:
            vf = el.summary.get("vF_m_s")
            r.highlights.append("Material 2D semimetálico con cono de Dirac en " + ", ".join(dirac) +
                                (f" (v_F = {vf:.2e} m/s)" if vf else "") + ".")
        if stable:
            r.highlights.append(f"Dinámicamente estable y electrónicamente: {cls}.")
        elif unstable:
            r.warnings.append(f"Electrónicamente «{cls}», pero con frecuencias imaginarias: la estructura no es un mínimo "
                              "estable, así que las propiedades electrónicas corresponden a una estructura inestable.")
    # R2: posibles causas de una inestabilidad
    if unstable:
        if ph.summary.get("asr_max", 0.0) > th.asr_tol:
            r.warnings.append(f"Los modos acústicos en Γ valen hasta {ph.summary['asr_max']:.2f} THz (> {th.asr_tol:g}): la "
                              "inestabilidad puede ser ruido numérico. Active «ASR» en la pestaña Fonones y compare.")
        if two_d and (st.summary.get("vacuum") or 0.0) < th.min_vacuum_2d:
            r.warnings.append(f"El vacío es de {st.summary['vacuum']:.1f} Å (< {th.min_vacuum_2d:g} Å): la interacción entre "
                              "capas puede alterar las ramas de baja frecuencia (ZA).")
        if ph.summary.get("min_cell_length", 99.0) < th.min_phonon_cell:
            r.warnings.append("La supercelda de fonones es pequeña: pruebe una mayor antes de concluir que hay inestabilidad.")
    # R3: la estructura calculada frente al mínimo de la serie de red
    a0 = st.params.get("a0")
    if a0 is not None and s.structure is not None and not getattr(s, "structure_source", None):
        a = s.structure.lattice.a
        if abs(a - a0) > th.lattice_agreement:
            r.suggestions.append(f"La estructura cargada tiene a = {a:.4f} Å y el mínimo de la serie es a₀ = {a0:.4f} Å: "
                                 "para los siguientes cálculos (bandas, fonones) use a₀.")
    # R4: magnetismo
    mag = el.summary.get("magnetization")
    if mag is not None and abs(mag) > 0.1:
        r.highlights.append(f"Momento magnético total {mag:+.2f} µB: compare configuraciones (FM/AFM) con las energías de OSZICAR.")
    # vista sugerida (una regla simple y visible)
    tab, why = None, None
    if unstable and ph.available:
        tab, why = "fonones", "hay frecuencias imaginarias"
    elif el.available and "gap" in el.summary:
        tab, why = "electronico", "hay bandas: brecha y puntos de Dirac"
    elif st.available:
        tab, why = "estructura", "estructura y series de convergencia"
    if tab:
        r.params["tab"], r.params["reason"] = tab, why
    r.available = bool(r.highlights or r.warnings or r.suggestions)
    return r


def analyze(session: Session, thresholds: Thresholds | None = None) -> Report:
    """Genera el reporte de todos los módulos a partir de una sesión cargada."""
    th = thresholds or Thresholds()
    mods = {"estructura": _structure(session, th), "fonones": _phonons(session, th),
            "electronico": _electronic(session, th), "densidad": _density(session, th),
            "consistencia": None}
    mods["consistencia"] = _consistency(session, th, mods["estructura"].params.get("a0"))
    mods["integracion"] = _integration(session, th, mods)
    return Report(mods)
