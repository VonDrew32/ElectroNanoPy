# Pendientes y sugerencias

Lista viva de lo que quedó anotado para añadir más adelante. Marque con `[x]` lo que se complete.

## Módulo estructural (OE1, semanas 3 a 5)

- [x] **Validar el parseo con archivos reales:** 56 POSCAR/CONTCAR de grafeno (VASP, clúster) leídos sin errores (v0.6). Falta SrMnSb2.
- [ ] **Validar y afinar la detección de enlaces con SrMnSb2** (semana 4). La regla actual es
      `d <= tolerancia * (r_i + r_j)` con radios atómicos de pymatgen y tolerancia 1.15; solo se verificó con grafeno.
      Los vecinos Sr-Sb y Mn-Sb tienen distancias muy distintas.
- [ ] **Exponer la tolerancia de enlaces como parámetro de la interfaz** (OE4).
- [x] **Leer OUTCAR/OSZICAR:** energía (E0) y esfuerzo por cálculo en las series (v0.6); magnetización con parser propio.
- [ ] **Evaluar otra librería de visualización 3D** (PyVista o pyqtgraph) (semana 5). El esqueleto usa
      matplotlib 3D, que es limitado en interactividad. Decidirlo junto con la elección de la GUI.
- [ ] **Mantener una prueba por cada función nueva.**
- [ ] **Control de versiones:** `git init` y no subir archivos pesados (CHGCAR, WAVECAR; ya están en `.gitignore`).

## Fonones (OE2, semanas 6 a 8)

Estado (v0.8): módulo funcional según el Anexo 1 (importa resultados de Phonopy, calcula bandas y DOS e identifica las
frecuencias imaginarias como criterio de estabilidad). Probado con modelos sintéticos y con NaCl real; **sin fonones
reales de un material 2D**. Añadido en v0.8: diagnóstico de calidad (regla de suma acústica en Γ, tamaño de la supercelda),
opción de imponer la ASR (copia) y selector de camino (automático, 2D, 3D) en la pestaña Fonones.

Se intentó validar con el ejemplo de grafeno de Phonopy (`example/Graphene-siesta`, SIESTA, base mínima, 3×3×1); no sirve:
el propio ejemplo dice que no está convergido, Phonopy 4 interpreta mal su celda (grupo Cmmm en vez de P6/mmm) y el
desplazamiento usado no está documentado. Con él los 3 acústicos en Γ salen 0.54/3.25/3.25 THz y f_max = 34.5 THz
(grafeno real: ≈ 47 THz), y los diagnósticos nuevos lo señalan. Queda esperar los cálculos de VASP.

- [ ] **Probar con datos reales de VASP + Phonopy.** El módulo se probó con un modelo sintético de resortes,
      no con fonones reales. Validar con grafeno (el modo flexural ZA suele mostrar frecuencias imaginarias
      pequeñas cerca de Γ por ruido numérico).
- [ ] **Definir con los asesores la tolerancia de estabilidad** (por defecto 0.05 THz fuera de Γ) y el radio de
      exclusión alrededor de Γ para los modos acústicos (por defecto 0.1 1/Å).
- [ ] **Confirmar el método de fonones** (se asumió desplazamientos finitos con Phonopy). Si usan otro método
      (por ejemplo DFPT o dinámica molecular ab initio), hay que ajustar el OE2.
- [x] **Advertir sobre calidad de la entrada:** regla de suma acústica (`Thresholds.asr_tol`, 0.3 THz) y vector más corto de la supercelda (`min_phonon_cell`, 10 Å). Ambos umbrales son propuestos por nosotros: consultar con los asesores.
- [ ] **Camino k para 2D y 3D:** revisar que el camino automático sea el esperado en SrMnSb2 (si es 3D, usa seekpath;
      si se trata como lámina, usa el camino 2D).

## Electrónico (OE3, semanas 9 a 11)

Estado: bandas, DOS, brecha, puntos de Dirac, velocidad de Fermi, pesos orbitales, magnetismo (OUTCAR), densidad de carga y gráficas,
probados con grafeno tight-binding sintético (no con un vasprun.xml real).

- [x] **Leer un vasprun.xml real de grafeno:** 28 vasprun.xml reales leídos (v0.6); E_F desde la DOS o el OUTCAR cuando falta.
- [ ] **Energía de Fermi en el cálculo de bandas.** En el grafeno el punto de Dirac está justo en E_F; si el cálculo de bandas
      es no autoconsistente, conviene tomar E_F del SCF (`efermi="scf/vasprun.xml"`). Verificar que el cono quede en 0 eV.
- [ ] **Umbrales del análisis** (valores por defecto): brecha nula si < 0.01 eV; punto de Dirac si dos bandas
      se separan < 0.05 eV a menos de 0.1 eV de E_F; velocidad de Fermi con ventana de 0.1 1/Å; DOS en E_F < 0.1 estados/eV
      para «semimetal». Confirmar con los asesores.
- [ ] **Densidad de carga (CHGCAR): probar con uno real.** Implementada (electrones totales, promedio planar, cortes 2D,
      interpolación, isosuperficies). Validada con el CHGCAR real de `2-Encuts/600`: integra 8.000 electrones (v0.6).
      Un CHGCAR con 15-20 Å de vacío puede pesar decenas de MB.
- [ ] **Proyecciones (fat bands):** el vasprun real de `2-Encuts/600` (LORBIT = 11) contiene proyecciones; falta probarlas sobre
      bandas en línea (Γ–M–K–Γ).
- [ ] **Magnetismo:** con un OUTCAR real no polarizado (grafeno, ISPIN=1) devuelve «sin polarización» (verificado). Falta un
      OUTCAR polarizado real (SrMnSb2 será el caso magnético).
- [ ] **DOS:** revisar suavizado (ISMEAR/SIGMA) y malla k del cálculo, que afectan la DOS en E_F.

## Decisión de los asesores (2026-11-19)

Los umbrales de `report.Thresholds` se **mantienen como están** (ASR 0.3 THz, supercelda ≥ 10 Å, convergencia 5 meV, acuerdo de red 0.005 Å, `gap_tol`, `dos_tol`, `encut_tol`). Por ahora las pruebas se hacen solo con grafeno. Las preguntas de abajo quedan como referencia y ya no bloquean el desarrollo.

## Preguntas para los asesores (umbrales y criterios)

1. **Estabilidad dinámica:** ¿qué frecuencia mínima toleran como «ruido numérico»? (Por defecto: imaginaria si f < -0.05 THz.)
2. **Modos acústicos cerca de Γ:** ¿se aceptan pequeñas frecuencias imaginarias del modo flexural (ZA) en el grafeno, y hasta
   qué radio alrededor de Γ? (Por defecto: se ignoran a menos de 0.1 1/Å.)
3. **Calidad de la entrada:** ¿qué tamaño de supercelda y qué criterios de convergencia (fuerzas, energía, malla k, ENCUT) debe
   advertir la plataforma?
4. **Método de fonones:** ¿desplazamientos finitos con Phonopy, DFPT u otro?
5. **Unidades:** ¿THz o cm^-1 por defecto en las gráficas?
6. **Brecha y funcional:** ¿qué funcional (PBE, HSE06…) usan, y debe la plataforma mostrar una advertencia de que PBE subestima la brecha?
7. **Punto de Dirac:** ¿qué tolerancia aceptan para decir que el cono está en E_F y qué velocidad de Fermi esperan (referencia de la literatura)?
8. **Cálculos de validación:** ¿quién genera los de grafeno y SrMnSb2 (estructura, bandas, DOS, fonones) y para cuándo?

9. **Parámetro de red:** ¿qué acuerdo (Å) piden entre a₀ de Birch–Murnaghan y a con σ = 0? (Por defecto: 0.005 Å, `Thresholds.lattice_agreement`.)
11. **Fonones (calidad):** ¿qué tolerancia aceptan para la regla de suma acústica en Γ (propuesta 0.3 THz ≈ 10 cm⁻¹) y qué tamaño mínimo de supercelda (propuesto 10 Å)?
10. **Convergencia (series de cálculos):** la plataforma da por convergida una serie de malla k, ENCUT o vacío si la energía
   varía ≤ 5 meV/celda entre los tres últimos puntos (`Thresholds.conv_tol_mev`; valor propuesto por nosotros, sin respaldo
   aún). ¿Qué criterio usan (meV/átomo, 1 meV/átomo…)?

## Hallazgos con los datos reales de grafeno (carpeta `Graphene`, pruebas de convergencia)

Hecho (v0.6): series de cálculos (`project.Series`), un único cálculo coherente por defecto, tabla y gráfica por serie
(`structural.series`), E_F desde OUTCAR si el vasprun no la trae, notas resumidas, lectura de CHGCAR real (8.000 e⁻).

Por hacer o vigilar:
- **Gap en el punto de Dirac:** con `ISMEAR=0` las bandas π/π* tienen ocupación ≈0.5 en K; un gap calculado por
  ocupaciones da 1.79 eV falso. Verificar `band_gap` con bandas reales (aún no hay cálculo en línea Γ–M–K–Γ).
- **Birch–Murnaghan** (v0.7): implementado en `structural/eos.py`; con los datos del usuario a₀ = 2.4683 Å (error 0.05 meV). B₀ solo se da en 3D; en la monocapa con c fijo no es físico. Falta probarlo con SrMnSb2 (3D).
- **E_vs_*.csv** con números mal formados (`-18.463.627`): la plataforma no los lee, usa OSZICAR/vasprun.
- **INCAR de plantilla** en la raíz de cada serie (`incar = """…`): se ignora; solo se leen los de cada cálculo.
- **ENCUT no monótono** (18 meV entre 400 y 700 eV) y malla k con oscilaciones: revisar que los cálculos compartan el resto de
  parámetros.
- **Faltan cálculos** para OE2/OE3: bandas en línea, DOS dedicada, fonones.

## Interfaz y reporte (OE4, semanas 12 a 14)

Estado: ventana PySide6 con 4 pestañas y opciones de visualización, detección de archivos, carga tolerante a fallos,
resumen con resaltados/advertencias/sugerencias. Probada sin pantalla y con datos sintéticos; **no** con una carpeta
real de VASP.

- [x] **Probar con la carpeta real de grafeno:** detección corregida (series, un cálculo coherente, profundidad 4) en v0.6. Probado en Windows por el usuario.
- [ ] **Revisar los umbrales del reporte** (`report.Thresholds`): vacío mínimo 10 Å, supercelda mínima de 3 celdas,
      diferencia de E_F 0.05 eV, desplazamiento del punto de Dirac 0.05 eV, tolerancia de NELECT 0.01.
- [ ] **Uso en un equipo sin pantalla o remoto** (clúster): la interfaz necesita entorno gráfico; el análisis (`project`,
      `report`) funciona sin él.
- [ ] **Empaquetado** (opcional): ejecutable con PyInstaller para quien no tenga Python.
- [ ] **Exportar** figuras y el resumen (PNG/SVG/PDF, HTML o JSON) desde la ventana.
- [ ] **Mejoras de la vista 3D** de la estructura (ejes iguales, rotación más fluida): ver la evaluación de PyVista/pyqtgraph.

## Revisión y validación (v0.11.0)

Ver `docs/VALIDACION.md`. Con datos de SrMnSb2 queda por confirmar (v0.12 ya filtra las degeneraciones generales de `dirac_points` y limita la tolerancia entre metales a 1.0, probado solo con datos sintéticos): puntos de Dirac en 3D y coordinación de Sr/Mn/Sb.
- [ ] Probar `packaging/electronanopy.spec` en Windows (solo se pudo ensayar en Linux).

## Validación de OE3 con datos reales (v0.9.0)

Hecha con el `vasprun.xml` de `2-Encuts/600` (malla uniforme, sin bandas en línea), con `examples/validar_oe3_malla_uniforme.py`:

- [x] Dirac en K, brecha 0, ħv_F ≈ 5.07 eV·Å (sesgo de muestreo −6.7 %; corregido ≈ 5.4 eV·Å ≈ 8.3e5 m/s), peso pz = 1 en K, 8.000 e⁻.
- [x] Tres errores corregidos: brecha por conteo con NELECT, DOS con paso grueso (NEDOS), pesos orbitales normalizados.
- [ ] Falta (necesita cálculos): bandas en línea Γ–M–K–Γ con LORBIT = 11, DOS con NEDOS ≥ 2000, un OUTCAR con espín, SrMnSb2.

## Pruebas con datos públicos reales (v0.13.0)

Ver `docs/VALIDACION.md`, sección 4. Hecho: lectura de OUTCAR de un átomo y no colineal (SOC), brecha aproximada en la malla k,
enlaces solo entre metal y no metal en compuestos, estabilidad sin camino seekpath, estructura de respaldo.

- [ ] Con las bandas en línea del grafeno reemplazar la brecha de la malla (ya se prefiere automáticamente).
- [ ] Camino q para fonones 3D cuando seekpath cambia la celda (hoy solo se evalúa la malla).
- [x] CHGCAR grande: probado hasta 590 MB (ver `docs/VALIDACION.md`); isosuperficie limitada a 1 M de caras.

## Fonones del grafeno con VASP (v0.14.0)

Recibidos y validados (`docs/VALIDACION.md`, sección 5). Queda por decidir con los profesores:
- [ ] Repetir con supercelda 6×6×1 (K y M conmensurables) para confirmar que el TO en K baja a ≈ 1290 cm⁻¹.
- [ ] Bandas Γ–M–K–Γ de VASP (LORBIT = 11) y DOS con NEDOS ≥ 2000 para cerrar OE3/OE5 del grafeno.

## Versión 1.0.0 (cambios respecto a la 0.14)

- [x] **Camino de bandas de fonones en la base del usuario:** `phonons/paths.py` convierte el camino de seekpath a la base de
      la celda del usuario (matriz de rotación de spglib, comprobada como entera/unimodular); si no hay mapeo, se evalúa en la
      malla q. Verificado con `tools/verificar_mapeo_seekpath.py` y pruebas de regresión.
- [x] **Ejecutable de consola** `ElectroNanoPy-consola` junto al de ventana (`packaging/electronanopy.spec`) para ver errores en Windows.
- [x] **CHGCAR desde otra serie:** si la carpeta principal no tiene CHGCAR, se toma el de una serie (p. ej. `2-Encuts`) y se avisa
      en las notas; así la pestaña Densidad funciona con la carpeta `Graphene`.
- [x] Pruebas: 178 aprobadas. Número de versión 1.0.0 (ajustes posteriores: 1.0.x).
- [ ] **Prueba del ejecutable en un equipo Windows** (a cargo del autor; ver `packaging/CONSTRUIR_WINDOWS.md`).
- [x] **OE5:** validado con grafeno; los asesores pidieron concentrar el desarrollo en este material.
- [ ] Trabajo futuro: validar con SrMnSb2 y, para el grafeno, bandas en línea y DOS de alta resolución.

## Siguiente en el plan

- [x] **OE4:** integración en una aplicación de escritorio (primera versión hecha; ver sección anterior).
- [ ] **OE5 (semanas 15 y 16):** validar con grafeno (y luego SrMnSb2) contra la literatura; requiere los cálculos de VASP.

## Decisiones abiertas (afectan al Anexo 1 y a OE3/OE4)

- [x] **Qué significa "visualización inteligente":** reglas automáticas (detección de brecha, punto de Dirac, frecuencias
      imaginarias; comprobaciones de consistencia; sugerencia de gráficos), sin machine learning (decidido).
- [x] **Framework de la interfaz:** PySide6 (decidido; actualizar el informe de avance, que menciona PyQt).
- [ ] **Quién genera los cálculos de validación** (grafeno y SrMnSb2) y en qué fecha estarán listos para OE5.
- [ ] **Informe de avance:** actualizar herramientas, pesos y porcentajes con lo que realmente se hizo.

## Posibles mejoras (no comprometidas en el Anexo 1)

- [x] **Ajuste Birch–Murnaghan** del barrido de red (v0.7, `structural/eos.py`); lee energías de OSZICAR, no de tablas.

- [x] Exportación de figuras (PNG/SVG/PDF), reporte JSON/HTML e informe con figuras (v0.10).
- [x] Estilos de exportación de figuras: pantalla, publicación 1/2 columnas y póster (v0.12). Pendiente: paleta configurable.
- [x] Reglas de consistencia entre cálculos (ENCUT, funcional, POTCAR, ISPIN; celda de fonones vs. estructura) (v0.10).
- [ ] Consistencia pendiente de probar con datos reales que tengan varios cálculos (SCF + bandas + fonones).
