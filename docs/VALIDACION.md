# Validación con grafeno (v0.11)

Qué se ha comprobado con datos reales y contra la literatura, y qué falta. Todo se puede reproducir con los scripts de `examples/`.

## 1. Revisión del código

Una revisión independiente (sin conocer el código) ejecutó la suite y probó cada módulo con los datos reales de grafeno. **No halló errores numéricos con grafeno.** Sí halló errores en casos que aparecerán con SrMnSb2, ya corregidos en v0.11 y cubiertos por `tests/test_auditoria.py`:

| Hallazgo | Corrección |
|---|---|
| Las pistas `scf`/`dos`/`band` se buscaban en la ruta absoluta (un directorio superior llamado `...dos...` cambiaba la clasificación) | Solo se mira la ruta relativa a la carpeta del proyecto; el orden de los archivos es determinista |
| Brecha con ISPIN = 2 dependía de dónde cayera E_F | Conteo de bandas por canal de espín con el momento total (N↑ = (N+M)/2, N↓ = (N−M)/2) |
| Brecha con espín-órbita contaba NELECT/2 bandas | Con `LSORBIT` se cuentan NELECT bandas |
| OUTCAR con ISPIN = 2 sin LORBIT se leía como «sin espín» | Se lee `number of electron … magnetization` y ISPIN |
| Pesos orbitales con LORBIT = 10 leían otra columna | Se valida el número de orbitales y se avisa |
| Series con carpetas `6x6x1`, `12x12x1` | Valor numérico = primer entero; orden natural |
| `Reduced_cell_tests`, `supercell`… se tomaban como serie de red | Tipo de serie por palabras completas |
| B₀ en 3D con solo `a` variando (n = 2) | B₀ solo si la celda escala por igual (|n − 3| ≤ 0.15) |
| Corte de densidad en 0.5 con una lámina en el borde de la celda | Media circular de z |
| Isosuperficie al 40 % del máximo (esferas en los núcleos) | 10 % del máximo |
| Un `vasprun` de relajación se usaba como SCF antes que el SCF | Las relajaciones solo son el último recurso |

En v0.12 también se filtraron las degeneraciones generales en `dirac_points` y se limitó la tolerancia de enlace entre metales a 1.0 (`tests/test_integracion.py`). Con datos sintéticos funcionan; falta confirmarlo con SrMnSb2 real.

## 2. Comparación con la literatura

`python examples/comparar_literatura_grafeno.py CARPETA_Graphene figura.png`

| Magnitud | ElectroNanoPy (tus cálculos VASP, PBE) | Literatura |
|---|---|---|
| a (mínimo de Birch–Murnaghan) | 2.4683 Å (σ_xx = 0 en 2.4674 Å) | 2.459 Å, experimento (Lindsay y Broido, arXiv:1003.2236). PBE suele sobreestimar la red ~0.3 %: aquí +0.38 % |
| Enlace C–C | 1.425 Å | 1.42 Å |
| C₁₁ + C₁₂ | 414 N/m (de E(a)); 416 N/m (de σ_xx(a)) | — |
| Y₂D = (C₁₁ + C₁₂)(1 − ν) con ν = 0.149 | 353 y 354 N/m | 345 N/m DFT; 340 ± 50 N/m experimento (arXiv:1611.01555) |
| Cono de Dirac | En K, separación 0 meV, E_D = E_F | — |
| ħv_F | 5.44 eV·Å = 8.3×10⁵ m/s | DFT-LDA 8.9×10⁵, GGA(PW91) 8.8×10⁵ (arXiv:1006.3443); LDA 8.3×10⁵ (arXiv:1303.3885); experimento 1.0–1.1×10⁶ m/s |
| Electrones integrados del CHGCAR | 8.000 | 8 (4 por átomo) |

Notas:
- El módulo elástico sale de dos medidas independientes (energía y esfuerzo del OUTCAR) y concuerdan entre sí al 0.3 %; ambas dentro del 3 % del valor DFT de la literatura.
- ħv_F se estima en la malla uniforme 18×18 (no hay bandas en línea): promedio de las direcciones K→Γ y K→M (la deformación trigonal cambia de signo entre ambas) y extrapolación de Richardson. El mismo estimador aplicado a un modelo tight-binding da el valor exacto (sesgo 0.0 %). El resultado cae en el rango de la DFT local/GGA publicada (8.3–8.9×10⁵ m/s); el experimento es mayor porque la DFT semilocal no incluye la renormalización por interacción electrón–electrón.
- No se encontró un valor de v_F específico de PBE en las fuentes consultadas; la comparación es con LDA/GGA-PW91.

## 3. Fonones (OE2) sin cálculo VASP

Todavía no hay fuerzas de VASP. Para ejercitar toda la cadena (carga → constantes de fuerza → bandas → estabilidad) se generaron fuerzas con el Tersoff optimizado de Lindsay y Broido (LAMMPS) y Phonopy en una supercelda 5×5×1 (`examples/fonones_grafeno_tersoff.py`; necesita `lammps` y `mpich` de PyPI). **Esto valida la plomería de la aplicación, no los números de DFT**: el potencial es empírico.

| Comprobación | Resultado |
|---|---|
| a de equilibrio del potencial | 2.492 Å, igual que el publicado para ese conjunto de parámetros |
| Modos acústicos en Γ | 0, 0, 0 THz (ASR cumplida sin corrección) |
| Estabilidad | Estable, f_min = −0.004 THz |
| Rama ZA | f ∝ q² cerca de Γ (0.0215, 0.0987 y 0.405 THz para q = 0.01, 0.02, 0.04) |
| Polarización (verificada con los autovectores) | ZO en Γ fuera del plano; pareja en Γ (50.6 THz = 1689 cm⁻¹) en el plano; pareja ZA/ZO degenerada en K (19.5 THz) fuera del plano |
| Modo en el plano en Γ | 1689 cm⁻¹ frente a ~1580 cm⁻¹ (banda G): el potencial sobreestima ~7 % |
| Rama ZO en Γ | 1303 cm⁻¹: el artículo del potencial advierte que ajusta peor la rama ZO |

Para comparar con la literatura a nivel DFT hacen falta las fuerzas reales de `05_fonones/disp-001`.

## 4. Pruebas con salidas reales públicas (v0.13)

Mientras no hay cálculos propios de otros materiales se probó la aplicación con salidas reales de VASP y Phonopy
que son públicas (repositorios `pymatgen-test-files` y `phonopy`). Sirven para detectar fallos de lectura y de
interpretación, no para validar valores físicos (son cálculos antiguos con parámetros distintos a los del proyecto).

| Caso | Rasgo que ejercita | Resultado |
|---|---|---|
| Si, Al, LiF, Si indirecta (SCF, malla uniforme) | un átomo, metal, aislante | Al: OUTCAR de un átomo rompía la lectura; LiF: Li–Li contado como enlace |
| LiFePO₄ (112 átomos, ISPIN = 2, LDAU) | espín, estructura grande | brecha 0.97 eV en la malla, sin errores |
| SiO₂ (ISPIN = 2) | espín sin momento neto | brecha 5.7 eV en la malla |
| NiO (LSORBIT, OUTCAR + CHGCAR no colineal) | bloques (x), (y), (z) | solo se leía (x): daba 0 µB; ahora ±1.67 µB en z y AFM detectado |
| Fe₃O₄ (CHGCAR con espín) | densidad con espín | 96 e⁻ integrados, vista de magnetización sugerida |
| CdTe 65 átomos (SOC) | SOC con defecto | metal (DOS en E_F alta), sin errores |
| Phonopy: NaCl, Si, Zr₃N₄, TiPN₃ (DFT) | fonones reales | Si: f_max 15.1 THz (lit. ≈ 15.3); NaCl con BORN (LO) estable; TiPN₃: seekpath cambia la celda |

Correcciones hechas a partir de estas pruebas:

- Un OUTCAR con un solo átomo (sin línea de suma) ya se lee.
- Cálculos no colineales (SOC): se leen los tres bloques y se informa el momento local de cada átomo y si se compensan.
- No se muestra «−0.000 µB».
- Entre metales no se cuentan enlaces si el compuesto tiene algún no metal (LiF: coordinación 6, no [6, 18]).
- **Brecha aproximada en la malla k:** si no hay bandas en línea, se calcula la brecha y la clasificación con los
  autovalores del SCF/DOS y se indica que es aproximada. Con el `vasprun.xml` de grafeno de `2-Encuts/600` da
  brecha 0 y «semimetal de Dirac».
- Si seekpath cambia la celda y no hay camino de alta simetría, la estabilidad se evalúa en la malla q (antes fallaba).
- Sin CONTCAR/POSCAR, la estructura se toma del CHGCAR, las bandas, la malla k o la celda de Phonopy, y se avisa.
- Con un metal de NELECT impar y sin momento neto el mensaje ya explica el motivo.

Lo que estas pruebas no pueden confirmar: números del grafeno con DFT propio, SrMnSb₂, rendimiento con CHGCAR grandes.

### Rendimiento con CHGCAR grande

`examples/prueba_rendimiento_chgcar.py N` genera un CHGCAR sintético (8 átomos tipo SrMnSb₂) y mide los pasos. Resultados en este equipo (8 GB de RAM):

| Malla | Puntos | Archivo | Sesión + informe | Pico de memoria | Isosuperficie + dibujo |
|---|---|---|---|---|---|
| 160³ | 4.1 M | 74 MB | < 1 s | 0.2 GB | 1 s |
| 256³ | 16.8 M | 302 MB | 1.9 s | 0.4 GB | 3 s (0.8 GB) |
| 320³ | 32.8 M | 590 MB | 2.6 s | 0.7 GB | 5 s (1.1 GB) |

Un nivel de isosuperficie muy bajo generaba millones de triángulos (6.6 M caras: 27 s y 4.7 GB al dibujar); ahora `isosurface` limita a
1 M de caras reduciendo la malla. Los CHGCAR con espín (×2) y no colineales (×4) necesitan proporcionalmente más memoria; no se
probaron a esa escala.

## 5. Fonones del grafeno calculados con VASP (OE2 / OE5, v0.14)

Carpeta `Step4-Phonopy-Bands-DOS` del usuario: grafeno, supercelda 5×5×1 (50 átomos), a = 2.4683 Å (igual al a₀ de Birch–Murnaghan de su
serie de red), vacío de 20.1 Å, PBE (PAW_PBE C), ENCUT 500 eV, k 4×4×1 centrada en Γ, EDIFF 1e-8, ISYM = 0. Phonopy 2.12 pidió **un único
desplazamiento** (grupo P6/mmm); el cálculo terminó y el SCF convergió. La carpeta contiene las fuerzas (`vasprun.xml`) pero no las
constantes de fuerza: `phonopy -f vasprun.xml` genera `FORCE_SETS` (en Phonopy ≥ 4, `phonopy-init -f`). Las fuerzas netas suman 0.

Frecuencias (cm⁻¹) frente a Maultzsch et al., PRL 92, 075501 (2004), Tabla 1, columna DFT-GGA (arXiv:cond-mat/0404637; cálculos de grafito,
cuyas ramas en el plano son las del grafeno):

| Punto | Rama | Calculado | Lit. GGA | Dif. | Experimento |
|---|---|---|---|---|---|
| Γ | ZO | 871.1 | 884 | −1.5 % | 868 |
| Γ | TO = LO | 1561.7 | 1569 | −0.5 % | 1582 |
| M | ZA | 467.9 | 476 | −1.7 % | – |
| M | TA | 619.6 | 634 | −2.3 % | – |
| M | ZO | 630.9 | 640 | −1.4 % | – |
| M | LA | 1324.8 | 1338 | −1.0 % | – |
| M | LO | 1339.4 | 1346 | −0.5 % | 1323 |
| M | TO | 1393.3 | 1396 | −0.2 % | 1390 |
| K | ZO = ZA | 531.5 | 539 | −1.4 % | – |
| K | TA | 997.1 | 1004 | −0.7 % | – |
| K | LO = LA | 1206.8 | 1221 | −1.2 % | 1194 |
| K | TO | 1345.9 | 1289 | **+4.4 %** | 1265 |

- Acústicos en Γ ≈ 0; no hay frecuencias imaginarias; la ZA es cuadrática cerca de Γ; f_max = 47.8 THz (1595 cm⁻¹).
- Todas las ramas coinciden con la literatura dentro de −2.3 %…−0.2 % (un corrimiento sistemático de ≈ −10 cm⁻¹), **salvo el TO en K**.
- El TO en K está sobreestimado porque K = (1/3, 1/3) no es conmensurable con una supercelda 5×5 (las constantes de fuerza se interpolan) y ese
  modo tiene una anomalía de Kohn (constantes de fuerza de largo alcance). Con 6×6×1 (o 9×9×1) K y M quedan dentro de la malla de la supercelda;
  no se recalculó. El mismo efecto aplana el TO/LO cerca de Γ.
- Reproducible con `examples/comparar_fonones_grafeno.py CARPETA` (figura bandas + DOS con la literatura) y con las pruebas
  `tests/test_datos_publicos.py` (datos en `tests/data/grafeno_5x5`).

Correcciones de la aplicación hechas con este caso:

- El `CONTCAR`/`POSCAR`/`OUTCAR` de una carpeta de Phonopy es de la supercelda desplazada: ya no se toma como «estructura» (antes salía «50 átomos, P1»
  y un aviso falso de que la celda de los fonones no coincidía). Se usa la celda unidad de los fonones, o un CONTCAR propio si existe.
- Esa carpeta ya no impide usar el SCF/DOS de las series de convergencia.
- Si la estructura sale de un punto de la serie de red, los fonones se comparan con el a₀ de Birch–Murnaghan («coincide con a₀ = 2.4683 Å»).
- NELECT se compara por átomo entre celda unidad y supercelda (200 vs 8 era un falso positivo); el `vasprun.xml` junto a `phonopy_disp.yaml` se
  reconoce como el cálculo de desplazamiento y ya no se rotula «Relajación».
- «-0.000 THz» se imprime como 0.000.

Observaciones sobre los archivos de entrada (no afectan a la aplicación): `POSCAR-unitcell` tiene c = 8.14 Å y la celda real c = 20.10 Å; `band.conf` usa
`DIM = 3 3 1` y `mesh.conf` `ATOM_NAME = Fe Se S` y `DIM = 5 5 1` (el primero debe coincidir con `phonopy_disp.yaml` y el segundo con el material); en
`band.conf` conviene escribir 1/3 en lugar de 0.333.

## 6. Qué falta (necesita cálculos)


- Bandas en línea Γ–M–K–Γ con LORBIT = 11 y DOS con NEDOS ≥ 2000.
- Fonones de grafeno con VASP (`05_fonones`).
- Un OUTCAR con espín y todo SrMnSb2.
