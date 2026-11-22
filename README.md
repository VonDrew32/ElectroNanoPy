# ElectroNanoPy

Plataforma para analizar y visualizar resultados de VASP. Incluye:
- **Módulo estructural (OE1):** lectura de POSCAR/CONTCAR, descriptores y visualización 3D.
- **Módulo de fonones (OE2):** lectura de cálculos de Phonopy, bandas, DOS y criterio de estabilidad dinámica.
- **Interfaz de escritorio (OE4):** ventana en PySide6 que integra los módulos y un resumen automático («visualización inteligente»).
- **Módulo electrónico (OE3):** bandas, DOS, brecha, puntos de Dirac, velocidad de Fermi, pesos orbitales, magnetismo y densidad de carga.

## Instalación (Python >= 3.10)
```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
pytest                           # deben pasar 178 pruebas
```

## Uso
```bash
python -m electronanopy.structural RUTA/A/CONTCAR --plot --supercell 3 3 1
python examples/demo_structural.py grafeno.png
```
```python
from electronanopy.structural import load_structure, summary, plot_structure
s = load_structure("CONTCAR"); print(summary(s)); plot_structure(s, (2, 2, 1))
```

## Estructura
- `src/electronanopy/structural/io.py`: lectura y errores claros.
- `src/electronanopy/structural/descriptors.py`: red, grupo espacial, vacío, enlaces, coordinación, ángulos.
- `src/electronanopy/structural/plotting.py`: figura 3D (matplotlib).
- `tests/`: pruebas con grafeno (valores conocidos).

## Módulo de fonones
Entrada: una carpeta con `phonopy_params.yaml`, o con `phonopy_disp.yaml` + `FORCE_SETS`
(se crea con `phonopy-init -f disp-*/vasprun.xml` en Phonopy ≥ 4, o `phonopy -f ...` en 2 y 3). Si hay un `BORN` se usa la corrección no analítica.
```python
from electronanopy.phonons import (load_phonopy, compute_bands, compute_dos,
                                   check_stability, plot_bands_dos)
ph = load_phonopy("carpeta_phonopy")
bands, dos = compute_bands(ph), compute_dos(ph)
rep = check_stability(ph, bands=bands)       # rep.stable, rep.imaginary, rep.message
plot_bands_dos(bands, dos, rep).savefig("fonones.png")
```
- Camino k: 2D (Γ–M–K–Γ en hexagonal, Γ–X–M–Γ en cuadrada) si la celda tiene vacío; seekpath si es 3D.
- Estabilidad: frecuencia imaginaria si f < -0.05 THz; los 3 modos acústicos se ignoran a menos de
  0.1 1/Å de Γ (ruido numérico). Ambos umbrales son parámetros (`tol`, `gamma_radius`).
- `electronanopy.phonons._synthetic` genera modelos **sintéticos** de resortes para pruebas y demos;
  no representan ningún material real.

Demo: `python examples/demo_phonons.py fonones.png`.

## Módulo electrónico
Entradas: `vasprun.xml` de un cálculo de bandas en modo línea (con su `KPOINTS`) y de un cálculo de DOS; `OUTCAR`
para los momentos magnéticos. Para proyecciones orbitales use `LORBIT = 11` en VASP.
```python
from electronanopy.electronic import (load_band_structure, load_dos, band_gap, dirac_points,
                                      fermi_velocity, classify, plot_bands_dos)
bs = load_band_structure("bandas/vasprun.xml", efermi="scf/vasprun.xml")   # E_F del cálculo SCF
tdos, cdos = load_dos("dos/vasprun.xml")
print(band_gap(bs).message, classify(bs, tdos))
dp = dirac_points(bs)                         # grafeno: punto de Dirac en K
print(fermi_velocity(bs, dp[0]))              # ħv_F en eV·Å y m/s
plot_bands_dos(bs, tdos).savefig("electronico.png")
```
- Todas las energías se grafican respecto a E_F. La brecha es la del cálculo de origen (no se corrige).
- La brecha se decide **contando bandas** con NELECT/2 (no depende de dónde caiga E_F, que en grafeno queda en el punto de Dirac).
- Para que el DOS(E_F) sea fiable, el paso de energía del DOS debe ser ≤ 0.05 eV (`NEDOS ≥ 2000`); si es más grueso, la clasificación lo ignora y avisa.
- Los pesos orbitales se normalizan (VASP solo proyecta dentro de las esferas PAW).
- `examples/validar_oe3_malla_uniforme.py` valida el módulo con una malla uniforme real (sin bandas en línea): reconstruye Γ–M–K–Γ y compara v_F con tight-binding.
- `electronic._synthetic` genera grafeno **tight-binding** (no DFT) para pruebas y demos.

Demo: `python examples/demo_electronic.py electronico.png`.

### Densidad de carga (CHGCAR)
Calcule el SCF con `LCHARG = .TRUE.`. VASP guarda ρ·V_celda; el módulo ya devuelve e/Å^3.
```python
from electronanopy.electronic import (load_chgcar, total_electrons, planar_average,
                                      plot_density_slice, plot_isosurface, density)
chg = load_chgcar("02_scf/CHGCAR")
print(total_electrons(chg))                    # debe igualar los electrones de valencia (grafeno: 8)
plot_density_slice(chg, axis=2, frac=0.5)      # corte en el plano del grafeno
plot_isosurface(chg, 0.4 * density(chg).max())
```
`planar_average`, `density_at` (interpolación periódica) y `net_magnetization` (cálculos con espín) también están disponibles.
Demo: `python examples/demo_density.py densidad.png` (densidad **sintética**, no DFT).

Los pendientes, las preguntas para los asesores y las decisiones abiertas están en `docs/PENDIENTES.md`.

- `electronanopy CARPETA` abre la aplicación con esa carpeta; `electronanopy --resumen CARPETA` imprime el reporte en la terminal sin ventanas.

## Integración entre módulos y estilos de figura (OE4, v0.12)
- El módulo «Integración entre módulos» del Resumen relaciona resultados de distintos módulos con reglas explícitas
  (semimetal 2D con cono de Dirac, estabilidad dinámica + carácter electrónico, posibles causas de una inestabilidad:
  ASR, vacío o supercelda; estructura cargada frente a a₀ de la serie; magnetismo) y propone qué pestaña abrir primero.
- Archivo → *Estilo de exportación de figuras*: pantalla, publicación (1 o 2 columnas) o póster (35 cm de ancho, letra de 24 pt, líneas gruesas).
- Puntos de Dirac: una degeneración en todo el camino (Kramers/espín) ya no se informa como punto de Dirac.
- Enlaces: entre dos metales la tolerancia se limita a 1.0 (no cuenta Sr–Sr o Mn–Mn de segunda capa).
- La vista 3D marca como átomos «fantasma» los vecinos periódicos que quedan fuera de la celda.

## Fonones de grafeno con VASP (v0.14)
- Con el cálculo 5×5×1 del usuario (PBE, ENCUT 500): estable, f_max 47.8 THz y todas las ramas en Γ, M y K dentro de ≈ 2 % de la literatura DFT-GGA (excepto el TO en K, +4 %, por la supercelda no conmensurable); ver `docs/VALIDACION.md`, sección 5, y `examples/comparar_fonones_grafeno.py`.

## Probado con datos públicos reales (v0.13)
- Se probó con salidas reales de VASP/Phonopy de repositorios públicos (Si, Al, LiF, LiFePO₄, SiO₂, NiO con SOC, Fe₃O₄, NaCl, Zr₃N₄…); ver `docs/VALIDACION.md`, sección 4.
- Sin bandas en línea se informa la **brecha aproximada en la malla k** del SCF/DOS (con aviso).
- OUTCAR no colineal (SOC): momentos locales por átomo y detección de ordenamiento antiferromagnético.
- Si seekpath no sirve para la celda de los fonones, la estabilidad se evalúa en la malla q; sin CONTCAR/POSCAR la estructura se toma de otro módulo.

## Validación y documentación
- `docs/VALIDACION.md`: revisión del código, comparación con la literatura (red, módulo elástico 2D, v_F) y fonones de prueba; `examples/comparar_literatura_grafeno.py` la reproduce con tus datos.
- `docs/MANUAL.md`: manual de usuario corto. `packaging/electronanopy.spec`: ejecutable de Windows con PyInstaller.

## Consistencia y exportación (OE4)
- `electronanopy.consistency` lee del `vasprun.xml` de cada cálculo (relajación, SCF, bandas, DOS y un `disp-*` de fonones)
  el corte de energía, el funcional, los pseudopotenciales, ISPIN, NELECT, LSORBIT y LDAU, y avisa si difieren
  (la malla k puede cambiar). También compara la celda de los fonones con la estructura cargada.
  Aparece como el módulo «Consistencia entre cálculos» del Resumen; el umbral es `Thresholds.encut_tol` (1 eV).
- `electronanopy.export`: figuras en PNG/SVG/PDF y reporte en JSON o HTML (con figuras incrustadas).
  En la ventana: menú **Archivo** → *Exportar figura actual* (Ctrl+E), *Exportar reporte* o *Exportar informe con figuras*.

## Interfaz de escritorio (OE4)
```bash
python -m electronanopy.app RUTA/DEL/PROYECTO     # o: electronanopy RUTA/DEL/PROYECTO
```
Abra la carpeta del proyecto (Ctrl+O). La plataforma detecta los archivos (admite la disposición
`01_relax / 02_scf / 03_bandas / 04_dos / 05_fonones` o una carpeta plana; ignora las carpetas `disp-*`), carga cada
módulo por separado (si uno falla, los demás siguen) y muestra:
- **Pestañas** Estructura, Fonones, Electrónico y Densidad, con opciones de visualización (supercelda, tolerancia de
  enlaces, unidades, ventana de energía, peso orbital, vista de la densidad, nivel de isosuperficie, mapa de color).
- **Resumen automático** con resaltados (brecha, punto de Dirac, frecuencias imaginarias, momentos magnéticos),
  **advertencias** de consistencia (E_F del SCF vs. bandas, electrones integrados vs. NELECT, tamaño de la supercelda,
  vacío en 2D) y **sugerencias** de qué graficar.

**Series de cálculos.** Si la carpeta contiene barridos (malla k, ENCUT, parámetro de red, vacío) en subcarpetas
numéricas hermanas (p. ej. `3-Lattice_Parameter/2.346 … 2.566`), la plataforma las agrupa en una serie, usa un único
cálculo coherente (el de menor energía en la serie de red; el último en las demás) y en la pestaña Estructura permite
elegir la serie y el cálculo, y ver la energía relativa frente al parámetro con la tabla (a, c, enlace, vacío, E). El
reporte avisa si una serie no está convergida (≤ 5 meV entre los tres últimos puntos) o no es monótona. En la serie de red ajusta Birch–Murnaghan
(a₀, V₀, E₀ y, en 3D, B₀) y lo compara con el a donde σ_xx = 0.

**Fonones.** Además de bandas, DOS y estabilidad, el reporte comprueba la calidad de la entrada: los 3 modos acústicos en Γ
deben ser ≈ 0 (regla de suma acústica; umbral 0.3 THz) y la supercelda debe tener un vector ≥ 10 Å. En la pestaña
Fonones se puede elegir el camino (automático, 2D o 3D) e imponer la regla de suma acústica sobre una copia.

«Visualización inteligente» = Umbrales en
`electronanopy.report.Thresholds`.

Capas: `project.py` (detección y carga, sin interfaz) → `report.py` (reglas) → `app/` (PySide6). Las dos primeras se
pueden usar sin Qt.

## Video de demostración

https://youtu.be/oD1LTm20JU4
