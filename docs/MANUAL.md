# Manual de usuario de ElectroNanoPy

## 1. Qué necesita
Carpetas con los resultados de VASP (y de Phonopy para fonones). La aplicación **solo lee**: no ejecuta cálculos ni toca el clúster.

| Pestaña | Archivos que usa |
|---|---|
| Estructura | `CONTCAR` (o `POSCAR`); para series: carpetas hermanas con `CONTCAR`, `OSZICAR`, `OUTCAR` |
| Fonones | `phonopy_params.yaml`, o `phonopy_disp.yaml` + `FORCE_SETS` |
| Electrónico | `vasprun.xml` de bandas (KPOINTS en modo línea) y de SCF/DOS; `OUTCAR` para el magnetismo. Use `LORBIT = 11` para los pesos orbitales y `NEDOS ≥ 2000` en la DOS |
| Densidad | `CHGCAR` del SCF (`LCHARG = .TRUE.`) |

## 2. Qué esperan los análisis de cada cálculo
La aplicación no genera entradas; esto es lo que necesita encontrar para que cada resultado sea fiable.

| Módulo | Qué debe tener el cálculo | Si falta |
|---|---|---|
| Estructura | `CONTCAR` de una relajación convergida. Para una serie: al menos tres carpetas hermanas con valor numérico en el nombre (corte de energía, parámetro de red, malla k, vacío) y `OSZICAR` con la energía | Sin `CONTCAR` la estructura se toma de otro módulo y se avisa |
| Fonones | Celda relajada, la misma de la estructura; supercelda de al menos 3×3 (2D: 3×3×1); `phonopy_params.yaml` o `phonopy_disp.yaml` + `FORCE_SETS`; `BORN` si el material es polar; en 2D, vacío suficiente | Se avisa si la supercelda es pequeña o la celda difiere de la estructura |
| Bandas | `vasprun.xml` con `KPOINTS` en modo línea (grafeno: Γ–M–K–Γ) y `LORBIT = 11` para los pesos orbitales | Sin bandas en línea solo hay una brecha **aproximada** en la malla k del SCF/DOS |
| DOS | `vasprun.xml` del SCF con malla centrada en Γ (en grafeno, múltiplo de 3 para incluir K) y `NEDOS ≥ 2000` | Con paso de energía > 0.05 eV se avisa que la DOS en E_F no es fiable |
| Magnetismo / SOC | `ISPIN = 2` (o `LSORBIT = .TRUE.`) y el `OUTCAR` del mismo cálculo | Sin `OUTCAR` no se informan momentos; con SOC se leen los tres componentes |
| Densidad | `CHGCAR` (`LCHARG = .TRUE.`); con `ISPIN = 2` trae también la magnetización | Se avisa si los electrones integrados no coinciden con `NELECT` |
| Consistencia | Mismo `ENCUT`, funcional y pseudopotenciales en los cálculos que se comparan | Se avisan las diferencias (la malla k puede cambiar) |

Rendimiento: se probó con un `CHGCAR` sintético de 590 MB (33 millones de puntos, 320³). Abrir la carpeta y generar el informe tarda
unos 3 s y usa unos 0.7 GB de memoria; dibujar la isosuperficie tarda unos 5 s y llega a 1.1 GB. Si el nivel produce más de un millón
de triángulos, la malla se reduce automáticamente. Un `CHGCAR` con espín ocupa el doble. Reproducible con
`examples/prueba_rendimiento_chgcar.py N`.

## 3. Abrir un proyecto
1. Active el entorno (`.venv\Scripts\Activate.ps1`) y ejecute `electronanopy`.
2. **Archivo → Abrir carpeta** y elija la carpeta del proyecto: se examinan hasta 4 niveles.
3. El panel «Archivos detectados» muestra qué se encontró y con qué archivo se trabaja; el «Resumen» muestra los resultados y advertencias. Si cierra un panel, **Ver** lo reabre.

Las pistas del nombre de la carpeta (dentro del proyecto) deciden qué es cada cálculo: `relax`, `scf`, `band`, `dos`, `phonon`/`fonon`. Tres o más carpetas hermanas con nombre numérico (`2-Encuts/400…700`, `3-Lattice_Parameter/2.40…`) se leen como una **serie** de convergencia.

## 4. Cada pestaña
- **Estructura**: vista (superior/lateral/isométrica), distancias de enlace, vacío en 2D y grupo espacial. Con «Serie» se ve la convergencia o, en la serie de red, la curva de Birch–Murnaghan: *a₀* es el parámetro para el siguiente cálculo; B₀ solo se da en 3D con la celda escalada por igual.
- **Fonones**: bandas y DOS, frecuencias imaginarias (estabilidad) y calidad (ASR en Γ, tamaño de la supercelda). «ASR» aplica la regla de suma acústica a una copia.
- **Electrónico**: bandas (respecto a E_F) y DOS, brecha (si no hay bandas en línea, una estimación en la malla k del SCF/DOS, marcada como aproximada), puntos de Dirac, ħv_F, proyección orbital y magnetismo.
- **Densidad**: corte 2D, promedio planar e isosuperficie de ρ(r); con espín, la magnetización («diff»).

## 5. Integración entre módulos
El Resumen termina con «Integración entre módulos»: reglas que combinan los resultados (por ejemplo «semimetal 2D con cono de Dirac», o «frecuencias imaginarias + modos acústicos en Γ distintos de cero: puede ser ruido, active ASR»). Al cargar un proyecto se abre la pestaña que esas reglas consideran más relevante (fonones si hay inestabilidad; electrónico si hay bandas; estructura en otro caso); el motivo aparece en la barra de estado.

## 6. Exportar
**Archivo → Estilo de exportación de figuras** elige el tamaño: *Pantalla*, *Publicación 1 columna* (8.5 cm), *Publicación 2 columnas* (17 cm) o *Póster* (35 cm de ancho, letra de 24 pt, líneas gruesas, para figuras de un póster). **Archivo → Exportar figura actual** (Ctrl+E: PNG, SVG o PDF) usa ese estilo; **Exportar reporte** (JSON o HTML) o **Exportar informe con figuras** (un HTML autocontenido).

## 7. Cómo leer las advertencias
Cada advertencia dice qué se vio y cuándo es importante. Las de «Consistencia entre cálculos» indican que el SCF, las bandas o los fonones se calcularon con distinto corte de energía, funcional o pseudopotenciales; en ese caso los resultados no son comparables. Los umbrales son valores de partida y están pendientes de acordar con los asesores (`docs/PENDIENTES.md`).

## 8. Programa de instalación (opcional)
`packaging/electronanopy.spec` construye una carpeta con `ElectroNanoPy.exe` (hay que construir en Windows): ver el encabezado del archivo.
