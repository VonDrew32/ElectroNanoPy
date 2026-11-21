# Construir y probar ElectroNanoPy.exe en Windows

Un `.exe` solo se puede construir en Windows. Estos pasos toman unos 15–20 minutos (la mayor parte es la descarga de paquetes y la compilación).

## 1. Preparar
1. Instale **Python 3.12 o 3.13 de 64 bits** desde python.org. En el instalador marque **«Add python.exe to PATH»**.
2. Descomprima `electronanopy_v1.0.0.zip`, por ejemplo en `C:\electronanopy`.
3. Abra **PowerShell** en esa carpeta (`C:\electronanopy\electronanopy`) y ejecute, línea por línea:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1          # si falla por permisos: Set-ExecutionPolicy -Scope Process Bypass
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pip install pyinstaller
```

## 2. Comprobar que funciona con Python (antes de empaquetar)
```powershell
python -m pytest -q                 # esperado: 177 passed
electronanopy tests\data\grafeno_5x5            # abre la ventana con los fonones del grafeno
electronanopy --resumen tests\data\grafeno_5x5  # imprime el informe en la terminal
```

## 3. Construir
```powershell
pyinstaller packaging/electronanopy.spec --noconfirm
```
Queda en `dist\ElectroNanoPy\`: **ElectroNanoPy.exe** (abre la ventana) y **ElectroNanoPy-consola.exe** (imprime en la terminal).
Se distribuye la carpeta completa (puede comprimirla en un `.zip`).

## 4. Lista de verificación del `.exe`
Haga esta prueba **en el equipo donde se construyó** y, si puede, en otro Windows **sin Python instalado**.

| # | Prueba | Resultado esperado |
|---|---|---|
| 1 | Doble clic en `ElectroNanoPy.exe` | Se abre la ventana en menos de ~15 s, sin ventana negra |
| 2 | Archivo → Abrir carpeta → `tests\data\grafeno_5x5` (copiada junto al `.exe`) | Pestaña Fonones con bandas y DOS; «estable» |
| 3 | Abrir su carpeta de grafeno (`Graphene`) | Estructura, series, Electrónico (brecha aproximada) y Densidad se dibujan |
| 4 | Archivo → Exportar figura (Ctrl+E) en estilo «Publicación» | Se crea el PNG |
| 5 | Archivo → Exportar informe completo | Se crea el HTML/JSON y se abre bien |
| 6 | En PowerShell: `.\ElectroNanoPy-consola.exe --resumen tests\data\grafeno_5x5` | Imprime el informe |
| 7 | Cerrar y abrir de nuevo | No deja procesos colgados en el Administrador de tareas |

## 5. Si algo falla
- **No abre / se cierra solo:** ejecute `ElectroNanoPy-consola.exe` desde PowerShell; el error aparece en la terminal. Envíe ese texto.
- **«No module named …»:** falta un módulo en el empaquetado; envíe el nombre y se agrega a `hiddenimports` en `electronanopy.spec`.
- **El antivirus lo bloquea:** es habitual con ejecutables nuevos sin firmar; agréguelo a las excepciones o use la instalación con Python (paso 2).

## Alternativa: construir en GitHub (sin instalar nada)
Si el proyecto está en GitHub, `.github/workflows/windows.yml` construye el `.exe` en un Windows de GitHub y lo deja como archivo descargable (pestaña *Actions* → ejecución → *Artifacts*). Se lanza solo, o a mano con *Run workflow*. No se ha ejecutado todavía.
