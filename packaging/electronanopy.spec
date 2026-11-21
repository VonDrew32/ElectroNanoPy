# -*- mode: python ; coding: utf-8 -*-
# Construir (en Windows, desde la carpeta del proyecto, con el entorno virtual activado):
#     python -m pip install pyinstaller
#     pyinstaller packaging/electronanopy.spec --noconfirm
# El resultado queda en dist/ElectroNanoPy/ (carpeta con ElectroNanoPy.exe, que abre la ventana, y ElectroNanoPy-consola.exe,
# que imprime en la terminal: ElectroNanoPy-consola.exe --resumen CARPETA). Hay que construir en cada sistema operativo.
# Guía paso a paso: packaging/CONSTRUIR_WINDOWS.md
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

datas, hidden = [], []
for pkg in ("pymatgen", "phonopy", "seekpath", "spglib", "matplotlib", "skimage", "scipy"):
    datas += collect_data_files(pkg)
for pkg in ("pymatgen.io.vasp", "pymatgen.electronic_structure", "pymatgen.symmetry", "phonopy", "skimage.measure",
            "scipy.special", "scipy.optimize"):
    hidden += collect_submodules(pkg)

a = Analysis(["lanzador.py"], pathex=["../src"], datas=datas, hiddenimports=hidden,
             excludes=["tkinter", "pytest", "IPython", "notebook", "PyQt5", "PyQt6", "PySide2"])
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="ElectroNanoPy", console=False)
exe_consola = EXE(pyz, a.scripts, [], exclude_binaries=True, name="ElectroNanoPy-consola", console=True)
coll = COLLECT(exe, exe_consola, a.binaries, a.datas, name="ElectroNanoPy")
