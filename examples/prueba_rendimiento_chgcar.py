"""Mide tiempo y memoria al leer y analizar un CHGCAR grande (sintético, con el formato de VASP).

Uso:  python examples/prueba_rendimiento_chgcar.py N [carpeta]      (malla N×N×N; N = 256 → ~300 MB de texto)
Escribe ``carpeta/CHGCAR`` (por defecto un directorio temporal), lo lee con ``load_chgcar`` y con la sesión completa
y muestra segundos y pico de memoria (MB) de cada paso.
"""
import resource
import sys
import tempfile
import time
from pathlib import Path

import numpy as np


def escribir_chgcar(ruta: Path, n: int) -> None:
    a, c = 4.4, 9.0   # vol vale a·a·c                                               # celda tetragonal de 8 átomos (tipo SrMnSb2)
    pos = [[0, 0, .33], [.5, .5, .67], [0, 0, .0], [.5, .5, 1.0 - .0], [0, .5, .1], [.5, 0, .9], [0, .5, .5], [.5, 0, .5]]
    with open(ruta, "w") as f:
        f.write("SrMnSb2 sintetico\n1.0\n")
        f.write(f"{a} 0 0\n0 {a} 0\n0 0 {c}\nSr Mn Sb\n2 2 4\nDirect\n")
        for p in pos:
            f.write("  " + " ".join(f"{x % 1.0:.6f}" for x in p) + "\n")
        f.write(f"\n  {n}  {n}  {n}\n")
        x = np.linspace(0, 1, n, endpoint=False)
        rng = np.random.default_rng(0)
        for k in range(n):                                        # una capa z (índice más lento = z) cada vez
            # CHGCAR guarda ρ·V_celda; tres máximos gaussianos por capa más un ruido pequeño, independiente de N
            capa = sum(np.exp(-(((x[:, None] - cx + .5) % 1 - .5) ** 2 + ((x[None, :] - cy + .5) % 1 - .5) ** 2) * 60)
                       for cx, cy in ((0, 0), (.5, .5), (0, .5)))
            capa = (capa * 30 * (1 + .3 * np.cos(2 * np.pi * x[k])) + rng.random((n, n)) * 1e-2).ravel()
            f.write("\n".join(" ".join(f"{v:.11E}" for v in capa[i:i + 5]) for i in range(0, capa.size, 5)) + "\n")
        f.write(f"augmentation occupancies   1  0\n")              # relleno al final como en VASP


def pico_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 128
    carpeta = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(tempfile.mkdtemp())
    carpeta.mkdir(parents=True, exist_ok=True)
    ruta = carpeta / "CHGCAR"
    t = time.time()
    escribir_chgcar(ruta, n)
    print(f"CHGCAR {n}^3 = {n ** 3 / 1e6:.1f} M puntos, {ruta.stat().st_size / 1e6:.0f} MB (escrito en {time.time() - t:.0f} s)")
    from electronanopy.electronic import load_chgcar, total_electrons, planar_average
    from electronanopy.electronic.density import slice_2d, isosurface
    t = time.time()
    chg = load_chgcar(ruta)
    print(f"load_chgcar      {time.time() - t:6.1f} s   pico {pico_mb():7.0f} MB")
    for nombre, fn in (("total_electrons", lambda: total_electrons(chg)),
                       ("planar_average", lambda: planar_average(chg, 2)),
                       ("slice_2d", lambda: slice_2d(chg, 2, .5)),
                       ("isosurface", lambda: isosurface(chg, 0.1 * float(chg.data["total"].max()) / chg.structure.volume))):
        t = time.time()
        fn()
        print(f"{nombre:16s} {time.time() - t:6.1f} s   pico {pico_mb():7.0f} MB")
