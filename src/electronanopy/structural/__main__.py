"""Uso: python -m electronanopy.structural RUTA/POSCAR [--plot]"""
import argparse, json
import matplotlib.pyplot as plt
from . import load_structure, summary, plot_structure


def main():
    ap = argparse.ArgumentParser(description="Resumen estructural de un POSCAR/CONTCAR de VASP")
    ap.add_argument("path"); ap.add_argument("--plot", action="store_true")
    ap.add_argument("--supercell", type=int, nargs=3, default=(1, 1, 1))
    args = ap.parse_args()
    s = load_structure(args.path)
    print(json.dumps(summary(s), indent=2, ensure_ascii=False))
    if args.plot:
        plot_structure(s, tuple(args.supercell)); plt.show()


if __name__ == "__main__":
    main()
