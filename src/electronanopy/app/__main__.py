"""Uso: python -m electronanopy.app [CARPETA_DEL_PROYECTO]"""
import sys
from .main_window import run

if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else None)
