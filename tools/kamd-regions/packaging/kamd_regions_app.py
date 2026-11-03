# Entry point of the frozen desktop app (PyInstaller).
import multiprocessing
import sys

from kamd_regions.desktop import main

if __name__ == "__main__":
    multiprocessing.freeze_support()
    sys.exit(main())
