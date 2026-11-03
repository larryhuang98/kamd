# PyInstaller spec for the KAMD Regions desktop app.
#   cd tools/kamd-regions/packaging && pyinstaller --noconfirm kamd-regions.spec
# Output: dist/KAMD Regions/ (Linux, Windows) or dist/KAMD Regions.app (macOS).
import sys
from PyInstaller.utils.hooks import collect_all, collect_submodules

datas, binaries, hiddenimports = [], [], []
for pkg in ("kamd_regions", "mdtraj", "parmed", "igraph", "leidenalg", "webview"):
    d, b, h = collect_all(pkg)
    datas += d
    binaries += b
    hiddenimports += h
hiddenimports += collect_submodules("uvicorn") + collect_submodules("kamd_openmm")
hiddenimports += ["scipy.sparse.csgraph", "scipy.spatial", "scipy.linalg"]
if sys.platform.startswith("linux"):
    hiddenimports += ["qtpy", "PySide6.QtWebEngineWidgets", "PySide6.QtWebChannel"]
    excludes = []
else:  # macOS and Windows use the system web view; keep Qt out of the bundle
    excludes = ["PySide6", "PyQt5", "PyQt6", "qtpy"]
excludes += ["tkinter", "matplotlib", "IPython", "pytest", "openmm", "MDAnalysis"]

a = Analysis(["kamd_regions_app.py"], pathex=["..", "../../../openmm"], datas=datas, binaries=binaries,
             hiddenimports=hiddenimports, excludes=excludes, noarchive=False)

# Qt parts the app never uses (Linux bundle); Qt WebEngine itself must stay
UNUSED = ("Quick3D", "QuickControls2Imagine", "QuickControls2Material", "QuickControls2Universal",
          "QuickControls2Fusion", "QuickControls2FluentWinUI3", "QuickControls2Basic", "Qt63D", "Qt3D",
          "Multimedia", "Charts", "Graphs", "DataVisualization", "Sql", "Designer", "Location", "Sensors",
          "Bluetooth", "Nfc", "SerialPort", "SerialBus", "SpatialAudio", "TextToSpeech", "Scxml",
          "RemoteObjects", "VirtualKeyboard", "Lottie", "QuickTimeline", "ShaderTools", "Help", "Test",
          "assetdownloader", "Quick/Dialogs", "VectorImage")


def keep(entry):
    path = entry[0].replace("\\", "/")
    if "/translations/" in path and "qtwebengine_locales" not in path:
        return False
    if "qtwebengine_locales/" in path and not path.endswith(("en-US.pak", "en-GB.pak")):
        return False
    return not any(u in path for u in UNUSED)


a.binaries = [b for b in a.binaries if keep(b)]
a.datas = [d for d in a.datas if keep(d)]
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="kamd-regions-app", console=False,
          argv_emulation=sys.platform == "darwin")
coll = COLLECT(exe, a.binaries, a.datas, name="KAMD Regions")
if sys.platform == "darwin":
    app = BUNDLE(coll, name="KAMD Regions.app", bundle_identifier="edu.uci.kamd-regions",
                 info_plist={"NSHighResolutionCapable": True, "CFBundleShortVersionString": "0.1.0"})
