# -*- mode: python ; coding: utf-8 -*-
# OligoForge desktop build.
#   pip install -r requirements.txt -r requirements-build.txt
#   pyinstaller oligoforge.spec
# Output: dist/OligoForge  (on Windows: dist/OligoForge.exe) — one double-click file.
from PyInstaller.utils.hooks import collect_submodules, collect_data_files

hiddenimports = (
    collect_submodules("uvicorn")
    + collect_submodules("primer3")
    + ["app", "anyio", "click", "h11", "httptools", "websockets"]
    + ["Bio.Entrez", "Bio.SeqIO", "Bio.Seq", "Bio.SeqRecord"]
)

datas = [("static", "static")]
datas += collect_data_files("primer3")          # primer3 thermodynamic parameter files
datas += collect_data_files("certifi")          # CA bundle for NCBI HTTPS

a = Analysis(
    ["launcher.py"],
    pathex=["."],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=["tkinter", "matplotlib", "pytest", "pandas", "webview", "IPython"],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="OligoForge",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    runtime_tmpdir=None,
    console=True,                 # visible "OligoForge is running" window; close it to stop
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
