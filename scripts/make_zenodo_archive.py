#!/usr/bin/env python3
"""Create a reproducible Zenodo manual-upload ZIP from this checkout, excluding .git and local generated outputs."""
from pathlib import Path
import sys,zipfile
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT.parent/"abs-css-reproducibility-v1.0.0-zenodo.zip"
SKIP={'.git','.venv','venv','outputs','__pycache__'}
if OUT.resolve().is_relative_to(ROOT.resolve()): raise RuntimeError('Archive output must be outside the repository')
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and not set(p.relative_to(ROOT).parts)&SKIP)
with zipfile.ZipFile(OUT,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9,allowZip64=True) as z:
 for p in files:
  rel=(ROOT.name+'/'+p.relative_to(ROOT).as_posix()); info=zipfile.ZipInfo(rel,date_time=(2026,10,10,0,0,0)); info.compress_type=zipfile.ZIP_DEFLATED
  info.external_attr=(0o100644<<16);z.writestr(info,p.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
print(OUT)
