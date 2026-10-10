#!/usr/bin/env python3
"""Offline release audit: prohibit publication/manuscript files and validate hashes.

Run BEFORE reproduce.py (which deliberately regenerates some tracked results).
"""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
BAD_EXT={'.tex','.bib','.pdf','.png','.jpeg','.jpg','.doc','.docx','.ppt','.pptx','.rtf','.odt','.html','.htm'}
BAD_PATH_SEGMENTS={'manuscript','editorial','submission','papers','paper','articles','article'}
IGNORE={'.git','.venv','venv','outputs','__pycache__'}

def main():
    paths=[p for p in ROOT.rglob('*') if p.is_file() and not set(p.relative_to(ROOT).parts)&IGNORE]
    rejected=[p.relative_to(ROOT).as_posix() for p in paths if p.suffix.lower() in BAD_EXT or set(p.relative_to(ROOT).parts[:-1])&BAD_PATH_SEGMENTS]
    if rejected:raise RuntimeError('Forbidden files: '+', '.join(rejected[:15]))
    for name in ('CITATION.cff','.zenodo.json','RELEASE.json','SOURCE_PROVENANCE.json','ORIGINAL_ARCHIVE_SHA256.json','EXECUTION_LEDGER.json'):
        if not (ROOT/name).is_file():raise RuntimeError('Required metadata missing: '+name)
    metadata=json.loads((ROOT/'.zenodo.json').read_text())
    if metadata.get('upload_type')!='software' or metadata.get('version')!='1.0.0':raise RuntimeError('Unexpected Zenodo record metadata')
    ledger=json.loads((ROOT/'EXECUTION_LEDGER.json').read_text())
    checked=0
    for rec in ledger['new_production_outputs']:
        p=ROOT/rec['file']
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=rec['sha256']:
            raise RuntimeError('New calculation ledger mismatch: '+rec['file'])
        checked+=1
    source=json.loads((ROOT/'SOURCE_PROVENANCE.json').read_text())
    if any(not (ROOT/p).exists() for p in ('repro/code/reproduce.py','repro/inputs/direct_mesh.csv','repro/results/spectral/direct_all.csv')):
        raise RuntimeError('Missing essential reproducibility inputs')
    expected_manifest={}
    for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
        sha,name=line.split('  ',1);expected_manifest[name]=sha
    matched=0
    for p in paths:
        name=p.relative_to(ROOT).as_posix()
        if name in ('SHA256SUMS.txt',):continue
        if name not in expected_manifest:raise RuntimeError('File not in release manifest: '+name)
        if hashlib.sha256(p.read_bytes()).hexdigest()!=expected_manifest[name]:
            raise RuntimeError('Bad hash: '+name)
        matched+=1
    if matched!=len(expected_manifest):raise RuntimeError(f'Manifest has {len(expected_manifest)-matched} missing files')
    print(f'PASS: {matched} release files match SHA256SUMS.txt; {checked} new-output ledger entries verified; no paper/figure/PDF/LaTeX files.')
if __name__=='__main__':main()
