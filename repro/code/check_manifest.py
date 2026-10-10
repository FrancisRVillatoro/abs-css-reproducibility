#!/usr/bin/env python3
"""Check the distributed content manifest before regenerating outputs."""
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parents[2]
def main():
    failures=[];count=0
    for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
        digest,name=line.split('  ',1);p=ROOT/name;count+=1
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:failures.append(name)
    if failures:raise RuntimeError('Checksum mismatch: '+', '.join(failures))
    print(f'PASS: {count} files match the distributed SHA-256 manifest.')
if __name__=='__main__':main()
