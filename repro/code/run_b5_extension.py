#!/usr/bin/env python3
"""Additional fixed-viscosity B5 levels; original production equations unchanged."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import argparse, json, sys, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'evidence/absCSS_B5_css_rh_v1/B5_css_rh_v1/code'
sys.path.insert(0,str(BASE))
from B5_css_rh import run_chase

def job(args):
    model,eta,m,outdir=args
    outfile=Path(outdir)/f'{model}_eta{eta:g}_m{m}.json'
    if outfile.exists() and os.environ.get('ABSCSS_FORCE_PDE')!='1':
        old=json.loads(outfile.read_text())
        if old.get('status')=='OK': return outfile.name,old['status'],old.get('runtime_s'), 'reused'
    data=run_chase(model,1.5,m,T=1.2,dt_ratio=.05,eta=eta,phase=.5,
        periods=12,cfast=3.,cslow=1.,separation=8.,opposite=True,
        tol=1e-10,maxit=20,c0=1.0,snap_times=[0,1.2])
    data['revision_purpose']='Independent extra spatial level at fixed eta'
    data['newton_tolerance']=1e-10
    data['newton_iteration_cap']=20
    temp=outfile.with_suffix('.tmp');temp.write_text(json.dumps(data));temp.replace(outfile)
    return outfile.name,data['status'],data.get('runtime_s'), 'executed'

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=3);a=ap.parse_args()
    out=ROOT/'repro/results/b5';out.mkdir(parents=True,exist_ok=True)
    specs=[(0.003,190),(0.001,225),(0.0003,300),(0.0003,450)]
    tasks=[(mod,eta,m,str(out)) for eta,m in specs for mod in ('CSS','RH')]
    start=time.time(); records=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(job,t) for t in tasks]):
            rec=f.result(); records.append(rec);print(rec,flush=True)
    (out/'execution_summary.json').write_text(json.dumps({'runs':records,'wall_s':time.time()-start},indent=2))
