#!/usr/bin/env python3
"""RH B5 refinement with algebraically equivalent time-scaled stage equations.
No unconverged steps are accepted. Store the scaling and stopping threshold.
"""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[k]='1'
import sys,json,time,argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.sparse import diags
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'evidence/absCSS_B5_css_rh_v1/B5_css_rh_v1/code'))
import B5_css_rh as b5
Original=b5.M6RHSparse
class ScaledRH(Original):
    def residual_jac(self,Z,W,dt,need_jac=True):
        v=super().residual_jac(Z,W,dt,need_jac)
        weights=np.repeat([dt/4.,2.*dt/3.],self.o.N)
        if need_jac:
            E,J=v;return weights*E,diags(weights)@J
        return weights*v
b5.M6RHSparse=ScaledRH

def work(spec):
    eta,m,tol,ratio=spec
    out=ROOT/'repro/results/b5'
    path=out/f'RH_scaled_eta{eta:g}_m{m}_tol{tol:g}_dt{ratio:g}.json'
    if path.exists() and os.environ.get('ABSCSS_FORCE_PDE')!='1':
        old=json.loads(path.read_text())
        if old.get('status')=='OK':return path.name,old['status'],old['runtime_s'],'reused'
    stats={'max_accepted_scaled_residual':0.,'failed_residual_history':[]}
    orig_step=b5.newton_step
    def step(*args,**kw):
        ans=orig_step(*args,**kw);uu,ok,hist=ans
        if ok and hist: stats['max_accepted_scaled_residual']=max(stats['max_accepted_scaled_residual'],hist[-1])
        elif not ok:stats['failed_residual_history']=hist
        return ans
    b5.newton_step=step
    data=b5.run_chase('RH',1.5,m,T=1.2,dt_ratio=ratio,eta=eta,periods=12,cfast=3.,cslow=1.,separation=8.,tol=tol,maxit=30,c0=1.,snap_times=[0,1.2])
    b5.newton_step=orig_step
    data.update(stats);data['newton_tolerance']=tol;data['newton_iteration_cap']=30
    data['residual_scaling']='E1 -> (dt/4) E1; E2 -> (2dt/3) E2; same row scaling in Jacobian'
    data['revision_purpose']='Additional fixed-viscosity spatial and solver-tolerance controls'
    path.write_text(json.dumps(data));return path.name,data['status'],data['runtime_s'],'executed'
if __name__=='__main__':
    specs=[(.003,190,1e-11,.05),(.001,225,1e-11,.05),(.0003,300,1e-11,.05),(.0003,450,1e-11,.05),(.001,125,1e-11,.05),(.001,225,1e-10,.05)]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for f in as_completed([pool.submit(work,s) for s in specs]):print(f.result(),flush=True)
