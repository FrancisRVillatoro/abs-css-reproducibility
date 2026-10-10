#!/usr/bin/env python3
import sys, math, json, time
from pathlib import Path
import numpy as np
ROOT=Path('/mnt/data/absCSS_B4_work')
sys.path.insert(0,str(ROOT/'B2src/B2_dipole_v1/code'))
sys.path.insert(0,str(ROOT/'B3src/B3_transverse_v1/code'))
from B1_pade_sparse_newton import PadeSparseOps, M1Sparse, M2Sparse, compact_u, compact_v
from B2_dipole_solver import newton_stage_adaptive, dipole_u, dipole_v, grid_from_m
from B3_regularized_solver import M34Regularized

def rel_l2(a,b): return float(np.linalg.norm(a-b)/np.linalg.norm(b))

def make_initial(kind,x,p,alpha,phase,dx):
    if kind=='compacton':
        x0=phase*dx
        return compact_u(x,p,alpha,x0=x0), compact_v(x,p,alpha,x0=x0)
    if kind=='dipole':
        xmid=phase*dx
        return dipole_u(x,p,alpha,xmid=xmid), dipole_v(x,p,alpha,xmid=xmid)
    raise ValueError(kind)

def run(kind,model,p,m,T=.2,dt_ratio=.1,eps=.06,eta=0.0,mu2=.1,phase=.5,tol=1e-10,refresh=2,maxit=20):
    N,dx,L=grid_from_m(m); x=-L/2+dx*np.arange(N); alpha=p*p/((p+1)*(p+2)); op=PadeSparseOps(N,dx)
    u0,v0=make_initial(kind,x,p,alpha,phase,dx)
    if model in ('M1','M5'):
        U=u0.copy(); M=M1Sparse(op,p,alpha,1.0,eta if model=='M5' else 0.0); recon=lambda z:z
    elif model in ('M2','M2reg'):
        U=v0.copy(); M=M2Sparse(op,p,alpha,1.0,0.0 if model=='M2' else mu2*dx*dx); recon=lambda z:np.sign(z)*np.abs(z)**(2/p)
    elif model in ('M3','M4'):
        U=u0.copy(); M=M34Regularized(op,p,alpha,1 if model=='M3' else 2,eps,c0=1.0); recon=lambda z:z
    else: raise ValueError(model)
    steps=max(1,int(math.ceil(T/(dt_ratio*dx)))); dt=T/steps
    nit=[]; t0=time.time(); fail=None
    for n in range(steps):
        U,ok,hist=newton_stage_adaptive(U,dt,M,tol=tol,maxit=maxit,refresh=refresh,line_search=10)
        nit.append(len(hist))
        if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U))>1e8:
            fail=n+1; break
    u=recon(U); tr=(fail-1)*dt if fail else T
    out=dict(status='FAIL' if fail else 'OK',kind=kind,model=model,p=p,m=m,N=N,L=L,dx=dx,T=T,T_reached=tr,dt=dt,dt_ratio=dt/dx,
             eps=eps if model in ('M3','M4') else None,eta=eta if model=='M5' else None,mu2=mu2 if model=='M2reg' else None,
             phase=phase,mean_newton=float(np.mean(nit)) if nit else None,max_newton=max(nit) if nit else None,
             shape_relL2=rel_l2(u,u0), x=x.tolist(),u=u.tolist(),u0=u0.tolist())
    return out

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--kind',choices=['compacton','dipole'],required=True); ap.add_argument('--model',required=True); ap.add_argument('--p',type=float,default=1.5); ap.add_argument('--m',type=int,required=True)
    ap.add_argument('--T',type=float,default=.2); ap.add_argument('--dt-ratio',type=float,default=.1); ap.add_argument('--eps',type=float,default=.06); ap.add_argument('--eta',type=float,default=0.0); ap.add_argument('--mu2',type=float,default=.1); ap.add_argument('--phase',type=float,default=.5); ap.add_argument('--tol',type=float,default=1e-10); ap.add_argument('--refresh',type=int,default=2); ap.add_argument('--maxit',type=int,default=20); ap.add_argument('--out')
    a=ap.parse_args(); r=run(a.kind,a.model,a.p,a.m,a.T,a.dt_ratio,a.eps,a.eta,a.mu2,a.phase,a.tol,a.refresh,a.maxit)
    if a.out: Path(a.out).parent.mkdir(parents=True,exist_ok=True); Path(a.out).write_text(json.dumps(r))
    print(json.dumps({k:v for k,v in r.items() if k not in ('x','u','u0')},indent=2))
