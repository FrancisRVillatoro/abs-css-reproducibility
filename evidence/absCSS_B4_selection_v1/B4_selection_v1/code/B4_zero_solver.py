#!/usr/bin/env python3
import sys, math, json, time
from pathlib import Path
import numpy as np
ROOT=Path('/mnt/data/absCSS_B4_work')
sys.path.insert(0,str(ROOT/'B2src/B2_dipole_v1/code'))
from B1_pade_sparse_newton import PadeSparseOps, M1Sparse, M2Sparse
from B2_dipole_solver import newton_stage_adaptive

def sine_u(x,L,a=1.0,shift=0.0):
    k=2*np.pi/L
    return (a/k)*np.sin(k*(x-shift))

def periodic_interp(x_src,u_src,x_dst,L):
    # periodic linear interpolation sufficient for diagnostics on nested-ish grids
    xs=np.r_[x_src, x_src[0]+L]
    us=np.r_[u_src, u_src[0]]
    xd=((x_dst-x_src[0])%L)+x_src[0]
    return np.interp(xd,xs,us)

def run(model,p,N,L=8.0,T=0.01,dt_ratio=0.1,eta=0.0,mu2=0.0,tol=1e-10,shift=0.0,refresh=2,maxit=20):
    dx=L/N; x=-L/2+dx*np.arange(N); alpha=p*p/((p+1)*(p+2)); op=PadeSparseOps(N,dx)
    u0=sine_u(x,L,shift=shift)
    if model in ('M1','M5'):
        U=u0.copy(); M=M1Sparse(op,p,alpha,0.0,eta); recon=lambda z:z
    elif model in ('M2','M2reg'):
        U=np.sign(u0)*np.abs(u0)**(p/2); e=0.0 if model=='M2' else mu2*dx*dx
        M=M2Sparse(op,p,alpha,0.0,e); recon=lambda v:np.sign(v)*np.abs(v)**(2/p)
    else: raise ValueError(model)
    steps=max(1,int(math.ceil(T/(dt_ratio*dx)))); dt=T/steps
    nit=[]; t0=time.time(); status='OK'; fail=None
    for n in range(steps):
        U,ok,hist=newton_stage_adaptive(U,dt,M,tol=tol,maxit=maxit,refresh=refresh,line_search=12)
        nit.append(len(hist))
        if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U))>1e8:
            status='FAIL'; fail=n+1; break
    u=recon(U); tr=(fail-1)*dt if fail else T
    # physical diagnostics
    C=dx*np.sum(u); P=0.5*dx*np.sum(u*u)
    return dict(status=status,model=model,p=p,N=N,L=L,dx=dx,T=T,T_reached=tr,dt=dt,steps=steps,
                eta=eta,mu2=mu2,mean_newton=float(np.mean(nit)) if nit else None, max_newton=max(nit) if nit else None,
                C=float(C),P=float(P),shift=shift,x=x.tolist(),u=u.tolist(),u0=u0.tolist())

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--model',required=True); ap.add_argument('--p',type=float,default=1.5); ap.add_argument('--N',type=int,required=True)
    ap.add_argument('--L',type=float,default=8.0); ap.add_argument('--T',type=float,default=.01); ap.add_argument('--dt-ratio',type=float,default=.1)
    ap.add_argument('--eta',type=float,default=0.0); ap.add_argument('--mu2',type=float,default=0.0); ap.add_argument('--tol',type=float,default=1e-10); ap.add_argument('--shift',type=float,default=0.0); ap.add_argument('--refresh',type=int,default=2); ap.add_argument('--maxit',type=int,default=20); ap.add_argument('--out')
    a=ap.parse_args(); r=run(a.model,a.p,a.N,a.L,a.T,a.dt_ratio,a.eta,a.mu2,a.tol,a.shift,a.refresh,a.maxit)
    if a.out:
        Path(a.out).parent.mkdir(parents=True,exist_ok=True); Path(a.out).write_text(json.dumps(r))
    print(json.dumps({k:v for k,v in r.items() if k not in ('x','u','u0')},indent=2))
