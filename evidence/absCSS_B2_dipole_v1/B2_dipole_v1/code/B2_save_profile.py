#!/usr/bin/env python3
import argparse, json, math, time
from pathlib import Path
import numpy as np
from B2_dipole_solver import grid_from_m,dipole_u,dipole_v,newton_stage_adaptive,diagnostics
from B1_pade_sparse_newton import PadeSparseOps,M1Sparse,M2Sparse

def run(model_name,p,m,T=1.,dt_ratio=.1,phase=.5,tol=1e-9,refresh=2,mu=.1,out='case.npz'):
    N,dx,L=grid_from_m(m); alpha=p*p/((p+1)*(p+2)); x=-L/2+dx*np.arange(N); xmid=phase*dx
    op=PadeSparseOps(N,dx)
    if model_name=='M1':
        U=dipole_u(x,p,alpha,xmid); model=M1Sparse(op,p,alpha,1.,0.); recon=lambda z:z
    elif model_name=='M2':
        U=dipole_v(x,p,alpha,xmid); model=M2Sparse(op,p,alpha,1.,mu*dx*dx); recon=lambda v:np.sign(v)*np.abs(v)**(2/p)
    else: raise ValueError
    steps=max(1,int(math.ceil(T/(dt_ratio*dx)))); dt=T/steps; nits=[]; t0=time.time()
    for n in range(steps):
        U,ok,h=newton_stage_adaptive(U,dt,model,tol=tol,maxit=24,refresh=refresh); nits.append(len(h))
        if not ok: return dict(status='FAIL',step=n+1,time=n*dt,runtime=time.time()-t0,residual=h)
    u=recon(U); ref=dipole_u(x,p,alpha,xmid); d=diagnostics(u,x,op,p,alpha,xmid,L,dx,ref)
    np.savez_compressed(out,x=x,u=u,ref=ref,raw=U)
    return dict(status='OK',model=model_name,p=p,m=m,N=N,dx=dx,L=L,dt=dt,steps=steps,runtime=time.time()-t0,mean_newton=float(np.mean(nits)),**d)

if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('--model',required=True);a.add_argument('--p',type=float,required=True);a.add_argument('--m',type=int,required=True);a.add_argument('--out',required=True);a.add_argument('--json',required=True);a.add_argument('--tol',type=float,default=1e-9);a.add_argument('--refresh',type=int,default=2);a.add_argument('--mu',type=float,default=.1); z=a.parse_args()
    r=run(z.model,z.p,z.m,tol=z.tol,refresh=z.refresh,mu=z.mu,out=z.out); Path(z.json).write_text(json.dumps(r,indent=2)+'\n'); print(json.dumps(r,indent=2))
