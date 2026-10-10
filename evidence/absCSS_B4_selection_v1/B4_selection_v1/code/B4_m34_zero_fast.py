#!/usr/bin/env python3
import sys,math,json,time
from pathlib import Path
import numpy as np
root=Path('/mnt/data/absCSS_B4_work')
sys.path.insert(0,str(root/'B2src/B2_dipole_v1/code'))
sys.path.insert(0,str(root/'B3src/B3_transverse_v1/code'))
from B1_pade_sparse_newton import PadeSparseOps
from B2_dipole_solver import newton_stage_adaptive
from B3_regularized_solver import M34Regularized

def run(p,r,eps,points=16,L=8,T=.01,dt_ratio=.1):
    N=int(round(L/(eps/points))); N += N%2; dx=L/N; x=-L/2+dx*np.arange(N); a=p*p/((p+1)*(p+2)); op=PadeSparseOps(N,dx); M=M34Regularized(op,p,a,r,eps,0.0)
    k=2*np.pi/L; U=np.sin(k*x)/k; U0=U.copy(); steps=max(1,int(math.ceil(T/(dt_ratio*dx))));dt=T/steps; nit=[];t0=time.time()
    for n in range(steps):
      U,ok,h=newton_stage_adaptive(U,dt,M,tol=1e-10,maxit=20,refresh=2,line_search=10);nit.append(len(h))
      if not ok: return dict(status='FAIL',p=p,r=r,eps=eps,N=N,dx=dx,T_reached=n*dt,hist=h)
    return dict(status='OK',p=p,r=r,eps=eps,N=N,dx=dx,dt=dt,T=T,mean_newton=float(np.mean(nit)),runtime=time.time()-t0,x=x.tolist(),u=U.tolist(),u0=U0.tolist())
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--p',type=float,default=1.5);ap.add_argument('--r',type=int,required=True);ap.add_argument('--eps',type=float,required=True);ap.add_argument('--points',type=int,default=16);ap.add_argument('--out');a=ap.parse_args();d=run(a.p,a.r,a.eps,a.points);print(json.dumps({k:v for k,v in d.items() if k not in ('x','u','u0')},indent=2));
 if a.out: Path(a.out).write_text(json.dumps(d))
