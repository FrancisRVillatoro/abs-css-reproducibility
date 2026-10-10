#!/usr/bin/env python3
import sys,math,json,time,argparse
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
HERE=Path('/mnt/data/B7_work/B6_css_chase_v1/code'); sys.path.insert(0,str(HERE))
from B5_css_rh import compatible_grid,pulse_sum,css_invariants,newton_step
from B1_pade_sparse_newton import PadeSparseOps,M1Sparse

def zeros_periodic(u,x,h):
 z=[];N=len(u)
 for i in range(N):
  j=(i+1)%N
  if u[i]*u[j]<0: z.append(float(x[i]+h*(-u[i])/(u[j]-u[i])))
 return z

def diag(u,x,h):
 return dict(max=float(u.max()),xmax=float(x[u.argmax()]),min=float(u.min()),xmin=float(x[u.argmin()]),positive_mass=float(h*np.sum(np.maximum(u,0))),negative_mass_abs=float(h*np.sum(np.maximum(-u,0))),nzero=len(zeros_periodic(u,x,h)))

def separation_locked(m,target=8.0,delta=0.49296585513721425,parity=0):
 h=2*math.pi/m; y=target/h-delta; n=int(round(y))
 if n%2!=parity:
  c=[n-1,n+1]; n=min(c,key=lambda k:abs((k+delta)*h-target))
 return (n+delta)*h,n

def run(m,eta,delta=0.49296585513721425,parity=0,tau=2.,dt_ratio=.1,p=1.5):
 x,h,L,N=compatible_grid(p,m,6); alpha=p*p/((p+1)*(p+2)); op=PadeSparseOps(N,h)
 sep,n=separation_locked(m,8.,delta,parity); cfast,cslow=3.,1.; c0=1.; midpoint=-2+0.5*h
 U=pulse_sum(x,p,alpha,[dict(center=midpoint-sep/2,c=cfast,sign=1.),dict(center=midpoint+sep/2,c=cslow,sign=-1.)])
 M=M1Sparse(op,p,alpha,c0=c0,eps=eta); inv0=css_invariants(U,op,h,p,alpha)
 T=sep/2+tau; steps=int(math.ceil(T/(dt_ratio*h)));dt=T/steps; hist=[]; fail=None;t0=time.time()
 for j in range(steps):
  U,ok,hh=newton_step(U,dt,M,tol=1e-9,maxit=22,refresh=1);hist.append(len(hh))
  if not ok or not np.all(np.isfinite(U)): fail=j+1;break
 Tr=(fail-1)*dt if fail else T
 return dict(status='FAIL' if fail else 'OK',T=T,T_reached=Tr,tau=tau,m=m,N=N,h=h,eta=eta,sep=sep,nsep=n,delta=delta,parity=parity,dt=dt,runtime_s=time.time()-t0,mean_newton=float(np.mean(hist)),fail_step=fail,inv0=inv0,inv=css_invariants(U,op,h,p,alpha),final_diag=diag(U,x,h),x=x.tolist(),u=U.tolist())
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--m',type=int,required=True);ap.add_argument('--eta',type=float,required=True);ap.add_argument('--out',required=True);ap.add_argument('--delta',type=float,default=.49296585513721425);a=ap.parse_args();r=run(a.m,a.eta,a.delta);Path(a.out).write_text(json.dumps(r));print(r['status'],r['T_reached'],r['sep'],r['final_diag'],r['runtime_s'])
