#!/usr/bin/env python3
import sys,math,json,time,argparse
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
sys.path.insert(0,'/mnt/data/B6_work/B5_css_rh_v1/code')
from B5_css_rh import compatible_grid,pulse_sum,css_invariants,newton_step
from B1_pade_sparse_newton import PadeSparseOps,M1Sparse

def zeros_periodic(u,x,h):
    z=[]; N=len(u)
    for i in range(N):
        j=(i+1)%N
        if u[i]==0: z.append(float(x[i]))
        elif u[i]*u[j]<0:
            xx=float(x[i]+h*(-u[i])/(u[j]-u[i]))
            z.append(xx)
    return z

def diag(u,x,h,L):
    imx=find_peaks(u,prominence=.01)[0]; imn=find_peaks(-u,prominence=.01)[0]
    maxs=sorted([(float(u[i]),float(x[i])) for i in imx],reverse=True)[:6]
    mins=sorted([(float(u[i]),float(x[i])) for i in imn])[:6]
    return dict(max=float(u.max()),xmin=float(x[u.argmin()]),min=float(u.min()),xmax=float(x[u.argmax()]),
                positive_mass=float(h*np.sum(np.maximum(u,0))),negative_mass_abs=float(h*np.sum(np.maximum(-u,0))),
                maxima=maxs,minima=mins,nzero=len(zeros_periodic(u,x,h)),zeros=zeros_periodic(u,x,h)[:20])

def run(m,eta,T=6,dt_ratio=.05,periods=6,offset=-2,p=1.5,C=None):
    x,h,L,N=compatible_grid(p,m,periods); alpha=p*p/((p+1)*(p+2)); op=PadeSparseOps(N,h)
    cfast,cslow,sep=3.,1.,8.; c0=1.; midpoint=offset+0.5*h
    U=pulse_sum(x,p,alpha,[dict(center=midpoint-sep/2,c=cfast,sign=1.),dict(center=midpoint+sep/2,c=cslow,sign=-1.)])
    M=M1Sparse(op,p,alpha,c0=c0,eps=eta)
    inv0=css_invariants(U,op,h,p,alpha)
    steps=max(1,int(math.ceil(T/(dt_ratio*h))));dt=T/steps
    targets=[0,.5,.8,1,1.2,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6];targets=[t for t in targets if t<=T]
    snaps=[dict(t=0.,diag=diag(U,x,h,L),inv=inv0)];it=1; hist=[];fail=None;t0=time.time()
    for n in range(steps):
        U,ok,hh=newton_step(U,dt,M,tol=1e-9,maxit=22,refresh=1);hist.append(len(hh));t=(n+1)*dt
        if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U))>1e8:
            fail=n+1;break
        while it<len(targets) and t+0.5*dt>=targets[it]:
            snaps.append(dict(t=float(t),diag=diag(U,x,h,L),inv=css_invariants(U,op,h,p,alpha)));it+=1
    Tr=(fail-1)*dt if fail else T
    inv=css_invariants(U,op,h,p,alpha)
    return dict(status='FAIL' if fail else 'OK',T=T,T_reached=Tr,p=p,m=m,N=N,L=L,h=h,eta=eta,C=C,dt=dt,dt_ratio=dt/h,periods=periods,offset=offset,
                t_touch=(8-2*math.pi)/2,t_centers=4.,runtime_s=time.time()-t0,mean_newton=float(np.mean(hist)),max_newton=max(hist),fail_step=fail,
                inv0=inv0,inv=inv,snapshots=snaps,final_diag=diag(U,x,h,L),x=x.tolist(),u=U.tolist())

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--m',type=int,required=True);ap.add_argument('--eta',type=float);ap.add_argument('--C',type=float);ap.add_argument('--T',type=float,default=6);ap.add_argument('--dt-ratio',type=float,default=.05);ap.add_argument('--out',required=True);a=ap.parse_args()
    h=2*math.pi/a.m; eta=a.eta if a.eta is not None else a.C*h*h
    r=run(a.m,eta,a.T,dt_ratio=a.dt_ratio,C=a.C);Path(a.out).write_text(json.dumps(r));print(r['status'],r['T_reached'],r['runtime_s'],r['N'],eta)
    for s in r['snapshots']:
        d=s['diag'];print(f"t={s['t']:.3f} max={d['max']:.4f}@{d['xmax']:.3f} min={d['min']:.4f}@{d['xmin']:.3f} M+={d['positive_mass']:.3f} M-={d['negative_mass_abs']:.3f} nz={d['nzero']}")
