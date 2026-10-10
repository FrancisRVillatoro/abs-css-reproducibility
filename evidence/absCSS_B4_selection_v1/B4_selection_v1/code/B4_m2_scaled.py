#!/usr/bin/env python3
import sys, math, json, time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import splu
ROOT=Path('/mnt/data/absCSS_B4_work')
sys.path.insert(0,str(ROOT/'B2src/B2_dipole_v1/code'))
from B1_pade_sparse_newton import PadeSparseOps, M2Sparse

def sine_u(x,L,a=1.0,shift=0.0):
    k=2*np.pi/L; return (a/k)*np.sin(k*(x-shift))

def scaled_norm(E,N,h):
    # derivative-definition blocks naturally scale like h^{-j}; scale them back.
    weights=(h,h*h,h**3,h**4,1.0,h,h*h,h**3,h**4,1.0)
    return max(float(np.max(np.abs(E[j*N:(j+1)*N])))*weights[j] for j in range(10))

def newton_scaled(U,dt,M,h,tol=1e-10,maxit=30,refresh=1,line_search=14):
    FU=M.F(U); W=U-(dt/3)*FU; Z=M.initial_Z(U,W,dt)
    N=M.o.N; nb=10; perm=np.arange(nb*N).reshape(nb,N).T.ravel(); hist=[]; lu=None
    for it in range(maxit):
        if lu is None or it%refresh==0:
            E,J=M.residual_jac(Z,W,dt,True)
            try: lu=splu(J[perm,:][:,perm].tocsc(),permc_spec='NATURAL')
            except Exception: return M.extract_X(Z),False,hist
        else: E=M.residual_jac(Z,W,dt,False)
        rr=scaled_norm(E,N,h); hist.append(rr)
        if rr<tol: return M.extract_X(Z),True,hist
        try: dp=lu.solve(-E[perm])
        except Exception: return M.extract_X(Z),False,hist
        d=np.empty_like(dp); d[perm]=dp
        accepted=False; lam=1.0
        for _ in range(line_search):
            Zt=Z+lam*d; Et=M.residual_jac(Zt,W,dt,False); rt=scaled_norm(Et,N,h)
            if rt < rr*(1-1e-12) or rt<tol:
                Z=Zt; accepted=True; break
            lam*=0.5
        if not accepted: return M.extract_X(Z),False,hist
    return M.extract_X(Z),False,hist

def run(p,N,L=8,T=.01,dt_ratio=.025,eta_v=1e-5,tol=1e-9,maxit=30):
    h=L/N; x=-L/2+h*np.arange(N); alpha=p*p/((p+1)*(p+2)); op=PadeSparseOps(N,h)
    u0=sine_u(x,L); V=np.sign(u0)*np.abs(u0)**(p/2); M=M2Sparse(op,p,alpha,0.0,eta_v)
    steps=max(1,int(math.ceil(T/(dt_ratio*h)))); dt=T/steps; nit=[]; status='OK'; fail=None; t0=time.time()
    for n in range(steps):
        V,ok,hist=newton_scaled(V,dt,M,h,tol=tol,maxit=maxit)
        nit.append(len(hist))
        if not ok or not np.all(np.isfinite(V)) or np.max(np.abs(V))>1e8:
            status='FAIL'; fail=n+1; break
    u=np.sign(V)*np.abs(V)**(2/p); tr=(fail-1)*dt if fail else T
    return dict(status=status,p=p,N=N,L=L,dx=h,T=T,T_reached=tr,dt=dt,dt_ratio=dt/h,steps=steps,
                eta_v=eta_v,mean_newton=float(np.mean(nit)),max_newton=max(nit),runtime_s=time.time()-t0,
                C=float(h*np.sum(u)),P=float(.5*h*np.sum(u*u)),x=x.tolist(),u=u.tolist(),u0=u0.tolist())

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--p',type=float,default=1.5); ap.add_argument('--N',type=int,required=True); ap.add_argument('--L',type=float,default=8); ap.add_argument('--T',type=float,default=.01); ap.add_argument('--dt-ratio',type=float,default=.025); ap.add_argument('--eta-v',type=float,default=1e-5); ap.add_argument('--tol',type=float,default=1e-9); ap.add_argument('--maxit',type=int,default=30); ap.add_argument('--out')
    a=ap.parse_args(); r=run(a.p,a.N,a.L,a.T,a.dt_ratio,a.eta_v,a.tol,a.maxit)
    if a.out: Path(a.out).write_text(json.dumps(r))
    print(json.dumps({k:v for k,v in r.items() if k not in ('x','u','u0')},indent=2))
