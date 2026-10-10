#!/usr/bin/env python3
import math,time
import numpy as np
from scipy.sparse.linalg import splu
from B1_pade_sparse_newton import PadeSparseOps,M1Sparse,compact_u,invariants_u
from B1_radiation_diagnostic import exterior_metrics

def newton_stage_adaptive(U,dt,model,tol=1e-10,maxit=16,refresh=4,line_search=8):
    FU=model.F(U); W=U-(dt/3.0)*FU; Z=model.initial_Z(U,W,dt)
    N=model.o.N; nb=len(Z)//N; perm=np.arange(nb*N).reshape(nb,N).T.ravel()
    hist=[]; lu=None
    for it in range(maxit):
        if lu is None or it%refresh==0:
            E,J=model.residual_jac(Z,W,dt,True)
            lu=splu(J[perm,:][:,perm].tocsc(),permc_spec='NATURAL')
        else:
            E=model.residual_jac(Z,W,dt,False)
        rr=float(np.linalg.norm(E,np.inf)); hist.append(rr)
        if rr<tol: return model.extract_X(Z),True,hist
        dp=lu.solve(-E[perm]); d=np.empty_like(dp); d[perm]=dp
        lam=1.; accepted=False
        for _ in range(line_search):
            Zt=Z+lam*d
            rt=float(np.linalg.norm(model.residual_jac(Zt,W,dt,False),np.inf))
            if rt<rr: Z=Zt; accepted=True; break
            lam*=0.5
        if not accepted:
            E,J=model.residual_jac(Z,W,dt,True)
            lu=splu(J[perm,:][:,perm].tocsc(),permc_spec='NATURAL')
            dp=lu.solve(-E[perm]); d=np.empty_like(dp); d[perm]=dp; lam=1.
            for _ in range(line_search):
                Zt=Z+lam*d
                rt=float(np.linalg.norm(model.residual_jac(Zt,W,dt,False),np.inf))
                if rt<rr: Z=Zt; accepted=True; break
                lam*=0.5
            if not accepted: return model.extract_X(Z),False,hist
    return model.extract_X(Z),False,hist

def run_case_fast(N,p,T,method,L,dt_ratio,eta_coeff,x0,tol=1e-10,maxit=16,refresh=2):
    dx=L/N; alpha=p*p/((p+1)*(p+2)); eps=0.0 if method=='M1' else eta_coeff*dx**3
    x=-L/2+dx*np.arange(N); op=PadeSparseOps(N,dx); U=compact_u(x,p,alpha,x0=x0); U0=U.copy(); model=M1Sparse(op,p,alpha,1.0,eps)
    steps=max(1,int(np.ceil(T/(dt_ratio*dx)))); dt=T/steps; nits=[]; t0=time.time()
    for n in range(steps):
        U,ok,hist=newton_stage_adaptive(U,dt,model,tol=tol,maxit=maxit,refresh=refresh)
        nits.append(len(hist))
        if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U))>1e6:
            return {'status':'FAIL','method':method,'N':N,'p':p,'T':T,'L':L,'dx':dx,'dt':dt,'dt_ratio':dt/dx,'eta':eps,'step':n+1,'residual':hist,'runtime_s':time.time()-t0,'mean_newton':float(np.mean(nits))}
    met,_=exterior_metrics(U,x,p,alpha,x0,L,dx); z0=invariants_u(U0,op,dx,p,alpha); z=invariants_u(U,op,dx,p,alpha)
    return {'status':'OK','method':method,'N':N,'p':p,'T':T,'L':L,'dx':dx,'dt':dt,'dt_ratio':dt/dx,'steps':steps,'eta':eps,'eta_coeff_dx3':eta_coeff if method=='M5' else 0.0,'runtime_s':time.time()-t0,'mean_newton':float(np.mean(nits)),'max_newton':int(max(nits)),'drift_C':float((z[0]-z0[0])/abs(z0[0])),'drift_P':float((z[1]-z0[1])/abs(z0[1])),'drift_H':float((z[2]-z0[2])/max(abs(z0[2]),1e-300)),'negative_mass_total':float(z[3]),**met}
