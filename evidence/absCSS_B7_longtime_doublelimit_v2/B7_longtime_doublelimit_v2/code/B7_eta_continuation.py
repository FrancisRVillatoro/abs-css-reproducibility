#!/usr/bin/env python3
import sys, math, time, json
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import splu
sys.path.insert(0,'/mnt/data/B7_work/B6_css_chase_v1/code')
from B5_css_rh import compatible_grid,pulse_sum,css_invariants
from B1_pade_sparse_newton import PadeSparseOps,M1Sparse
sys.path.insert(0,'/mnt/data/B7_work/B7_longtime_doublelimit_v1/code')
from B7_phase_locked import separation_locked, diag

def newton_step_pred(U,dt,model,Xpred=None,tol=1e-9,maxit=22,line_search=10):
    FU=model.F(U); W=U-(dt/3.0)*FU
    if Xpred is None:
        Z=model.initial_Z(U,W,dt)
    else:
        X=Xpred.copy(); V=X+(dt/4.0)*model.F(X)
        Z=np.concatenate([X,model.q(X),V,model.q(V)])
    N=model.o.N; nb=len(Z)//N; perm=np.arange(nb*N).reshape(nb,N).T.ravel()
    hist=[]
    for it in range(maxit):
        E,J=model.residual_jac(Z,W,dt,True)
        rr=float(np.linalg.norm(E,np.inf)); hist.append(rr)
        if rr<tol: return model.extract_X(Z),True,hist
        try: lu=splu(J[perm,:][:,perm].tocsc(),permc_spec='NATURAL'); dp=lu.solve(-E[perm])
        except Exception: return model.extract_X(Z),False,hist
        d=np.empty_like(dp); d[perm]=dp
        lam=1.; ok=False
        for _ in range(line_search):
            zt=Z+lam*d; rt=float(np.linalg.norm(model.residual_jac(zt,W,dt,False),np.inf))
            if rt<rr: Z=zt; ok=True; break
            lam*=.5
        if not ok: return model.extract_X(Z),False,hist
    return model.extract_X(Z),False,hist

def setup(m,eta,dt_ratio=.1,delta=.49296585513721425,p=1.5,tau=2.):
    x,h,L,N=compatible_grid(p,m,6); alpha=p*p/((p+1)*(p+2)); op=PadeSparseOps(N,h)
    sep,n=separation_locked(m,8.,delta,0); midpoint=-2+0.5*h
    U=pulse_sum(x,p,alpha,[dict(center=midpoint-sep/2,c=3.,sign=1.),dict(center=midpoint+sep/2,c=1.,sign=-1.)])
    T=sep/2+tau; steps=int(math.ceil(T/(dt_ratio*h))); dt=T/steps
    return x,h,L,N,alpha,op,sep,U,T,steps,dt

def run(m,eta,dt_ratio=.1,pred_states=None,out=None):
    x,h,L,N,alpha,op,sep,U,T,steps,dt=setup(m,eta,dt_ratio)
    M=M1Sparse(op,1.5,alpha,c0=1.,eps=eta); inv0=css_invariants(U,op,h,1.5,alpha)
    states=[]; hist=[]; fail=None; t0=time.time()
    if pred_states is not None and len(pred_states)!=steps: raise ValueError((len(pred_states),steps))
    for j in range(steps):
        xp=None if pred_states is None else pred_states[j]
        U,ok,hh=newton_step_pred(U,dt,M,xp); hist.append(len(hh))
        if not ok or not np.all(np.isfinite(U)):
            fail=j+1; break
        states.append(U.copy())
    Tr=(fail-1)*dt if fail else T
    result=dict(status='FAIL' if fail else 'OK',eta=eta,m=m,N=N,h=h,sep=sep,T=T,T_reached=Tr,dt=dt,steps=steps,fail_step=fail,runtime_s=time.time()-t0,mean_newton=float(np.mean(hist)),final_diag=diag(U,x,h),inv0=inv0,inv=css_invariants(U,op,h,1.5,alpha),x=x.tolist(),u=U.tolist())
    if out: Path(out).write_text(json.dumps(result))
    return result,states

if __name__=='__main__':
    m=115; dr=.1
    seq=[.01,.0075,.005,.004,.003]
    pred=None
    for eta in seq:
        r,states=run(m,eta,dr,pred,out=f'/mnt/data/B7_work/cont_m{m}_eta{eta:.4g}.json')
        print('ETA',eta,r['status'],r['T_reached'],'runtime',r['runtime_s'],'nit',r['mean_newton'],r['final_diag'],flush=True)
        if r['status']!='OK': break
        pred=states
