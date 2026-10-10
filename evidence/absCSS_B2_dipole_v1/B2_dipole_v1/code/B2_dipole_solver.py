#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import splu

from B1_pade_sparse_newton import (
    PadeSparseOps, M1Sparse, M2Sparse, compact_u, compact_v, phase_shape, invariants_u
)


def newton_stage_adaptive(U,dt,model,tol=1e-10,maxit=18,refresh=2,line_search=8):
    FU=model.F(U); W=U-(dt/3.0)*FU; Z=model.initial_Z(U,W,dt)
    N=model.o.N; nb=len(Z)//N; perm=np.arange(nb*N).reshape(nb,N).T.ravel()
    hist=[]; lu=None
    for it in range(maxit):
        if lu is None or it % refresh == 0:
            E,J=model.residual_jac(Z,W,dt,True)
            try:
                lu=splu(J[perm,:][:,perm].tocsc(),permc_spec='NATURAL')
            except Exception:
                return model.extract_X(Z),False,hist
        else:
            E=model.residual_jac(Z,W,dt,False)
        rr=float(np.linalg.norm(E,np.inf)); hist.append(rr)
        if rr < tol:
            return model.extract_X(Z),True,hist
        try:
            dp=lu.solve(-E[perm])
        except Exception:
            return model.extract_X(Z),False,hist
        if not np.all(np.isfinite(dp)):
            return model.extract_X(Z),False,hist
        d=np.empty_like(dp); d[perm]=dp
        lam=1.0; accepted=False
        for _ in range(line_search):
            Zt=Z+lam*d
            rt=float(np.linalg.norm(model.residual_jac(Zt,W,dt,False),np.inf))
            if rt < rr:
                Z=Zt; accepted=True; break
            lam*=0.5
        if not accepted:
            # refresh once before declaring failure
            E,J=model.residual_jac(Z,W,dt,True)
            try:
                lu=splu(J[perm,:][:,perm].tocsc(),permc_spec='NATURAL')
                dp=lu.solve(-E[perm])
            except Exception:
                return model.extract_X(Z),False,hist
            d=np.empty_like(dp); d[perm]=dp
            lam=1.0
            for _ in range(line_search):
                Zt=Z+lam*d
                rt=float(np.linalg.norm(model.residual_jac(Zt,W,dt,False),np.inf))
                if rt < rr:
                    Z=Zt; accepted=True; break
                lam*=0.5
            if not accepted:
                return model.extract_X(Z),False,hist
    return model.extract_X(Z),False,hist


def dipole_u(x,p,alpha,xmid=0.0,c=1.0):
    # Normalisation gives half-width Lc=pi, but compute through compact_u shifts.
    out=np.zeros_like(x,dtype=float)
    # Support positive lobe [xmid-2pi,xmid], centre xmid-pi.
    # Support negative lobe [xmid,xmid+2pi], centre xmid+pi.
    z=x-xmid
    left=(z>=-2*math.pi-1e-14)&(z<=0+1e-14)
    right=(z>=0-1e-14)&(z<=2*math.pi+1e-14)
    if np.any(left): out[left]=compact_u(x[left],p,alpha,x0=xmid-math.pi,c=c)
    if np.any(right): out[right]=-compact_u(x[right],p,alpha,x0=xmid+math.pi,c=c)
    # force exact internal junction zero if represented twice
    out[np.abs(z)<1e-13]=0.0
    return out


def dipole_v(x,p,alpha,xmid=0.0,c=1.0):
    out=np.zeros_like(x,dtype=float)
    z=x-xmid
    left=(z>=-2*math.pi-1e-14)&(z<=0+1e-14)
    right=(z>=0-1e-14)&(z<=2*math.pi+1e-14)
    if np.any(left): out[left]=compact_v(x[left],p,alpha,x0=xmid-math.pi,c=c)
    if np.any(right): out[right]=-compact_v(x[right],p,alpha,x0=xmid+math.pi,c=c)
    out[np.abs(z)<1e-13]=0.0
    return out


def periodic_distance(x,center,L):
    return (x-center+0.5*L)%L-0.5*L


def shifted_dipole(x,p,alpha,center,c=1.0):
    # x lies on periodic domain but support remains far from boundaries in our runs.
    return dipole_u(x,p,alpha,xmid=center,c=c)


def zero_crossing_near(u,x,expected,window=math.pi/2):
    mask=np.where(np.abs(x-expected)<=window)[0]
    if len(mask)<2: return float('nan')
    best=None
    for ia,ib in zip(mask[:-1],mask[1:]):
        ya,yb=u[ia],u[ib]
        if ya==0.0:
            root=x[ia]
        elif yb==0.0:
            root=x[ib]
        elif ya*yb<0:
            root=x[ia]-ya*(x[ib]-x[ia])/(yb-ya)
        else:
            continue
        score=abs(root-expected)
        if best is None or score<best[0]: best=(score,root)
    if best is not None: return float(best[1])
    j=mask[np.argmin(np.abs(u[mask]))]
    return float(x[j])


def analytic_dipole_invariants(p,c=1.0):
    alpha=p*p/((p+1)*(p+2))
    A=(c*(p+1)*(p+2)/2.0)**(1.0/p)
    k=0.5
    I=lambda power: math.sqrt(math.pi)*math.gamma((power+1)/2)/math.gamma((power+2)/2)
    Pc=A*A/(2*k)*I(4.0/p)
    Hc=-2*c*Pc/(p+2)
    return dict(C=0.0,P=2*Pc,H=2*Hc,alpha=alpha,A=A)


def diagnostics(u,x,op,p,alpha,xmid,L,dx,ref0=None):
    if ref0 is None: ref0=dipole_u(x,p,alpha,xmid=xmid)
    phase,shape=phase_shape(u,ref0,dx)
    center=xmid+phase
    ref=shifted_dipole(x,p,alpha,center)
    # Recompute direct residual to avoid any FFT interpolation artefact in diagnostics.
    residual=u-ref
    d=periodic_distance(x,center,L)
    A=float(np.max(np.abs(ref0)))
    internal=np.abs(d)<=0.5
    exterior=np.abs(d)>(2*math.pi+0.5)
    thr=1e-3*A
    bad=internal & (np.abs(residual)>thr)
    if np.any(bad):
        coords=d[bad]; width=float(coords.max()-coords.min()+dx)
    else: width=0.0
    zero=zero_crossing_near(u,x,center)
    C,P,H,neg=invariants_u(u,op,dx,p,alpha)
    ex=analytic_dipole_invariants(p)
    return {
        'phase':float(phase),'shape_relL2':float(np.linalg.norm(residual)/np.linalg.norm(ref)),
        'internal_zero':zero,'zero_minus_aligned_exact':float(zero-center),
        'C_abs':float(abs(C)),'C_signed':float(C),
        'P_rel_exact':float((P-ex['P'])/ex['P']),
        'H_rel_exact':float((H-ex['H'])/abs(ex['H'])),
        'P':float(P),'H':float(H),'negative_mass':float(neg),
        'internal_max_abs_error':float(np.max(np.abs(residual[internal]))),
        'internal_L1_abs_error':float(dx*np.sum(np.abs(residual[internal]))),
        'internal_width_err_gt_1e3A':width,
        'exterior_max_abs':float(np.max(np.abs(u[exterior]))) if np.any(exterior) else 0.0,
        'exterior_L1_abs':float(dx*np.sum(np.abs(u[exterior]))) if np.any(exterior) else 0.0,
        'max_abs':float(np.max(np.abs(u))),
    }


def grid_from_m(m,Ltarget=40.0):
    h=2*math.pi/m
    n0=int(round(Ltarget/h))
    # even N -> x=0 is exactly a grid node on [-L/2,L/2)
    N=n0 if n0%2==0 else (n0+1 if abs((n0+1)*h-Ltarget)<=abs((n0-1)*h-Ltarget) else n0-1)
    L=N*h
    return N,h,L


def run_case(model_name,p,m,T=1.0,dt_ratio=0.1,c0=1.0,mu_m2=0.1,eta_coeff_m5=0.05,
             tol=1e-10,maxit=18,refresh=2,track=False,interface_phase=0.0):
    N,dx,L=grid_from_m(m)
    alpha=p*p/((p+1)*(p+2)); x=-L/2+dx*np.arange(N); xmid=interface_phase*dx
    op=PadeSparseOps(N,dx)
    if model_name in ('M1','M5'):
        U=dipole_u(x,p,alpha,xmid)
        eps=0.0 if model_name=='M1' else eta_coeff_m5*dx**3
        model=M1Sparse(op,p,alpha,c0,eps)
        reconstruct=lambda z:z
    elif model_name in ('M2','M2_0'):
        U=dipole_v(x,p,alpha,xmid)
        eps=0.0 if model_name=='M2_0' else mu_m2*dx**2
        model=M2Sparse(op,p,alpha,c0,eps)
        reconstruct=lambda v:np.sign(v)*np.abs(v)**(2.0/p)
    else: raise ValueError(model_name)
    U0=U.copy(); u0=reconstruct(U0); ref0=dipole_u(x,p,alpha,xmid)
    inv0=invariants_u(u0,op,dx,p,alpha)
    steps=max(1,int(math.ceil(T/(dt_ratio*dx)))); dt=T/steps
    nit=[]; t0=time.time(); track_rows=[]
    if track:
        d0=diagnostics(u0,x,op,p,alpha,xmid,L,dx,ref0)
        track_rows.append((0.0,d0['internal_zero'],d0['shape_relL2']))
    for n in range(steps):
        U,ok,hist=newton_stage_adaptive(U,dt,model,tol,maxit,refresh)
        nit.append(len(hist))
        if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U))>1e6:
            return {'status':'FAIL','model':model_name,'p':p,'m':m,'N':N,'dx':dx,'L':L,
                    'T':T,'dt':dt,'dt_ratio':dt/dx,'eps':eps,'interface_phase':interface_phase,'step':n+1,'time_reached':n*dt,
                    'residual':hist,'runtime_s':time.time()-t0,'mean_newton':float(np.mean(nit))}
        if track and ((n+1)%max(1,steps//10)==0 or n+1==steps):
            uu=reconstruct(U); d=diagnostics(uu,x,op,p,alpha,xmid,L,dx,ref0)
            track_rows.append(((n+1)*dt,d['internal_zero'],d['shape_relL2']))
    u=reconstruct(U); d=diagnostics(u,x,op,p,alpha,xmid,L,dx,ref0)
    inv=invariants_u(u,op,dx,p,alpha)
    d.update({
        'status':'OK','model':model_name,'p':p,'m':m,'N':N,'dx':dx,'L':L,'T':T,
        'dt':dt,'dt_ratio':dt/dx,'steps':steps,'eps':eps,'interface_phase':interface_phase,'runtime_s':time.time()-t0,
        'mean_newton':float(np.mean(nit)),'max_newton':int(max(nit)),
        'C_change_abs':float(abs(inv[0]-inv0[0])),
        'P_drift_rel':float((inv[1]-inv0[1])/abs(inv0[1])),
        'H_drift_rel':float((inv[2]-inv0[2])/max(abs(inv0[2]),1e-300)),
        'track':track_rows,
    })
    return d


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model',choices=['M1','M5','M2','M2_0'],required=True)
    ap.add_argument('--p',type=float,required=True); ap.add_argument('--m',type=int,required=True)
    ap.add_argument('--T',type=float,default=1.0); ap.add_argument('--dt-ratio',type=float,default=0.1)
    ap.add_argument('--c0',type=float,default=1.0); ap.add_argument('--mu-m2',type=float,default=0.1)
    ap.add_argument('--eta-coeff-m5',type=float,default=0.05); ap.add_argument('--tol',type=float,default=1e-10)
    ap.add_argument('--maxit',type=int,default=18); ap.add_argument('--refresh',type=int,default=2); ap.add_argument('--interface-phase',type=float,default=0.0)
    ap.add_argument('--track',action='store_true'); ap.add_argument('--out')
    a=ap.parse_args(); r=run_case(a.model,a.p,a.m,a.T,a.dt_ratio,a.c0,a.mu_m2,a.eta_coeff_m5,a.tol,a.maxit,a.refresh,a.track,a.interface_phase)
    s=json.dumps(r,indent=2); print(s)
    if a.out: Path(a.out).write_text(s+'\n')
if __name__=='__main__': main()
