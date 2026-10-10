#!/usr/bin/env python3
from __future__ import annotations
import math, json, time, sys
from pathlib import Path
import numpy as np
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import splu

# load validated CSS/B2 code shipped in the same directory
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from B1_pade_sparse_newton import PadeSparseOps, M1Sparse, compact_u, phase_shape, invariants_u
from B2_dipole_solver import dipole_u


def diag(v): return diags(v,0,format='csr')

class M6RHSparse:
    """Matched absolute Rosenau-Hyman comparator.
    F is the left-hand spatial operator in u_t + F(u)=0, optionally in
    frame speed c0: F=-c0*u_x+a Dx(g)+b Dx^3(g)+eta u_xxxx.
    Exact R12 stages are multiplied by one Padé A.
    """
    def __init__(self,ops,p,alpha,c0=0.0,eta=0.0):
        self.o=ops; self.p=float(p); self.alpha=float(alpha); self.c0=float(c0); self.eta=float(eta)
        self.q=self.p+1.0
        self.a=4.0/(self.p+2.0)**2
        self.b=4.0*self.alpha/((self.p+1.0)*(self.p+2.0))
        self.kappa=math.sqrt(self.a/self.b)

    def g(self,u): return np.abs(u)**self.p*u
    def gp(self,u): return (self.p+1.0)*np.abs(u)**self.p
    def R(self,u):
        o=self.o; g=self.g(u)
        return -self.c0*(o.B@u)+self.a*(o.B@g)+self.b*(o.C@g)+self.eta*(o.D@u)
    def Rjac(self,u):
        o=self.o; G=diag(self.gp(u))
        return (-self.c0*o.B + self.a*(o.B@G)+self.b*(o.C@G)+self.eta*o.D).tocsr()
    def F(self,u):
        o=self.o; g=self.g(u)
        return -self.c0*o.D1(u)+self.a*o.D1(g)+self.b*o.D3(g)+self.eta*o.D4(u)
    def residual_jac(self,Z,W,dt,need_jac=True):
        N=self.o.N; o=self.o; X=Z[:N]; V=Z[N:]
        k1=4.0/dt; k2=3.0/(2.0*dt)
        E1=k1*(o.A@(V-X))-self.R(X)
        E2=k2*(o.A@(W-X))-self.R(V)
        E=np.concatenate([E1,E2])
        if not need_jac: return E
        RX=self.Rjac(X); RV=self.Rjac(V)
        J=bmat([[-k1*o.A-RX, k1*o.A],[-k2*o.A,-RV]],format='csr')
        return E,J
    def initial_Z(self,U,W,dt):
        V=U+(dt/4.0)*self.F(U)
        return np.concatenate([U.copy(),V])
    def extract_X(self,Z): return Z[:self.o.N]


def newton_step(U,dt,model,tol=1e-10,maxit=18,refresh=1,line_search=10):
    FU=model.F(U); W=U-(dt/3.0)*FU; Z=model.initial_Z(U,W,dt)
    N=model.o.N; nb=len(Z)//N; perm=np.arange(nb*N).reshape(nb,N).T.ravel()
    hist=[]; lu=None
    for it in range(maxit):
        if lu is None or it%refresh==0:
            E,J=model.residual_jac(Z,W,dt,True)
            try: lu=splu(J[perm,:][:,perm].tocsc(),permc_spec='NATURAL')
            except Exception: return model.extract_X(Z),False,hist
        else: E=model.residual_jac(Z,W,dt,False)
        rr=float(np.linalg.norm(E,np.inf)); hist.append(rr)
        if rr<tol: return model.extract_X(Z),True,hist
        try: dp=lu.solve(-E[perm])
        except Exception: return model.extract_X(Z),False,hist
        if not np.all(np.isfinite(dp)): return model.extract_X(Z),False,hist
        d=np.empty_like(dp); d[perm]=dp
        lam=1.; ok=False
        for _ in range(line_search):
            zt=Z+lam*d; rt=float(np.linalg.norm(model.residual_jac(zt,W,dt,False),np.inf))
            if rt<rr: Z=zt; ok=True; break
            lam*=.5
        if not ok: return model.extract_X(Z),False,hist
    return model.extract_X(Z),False,hist


def compatible_grid(p,m,periods=8):
    # L=periods * 2pi/kappa = periods*2pi*p/(p+1), h=2pi/m.
    # Caller chooses m such that N=periods*m*p/(p+1) is integer.
    h=2*math.pi/m; L=periods*2*math.pi*p/(p+1.0); N=int(round(L/h))
    err=abs(N*h-L)
    if err>1e-11: raise ValueError((p,m,N,L,h,err))
    x=-L/2+h*np.arange(N)
    return x,h,L,N


def exact_lobe(x,p,alpha,center,c=1.0,sign=1.0): return sign*compact_u(x,p,alpha,x0=center,c=c)

def exact_dipole(x,p,alpha,xmid,c=1.0): return dipole_u(x,p,alpha,xmid=xmid,c=c)


def rh_invariants(u,x,dx,p,alpha,op,t=0.0,c0=0.0):
    a=4/(p+2.)**2; b=4*alpha/((p+1)*(p+2)); kap=math.sqrt(a/b)
    C0=dx*np.sum(u); Ccy=dx*np.sum(u*np.cos(kap*x)); Csy=dx*np.sum(u*np.sin(kap*x))
    th=kap*c0*t; Cc=math.cos(th)*Ccy-math.sin(th)*Csy; Cs=math.cos(th)*Csy+math.sin(th)*Ccy
    H=-dx*np.sum(np.abs(u)**(p+2))/(p+2)
    P=.5*dx*np.sum(u*u)
    ux=op.D1(u)
    rhs=b*(p+1)*p/2*dx*np.sum(np.sign(u)*np.abs(u)**(p-1)*ux**3)
    dPdt=-dx*np.sum(u*(-c0*op.D1(u)+a*op.D1(np.abs(u)**p*u)+b*op.D3(np.abs(u)**p*u)))
    return dict(C0=float(C0),Cc=float(Cc),Cs=float(Cs),Cc_frame=float(Ccy),Cs_frame=float(Csy),H=float(H),P=float(P),Pdot_rhs=float(rhs),Pdot_semidisc=float(dPdt),Pdot_resid=float(dPdt-rhs))


def css_invariants(u,op,dx,p,alpha):
    C,P,H,neg=invariants_u(u,op,dx,p,alpha); return dict(C=float(C),P=float(P),H=float(H),negative_mass=float(neg))


def run_profile(kind,model,p,m,T=1.0,dt_ratio=.05,phase=.5,c=1.0,eta=0.0,periods=8,tol=1e-10,maxit=18,save_snap=False,c0=1.0):
    x,h,L,N=compatible_grid(p,m,periods); alpha=p*p/((p+1)*(p+2)); op=PadeSparseOps(N,h)
    # center/junction at phase*h around zero; exact support is phase locked since 2pi/h=m integer
    center=phase*h
    if kind=='lobe': U=exact_lobe(x,p,alpha,center,c); ref0=U.copy()
    elif kind=='dipole': U=exact_dipole(x,p,alpha,center,c); ref0=U.copy()
    else: raise ValueError(kind)
    if model=='CSS': M=M1Sparse(op,p,alpha,c0=c0,eps=eta)
    elif model=='RH': M=M6RHSparse(op,p,alpha,c0=c0,eta=eta)
    else: raise ValueError(model)
    inv0=css_invariants(U,op,h,p,alpha) if model=='CSS' else rh_invariants(U,x,h,p,alpha,op,t=0.0,c0=c0)
    steps=max(1,int(math.ceil(T/(dt_ratio*h)))); dt=T/steps; hist=[]; snaps=[]; t0=time.time(); fail=None
    snap_every=max(1,steps//20)
    for n in range(steps):
        U,ok,hh=newton_step(U,dt,M,tol=tol,maxit=maxit,refresh=1); hist.append(len(hh))
        if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U))>1e8:
            fail=n+1; break
        if save_snap and ((n+1)%snap_every==0 or n+1==steps): snaps.append((n+1)*dt,U.copy())
    Tr=(fail-1)*dt if fail else T
    # compare to initial after phase alignment; exact shift in moving coordinates is (c-c0)T
    ph,shape=phase_shape(U,ref0,h)
    expected=(c-c0)*Tr
    phase_err=min(abs(ph-expected),abs(ph-expected+L),abs(ph-expected-L))
    inv=css_invariants(U,op,h,p,alpha) if model=='CSS' else rh_invariants(U,x,h,p,alpha,op,t=Tr,c0=c0)
    out=dict(status='FAIL' if fail else 'OK',kind=kind,model=model,p=p,m=m,N=N,L=L,h=h,T=T,T_reached=Tr,dt=dt,dt_ratio=dt/h,
             c=c,c0=c0,eta=eta,phase_shift=ph,phase_error=phase_err,shape_relL2=shape,mean_newton=float(np.mean(hist)) if hist else None,max_newton=max(hist) if hist else None,runtime_s=time.time()-t0,
             inv0=inv0,inv=inv)
    if fail: out['fail_step']=fail
    if save_snap: out['snapshots']=[dict(t=float(t),u=u.tolist()) for t,u in snaps]
    out['x']=x.tolist(); out['u']=U.tolist(); out['u0']=ref0.tolist()
    return out


def pulse_sum(x,p,alpha,pulses):
    y=np.zeros_like(x)
    for d in pulses:
        y += exact_lobe(x,p,alpha,d['center'],d['c'],d.get('sign',1.0))
    return y


def periodic_dist(x,c,L): return (x-c+L/2)%L-L/2

def components_threshold(u,x,L,thr):
    mask=np.abs(u)>thr
    # components on periodic grid, return sign, max, center weighted abs, width points
    N=len(u); visited=np.zeros(N,bool); comps=[]
    for i in range(N):
        if not mask[i] or visited[i]: continue
        inds=[]; stack=[i]; visited[i]=True
        while stack:
            j=stack.pop(); inds.append(j)
            for k in ((j-1)%N,(j+1)%N):
                if mask[k] and not visited[k]: visited[k]=True; stack.append(k)
        arr=np.array(inds); vals=u[arr]; weights=np.abs(vals)
        # circular center via phases
        ang=2*np.pi*(x[arr]+L/2)/L; z=np.sum(weights*np.exp(1j*ang)); a=np.angle(z); cen=(a%(2*np.pi))*L/(2*np.pi)-L/2
        jpk=arr[np.argmax(np.abs(vals))]
        comps.append(dict(n=int(len(arr)),center=float(cen),peak=float(u[jpk]),maxabs=float(abs(u[jpk])),sign=int(np.sign(u[jpk]))))
    comps.sort(key=lambda d:d['maxabs'],reverse=True)
    return comps


def run_chase(model,p,m,T=6.0,dt_ratio=.05,eta=0.0,phase=.5,periods=12,cfast=2.0,cslow=1.0,separation=4.0,opposite=True,tol=1e-10,maxit=20,snap_times=None,c0=1.0):
    # need compatible grid with 'periods' harmonic periods
    x,h,L,N=compatible_grid(p,m,periods); alpha=p*p/((p+1)*(p+2)); op=PadeSparseOps(N,h)
    # put faster pulse behind (left); slower ahead (right)
    midpoint=phase*h
    c1=midpoint-separation/2; c2=midpoint+separation/2
    pulses=[dict(center=c1,c=cfast,sign=1.0),dict(center=c2,c=cslow,sign=-1.0 if opposite else 1.0)]
    U=pulse_sum(x,p,alpha,pulses)
    U0=U.copy()
    if model=='CSS': M=M1Sparse(op,p,alpha,c0=c0,eps=eta)
    elif model=='RH': M=M6RHSparse(op,p,alpha,c0=c0,eta=eta)
    else: raise ValueError
    inv0=css_invariants(U,op,h,p,alpha) if model=='CSS' else rh_invariants(U,x,h,p,alpha,op,t=0.0,c0=c0)
    steps=max(1,int(math.ceil(T/(dt_ratio*h)))); dt=T/steps
    if snap_times is None: snap_times=[0,2,3.5,4,4.5,5,6]
    targets=sorted(set([float(t) for t in snap_times if 0<=t<=T])); snaps=[]; itarget=0
    if targets and targets[0]==0: snaps.append((0.0,U.copy())); itarget=1
    hist=[]; t0=time.time(); fail=None
    for n in range(steps):
        U,ok,hh=newton_step(U,dt,M,tol=tol,maxit=maxit,refresh=1); hist.append(len(hh)); t=(n+1)*dt
        if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U))>1e8: fail=n+1; break
        while itarget<len(targets) and t+0.5*dt>=targets[itarget]: snaps.append((t,U.copy())); itarget+=1
    Tr=(fail-1)*dt if fail else T
    inv=css_invariants(U,op,h,p,alpha) if model=='CSS' else rh_invariants(U,x,h,p,alpha,op,t=Tr,c0=c0)
    A1=(cfast*(p+1)*(p+2)/2)**(1/p); A2=(cslow*(p+1)*(p+2)/2)**(1/p); thr=.02*min(A1,A2)
    snapdiag=[]
    for t,u in snaps:
        comps=components_threshold(u,x,L,thr)
        # independent-pulse reference before collision (formal superposition)
        pref=[dict(center=c1+(cfast-c0)*t,c=cfast,sign=1.0),dict(center=c2+(cslow-c0)*t,c=cslow,sign=-1.0 if opposite else 1.0)]
        ref=pulse_sum(x,p,alpha,pref)
        res=np.linalg.norm(u-ref)/max(np.linalg.norm(ref),1e-300)
        # radiation outside +- (pi+0.5) around two strongest components
        if len(comps)>=2:
            mask=np.ones(N,dtype=bool)
            for cc in comps[:2]: mask &= (np.abs(periodic_dist(x,cc['center'],L))>(math.pi+0.5))
            rad_l2=float(np.sqrt(h*np.sum(u[mask]**2)))
            rad_max=float(np.max(np.abs(u[mask]))) if np.any(mask) else 0.
        else: rad_l2=float('nan'); rad_max=float('nan')
        snapdiag.append(dict(t=float(t),ncomp=len(comps),components=comps[:8],independent_relL2=float(res),radiation_L2=rad_l2,radiation_max=rad_max))
    return dict(status='FAIL' if fail else 'OK',model=model,p=p,m=m,N=N,L=L,h=h,T=T,T_reached=Tr,dt=dt,dt_ratio=dt/h,eta=eta,c0=c0,cfast=cfast,cslow=cslow,separation=separation,opposite=opposite,
                mean_newton=float(np.mean(hist)) if hist else None,max_newton=max(hist) if hist else None,runtime_s=time.time()-t0,inv0=inv0,inv=inv,snapshots=snapdiag,x=x.tolist(),u=U.tolist(),u0=U0.tolist(), fail_step=fail)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='cmd',required=True)
    q=sub.add_parser('profile'); q.add_argument('--kind',choices=['lobe','dipole'],required=True); q.add_argument('--model',choices=['CSS','RH'],required=True); q.add_argument('--p',type=float,required=True); q.add_argument('--m',type=int,required=True); q.add_argument('--T',type=float,default=1); q.add_argument('--dt-ratio',type=float,default=.05); q.add_argument('--phase',type=float,default=.5); q.add_argument('--c',type=float,default=1); q.add_argument('--eta',type=float,default=0); q.add_argument('--periods',type=int,default=8); q.add_argument('--out')
    q=sub.add_parser('chase'); q.add_argument('--model',choices=['CSS','RH'],required=True); q.add_argument('--p',type=float,required=True); q.add_argument('--m',type=int,required=True); q.add_argument('--T',type=float,default=6); q.add_argument('--dt-ratio',type=float,default=.05); q.add_argument('--phase',type=float,default=.5); q.add_argument('--eta',type=float,default=0); q.add_argument('--periods',type=int,default=12); q.add_argument('--cfast',type=float,default=2); q.add_argument('--cslow',type=float,default=1); q.add_argument('--separation',type=float,default=4); q.add_argument('--same-sign',action='store_true'); q.add_argument('--out')
    a=ap.parse_args()
    if a.cmd=='profile': r=run_profile(a.kind,a.model,a.p,a.m,a.T,a.dt_ratio,a.phase,a.c,a.eta,a.periods)
    else: r=run_chase(a.model,a.p,a.m,a.T,a.dt_ratio,a.eta,a.phase,a.periods,a.cfast,a.cslow,a.separation,not a.same_sign)
    out={k:v for k,v in r.items() if k not in ('x','u','u0')}
    print(json.dumps(out,indent=2))
    if a.out: Path(a.out).parent.mkdir(parents=True,exist_ok=True); Path(a.out).write_text(json.dumps(r))
