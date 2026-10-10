#!/usr/bin/env python3
from __future__ import annotations
import math, json, argparse
from pathlib import Path
import numpy as np
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import spsolve
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from B1_pade_sparse_newton import PadeSparseOps


def diag(v): return diags(v,0,format='csr')

class M34Regularized:
    """Hamiltonian regularization r=1 or 2 of absCSS on l=p+2.
    Equation form used by R12: u_t + F=0 with F=-c0*u_x + D_x[Q+alpha A' u_x^2+2alpha A u_xx].
    Q=G' regularized consistently with H_epsilon,r.
    """
    def __init__(self, ops, p, alpha, reg_r, epsilon, c0=0.0):
        self.o=ops; self.p=float(p); self.l=float(p+2); self.alpha=float(alpha)
        self.rr=int(reg_r); self.eps=float(epsilon); self.c0=float(c0)
        if self.rr not in (1,2): raise ValueError('r must be 1 or 2')
        if self.eps<=0: raise ValueError('epsilon>0')

    def funcs(self,u):
        p,l,r,e = self.p,self.l,self.rr,self.eps
        S=u**(2*r)+e**(2*r)
        # A and first/second derivatives
        aexp=p/(2*r)
        A=S**aexp
        Ap=p*u**(2*r-1)*S**(aexp-1)
        # derivative of Ap
        App=p*((2*r-1)*u**(2*r-2)*S**(aexp-1)
               +(p-2*r)*u**(4*r-2)*S**(aexp-2))
        # Q=G'(u), G=((S)^(l/2r)-e^l)/(l(l-1))
        gexp=l/(2*r)
        Q=u**(2*r-1)*S**(gexp-1)/(l-1)
        Qp=((2*r-1)*u**(2*r-2)*S**(gexp-1)
            +(l-2*r)*u**(4*r-2)*S**(gexp-2))/(l-1)
        return A,Ap,App,Q,Qp

    def flux(self,u,q,uxx):
        A,Ap,_,Q,_=self.funcs(u)
        return Q + self.alpha*Ap*q*q + 2*self.alpha*A*uxx

    def flux_jac(self,u,q,uxx):
        A,Ap,App,Q,Qp=self.funcs(u)
        fu=Qp+self.alpha*App*q*q+2*self.alpha*Ap*uxx
        fq=2*self.alpha*Ap*q
        fr=2*self.alpha*A
        return fu,fq,fr

    def aux(self,u):
        q=self.o.D1(u); r2=self.o.D2(u); return q,r2

    def F(self,u):
        q,r2=self.aux(u)
        return -self.c0*q + self.o.D1(self.flux(u,q,r2))

    def initial_Z(self,U,W,dt):
        V=U+(dt/4)*self.F(U)
        qX,rX=self.aux(U); qV,rV=self.aux(V)
        return np.concatenate([U,qX,rX,V,qV,rV])

    def residual_jac(self,Z,W,dt,need_jac=True):
        N=self.o.N; o=self.o
        X=Z[:N]; qX=Z[N:2*N]; rX=Z[2*N:3*N]
        V=Z[3*N:4*N]; qV=Z[4*N:5*N]; rV=Z[5*N:6*N]
        PhiX=self.flux(X,qX,rX); PhiV=self.flux(V,qV,rV)
        k1=4/dt; k2=3/(2*dt)
        # F=-c0 D1 u - D1 Phi; multiply stage equations by A
        E=np.concatenate([
            o.A@qX-o.B@X,
            o.A@rX-o.B@qX,
            k1*(o.A@(V-X)) + self.c0*(o.B@X) - o.B@PhiX,
            o.A@qV-o.B@V,
            o.A@rV-o.B@qV,
            k2*(o.A@(W-X)) + self.c0*(o.B@V) - o.B@PhiV,
        ])
        if not need_jac: return E
        I=diags(np.ones(N),0,format='csr'); Z0=None
        fuX,fqX,frX=self.flux_jac(X,qX,rX)
        fuV,fqV,frV=self.flux_jac(V,qV,rV)
        J=bmat([
            [-o.B, o.A, Z0, Z0, Z0, Z0],
            [ Z0,-o.B, o.A, Z0, Z0, Z0],
            [-k1*o.A+self.c0*o.B-o.B@diag(fuX), -o.B@diag(fqX), -o.B@diag(frX), k1*o.A, Z0, Z0],
            [ Z0,Z0,Z0,-o.B,o.A,Z0],
            [ Z0,Z0,Z0,Z0,-o.B,o.A],
            [-k2*o.A,Z0,Z0, self.c0*o.B-o.B@diag(fuV),-o.B@diag(fqV),-o.B@diag(frV)],
        ],format='csr')
        return E,J

    def extract_X(self,Z): return Z[:self.o.N]


def newton_stage(U,dt,model,tol=1e-10,maxit=12,line_search=10):
    FU=model.F(U); W=U-(dt/3)*FU; Z=model.initial_Z(U,W,dt)
    hist=[]
    for it in range(maxit):
        E,J=model.residual_jac(Z,W,dt,True); rn=float(np.linalg.norm(E,np.inf)); hist.append(rn)
        if rn<tol: return model.extract_X(Z),True,hist
        try: d=spsolve(J,-E,permc_spec='COLAMD',use_umfpack=False)
        except Exception: return model.extract_X(Z),False,hist
        if not np.all(np.isfinite(d)): return model.extract_X(Z),False,hist
        lam=1.0; accepted=False
        for _ in range(line_search):
            Zt=Z+lam*d; rt=float(np.linalg.norm(model.residual_jac(Zt,W,dt,False),np.inf))
            if rt<rn: Z=Zt; accepted=True; break
            lam*=0.5
        if not accepted: return model.extract_X(Z),False,hist
    return model.extract_X(Z),False,hist


def u0_sine(x,a=1.0,L=2.0):
    k=2*np.pi/L
    return (a/k)*np.sin(k*x)


def smooth_switch01(t):
    t=np.asarray(t,float); out=np.zeros_like(t)
    out[t>=1]=1.0
    m=(t>0)&(t<1)
    tm=t[m]
    a=np.exp(-1.0/tm); b=np.exp(-1.0/(1.0-tm))
    out[m]=a/(a+b)
    return out


def u0_prepared_periodic(x,a=1.0,w=0.25,W=0.75,b=1.0,L=2.0):
    # For L=2: odd periodic profile, exactly a*x near zero and exactly
    # b*(1-|x|) with the appropriate sign near the periodic crossing.
    if abs(L-2.0)>1e-12:
        # scale coordinates to y in [-1,1) while preserving slope a at 0
        y=2*x/L; scale=L/2
    else:
        y=x; scale=1.0
    s=np.abs(y); left=a*scale*s; right=b*scale*(1.0-s)
    f=left.copy(); f[s>=W]=right[s>=W]
    m=(s>w)&(s<W); t=(s[m]-w)/(W-w); sw=smooth_switch01(t)
    f[m]=(1-sw)*left[m]+sw*right[m]
    return np.sign(y)*f

def u0_prepared(x,a=1.0,w=0.25,W=0.75):
    ax=np.abs(x); chi=np.ones_like(x)
    chi[ax>=W]=0.0
    m=(ax>w)&(ax<W)
    t=(W-ax[m])/(W-w)
    chi[m]=smooth_switch01(t)
    return a*x*chi

def periodic_delta(x,s,L): return (x-s+L/2)%L-L/2

def zero_candidates(x,u,L):
    N=len(x); zs=[]
    for i in range(N):
        j=(i+1)%N; ui,uj=u[i],u[j]
        if ui==0: z=x[i]
        elif ui*uj<0:
            xj=x[j] if j!=0 else x[0]+L
            z=x[i]-ui*(xj-x[i])/(uj-ui); z=((z+L/2)%L)-L/2
        else: continue
        if not any(abs(periodic_delta(z,z0,L))<1e-10 for z0 in zs): zs.append(float(z))
    return sorted(zs)

def locate_zero(x,u,guess,L):
    # sign-changing crossing closest to guess, periodic segments including boundary
    cand=[(abs(periodic_delta(z,guess,L)),z) for z in zero_candidates(x,u,L)]
    if not cand:
        i=int(np.argmin(np.abs(u))); return float(x[i])
    return float(min(cand,key=lambda t:t[0])[1])

def _old_locate_zero_unused(x,u,guess,L):
    N=len(x); cand=[]
    for i in range(N):
        j=(i+1)%N
        ui,uj=u[i],u[j]
        if ui==0: z=x[i]
        elif ui*uj<0:
            xj=x[j]
            if j==0: xj=x[0]+L
            z=x[i]-ui*(xj-x[i])/(uj-ui)
            # wrap
            z=((z+L/2)%L)-L/2
        else: continue
        d=abs(periodic_delta(z,guess,L)); cand.append((d,z))
    if not cand:
        i=int(np.argmin(np.abs(u))); return float(x[i])
    return float(min(cand,key=lambda t:t[0])[1])

def local_slope(x,u,s,L,R):
    d=periodic_delta(x,s,L); m=np.abs(d)<=R
    if np.sum(m)<5: return float('nan')
    cc=np.polyfit(d[m],u[m],3); return float(cc[-2])

def fit_exponent(x,u,s,L,R,h,p):
    d=np.abs(periodic_delta(x,s,L)); lo=max(2*h,0.05*R); m=(d>=lo)&(d<=R)&(np.abs(u)>1e-14)
    if np.sum(m)<8: return float('nan'),float('nan')
    ru=float(np.polyfit(np.log(d[m]),np.log(np.abs(u[m])),1)[0])
    v=np.sign(u[m])*np.abs(u[m])**(p/2)
    rv=float(np.polyfit(np.log(d[m]),np.log(np.abs(v)),1)[0])
    return ru,rv

def invariants(x,u,model):
    dx=model.o.dx; p=model.p; l=model.l; r=model.rr; e=model.eps
    q=model.o.D1(u); A,_,_,_,_=model.funcs(u)
    S=u**(2*r)+e**(2*r)
    G=(S**(l/(2*r))-e**l)/(l*(l-1))
    H=float(dx*np.sum(model.alpha*A*q*q-G))
    C=float(dx*np.sum(u)); P=float(0.5*dx*np.sum(u*u))
    return C,P,H

def run_case(p=1.5,reg_r=1,eps=0.1,points=24,L=2.0,a=1.0,T=0.1,dt_ratio=0.05,c0=0.0,out=None,sample_count=80,initial='prepared'):
    htarget=eps/(points*abs(a)); N=int(round(L/htarget)); N += N%2 # even
    dx=L/N; x=-L/2+dx*np.arange(N)
    alpha=p*p/((p+1)*(p+2)); ops=PadeSparseOps(N,dx); model=M34Regularized(ops,p,alpha,reg_r,eps,c0)
    u=(u0_prepared_periodic(x,a,L=L) if initial=='prepared2' else (u0_prepared(x,a) if initial=='prepared' else u0_sine(x,a,L))); uinit=u.copy(); C0,P0,H0=invariants(x,u,model)
    dt0=dt_ratio*dx; nsteps=max(1,int(math.ceil(T/dt0))); dt=T/nsteps
    # dense-at-start sampling: geometric targets + linear later
    tg=np.unique(np.concatenate(([0.0],np.geomspace(max(dt,1e-12),max(min(T,0.2*T),dt),max(10,sample_count//2)),np.linspace(0,T,max(10,sample_count//2)))))
    sample_steps=sorted(set(int(np.clip(round(t/dt),0,nsteps)) for t in tg if t<=T+1e-15)); sample_set=set(sample_steps)
    rows=[]; snaps_t=[]; snaps_u=[]; guess=0.0
    def record(t,u,guess):
        s=locate_zero(x,u,guess,L); slope=local_slope(x,u,s,L,max(4*dx,0.25*eps))
        zs=zero_candidates(x,u,L); rec={'t':float(t),'zero_x':s,'slope':slope,'zero_count':len(zs),'zero_candidates':zs}
        for mult in (2,4,8):
            R=min(mult*eps,0.35*L)
            ru,rv=fit_exponent(x,u,s,L,R,dx,p)
            rec[f'r_u_R{mult}eps']=ru; rec[f'r_v_R{mult}eps']=rv
        C,P,H=invariants(x,u,model); rec.update(C=C,P=P,H=H,dP_rel=(P-P0)/P0,dH_rel=(H-H0)/(abs(H0)+1e-300))
        return s,rec
    guess,rec=record(0,u,guess); rows.append(rec); snaps_t.append(0.0); snaps_u.append(u.copy())
    nit=[]; ok=True; fail_step=None
    for n in range(1,nsteps+1):
        un,good,hist=newton_stage(u,dt,model,tol=1e-10,maxit=14); nit.append(len(hist));
        if not good:
            ok=False; fail_step=n; break
        u=un; t=n*dt
        if n in sample_set:
            guess,rec=record(t,u,guess); rows.append(rec); snaps_t.append(t); snaps_u.append(u.copy())
    # always final accepted
    tf=(fail_step-1)*dt if fail_step else T
    if rows[-1]['t']<tf-0.5*dt:
        guess,rec=record(tf,u,guess); rows.append(rec); snaps_t.append(tf); snaps_u.append(u.copy())
    # initial exact zero velocity from grid F and analytic local formula
    F0=model.F(uinit); i0=int(np.argmin(np.abs(x)))
    ut0=float(-F0[i0]); slope0=a
    vgrid=-ut0/slope0
    if reg_r==1:
        vanal=eps**p/(p+1) + alpha*p*a*a*eps**(p-2)
    else:
        vanal=0.0
    result={'p':p,'r':reg_r,'epsilon':eps,'points_per_epsilon':eps/(abs(a)*dx),'N':N,'L':L,'dx':dx,'dt':dt,'dt_over_dx':dt/dx,'T_target':T,'T_reached':tf,'ok':ok,'fail_step':fail_step,'alpha':alpha,'a':a,'c0':c0,'initial':initial,'C0':C0,'P0':P0,'H0':H0,'initial_ut_at_zero_grid':ut0,'initial_zero_speed_grid':vgrid,'initial_zero_speed_local_exact':vanal,'mean_newton_iters':float(np.mean(nit)) if nit else None,'max_newton_iters':int(max(nit)) if nit else None,'timeseries':rows}
    if out:
        Path(out).parent.mkdir(parents=True,exist_ok=True); Path(out).write_text(json.dumps(result,indent=2)); np.savez_compressed(str(Path(out).with_suffix('.npz')),x=x,t=np.array(snaps_t),u=np.array(snaps_u))
    return result

def main():
    ap=argparse.ArgumentParser();
    ap.add_argument('--p',type=float,default=1.5); ap.add_argument('--initial',choices=['prepared2','prepared','sine'],default='prepared2'); ap.add_argument('--r',type=int,default=1); ap.add_argument('--epsilon',type=float,default=0.1); ap.add_argument('--points',type=float,default=24); ap.add_argument('--L',type=float,default=2.0); ap.add_argument('--a',type=float,default=1.0); ap.add_argument('--T',type=float,default=0.1); ap.add_argument('--dt-ratio',type=float,default=0.05); ap.add_argument('--out',type=str,default=None)
    z=ap.parse_args(); res=run_case(z.p,z.r,z.epsilon,z.points,z.L,z.a,z.T,z.dt_ratio,0,z.out,initial=z.initial); print(json.dumps({k:res[k] for k in res if k!='timeseries'},indent=2))
if __name__=='__main__': main()
