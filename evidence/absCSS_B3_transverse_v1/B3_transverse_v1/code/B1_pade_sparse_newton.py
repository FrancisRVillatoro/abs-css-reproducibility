#!/usr/bin/env python3
"""Scalable sparse Newton solver for absCSS benchmark B1.

Exact reformulation of the Padé/Petrov--Galerkin + R_{1,2} pilot.
It avoids dense A^{-1} Jacobians by introducing Padé auxiliary variables.
The nonlinear Jacobian updates are diagonal; A,B,C,D remain sparse cyclic
five-point operators.
"""
import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.sparse import bmat, csr_matrix, diags
from scipy.sparse.linalg import spsolve


def cyclic_stencil(N, offsets, coeffs):
    rows = []
    cols = []
    data = []
    ii = np.arange(N, dtype=np.int64)
    for off, coeff in zip(offsets, coeffs):
        rows.append(ii)
        cols.append((ii + off) % N)
        data.append(np.full(N, coeff, dtype=float))
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    data = np.concatenate(data)
    return csr_matrix((data, (rows, cols)), shape=(N, N))


class PadeSparseOps:
    def __init__(self, N, dx):
        self.N = N
        self.dx = dx
        # A symbol: (33 + 26 cos(theta) + cos(2 theta))/60
        self.A = cyclic_stencil(
            N, [-2, -1, 0, 1, 2],
            [1/120, 13/60, 33/60, 13/60, 1/120],
        )
        self.B = cyclic_stencil(
            N, [-2, -1, 1, 2],
            [-1/(24*dx), -10/(24*dx), 10/(24*dx), 1/(24*dx)],
        )
        self.C = cyclic_stencil(
            N, [-2, -1, 1, 2],
            [-1/(2*dx**3), 2/(2*dx**3), -2/(2*dx**3), 1/(2*dx**3)],
        )
        self.D = cyclic_stencil(
            N, [-2, -1, 0, 1, 2],
            [1/dx**4, -4/dx**4, 6/dx**4, -4/dx**4, 1/dx**4],
        )
        self.A2 = (self.A @ self.A).tocsr()
        self.AB = (self.A @ self.B).tocsr()
        self.B2 = (self.B @ self.B).tocsr()
        self.AD = (self.A @ self.D).tocsr()

        th = 2*np.pi*np.fft.fftfreq(N)
        asym = (33.0 + 26.0*np.cos(th) + np.cos(2*th))/60.0
        bsym = 1j*(10*np.sin(th) + np.sin(2*th))/(12*dx)
        csym = 1j*(np.sin(2*th) - 2*np.sin(th))/(dx**3)
        dsym = 16*np.sin(th/2)**4/(dx**4)
        self.s1 = bsym/asym
        self.s2 = self.s1*self.s1
        self.s3 = csym/asym
        self.s4 = dsym/asym

    @staticmethod
    def _app(sym, z):
        return np.fft.ifft(sym*np.fft.fft(z)).real

    def D1(self, z): return self._app(self.s1, z)
    def D2(self, z): return self._app(self.s2, z)
    def D3(self, z): return self._app(self.s3, z)
    def D4(self, z): return self._app(self.s4, z)


def compact_u(x, p, alpha, x0=-5.0, c=1.0):
    A = (c*(p+1)*(p+2)/2.0)**(1.0/p)
    k = p/(2*np.sqrt(alpha*(p+1)*(p+2)))
    z = x-x0
    y = np.zeros_like(x)
    m = np.abs(k*z) <= np.pi/2
    y[m] = A*np.cos(k*z[m])**(2.0/p)
    return y


def compact_v(x, p, alpha, x0=-5.0, c=1.0):
    B = np.sqrt(c*(p+1)*(p+2)/2.0)
    k = p/(2*np.sqrt(alpha*(p+1)*(p+2)))
    z = x-x0
    y = np.zeros_like(x)
    m = np.abs(k*z) <= np.pi/2
    y[m] = B*np.cos(k*z[m])
    return y


def diag(v):
    return diags(v, 0, format='csr')


class M1Sparse:
    """M1 after multiplying F by A^2 and introducing q=D1 u."""
    def __init__(self, ops, p, alpha, c0, eps):
        self.o = ops
        self.p = p
        self.a = alpha
        self.c0 = c0
        self.eps = eps

    def q(self, u):
        return self.o.D1(u)

    def R(self, u, q):
        o, p, a = self.o, self.p, self.a
        au = np.abs(u)
        f = au**p*u/(p+1.0)
        phi = np.sign(u)*au**(p-1.0)
        g = phi*q*q
        h = au**p*q
        return (-self.c0*(o.AB@u) + o.AB@f - a*p*(o.AB@g)
                + 2*a*(o.B2@h) + self.eps*(o.AD@u))

    def Rjac(self, u, q):
        o, p, a = self.o, self.p, self.a
        au = np.abs(u)
        sg = np.sign(u)
        fp = au**p
        phi = sg*au**(p-1.0)
        # Generalised derivative used by the dense pilot: 0 at the exact zero set.
        phip = np.zeros_like(u)
        m = au > 1e-13
        # At p=2, phi(u)=u is smooth and phi'(0)=1 exactly.  The generic
        # zero-set convention phip=0 is reserved for the nonsmooth/singular
        # cases p<2; using it at p=2 can create a spurious Newton defect.
        if abs(p-2.0) < 1e-14:
            phip.fill(1.0)
        elif p != 1.0:
            phip[m] = (p-1.0)*au[m]**(p-2.0)
        hp_u = p*sg*au**(p-1.0)*q
        Ru = (-self.c0*o.AB + o.AB@diag(fp)
              - a*p*(o.AB@diag(phip*q*q))
              + 2*a*(o.B2@diag(hp_u)) + self.eps*o.AD)
        Rq = (-a*p*(o.AB@diag(2*phi*q))
              + 2*a*(o.B2@diag(fp)))
        return Ru.tocsr(), Rq.tocsr()

    def F(self, u):
        o, p, a = self.o, self.p, self.a
        q = o.D1(u)
        au = np.abs(u)
        phi = np.sign(u)*au**(p-1.0)
        return (-self.c0*q + o.D1(au**p*u/(p+1.0))
                - a*p*o.D1(phi*q*q)
                + 2*a*o.D2(au**p*q) + self.eps*o.D4(u))

    def residual_jac(self, Z, W, dt, need_jac=True):
        N, o = self.o.N, self.o
        X = Z[0:N]; qX = Z[N:2*N]; V = Z[2*N:3*N]; qV = Z[3*N:4*N]
        k1 = 4.0/dt
        k2 = 3.0/(2.0*dt)
        E1 = o.A@qX - o.B@X
        E2 = k1*(o.A2@(V-X)) - self.R(X, qX)
        E3 = o.A@qV - o.B@V
        E4 = k2*(o.A2@(W-X)) - self.R(V, qV)
        E = np.concatenate([E1, E2, E3, E4])
        if not need_jac:
            return E
        RuX, RqX = self.Rjac(X, qX)
        RuV, RqV = self.Rjac(V, qV)
        J = bmat([
            [-o.B,                  o.A,               None,               None],
            [-k1*o.A2-RuX,          -RqX,               k1*o.A2,            None],
            [None,                  None,               -o.B,               o.A],
            [-k2*o.A2,              None,               -RuV,               -RqV],
        ], format='csr')
        return E, J

    def initial_Z(self, U, W, dt):
        FX = self.F(U)
        V = U + (dt/4.0)*FX
        return np.concatenate([U.copy(), self.q(U), V, self.q(V)])

    def extract_X(self, Z):
        return Z[:self.o.N]


class M2Sparse:
    """M2 with q=D1v, r=D2v, z=D3v, w=D4v as local Padé auxiliaries."""
    def __init__(self, ops, p, alpha, c0, eps):
        self.o = ops
        self.p = p
        self.a = alpha
        self.c0 = c0
        self.eps = eps
        self.aa = 1.0 + 6.0/p
        self.bb = 4.0/(p*p)

    def aux(self, v):
        o = self.o
        return o.D1(v), o.D2(v), o.D3(v), o.D4(v)

    def S(self, v, q, r, z, w):
        return (-self.c0*q + v*v*q
                + 2*self.a*(v*v*z + self.aa*v*q*r + self.bb*q**3)
                + self.eps*w)

    def Sjac_diag(self, v, q, r, z):
        a = self.a
        sv = 2*v*q + 2*a*(2*v*z + self.aa*q*r)
        sq = -self.c0 + v*v + 2*a*(self.aa*v*r + 3*self.bb*q*q)
        sr = 2*a*self.aa*v*q
        sz = 2*a*v*v
        return sv, sq, sr, sz

    def F(self, v):
        q, r, z, w = self.aux(v)
        return self.S(v, q, r, z, w)

    def residual_jac(self, Z, W, dt, need_jac=True):
        N, o = self.o.N, self.o
        # X,q,r,z,w,V,q,r,z,w
        X = Z[0*N:1*N]; qX = Z[1*N:2*N]; rX = Z[2*N:3*N]; zX = Z[3*N:4*N]; wX = Z[4*N:5*N]
        V = Z[5*N:6*N]; qV = Z[6*N:7*N]; rV = Z[7*N:8*N]; zV = Z[8*N:9*N]; wV = Z[9*N:10*N]
        SX = self.S(X,qX,rX,zX,wX)
        SV = self.S(V,qV,rV,zV,wV)
        E = np.concatenate([
            o.A@qX-o.B@X,
            o.A@rX-o.B@qX,
            o.A@zX-o.C@X,
            o.A@wX-o.D@X,
            V-X-(dt/4.0)*SX,
            o.A@qV-o.B@V,
            o.A@rV-o.B@qV,
            o.A@zV-o.C@V,
            o.A@wV-o.D@V,
            X-W+(2*dt/3.0)*SV,
        ])
        if not need_jac:
            return E
        svX,sqX,srX,szX = self.Sjac_diag(X,qX,rX,zX)
        svV,sqV,srV,szV = self.Sjac_diag(V,qV,rV,zV)
        I = diags(np.ones(N),0,format='csr')
        dSX = [diag(svX),diag(sqX),diag(srX),diag(szX),self.eps*I]
        dSV = [diag(svV),diag(sqV),diag(srV),diag(szV),self.eps*I]
        q = dt/4.0
        r = 2*dt/3.0
        Z0 = None
        J = bmat([
            [-o.B, o.A,   Z0,   Z0,   Z0,      Z0, Z0, Z0, Z0, Z0],
            [ Z0, -o.B,  o.A,   Z0,   Z0,      Z0, Z0, Z0, Z0, Z0],
            [-o.C, Z0,    Z0,  o.A,   Z0,      Z0, Z0, Z0, Z0, Z0],
            [-o.D, Z0,    Z0,   Z0,  o.A,      Z0, Z0, Z0, Z0, Z0],
            [-I-q*dSX[0], -q*dSX[1], -q*dSX[2], -q*dSX[3], -q*dSX[4], I, Z0, Z0, Z0, Z0],
            [ Z0, Z0, Z0, Z0, Z0,                 -o.B, o.A, Z0, Z0, Z0],
            [ Z0, Z0, Z0, Z0, Z0,                  Z0, -o.B, o.A, Z0, Z0],
            [ Z0, Z0, Z0, Z0, Z0,                 -o.C, Z0, Z0, o.A, Z0],
            [ Z0, Z0, Z0, Z0, Z0,                 -o.D, Z0, Z0, Z0, o.A],
            [ I, Z0, Z0, Z0, Z0, r*dSV[0], r*dSV[1], r*dSV[2], r*dSV[3], r*dSV[4]],
        ], format='csr')
        return E, J

    def initial_Z(self, U, W, dt):
        FX = self.F(U)
        V = U + (dt/4.0)*FX
        qX,rX,zX,wX = self.aux(U)
        qV,rV,zV,wV = self.aux(V)
        return np.concatenate([U,qX,rX,zX,wX,V,qV,rV,zV,wV])

    def extract_X(self, Z):
        return Z[:self.o.N]


def newton_stage(U, dt, model, tol=1e-10, maxit=12, line_search=8):
    FU = model.F(U)
    W = U - (dt/3.0)*FU
    Z = model.initial_Z(U,W,dt)
    hist = []
    nnz = []
    for it in range(maxit):
        E, J = model.residual_jac(Z,W,dt,True)
        r = float(np.linalg.norm(E,np.inf))
        hist.append(r)
        nnz.append(int(J.nnz))
        if r < tol:
            return model.extract_X(Z), True, hist, nnz
        try:
            d = spsolve(J, -E, permc_spec='COLAMD', use_umfpack=False)
        except Exception:
            return model.extract_X(Z), False, hist, nnz
        if not np.all(np.isfinite(d)):
            return model.extract_X(Z), False, hist, nnz
        lam = 1.0
        accepted = False
        for _ in range(line_search):
            Zt = Z + lam*d
            rt = float(np.linalg.norm(model.residual_jac(Zt,W,dt,False),np.inf))
            if rt < r:
                Z = Zt
                accepted = True
                break
            lam *= 0.5
        if not accepted:
            return model.extract_X(Z), False, hist, nnz
    return model.extract_X(Z), False, hist, nnz


def phase_shape(y, ref, dx):
    Y=np.fft.fft(y); R=np.fft.fft(ref); cc=np.fft.ifft(Y*np.conj(R)).real
    j=int(np.argmax(cc)); n=len(y)
    jm=(j-1)%n; jp=(j+1)%n; den=cc[jm]-2*cc[j]+cc[jp]
    delta=0.0 if den==0 else 0.5*(cc[jm]-cc[jp])/den
    delta=float(np.clip(delta,-0.5,0.5)); lag=j+delta
    if lag>n/2: lag-=n
    k=2*np.pi*np.fft.fftfreq(n,d=dx)
    aligned=np.fft.ifft(R*np.exp(-1j*k*(lag*dx))).real
    return float(lag*dx), float(np.linalg.norm(y-aligned)/np.linalg.norm(ref))


def invariants_u(u,op,dx,p,alpha):
    ux=op.D1(u); au=np.abs(u)
    return (dx*np.sum(u), 0.5*dx*np.sum(u*u),
            dx*np.sum(alpha*au**p*ux**2-au**(p+2)/((p+1)*(p+2))),
            dx*np.sum(np.maximum(-u,0)))


def run(model_name,N,p,T,L=40.0,c0=1.0,dt_ratio=0.1,mu=0.0,eps=None,x0=-5.0,
        newton_tol=1e-10,newton_max=12):
    dx=L/N
    alpha=p*p/((p+1)*(p+2))
    eps = mu*dx*dx if eps is None else eps
    x=-L/2+dx*np.arange(N)
    o=PadeSparseOps(N,dx)
    if model_name=='M1':
        U=compact_u(x,p,alpha,x0); ref=U.copy(); M=M1Sparse(o,p,alpha,c0,eps)
    else:
        U=compact_v(x,p,alpha,x0); ref=U.copy(); M=M2Sparse(o,p,alpha,c0,eps)
    U0=U.copy()
    steps=max(1,int(np.ceil(T/(dt_ratio*dx))))
    dt=T/steps
    nit=[]; nnz=[]; t0=time.time()
    for n in range(steps):
        U,ok,h,nz=newton_stage(U,dt,M,newton_tol,newton_max)
        nit.append(len(h)); nnz += nz
        if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U))>1e6:
            return dict(status='FAIL',model=model_name,N=N,p=p,T=T,dx=dx,dt=dt,dt_ratio=dt/dx,
                        mu=mu,eps=eps,step=n+1,residual=h,runtime_s=time.time()-t0,
                        mean_newton=float(np.mean(nit)),max_jac_nnz=max(nnz) if nnz else 0)
    ph,sh=phase_shape(U,ref,dx)
    R=dict(status='OK',model=model_name,N=N,p=p,T=T,dx=dx,dt=dt,dt_ratio=dt/dx,
           mu=mu,eps=eps,steps=steps,shape_relL2=sh,phase=ph,runtime_s=time.time()-t0,
           mean_newton=float(np.mean(nit)),max_newton=int(max(nit)),max_jac_nnz=int(max(nnz) if nnz else 0),
           max_abs=float(np.max(np.abs(U))))
    if model_name=='M1':
        z0=invariants_u(U0,o,dx,p,alpha); z=invariants_u(U,o,dx,p,alpha)
        R.update(drift_C=(z[0]-z0[0])/abs(z0[0]), drift_P=(z[1]-z0[1])/abs(z0[1]),
                 drift_H=(z[2]-z0[2])/max(abs(z0[2]),1e-300), negative_mass=z[3])
    else:
        u=np.sign(U)*np.abs(U)**(2/p); ur=compact_u(x,p,alpha,x0)
        _,e=phase_shape(u,ur,dx); R['shape_relL2_u']=e
    return R


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model',choices=['M1','M2'],required=True)
    ap.add_argument('--N',type=int,required=True)
    ap.add_argument('--p',type=float,default=1.5)
    ap.add_argument('--T',type=float,default=.2)
    ap.add_argument('--dt-ratio',type=float,default=.1)
    ap.add_argument('--mu',type=float,default=0.0)
    ap.add_argument('--eps',type=float,default=None)
    ap.add_argument('--newton-tol',type=float,default=1e-10)
    ap.add_argument('--newton-max',type=int,default=12)
    ap.add_argument('--out')
    a=ap.parse_args()
    res=run(a.model,a.N,a.p,a.T,dt_ratio=a.dt_ratio,mu=a.mu,eps=a.eps,
            newton_tol=a.newton_tol,newton_max=a.newton_max)
    print(json.dumps(res,indent=2))
    if a.out: Path(a.out).write_text(json.dumps(res,indent=2))

if __name__=='__main__':
    main()
