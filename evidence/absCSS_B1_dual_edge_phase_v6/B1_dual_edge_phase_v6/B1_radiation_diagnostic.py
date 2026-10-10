#!/usr/bin/env python3
"""B1 compacton radiation diagnostic for the scalable absCSS solver.

Measures numerical radiation after best-fit phase alignment of the exact compacton.
The compacton is stationary in the c0=1 comoving frame and has half-support pi.
Two exterior regions are reported:
  * stencil guard: outside aligned support plus 4*dx;
  * fixed guard:   outside aligned support plus 0.25 and 0.50.
The fixed guards distinguish emitted radiation from the shrinking edge stencil layer.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np

from B1_pade_sparse_newton import (
    PadeSparseOps, M1Sparse, compact_u, phase_shape, invariants_u
)
from scipy.sparse.linalg import spsolve



def newton_stage_interleaved(U, dt, model, tol=1e-10, maxit=12, line_search=8):
    """Algebraically identical Newton step, solved after node-wise permutation.

    Old ordering is [field_1(all nodes), field_2(all nodes), ...].  For M1 the
    four field/equation blocks are permuted to [field_1(i),...,field_4(i)]
    node by node.  This exposes the local cyclic stencil bandwidth and reduces
    sparse-LU fill without altering the nonlinear system.
    """
    FU = model.F(U)
    W = U - (dt/3.0)*FU
    Z = model.initial_Z(U,W,dt)
    N = model.o.N
    nblocks = len(Z)//N
    perm = np.arange(nblocks*N).reshape(nblocks,N).T.ravel()
    hist=[]
    for _it in range(maxit):
        E,J = model.residual_jac(Z,W,dt,True)
        rr=float(np.linalg.norm(E,np.inf)); hist.append(rr)
        if rr < tol:
            return model.extract_X(Z), True, hist
        Jp=J[perm,:][:,perm].tocsr(); Ep=E[perm]
        try:
            dp=spsolve(Jp,-Ep,permc_spec='NATURAL',use_umfpack=False)
        except Exception:
            return model.extract_X(Z), False, hist
        if not np.all(np.isfinite(dp)):
            return model.extract_X(Z), False, hist
        d=np.empty_like(dp); d[perm]=dp
        lam=1.0; accepted=False
        for _ in range(line_search):
            Zt=Z+lam*d
            rt=float(np.linalg.norm(model.residual_jac(Zt,W,dt,False),np.inf))
            if rt < rr:
                Z=Zt; accepted=True; break
            lam*=0.5
        if not accepted:
            return model.extract_X(Z), False, hist
    return model.extract_X(Z), False, hist

def periodic_distance(x, center, L):
    return (x - center + 0.5*L) % L - 0.5*L


def shifted_reference(x, p, alpha, center, c=1.0):
    return compact_u(x, p, alpha, x0=center, c=c)


def exterior_metrics(u, x, p, alpha, x0, L, dx):
    ref0 = compact_u(x, p, alpha, x0=x0)
    phase, shape = phase_shape(u, ref0, dx)
    center = x0 + phase
    ref = shifted_reference(x, p, alpha, center)
    dist = np.abs(periodic_distance(x, center, L))
    refnorm = np.linalg.norm(ref)
    amp = float(np.max(np.abs(ref)))

    out = {
        'phase': float(phase),
        'shape_relL2_rebuilt': float(np.linalg.norm(u-ref)/refnorm),
        'aligned_center': float(center),
    }
    guards = [('edge0', 0.0), ('stencil4', 4.0*dx), ('fixed025', 0.25), ('fixed050', 0.50)]
    for name, guard in guards:
        mask = dist > (math.pi + guard)
        vals = u[mask]
        out[f'{name}_npts'] = int(mask.sum())
        out[f'{name}_max_abs'] = float(np.max(np.abs(vals))) if vals.size else 0.0
        out[f'{name}_max_relA'] = out[f'{name}_max_abs']/amp if amp else 0.0
        out[f'{name}_L2_rel'] = float(np.linalg.norm(vals)/refnorm) if vals.size else 0.0
        out[f'{name}_L1_abs'] = float(dx*np.sum(np.abs(vals))) if vals.size else 0.0
        out[f'{name}_negative_mass'] = float(dx*np.sum(np.maximum(-vals,0.0))) if vals.size else 0.0
    return out, ref


def run_case(N, p, T, method, L=40.0, c0=1.0, dt_ratio=0.1,
             eta_coeff=0.05, x0=-5.0, newton_tol=1e-10, newton_max=12,
             save_profile=None):
    dx = L/N
    alpha = p*p/((p+1)*(p+2))
    if method == 'M1':
        eps = 0.0
    elif method == 'M5':
        eps = eta_coeff*dx**3
    else:
        raise ValueError(method)
    mu_equiv = eps/(dx*dx)
    x = -L/2 + dx*np.arange(N)
    op = PadeSparseOps(N, dx)
    U = compact_u(x, p, alpha, x0=x0)
    U0 = U.copy()
    model = M1Sparse(op, p, alpha, c0, eps)
    steps = max(1, int(np.ceil(T/(dt_ratio*dx))))
    dt = T/steps
    nit = []
    t0 = time.time()
    for n in range(steps):
        U, ok, hist = newton_stage_interleaved(U, dt, model, newton_tol, newton_max)
        nit.append(len(hist))
        if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U)) > 1e6:
            return {
                'status':'FAIL','method':method,'N':N,'p':p,'T':T,'dx':dx,'dt':dt,
                'dt_ratio':dt/dx,'eta':eps,'mu_equiv':mu_equiv,'step':n+1,
                'residual':hist,'runtime_s':time.time()-t0,
                'mean_newton':float(np.mean(nit)) if nit else 0.0,
            }
    met, ref = exterior_metrics(U, x, p, alpha, x0, L, dx)
    z0 = invariants_u(U0, op, dx, p, alpha)
    z = invariants_u(U, op, dx, p, alpha)
    result = {
        'status':'OK','method':method,'N':N,'p':p,'T':T,'L':L,'dx':dx,'dt':dt,
        'dt_ratio':dt/dx,'steps':steps,'eta':eps,'eta_coeff_dx3':eta_coeff if method=='M5' else 0.0,
        'mu_equiv':mu_equiv,'runtime_s':time.time()-t0,
        'mean_newton':float(np.mean(nit)),'max_newton':int(max(nit)),
        'drift_C':float((z[0]-z0[0])/abs(z0[0])),
        'drift_P':float((z[1]-z0[1])/abs(z0[1])),
        'drift_H':float((z[2]-z0[2])/max(abs(z0[2]),1e-300)),
        'negative_mass_total':float(z[3]),
        **met,
    }
    if save_profile:
        path = Path(save_profile)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['x','u_final','u_exact_phase_aligned','residual'])
            for xi, ui, ri in zip(x, U, ref):
                w.writerow([f'{xi:.17g}', f'{ui:.17g}', f'{ri:.17g}', f'{ui-ri:.17g}'])
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=int, required=True)
    ap.add_argument('--p', type=float, required=True)
    ap.add_argument('--T', type=float, default=1.0)
    ap.add_argument('--method', choices=['M1','M5'], required=True)
    ap.add_argument('--dt-ratio', type=float, default=0.1)
    ap.add_argument('--eta-coeff', type=float, default=0.05,
                    help='M5 uses eta = eta_coeff * dx^3')
    ap.add_argument('--newton-tol', type=float, default=1e-10)
    ap.add_argument('--newton-max', type=int, default=12)
    ap.add_argument('--out')
    ap.add_argument('--profile')
    a = ap.parse_args()
    r = run_case(a.N, a.p, a.T, a.method, dt_ratio=a.dt_ratio,
                 eta_coeff=a.eta_coeff, newton_tol=a.newton_tol,
                 newton_max=a.newton_max, save_profile=a.profile)
    s = json.dumps(r, indent=2)
    print(s)
    if a.out:
        Path(a.out).write_text(s+'\n')

if __name__ == '__main__':
    main()
