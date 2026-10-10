#!/usr/bin/env python3
"""Symbolic kernel and translation checks plus independent high-precision points.
Point evaluations are regression tests, NOT part of the proof of interval coverage.
"""
from __future__ import annotations
import json,sys,time,subprocess
from fractions import Fraction as F
from pathlib import Path
import sympy as s
import mpmath as mp
from spectral_verify import Problem,OUT
z,p,eps,a=s.symbols('z p eps a',positive=True)
checks={}
for r in [1,2]:
    f=(1+z**(2*r))**(p/(2*r))
    K=p*z**(2*r-2)*(1+z**(2*r))**(p/(2*r)-2)*((2*r-1)+(p-1)*z**(2*r))
    checks[f'kernel_r{r}']=s.simplify(s.diff(f,z,2)-K)==0
th=s.symbols('theta',real=True);q=-s.sin(th)
Lp=p*p/4*s.cos(th)**2*s.diff(q,th,3)-p*(p+6)/4*s.sin(th)*s.cos(th)*s.diff(q,th,2)+((p*p-8)*s.cos(th)**2+8)/4*s.diff(q,th)-(p+2)*(p+4)/4*s.sin(th)*s.cos(th)*q
checks['translation_all_p']=s.trigsimp(Lp)==0
mp.mp.dps=75
N=200
g=[mp.mpf(0)]*(N+1);g[0]=mp.mpf(1)
for m in range(1,N//2+1):g[2*m]=(-1)**m*2**(2*m)*mp.bernoulli(2*m)/mp.factorial(2*m)
h=[(1-k)*g[k] for k in range(N+1)]
def point(pp,oo):
    pp,oo=mp.mpf(str(pp)),mp.mpf(str(oo));A2=1+6/pp;Ah=8/pp**2;At=1-Ah;C=1+6/pp+8/pp**2;D=4/pp**2
    aa=[mp.mpc(1)]+[mp.mpc(0)]*N
    for n in range(1,N+1):
        val=mp.mpc(0)
        for k in range(2,n+1,2):
            m=n-k;val+=(A2*g[k]*m*(m-1)+Ah*h[k]*m+C*g[k-2])*aa[m]
        if n>=2:val+=At*(n-2)*aa[n-2]
        for k in range(1,n+1,2):val+=D*1j*oo*h[k-1]*aa[n-k]
        aa[n]=-val/(n*((n-1)**2+6*(n-1)/pp+8/pp**2))
    t=mp.pi/2;q=mp.fsum(aa[n]*t**n for n in range(N+1));qp=mp.fsum(n*aa[n]*t**(n-1) for n in range(1,N+1))
    return q,qp,mp.re(qp*mp.conj(q))
rows=[]
for pp,oo in [('1','2.025'),('1.415','2.05'),('1.0025','28'),('1.705','47.5'),('1.995','39.9')]:
    t=time.time();q,qp,f=point(pp,oo);pr=Problem(F(pp),F(pp),1200);code,v,st=pr.solve(F(oo),F(oo),1200)
    passed=code==0 and mp.mpf(v[-1].lo)<=f<=mp.mpf(v[-1].hi)
    rows.append(dict(p=pp,omega=oo,F1_mpmath=str(f),interval_F1=v[-1].tolist(),inside_interval=passed,seconds=time.time()-t))
    print(rows[-1],flush=True)
assert all(checks.values()) and all(r['inside_interval'] for r in rows)
(OUT/'symbolic_and_point_checks.json').write_text(json.dumps({'symbolic_checks':checks,'point_crosschecks':rows,'point_checks_are_not_coverage_proof':True},indent=2))
