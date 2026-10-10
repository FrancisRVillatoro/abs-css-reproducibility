#!/usr/bin/env python3
"""Independent complete spectral certificate, directed IEEE-754 + MPFR endpoints.
All parameter and frequency boxes have exact rational endpoints. The old
certificate is used only to choose a mesh, never as numerical evidence.
"""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import ctypes as ct,json,sys,math,time,csv,argparse,hashlib,platform
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as sp
from mpfr_bounds import I,rat,pi,version,binary
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
OUT=ROOT/'repro/results/spectral';OUT.mkdir(parents=True,exist_ok=True)
INPUT=ROOT/'repro/inputs'
class CI(ct.Structure):
    _fields_=[('lo',ct.c_double),('hi',ct.c_double)]
    def py(self):return I(self.lo,self.hi)
def ci(a):
    a=I.cast(a);return CI(a.lo,a.hi)
lib=ct.CDLL(str(Path(os.environ.get('ABSCSS_INTERVAL_LIBRARY',str(HERE/'interval_core.so')))))
lib.iv_environment.restype=ct.c_int
lib.iv_operations.argtypes=[ct.c_int,CI,CI,ct.POINTER(CI)]
lib.iv_problem.argtypes=[CI,ct.c_int,ct.POINTER(CI),CI];lib.iv_problem.restype=ct.c_void_p
lib.iv_free.argtypes=[ct.c_void_p]
lib.iv_solve.argtypes=[ct.c_void_p,CI,ct.c_int,CI,CI,CI,ct.POINTER(CI),ct.POINTER(ct.c_double)]
lib.iv_low.argtypes=[ct.c_void_p,CI,ct.c_int,CI,CI,CI,ct.POINTER(CI),ct.POINTER(CI),ct.POINTER(CI),ct.POINTER(ct.c_double)]
lib.iv_wkb_bound.argtypes=[CI,ct.c_double,ct.c_double,ct.POINTER(ct.c_double),ct.POINTER(ct.c_double)];lib.iv_wkb_bound.restype=ct.c_double

_CACHE={}
def coeffs(N):
    if N in _CACHE:return _CACHE[N]
    g=(CI*(N+1))();x0=(CI*(N+1))();x1=(CI*(N+1))()
    g[0]=ci(1);x0[0]=ci(1)
    fact=1
    for n in range(1,N+1):
        fact*=n
        if n%2==0:
            m=n//2; bn=sp.bernoulli(n);q=F((-1)**m*2**n*int(sp.numer(bn)),int(sp.denom(bn))*fact)
            g[n]=ci(rat(q));x0[n]=ci(rat(F((-1)**m,fact)))
        else:x1[n]=ci(rat(F(-(-1)**((n-1)//2),2*fact)))
    _CACHE[N]=(g,x0,x1);return _CACHE[N]
class Problem:
    def __init__(self,plo,phi,N=1200,delta=I(1)):
        self.P=I(rat(plo).lo,rat(phi).hi);self.N=N;self.g,self.x0,self.x1=coeffs(N)
        self.ptr=lib.iv_problem(ci(self.P),N,self.g,ci(delta))
        if not self.ptr:raise RuntimeError('unsupported floating-point environment')
    def __del__(self):
        if getattr(self,'ptr',None):lib.iv_free(self.ptr)
    def solve(self,olo,ohi,N,R=F(33,20),tstar=F(11,7),t=None):
        t=pi()/2 if t is None else I.cast(t)
        out=(CI*7)();stats=(ct.c_double*3)();om=I(rat(olo).lo,rat(ohi).hi)
        status=lib.iv_solve(self.ptr,ci(om),N,ci(rat(R)),ci(rat(tstar)),ci(t),out,stats)
        return status,[x.py() for x in out],list(stats)
    def low(self,zmax=4,N=520):
        out=(CI*7)();stats=(ct.c_double*3)()
        status=lib.iv_low(self.ptr,ci(I(0,rat(zmax).hi)),N,ci(rat(F(33,20))),ci(rat(F(11,7))),ci(pi()/2),self.x0,self.x1,out,stats)
        return status,[x.py() for x in out],list(stats)

# Rectangular complex interval arithmetic; only MPFR-directed endpoint primitives.
def cadd(a,b):return (a[0]+b[0],a[1]+b[1])
def cscale(a,b):return (a[0]*b,a[1]*b)
def cmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cabs(a):
    x,y=a[0].abs(),a[1].abs();return (x*x+y*y).sqrt()

def high(plo,phi,omega0,N=180,T=7,Rcone=F(49,100)):
    pp=I(rat(plo).lo,rat(phi).hi);delta=I(0,rat(F(1,omega0**2)).hi)
    pr=Problem(plo,phi,N,delta)
    status,vals,stats=pr.solve(1,1,N,R=F(14),tstar=F(T),t=T)
    if status:raise RuntimeError(('high recurrence',status,stats))
    q=(vals[0],vals[1]);X1=(vals[2],vals[3]);X2=(vals[4],vals[5])
    tm=rat(F(T,omega0))
    # IMPORTANT: the complete parameter range is 0 <= t <= T/omega0,
    # not only its upper endpoint. These enclosures include t=0.
    sinc=I((1-tm**2/6).lo,1.);cost=I((1-tm**2/2).lo,1.)
    f=T*sinc;a=2/pp-rat(F(1,3));kp=(4/(pp**2)).cbrt()
    Y0=cscale(q,f**a)
    Y1=cscale(cadd(cscale(q,a*cost*(f**(a-rat(F(1,3))))),cscale(X1,(f**(a+rat(F(2,3))))/T)),1/kp)
    qa=a*((a-rat(F(1,3)))*cost**2*f**(a-rat(F(2,3)))-delta*f**(a+rat(F(4,3))))
    Y2=cscale(cadd(cadd(cscale(q,qa),cscale(X1,(2*a+rat(F(2,3)))*cost*f**(a+rat(F(1,3)))/T)),cscale(X2,f**(a+rat(F(4,3)))/(T*T))),1/kp**2)
    sq3=I(3).sqrt();roots=[(sq3/2,I(-.5)),(I(0),I(1)),(-sq3/2,I(-.5))];zs=[]
    for mu in roots:
        co=(mu[0],-mu[1]);zs.append(cscale(cadd(cadd(Y0,cmul(co,Y1)),cmul(cmul(co,co),Y2)),rat(F(1,3))))
    zp=cabs(zs[0]);rm=(cabs(zs[1])+cabs(zs[2]))/zp
    # At the matching point, use an actual uniform lower bound for sinc.
    rb=(rat(phi)**2/(4*T*sinc)).cbrt().hi
    hb=(rat(phi)**2/(4*omega0)).cbrt().hi
    sqrt_sup=(ct.c_double*2048)(*[rat(1-F(j,2048)**6).sqrt().hi for j in range(2048)])
    rows=(ct.c_double*2048)();bmax=lib.iv_wkb_bound(ci(pp),rb,hb,sqrt_sup,rows)
    rc=rat(Rcone);rhs=-sq3*rc/2+I(0,bmax)/3*(1+rc)*(2+rc)
    logder=(sq3/2-sq3*rc)/(1+rc)**2
    good=rm.hi<rc.lo and rc.hi<.5 and rhs.hi<0 and logder.lo>0
    return dict(p_lo=str(plo),p_hi=str(phi),omega_min=omega0,N=N,T=T,S_upper=stats[0],tail_upper=stats[1],M_upper=stats[2],matching_ratio_upper=rm.hi,zplus_lower=zp.lo,b_upper=bmax,cone_R=str(Rcone),cone_rhs_upper=rhs.hi,logder_lower=logder.lo,pass_certificate=good,sinc_enclosure=sinc.tolist(),cos_enclosure=cost.tolist()),list(rows)

def prim_test(n=3000):
    if not lib.iv_environment():raise RuntimeError('IEEE environment unsupported')
    rng=np.random.default_rng(1032026);count=0
    vals=[0.,-0.,1.,-1.,2.**-1022,math.nextafter(0.,1.),1e-250,-1e-250,1e250,-1e250]
    for _ in range(n):vals.append(float(np.ldexp(rng.uniform(-1,1),int(rng.integers(-1000,1000)))))
    for j in range(n):
        a=sorted([vals[j],vals[(j*17+3)%len(vals)]]);b=sorted([vals[(j*19+5)%len(vals)],vals[(j*31+7)%len(vals)]])
        for op,name in [(0,'add'),(1,'sub'),(2,'mul')]:
            out=CI();lib.iv_operations(op,CI(*a),CI(*b),ct.byref(out))
            # Endpoints of exact real ranges occur at corner inputs.
            true=[]
            for x in a:
                for y in b:
                    xx,yy=F.from_float(x),F.from_float(y)
                    true.append(xx+yy if op==0 else xx-yy if op==1 else xx*yy)
            if math.isfinite(out.lo) and F.from_float(out.lo)>min(true):raise AssertionError((name,'lo',a,b,out.lo))
            if math.isfinite(out.hi) and F.from_float(out.hi)<max(true):raise AssertionError((name,'hi',a,b,out.hi))
            count+=1
        pos=sorted([abs(b[0])+2.**-1000,abs(b[1])+2.**-1000])
        if pos[0]>0 and all(math.isfinite(v) for v in pos):
            out=CI();lib.iv_operations(3,CI(*a),CI(*pos),ct.byref(out))
            true=[F.from_float(x)/F.from_float(y) for x in a for y in pos]
            if math.isfinite(out.lo) and F.from_float(out.lo)>min(true):raise AssertionError(('div','lo',a,pos))
            if math.isfinite(out.hi) and F.from_float(out.hi)<max(true):raise AssertionError(('div','hi',a,pos))
            count+=1
    rec={'pass':True,'exact_rational_primitive_tests':count,'mpfr_version':version(),'platform':platform.platform(),'python':sys.version,'precision_for_mpfr_primitives':192,'pi_enclosure':pi().tolist(),'notes':'Tests supplement the primitive enclosure proof; they do not replace it.'}
    (OUT/'primitive_tests.json').write_text(json.dumps(rec,indent=2));return rec

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--test',action='store_true');ap.add_argument('--high-smoke',action='store_true');args=ap.parse_args()
    if args.test:print(prim_test())
    if args.high_smoke:
        for lo,hi,om in [(F(1),F(201,200),28),(F(3,2),F(151,100),44),(F(199,100),F(2),40)]:
            t=time.time();r,b=high(lo,hi,om);print(r,'time',time.time()-t,flush=True)
