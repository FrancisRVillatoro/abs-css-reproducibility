"""Directed MPFR endpoints, via the public C ABI (no Python extension needed).
Each elementary result is rounded at 192 bits before directed conversion to
binary64. Fractions are converted using mpfr_set_str/div with correct rounding.
"""
from __future__ import annotations
import ctypes as C, ctypes.util, math
from fractions import Fraction
_lib=C.CDLL(ctypes.util.find_library('mpfr'))
# The supported platform is the LP64 MPFR ABI; checked by the runner.
class _S(C.Structure):
    _fields_=[('_mpfr_prec',C.c_long),('_mpfr_sign',C.c_int),('_mpfr_exp',C.c_long),('_mpfr_d',C.POINTER(C.c_ulong))]
P=C.POINTER(_S)
for name,args,res in [('mpfr_init2',[P,C.c_long],None),('mpfr_clear',[P],None),('mpfr_set_str',[P,C.c_char_p,C.c_int,C.c_int],C.c_int),('mpfr_set_d',[P,C.c_double,C.c_int],C.c_int),('mpfr_get_d',[P,C.c_int],C.c_double),('mpfr_get_version',[],C.c_char_p),('mpfr_const_pi',[P,C.c_int],C.c_int)]:
    f=getattr(_lib,name);f.argtypes=args;f.restype=res
for name in ['add','sub','mul','div','pow']:
    f=getattr(_lib,'mpfr_'+name);f.argtypes=[P,P,P,C.c_int];f.restype=C.c_int
for name in ['sqrt','cbrt','sin','cos','exp','log']:
    f=getattr(_lib,'mpfr_'+name);f.argtypes=[P,P,C.c_int];f.restype=C.c_int
RNDD=3;RNDU=2;RNDN=0
class M:
    def __init__(self,x=0,prec=192):
        self.s=_S();_lib.mpfr_init2(C.byref(self.s),prec)
        if isinstance(x,float):_lib.mpfr_set_d(C.byref(self.s),x,RNDN)
        else:
            if isinstance(x,int) and x.bit_length()>prec:raise ValueError('use sufficiently large precision for integer')
            _lib.mpfr_set_str(C.byref(self.s),str(x).encode(),10,RNDN)
    def __del__(self):
        if hasattr(self,'s'):_lib.mpfr_clear(C.byref(self.s))

def binary(name,a,b,upper=False):
    ma,mb,mc=M(float(a)),M(float(b)),M();rnd=RNDU if upper else RNDD
    getattr(_lib,'mpfr_'+name)(C.byref(mc.s),C.byref(ma.s),C.byref(mb.s),rnd)
    return _lib.mpfr_get_d(C.byref(mc.s),rnd)
def unary(name,a,upper=False):
    ma,mb=M(float(a)),M();rnd=RNDU if upper else RNDD
    getattr(_lib,'mpfr_'+name)(C.byref(mb.s),C.byref(ma.s),rnd)
    return _lib.mpfr_get_d(C.byref(mb.s),rnd)

def rat(q):
    q=Fraction(q)
    # This exact-rational conversion does not rely on conversion accuracy:
    # compare the candidate binary rational with the exact input.
    v=float(q)
    if not math.isfinite(v):raise OverflowError(q)
    f=Fraction.from_float(v)
    return I(v if f<=q else math.nextafter(v,-math.inf),v if f>=q else math.nextafter(v,math.inf))

class I:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):self.lo=float(lo);self.hi=float(lo if hi is None else hi)
    @staticmethod
    def cast(a):
        if isinstance(a,I):return a
        return rat(a) if isinstance(a,(Fraction,int,str)) else I(a)
    def __add__(self,o):
        o=I.cast(o);return I(binary('add',self.lo,o.lo),binary('add',self.hi,o.hi,True))
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,o):return self+-I.cast(o)
    def __rsub__(self,o):return I.cast(o)+-self
    def __mul__(self,o):
        o=I.cast(o);pairs=[(x,y) for x in (self.lo,self.hi) for y in (o.lo,o.hi)]
        return I(min(binary('mul',x,y) for x,y in pairs),max(binary('mul',x,y,True) for x,y in pairs))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=I.cast(o)
        if o.lo<=0<=o.hi:raise ZeroDivisionError(o)
        pairs=[(x,y) for x in (self.lo,self.hi) for y in (o.lo,o.hi)]
        return I(min(binary('div',x,y) for x,y in pairs),max(binary('div',x,y,True) for x,y in pairs))
    def __rtruediv__(self,o):return I.cast(o)/self
    def __pow__(self,e):
        if isinstance(e,int):
            if e<0:return 1/self**(-e)
            a=self;v=I(1)
            while e:
                if e&1:v=v*a
                e//=2
                if e:a=a*a
            return v
        e=I.cast(e)
        if self.lo<=0:raise ValueError('positive base required')
        pairs=[(x,y) for x in (self.lo,self.hi) for y in (e.lo,e.hi)]
        return I(min(binary('pow',x,y) for x,y in pairs),max(binary('pow',x,y,True) for x,y in pairs))
    def sqrt(self):
        if self.lo<0:raise ValueError(self)
        return I(unary('sqrt',self.lo),unary('sqrt',self.hi,True))
    def cbrt(self):return I(unary('cbrt',self.lo),unary('cbrt',self.hi,True))
    def abs(self):return I(0 if self.lo<=0<=self.hi else min(abs(self.lo),abs(self.hi)),max(abs(self.lo),abs(self.hi)))
    def __repr__(self):return f'I({self.lo!r},{self.hi!r})'
    def tolist(self):return [self.lo,self.hi]

def pi():
    a=M();b=M();_lib.mpfr_const_pi(C.byref(a.s),RNDD);_lib.mpfr_const_pi(C.byref(b.s),RNDU)
    return I(_lib.mpfr_get_d(C.byref(a.s),RNDD),_lib.mpfr_get_d(C.byref(b.s),RNDU))

def version():return _lib.mpfr_get_version().decode()
