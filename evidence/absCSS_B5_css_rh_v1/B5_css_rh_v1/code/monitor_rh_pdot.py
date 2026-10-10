import sys,math,json,numpy as np
sys.path.insert(0,'/mnt/data/B5_work')
from B5_css_rh import *

def run(m):
 p=1.5; T=1.2; dt_ratio=.05; x,h,L,N=compatible_grid(p,m,12);alpha=p*p/((p+1)*(p+2));op=PadeSparseOps(N,h);mid=.5*h
 U=pulse_sum(x,p,alpha,[dict(center=mid-4,c=3,sign=1),dict(center=mid+4,c=1,sign=-1)])
 M=M6RHSparse(op,p,alpha,c0=1,eta=0);steps=int(math.ceil(T/(dt_ratio*h)));dt=T/steps
 ts=[0.]; Ps=[]; rhs=[]
 q=rh_invariants(U,x,h,p,alpha,op,t=0,c0=1);Ps.append(q['P']);rhs.append(q['Pdot_rhs'])
 for n in range(steps):
  U,ok,hist=newton_step(U,dt,M,tol=1e-9,maxit=20,refresh=1)
  if not ok: raise RuntimeError((m,n,hist))
  t=(n+1)*dt; q=rh_invariants(U,x,h,p,alpha,op,t=t,c0=1);ts.append(t);Ps.append(q['P']);rhs.append(q['Pdot_rhs'])
 ts=np.array(ts);Ps=np.array(Ps);rhs=np.array(rhs)
 integ=np.trapz(rhs,ts); actual=Ps[-1]-Ps[0]
 # finite-difference residual L2 over interior time nodes
 fd=(Ps[2:]-Ps[:-2])/(ts[2:]-ts[:-2]); rr=rhs[1:-1]; rel=np.linalg.norm(fd-rr)/max(np.linalg.norm(rr),1e-300)
 return dict(m=m,N=N,h=h,dt=dt,steps=steps,P0=Ps[0],P1=Ps[-1],deltaP=actual,int_rhs=integ,integral_resid=actual-integ,integral_rel=(actual-integ)/max(abs(actual),1e-300),fd_rhs_relL2=rel,max_rhs=float(np.max(np.abs(rhs))))
for m in [125,155,190]:
 r=run(m);print(json.dumps(r));open(f'/mnt/data/B5_work/rh_pdot_m{m}.json','w').write(json.dumps(r,indent=2))
