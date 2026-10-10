#!/usr/bin/env python3
import argparse,json,math,time
from pathlib import Path
import numpy as np
from B1_pade_sparse_newton import PadeSparseOps,M1Sparse,compact_u,invariants_u
from B1_radiation_diagnostic import newton_stage_interleaved,exterior_metrics
from dual_edge_phase_study import geometry

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--m',type=int,required=True);ap.add_argument('--p',type=float,required=True);ap.add_argument('--psiL',type=float,required=True);ap.add_argument('--psiR',type=float,required=True);ap.add_argument('--method',choices=['M1','M5'],required=True);ap.add_argument('--chunk',type=int,default=25);ap.add_argument('--state',required=True);ap.add_argument('--eta-coeff',type=float,default=.05);a=ap.parse_args()
 g=geometry(a.m,a.psiL,a.psiR);N=g['N'];L=g['L'];dx=g['dx'];alpha=a.p*a.p/((a.p+1)*(a.p+2));eps=0.0 if a.method=='M1' else a.eta_coeff*dx**3
 x=-L/2+dx*np.arange(N);ops=PadeSparseOps(N,dx);model=M1Sparse(ops,a.p,alpha,1.0,eps);total=max(1,int(np.ceil(1.0/(.1*dx))));dt=1.0/total
 sp=Path(a.state)
 if sp.exists():
  z=np.load(sp,allow_pickle=False);U=z['U'];step=int(z['step'])
 else:
  U=compact_u(x,a.p,alpha,x0=g['x0']);step=0
 U0=compact_u(x,a.p,alpha,x0=g['x0']); t0=time.time(); nits=[]; target=min(total,step+a.chunk); status='OK'
 for n in range(step,target):
  U,ok,hist=newton_stage_interleaved(U,dt,model,1e-10,12);nits.append(len(hist))
  if not ok or not np.all(np.isfinite(U)) or np.max(np.abs(U))>1e6:
   status='FAIL'; target=n+1; break
 step=target
 np.savez_compressed(sp,U=U,step=np.array(step,dtype=np.int64))
 met,_=exterior_metrics(U,x,a.p,alpha,g['x0'],L,dx); z0=invariants_u(U0,ops,dx,a.p,alpha);zz=invariants_u(U,ops,dx,a.p,alpha)
 r={'status':status,'complete':bool(step>=total and status=='OK'),'step':step,'total_steps':total,'t':step*dt,'dt':dt,'N':N,'m':a.m,'p':a.p,'method':a.method,'L':L,'dx':dx,'psi_L':g['psi_L_actual'],'psi_R':g['psi_R_actual'],'x0':g['x0'],'runtime_chunk_s':time.time()-t0,'mean_newton_chunk':float(np.mean(nits)) if nits else 0.,'max_newton_chunk':int(max(nits)) if nits else 0,'drift_P':float((zz[1]-z0[1])/abs(z0[1])),'drift_H':float((zz[2]-z0[2])/max(abs(z0[2]),1e-300)),**met}
 print(json.dumps(r,indent=2))
 Path(str(sp)+'.json').write_text(json.dumps(r,indent=2)+'\n')
if __name__=='__main__': main()
