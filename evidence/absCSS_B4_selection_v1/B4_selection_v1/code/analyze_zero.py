import json, numpy as np, math, os
from pathlib import Path
L=8.0

def loadj(f):
 d=json.load(open(f)); return np.array(d['x']),np.array(d['u']),d

def loadnpz(f):
 z=np.load(f); return z['x'],z['u'][-1],{'T':float(z['t'][-1])}

def interp_periodic(xs,us,xd,L=8.0):
 xs=np.asarray(xs); us=np.asarray(us); xd=np.asarray(xd)
 order=np.argsort(xs); xs=xs[order]; us=us[order]
 xx=np.r_[xs,xs[0]+L]; uu=np.r_[us,us[0]]
 xw=((xd-xs[0])%L)+xs[0]
 return np.interp(xw,xx,uu)

def rel(a,b,dx): return math.sqrt(dx*np.sum((a-b)**2))/math.sqrt(dx*np.sum(b*b))
# common grid M1 N3200 as reference spatially
x1,u1,d1=loadj('/mnt/data/absCSS_B4_work/z_m1_3200.json'); dx=x1[1]-x1[0]
rows=[]
for eps,label1,label2 in [(.1,'/mnt/data/absCSS_B4_work/p15_e01_r1.npz','/mnt/data/absCSS_B4_work/p15_e01_r2.npz'),(.06,'/mnt/data/absCSS_B4_work/B3src/B3_transverse_v1/raw/L8_eps_p15/e006_r1.npz','/mnt/data/absCSS_B4_work/B3src/B3_transverse_v1/raw/L8_eps_p15/e006_r2.npz'),(.03,'/mnt/data/absCSS_B4_work/B3src/B3_transverse_v1/raw/L8_eps_p15/e003_r1.npz','/mnt/data/absCSS_B4_work/B3src/B3_transverse_v1/raw/L8_eps_p15/e003_r2.npz')]:
 x3,u3,_=loadnpz(label1); x4,u4,_=loadnpz(label2)
 a=interp_periodic(x3,u3,x1); b=interp_periodic(x4,u4,x1)
 rows.append((eps,rel(a,b,dx),rel(a,u1,dx),rel(b,u1,dx)))
print('eps,D34,D3-M1,D4-M1')
for r in rows: print(*r)
for col,name in [(1,'D34'),(2,'D3M1'),(3,'D4M1')]:
 e=np.array([r[0] for r in rows]); y=np.array([r[col] for r in rows]); q=np.polyfit(np.log(e),np.log(y),1)[0]; print(name,'slope',q)
# M5 vs M1 on N1600 common grid
x0,u0,_=loadj('/mnt/data/absCSS_B4_work/z_m5_N1600_eta0.json')
print('M1 N1600 vs M1 N3200',rel(interp_periodic(x0,u0,x1),u1,dx))
print('eta D(M5,M1_N1600)')
for eta in ['0.0001','3e-05','1e-05','3e-06','1e-06']:
 x,u,d=loadj('/mnt/data/absCSS_B4_work/z_m5_N1600_eta'+eta+'.json')
 # compare to eta0 same N grid
 xz,uz,_=loadj('/mnt/data/absCSS_B4_work/z_m5_N1600_eta0.json')
 dxz=xz[1]-xz[0]
 print(eta,rel(u,uz,dxz))
# spatial fixed eta 1e-5 vs N1600
xref,uref,_=loadj('/mnt/data/absCSS_B4_work/z_m5_eta1e-5_N1600.json'); dxr=xref[1]-xref[0]
print('spatial eta=1e-5 vs N1600')
for N in [400,800]:
 x,u,d=loadj(f'/mnt/data/absCSS_B4_work/z_m5_eta1e-5_N{N}.json')
 print(N, rel(interp_periodic(x,u,xref),uref,dxr))
# time conv N800 eta=1e-5
xa,ua,_=loadj('/mnt/data/absCSS_B4_work/z_m5_eta1e-5_N800.json'); xb,ub,_=loadj('/mnt/data/absCSS_B4_work/z_m5_eta1e-5_N800_dt005.json'); print('time eta1e-5 N800',rel(ub,ua,xa[1]-xa[0]))
