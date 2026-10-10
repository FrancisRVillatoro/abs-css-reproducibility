#!/usr/bin/env python3
import argparse, csv, json, math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from B1_radiation_diagnostic import run_case

L_TARGET=40.0
X0_REF=-5.0
PSIS=(0.0,0.25,0.5,0.75)

def frac(x):
    y=x-math.floor(x)
    if abs(y-1.0)<5e-12 or abs(y)<5e-12:
        return 0.0
    return y

def geometry(m, psiL, psiR, L_target=L_TARGET, x0_ref=X0_REF):
    # Exact phase relation: 2*pi/h = m + psiR - psiL.
    den=m + psiR - psiL
    if den <= 0:
        raise ValueError('nonpositive denominator')
    h=2.0*math.pi/den
    N=max(32, int(round(L_target/h)))
    L=N*h
    xmin=-0.5*L
    # Enforce left edge phase while keeping center as close as possible to x0_ref.
    kL=int(round((x0_ref-math.pi-xmin)/h - psiL))
    x0=xmin+math.pi+(kL+psiL)*h
    aL=frac((x0-math.pi-xmin)/h)
    aR=frac((x0+math.pi-xmin)/h)
    return {'m':m,'N':N,'L':L,'dx':h,'x0':x0,'kL':kL,
            'psi_L_actual':aL,'psi_R_actual':aR,
            'L_minus_target':L-L_target,'x0_minus_ref':x0-x0_ref}

def one_case(spec):
    m,p,psiL,psiR,method,T,dt_ratio,eta_coeff=spec
    g=geometry(m,psiL,psiR)
    r=run_case(N=g['N'],p=p,T=T,method=method,L=g['L'],dt_ratio=dt_ratio,
               eta_coeff=eta_coeff,x0=g['x0'],newton_tol=1e-10,newton_max=12)
    r.update(g)
    r['psi_L_target']=psiL; r['psi_R_target']=psiR
    return r

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ms',default='63,126,188,251')
    ap.add_argument('--ps',default='1.75,2.0')
    ap.add_argument('--methods',default='M1')
    ap.add_argument('--psis',default='0,0.25,0.5,0.75')
    ap.add_argument('--T',type=float,default=1.0)
    ap.add_argument('--dt-ratio',type=float,default=0.1)
    ap.add_argument('--eta-coeff',type=float,default=0.05)
    ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--outdir',default='dual_edge_results')
    a=ap.parse_args()
    ms=[int(x) for x in a.ms.split(',') if x]
    ps=[float(x) for x in a.ps.split(',') if x]
    methods=[x.strip() for x in a.methods.split(',') if x.strip()]
    psis=[float(x) for x in a.psis.split(',') if x]
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    raw=out/'raw'; raw.mkdir(exist_ok=True)
    specs=[(m,p,pl,pr,method,a.T,a.dt_ratio,a.eta_coeff)
           for method in methods for p in ps for pl in psis for pr in psis for m in ms]
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs={ex.submit(one_case,s):s for s in specs}
        for i,f in enumerate(as_completed(futs),1):
            s=futs[f]
            try: r=f.result()
            except Exception as e:
                r={'status':'EXCEPTION','m':s[0],'p':s[1],'psi_L_target':s[2],
                   'psi_R_target':s[3],'method':s[4],'error':repr(e)}
            results.append(r)
            tag=f"{r.get('method')}_p{r.get('p')}_L{r.get('psi_L_target')}_R{r.get('psi_R_target')}_m{r.get('m')}"
            (raw/f'{tag}.json').write_text(json.dumps(r,indent=2)+'\n')
            if r.get('status')=='OK':
                print(f"[{i}/{len(specs)}] {tag} N={r['N']} h={r['dx']:.5g} E={r['shape_relL2_rebuilt']:.3e} A={r['fixed050_max_abs']:.3e}",flush=True)
            else:
                print(f"[{i}/{len(specs)}] {tag} {r.get('status')} step={r.get('step','?')}",flush=True)
    results.sort(key=lambda r:(str(r.get('method')),float(r.get('p',0)),float(r.get('psi_L_target',0)),float(r.get('psi_R_target',0)),int(r.get('m',0))))
    (out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    keys=sorted({k for r in results for k in r})
    with (out/'results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(results)

if __name__=='__main__': main()
