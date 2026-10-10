#!/usr/bin/env python3
"""Recompute every finite spectral enclosure and its exact coverage mesh."""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import csv,json,time,math,argparse
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
from spectral_verify import *

def qround(s):return F(s)
def direct_group(item):
    ps,pe,rows=item;plo,phi=F(ps),F(pe)
    outdir=OUT/'direct';outdir.mkdir(exist_ok=True)
    name=f'p_{float(plo):.4f}_{float(phi):.4f}.csv';fp=outdir/name
    start=time.time();pr=Problem(plo,phi,N=1200);out=[]
    for olo,ohi in rows:
        olo,ohi=F(olo),F(ohi);N=400
        code,vals,stats=pr.solve(olo,ohi,N)
        if code==1:N=min(1200,int(math.ceil(stats[0]))+100);code,vals,stats=pr.solve(olo,ohi,N)
        while (code!=0 or vals[-1].lo<=0.) and N<1200:
            N=min(1200,N+150);code,vals,stats=pr.solve(olo,ohi,N)
        # Refinement if dependency widths, rather than Taylor remainder, prevent a sign.
        stack=[(olo,ohi,code,vals,stats,N)];certified=[]
        while stack:
            a,b,code,vals,stats,N=stack.pop()
            if code==0 and vals[-1].lo>0:
                certified.append({'p_lo':str(plo),'p_hi':str(phi),'omega_lo':str(a),'omega_hi':str(b),'F1_lower':vals[-1].lo,'F1_upper':vals[-1].hi,'tail_upper':stats[1],'M_upper':stats[2],'S_upper':stats[0],'N':N})
            else:
                if b-a<F(1,100000):raise RuntimeError(('unresolved direct',ps,pe,str(a),str(b),code,vals[-1],stats))
                mid=(a+b)/2
                for aa,bb in [(a,mid),(mid,b)]:
                    c,v,s=pr.solve(aa,bb,N);stack.append((aa,bb,c,v,s,N))
        out.extend(sorted(certified,key=lambda r:F(r['omega_lo'])))
    with fp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
    return name,len(out),min(float(x['F1_lower']) for x in out),time.time()-start

def high_one(spec):
    plo,phi,om=spec;rec,rows=high(F(plo),F(phi),om)
    if not rec['pass_certificate']:raise RuntimeError(('unresolved high',rec))
    out=OUT/'high';out.mkdir(exist_ok=True)
    name=f'p_{float(F(plo)):.4f}_{float(F(phi)):.4f}'
    (out/(name+'.json')).write_text(json.dumps(rec,indent=2))
    with (out/(name+'_remainder.csv')).open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['w_lo','w_hi','b_upper']);w.writerows((str(F(i,2048)),str(F(i+1,2048)),b) for i,b in enumerate(rows))
    return name,rec['matching_ratio_upper'],rec['cone_rhs_upper']

def low_all():
    low=[]
    with (INPUT/'low_mesh.csv').open() as f:
        for row in csv.DictReader(f):
            lo,hi=qround(row['p_lo']),qround(row['p_hi']);pr=Problem(lo,hi,N=520)
            code,v,s=pr.low(N=520)
            if code or not v[0].lo>0:raise RuntimeError(('unresolved low',lo,hi,code,v[0],s))
            low.append(dict(p_lo=str(lo),p_hi=str(hi),z_lo='0',z_hi='4',G_lower=v[0].lo,G_upper=v[0].hi,tail_upper=s[1],M_upper=s[2],S_upper=s[0],N=520))
    with (OUT/'low.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(low[0]));w.writeheader();w.writerows(low)
    return low

def coverage():
    direct=[]
    for fp in sorted((OUT/'direct').glob('*.csv')):
        with fp.open() as f:direct.extend(list(csv.DictReader(f)))
    dgroups={}
    for r in direct:dgroups.setdefault((F(r['p_lo']),F(r['p_hi'])),[]).append(r)
    highrows=[json.loads(p.read_text()) for p in (OUT/'high').glob('*.json')]
    hmap={(F(x['p_lo']),F(x['p_hi'])):x for x in highrows}
    assert set(hmap)==set(dgroups)
    bounds=sorted(dgroups);assert bounds[0][0]==1 and bounds[-1][1]==2
    assert all(a[1]==b[0] for a,b in zip(bounds,bounds[1:]))
    for k,rs in dgroups.items():
        rs.sort(key=lambda x:F(x['omega_lo']))
        assert F(rs[0]['omega_lo'])==2
        assert all(F(a['omega_hi'])==F(b['omega_lo']) for a,b in zip(rs,rs[1:]))
        assert F(rs[-1]['omega_hi'])>=hmap[k]['omega_min']
        assert all(float(r['F1_lower'])>0 and int(r['N'])>float(r['S_upper']) for r in rs)
    with (OUT/'low.csv').open() as f:low=list(csv.DictReader(f))
    low.sort(key=lambda r:F(r['p_lo']));assert F(low[0]['p_lo'])==1 and F(low[-1]['p_hi'])==2
    assert all(F(a['p_hi'])==F(b['p_lo']) for a,b in zip(low,low[1:]))
    assert all(float(r['G_lower'])>0 and int(r['N'])>float(r['S_upper']) for r in low)
    summary=dict(all_pass=True,backend='IEC 60559 binary64 with FE_DOWNWARD, no contraction/fast-math; MPFR 4.2.2 directed transcendental endpoints at 192 bits',direct_rectangles=len(direct),direct_p_boxes=len(bounds),low_boxes=len(low),high_boxes=len(highrows),F1_min=min(float(r['F1_lower']) for r in direct),direct_tail_max=max(float(r['tail_upper']) for r in direct),G_min=min(float(r['G_lower']) for r in low),low_tail_max=max(float(r['tail_upper']) for r in low),high_Rmatch_max=max(r['matching_ratio_upper'] for r in highrows),high_b_max=max(r['b_upper'] for r in highrows),high_R='49/100',high_cone_rhs_max=max(r['cone_rhs_upper'] for r in highrows),high_logder_min=min(r['logder_lower'] for r in highrows),high_tail_max=max(r['tail_upper'] for r in highrows),high_S_max=max(r['S_upper'] for r in highrows),exact_rational_coverage='[1,2] in p; [0,4] in z at low frequency, adjoining finite frequency cells beginning at omega=2, overlap with high frequency for every p-box',legacy_outputs_used_as_proof=False)
    (OUT/'certificate_summary.json').write_text(json.dumps(summary,indent=2))
    with (OUT/'direct_all.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(direct[0]));w.writeheader();w.writerows(direct)
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=3);ap.add_argument('--coverage-only',action='store_true');a=ap.parse_args()
    if a.coverage_only:print(json.dumps(coverage(),indent=2));raise SystemExit(0)
    start=time.time();print('primitives',prim_test(),flush=True)
    low=low_all();print('low PASS',len(low),min(r['G_lower'] for r in low),flush=True)
    groups={}
    with (INPUT/'direct_mesh.csv').open() as f:
        for r in csv.DictReader(f):
            lo,hi=qround(r['p_lo']),qround(r['p_hi']);groups.setdefault((str(lo),str(hi)),[]).append((str(qround(r['omega_lo'])),str(qround(r['omega_hi']))))
    # The existing direct bands already extend beyond all old WKB thresholds.
    hs=[]
    for (ps,pe),rs in groups.items():
        plo=F(ps)
        if plo<F(5,4):om=28
        elif plo<F(69,50):om=40
        elif plo<F(37,25):om=42
        elif plo<F(79,50):om=44
        elif plo<F(42,25):om=46
        elif plo<F(7,4):om=48
        else:om=40
        hs.append((ps,pe,om))
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(high_one,s) for s in hs]):print('high',f.result(),flush=True)
        for f in as_completed([pool.submit(direct_group,(ps,pe,rs)) for (ps,pe),rs in groups.items()]):print('direct',f.result(),flush=True)
    summary=coverage();summary['execution_wall_seconds']=time.time()-start
    (OUT/'certificate_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
