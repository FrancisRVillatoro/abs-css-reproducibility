#!/usr/bin/env python3
"""Generate reproducible numerical evidence from archived and new fields.
No PDE time-stepping is performed here. Failed runs are never plotted.
"""
from __future__ import annotations
import json,csv,sys,math
from pathlib import Path
import numpy as np
from scipy.signal import resample
ROOT=Path(__file__).resolve().parents[2]
B4=ROOT/'evidence/absCSS_B4_selection_v1/B4_selection_v1'
B5=ROOT/'evidence/absCSS_B5_css_rh_v1/B5_css_rh_v1'
B6=ROOT/'evidence/absCSS_B6_css_chase_v1/B6_css_chase_v1'
NEW=ROOT/'repro/results/b5';OUT=ROOT/'repro/results/analysis';OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(B5/'code'))
from compare_fields import align_diff

def load(path):
    d=json.loads(Path(path).read_text())
    if d.get('status') not in ('OK',):raise ValueError((path,d.get('status')))
    return d

def csvout(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def diff(a,b,n=3600):
    # normalization by first field, exactly as the B5 production comparator
    aa=resample(np.asarray(a['u']),n);bb=resample(np.asarray(b['u']),n)
    v,lag=align_diff(aa,bb);return v,lag

def periodic_interp(x,u,xd,L):
    x=np.asarray(x);u=np.asarray(u);order=np.argsort(x);x=x[order];u=u[order]
    return np.interp((np.asarray(xd)-x[0])%L+x[0],np.r_[x,x[0]+L],np.r_[u,u[0]])

def b5sets():
    raw=B5/'raw/chase'
    mapping={.003:[(95,raw/'e0.003_m95.json',raw/'RH_m95_T1.2_eta0.003.json'),(125,raw/'e0.003_m125.json',raw/'RH_m125_T1.2_eta0.003.json')],
        .001:[(125,raw/'CSS_dt0.05.json',raw/'RH_dt0.05.json'),(155,raw/'e0.001_m155.json',raw/'RH_m155_T1.2_eta0.001.json')],
        .0003:[(190,raw/'e0.0003_m190.json',raw/'RH_m190_T1.2_eta0.0003.json'),(225,raw/'e0.0003_m225.json',raw/'RH_m225_T1.2_eta0.0003.json')]}
    for eta,ms in [(.003,[190]),(.001,[225]),(.0003,[300,450])]:
        mapping[eta]+= [(m,NEW/f'CSS_eta{eta:g}_m{m}.json',NEW/f'RH_scaled_eta{eta:g}_m{m}_tol1e-11_dt0.05.json') for m in ms]
    return mapping

def main():
    allrows=[];summary=[];refs=[]
    for eta,triples in b5sets().items():
        prev=None;rr=[]
        for m,cf,rf in triples:
            c,r=load(cf),load(rf)
            assert c['p']==r['p']==1.5 and c['T_reached']==r['T_reached']==1.2
            assert math.isclose(c['eta'],eta) and math.isclose(r['eta'],eta)
            assert c['N']==r['N'] and np.allclose(c['x'],r['x'],rtol=0,atol=1e-12)
            D,lag=diff(c,r);rawD=np.linalg.norm(np.array(c['u'])-r['u'])/np.linalg.norm(c['u'])
            dc,dr=(diff(c,prev[0])[0],diff(r,prev[1])[0]) if prev else (math.nan,math.nan)
            row=dict(eta=eta,m=m,N=c['N'],h=c['h'],dt_ratio_CSS=c['dt_ratio'],dt_ratio_RH=r['dt_ratio'],D_CSS_RH=D,D_unaligned=rawD,CSS_dP_rel=(c['inv']['P']-c['inv0']['P'])/c['inv0']['P'],RH_dP_rel=(r['inv']['P']-r['inv0']['P'])/r['inv0']['P'],phase_shift=lag*c['L']/3600,CSS_grid_difference=dc,RH_grid_difference=dr,CSS_file=str(cf.relative_to(ROOT)),RH_file=str(rf.relative_to(ROOT)))
            rr.append(row);allrows.append(row);prev=(c,r)
        last=rr[-1];pen=rr[-2]
        summary.append(dict(eta=eta,m_values=';'.join(str(t[0]) for t in triples),finest_m=last['m'],D_finest=last['D_CSS_RH'],CSS_final_difference=last['CSS_grid_difference'],RH_final_difference=last['RH_grid_difference'],CSS_ratio_to_previous=last['CSS_grid_difference']/pen['CSS_grid_difference'],RH_ratio_to_previous=last['RH_grid_difference']/pen['RH_grid_difference'],D_over_sum_final_differences=last['D_CSS_RH']/(last['CSS_grid_difference']+last['RH_grid_difference'])))
        refs.extend([last['CSS_file'],last['RH_file']])
    csvout('b5_refinement_full.csv',allrows);csvout('b5_refinement_summary.csv',summary)
    # Normalized discrepancy against finest raw M2 field, without phase alignment.
    m2=[]
    for low,high in [(800,1200),(1200,1600),(800,1600)]:
        a=load(B4/f'raw_selected/z_M2scaled_eta1e-5_N{low}.json');b=load(B4/f'raw_selected/z_M2scaled_eta1e-5_N{high}.json')
        bi=periodic_interp(b['x'],b['u'],a['x'],8);v=float(np.linalg.norm(np.asarray(a['u'])-bi)/np.linalg.norm(bi))
        m2.append(dict(N_coarse=low,N_fine=high,relative_L2=v,coarse_dt_ratio=a.get('dt_ratio',a.get('dt',0)/(8/low)),fine_dt_ratio=b.get('dt_ratio',b.get('dt',0)/(8/high)),denominator='fine field interpolated onto coarse grid',alignment='none',interpolation='periodic piecewise linear, fine-to-coarse'))
    csvout('b4_m2_spatial_differences.csv',m2)
    # Scaled-Newton controls: same PDE, stages and meshes, only algebraic norm changes.
    ctrl=[]
    for label,f1,f2 in [
        ('RH original vs time-scaled stages, m125 eta.001',B5/'raw/chase/RH_dt0.05.json',NEW/'RH_scaled_eta0.001_m125_tol1e-11_dt0.05.json'),
        ('RH scaled tol 1e-10 vs 1e-11, m225 eta.001',NEW/'RH_scaled_eta0.001_m225_tol1e-10_dt0.05.json',NEW/'RH_scaled_eta0.001_m225_tol1e-11_dt0.05.json')]:
        a,b=load(f1),load(f2);d=float(np.linalg.norm(np.array(a['u'])-b['u'])/np.linalg.norm(b['u']))
        ctrl.append(dict(control=label,relative_L2_unaligned=d,phase_aligned_L2=diff(b,a)[0]))
    csvout('b5_solver_controls.csv',ctrl)
    # Native RH invariant drift, signed relative H and P; absolute changes in Casimirs.
    rh=[]
    for m in (125,155,190,250):
        d=load(B5/f'raw/chase/RH_m{m}_T1.2_tol1e-9.json');a=d['inv0'];b=d['inv']
        row={'m':m,'P_relative_change':(b['P']-a['P'])/a['P'],'H_relative_change':abs(b['H']-a['H'])/abs(a['H'])}
        row.update({f'{k}_absolute_change':abs(b[k]-a[k]) for k in ('C0','Cc','Cs')});rh.append(row)
    csvout('b5_native_invariants.csv',rh)
    # Preserve source tables verbatim under unique names used by the figures.
    for name in ('compacton_hamiltonian_eps.csv','compacton_hypervisc_eta.csv','dipole_hamiltonian_eps.csv','dipole_hypervisc_eta.csv','zero_p15_hamiltonian_eps.csv','M2_p15_time_convergence.csv','M2_p15_fixed_eta1e-5_space.csv'):
        (OUT/('b4_'+name)).write_bytes((B4/'data'/name).read_bytes())
    log={'B5':summary,'M2':m2,'solver_controls':ctrl,'reproduction_checks':{'M2_last_matches_archive':abs(m2[1]['relative_L2']-4.205695024511644e-5)<1e-12},'source_files_for_B5_profiles':refs,'note':'Grid differences and ratios are empirical resolution diagnostics, not rigorous discretisation error bounds.'}
    (OUT/'analysis_summary.json').write_text(json.dumps(log,indent=2));print(json.dumps(log,indent=2))
if __name__=='__main__':main()
