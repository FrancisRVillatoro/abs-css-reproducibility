#!/usr/bin/env python3
"""Regenerate scientific visualizations using source CSVs and verified raw fields.
Each panel is a separate matplotlib figure. No synthetic trajectory is plotted.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
from scipy.signal import resample
from postprocess_revision import ROOT,B4,B5,B6,OUT,NEW,load,b5sets
FIG=ROOT/'outputs/figures';FIG.mkdir(parents=True,exist_ok=True)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['ps.fonttype']=42
matplotlib.rcParams['font.size']=8.5
matplotlib.rcParams['axes.labelsize']=9
matplotlib.rcParams['legend.fontsize']=7.1
matplotlib.rcParams['xtick.labelsize']=8
matplotlib.rcParams['ytick.labelsize']=8

def finish(fig,name):
    fig.tight_layout(pad=.65)
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.035)
    fig.savefig(FIG/(name+'.png'),dpi=170,bbox_inches='tight',pad_inches=.035)
    plt.close(fig)

fig,ax=plt.subplots(figsize=(3.25,2.65))
for filename,label,mark in [('b4_compacton_hamiltonian_eps.csv','Compacton, $T=0.2$','o'),('b4_dipole_hamiltonian_eps.csv','Dipole, $T=0.2$','s'),('b4_zero_p15_hamiltonian_eps.csv','Transverse zero, $T=0.01$','^')]:
    df=pd.read_csv(OUT/filename).sort_values('eps')
    ax.loglog(df.eps,df.D34,marker=mark,markersize=4,label=label)
ax.xaxis.set_major_locator(FixedLocator([.03,.06,.1]));ax.xaxis.set_major_formatter(FixedFormatter(['0.03','0.06','0.10']));ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel(r'Hamiltonian smoothing $\varepsilon$');ax.set_ylabel(r'Mutual difference $D_{34}$')
ax.legend(loc='upper left',frameon=False);ax.set_title(r'(a) Two Hamiltonian paths, $p=3/2$',fontsize=9)
finish(fig,'fig1a_hamiltonian_agreement')

fig,ax=plt.subplots(figsize=(3.25,2.65));df=pd.read_csv(OUT/'b4_zero_p15_hamiltonian_eps.csv').sort_values('eps')
for col,label,mark in [('D34',r'$D_{34}$','o'),('D3_hyperref',r'$E_{3;5_*}$','s'),('D4_hyperref',r'$E_{4;5_*}$','^')]:
    ax.loglog(df.eps,df[col],marker=mark,markersize=4,label=label)
ax.xaxis.set_major_locator(FixedLocator([.03,.06,.1]));ax.xaxis.set_major_formatter(FixedFormatter(['0.03','0.06','0.10']));ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel(r'Hamiltonian smoothing $\varepsilon$');ax.set_ylabel('Relative field difference')
ax.legend(loc='upper left',frameon=False);ax.set_title('(b) Comparison with a finite-viscosity field',fontsize=9)
finish(fig,'fig1b_hyperviscous_comparison')

# The two profile panels use exactly the same finest CSS/RH pair.
cf=NEW/'CSS_eta0.0003_m450.json';rf=NEW/'RH_scaled_eta0.0003_m450_tol1e-11_dt0.05.json'
c,r=load(cf),load(rf);n=3600;u=resample(c['u'],n);v=resample(r['u'],n)
A=np.fft.fft(u);B=np.fft.fft(v);cc=np.fft.ifft(A*np.conj(B)).real;j=int(np.argmax(cc));den=cc[(j-1)%n]-2*cc[j]+cc[(j+1)%n]
delta=0. if den==0 else float(np.clip(.5*(cc[(j-1)%n]-cc[(j+1)%n])/den,-.5,.5))
lag=j+delta
if lag>n/2:lag-=n
vs=np.fft.ifft(B*np.exp(-1j*2*np.pi*np.fft.fftfreq(n)*lag)).real
x=-c['L']/2+c['L']*np.arange(n)/n
np.savetxt(OUT/'b5_profile_and_difference.csv',np.c_[x,u,v,vs,u-vs],delimiter=',',header='x,CSS,RH_unaligned,RH_aligned,CSS_minus_RH_aligned',comments='')
fig,ax=plt.subplots(figsize=(3.25,2.65));ax.plot(x,u,label='CSS');ax.plot(x,vs,'--',label='RH (phase aligned)')
ax.set_xlim(-6,8);ax.set_xlabel('Moving-frame coordinate');ax.set_ylabel('Field');ax.legend(frameon=False)
ax.set_title(r'(a) Matched profiles after contact, $T=1.2$',fontsize=9)
finish(fig,'fig2a_css_rh_profiles')
fig,ax=plt.subplots(figsize=(3.25,2.65));ax.plot(x,u-vs)
ax.set_xlim(-6,8);ax.set_xlabel('Moving-frame coordinate');ax.set_ylabel('CSS minus aligned RH')
ax.set_title(r'(b) Field difference, $\eta=3\times10^{-4}$',fontsize=9)
finish(fig,'fig2b_css_rh_difference')

fig,ax=plt.subplots(figsize=(5.55,3.0));df=pd.read_csv(OUT/'b5_refinement_full.csv')
for (eta,g),mark in zip(df.groupby('eta',sort=False),['o','s','^']):
    g=g.sort_values('h');ax.plot(g.h,g.D_CSS_RH,marker=mark,markersize=4,label=rf'$\eta={eta:g}$')
ax.set_xlabel(r'Spatial step $h$');ax.set_ylabel(r'Aligned inter-model difference $D_{\rm CSS,RH}$')
ax.set_ylim(0,.0076);ax.legend(frameon=False,ncol=3,loc='lower center');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
finish(fig,'fig3_b5_spatial_refinement')

# Supplement: actual archived B6 extrema, not fabricated snapshots.
d=load(B6/'data/B6_p15_m85_C15h2_T6_dt1.json')
t=np.array([s['t'] for s in d['snapshots']]);mx=np.array([s['diag']['max'] for s in d['snapshots']]);mn=np.array([s['diag']['min'] for s in d['snapshots']])
fig,ax=plt.subplots(figsize=(3.25,2.65));ax.plot(t,mx,'o-',markersize=3,label='Maximum');ax.plot(t,mn,'s-',markersize=3,label='Minimum')
ax.set_xlabel('Time');ax.set_ylabel('Field extremum');ax.legend(frameon=False,loc='center right');ax.set_title(r'(a) Survival of both signs, $C=1.5$',fontsize=9)
finish(fig,'figS1a_b6_extrema')
fig,ax=plt.subplots(figsize=(3.25,2.65))
for tag,C,ls in [('10',1.,'-'),('12',1.2,'--'),('15',1.5,':')]:
    d=load(B6/f'data/B6_p15_m75_C{tag}h2_T6_dt1.json');ax.plot(d['x'],d['u'],linestyle=ls,label=rf'$C={C:g}$')
ax.set_xlabel('Moving-frame coordinate');ax.set_ylabel('Field at $T=6$');ax.legend(frameon=False);ax.set_title(r'(b) Dependence on $\eta=C h^2$',fontsize=9)
finish(fig,'figS1b_b6_viscosity')
print('Generated',len(list(FIG.glob('*.pdf'))),'vector plot panels from actual data.')
