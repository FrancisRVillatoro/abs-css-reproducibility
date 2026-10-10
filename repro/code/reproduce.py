#!/usr/bin/env python3
"""Offline entry point for the public code-and-data reproducibility release.
Default: rebuild the directed interval library, verify every spectral box,
postprocess stored PDE fields and regenerate the scientific plots into outputs/.
No manuscript, paper, LaTeX supplement or publication files are generated.
--rerun-pde also repeats the additional B5 spatial runs (longer computation).
"""
from __future__ import annotations
import argparse,ctypes,json,os,platform,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];CODE=ROOT/'repro/code'

def run(args:list[str],log:Path,env=None)->None:
    print('RUN:', ' '.join(args),flush=True)
    log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w') as f:
        result=subprocess.run(args,cwd=CODE,stdout=f,stderr=subprocess.STDOUT,env=env)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}); see {log}')

def main()->None:
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--workers',type=int,default=3)
    a.add_argument('--skip-spectral',action='store_true',help='verify stored coverage instead of re-evaluating the recurrence')
    a.add_argument('--rerun-pde',action='store_true',help='recompute new B5 fields; replaces generated outputs but never original evidence')
    opt=a.parse_args()
    if platform.system()!='Linux' or platform.machine() not in ('x86_64','AMD64') or ctypes.sizeof(ctypes.c_long)!=8:
        raise RuntimeError('This audited build requires Linux x86-64 LP64. Port and audit the arithmetic backend before using another ABI.')
    if opt.workers<1:raise ValueError('workers must be positive')
    logs=ROOT/'repro/results/verification_logs';start=time.time()
    flags=['gcc','-std=c11','-O3','-fPIC','-shared','-frounding-math','-ffp-contract=off','-fno-fast-math','interval_core.c','-lm','-o','interval_core.so']
    run(flags,logs/'compile.log')
    env=dict(os.environ)
    env.pop('ABSCSS_INTERVAL_LIBRARY',None)
    for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):env[k]='1'
    run([sys.executable,'verify_original_symbolic.py'],logs/'symbolic.log',env)
    cmd=[sys.executable,'run_spectral_full.py','--workers',str(opt.workers)]
    if opt.skip_spectral:cmd.append('--coverage-only')
    run(cmd,logs/'spectral.log',env)
    run([sys.executable,'check_symbolic_and_points.py'],logs/'point_checks.log',env)
    if opt.rerun_pde:
        env['ABSCSS_FORCE_PDE']='1'
        run([sys.executable,'run_b5_extension.py','--workers',str(opt.workers)],logs/'b5_extension.log',env)
        run([sys.executable,'run_b5_rh_scaled.py'],logs/'b5_scaled.log',env)
        env.pop('ABSCSS_FORCE_PDE')
    run([sys.executable,'postprocess_revision.py'],logs/'postprocess.log',env)
    run([sys.executable,'generate_figures.py'],logs/'figures.log',env)
    cert=json.loads((ROOT/'repro/results/spectral/certificate_summary.json').read_text())
    data=json.loads((ROOT/'repro/results/analysis/analysis_summary.json').read_text())
    if not cert['all_pass'] or not all(data['reproduction_checks'].values()):raise RuntimeError('A validation check did not pass')
    result={'all_pass':True,'spectral_recomputed':not opt.skip_spectral,'additional_pde_runs_recomputed':opt.rerun_pde,'elapsed_seconds':time.time()-start,'compiler_flags':flags,'python':sys.version,'platform':platform.platform()}
    (logs/'verification_summary.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
