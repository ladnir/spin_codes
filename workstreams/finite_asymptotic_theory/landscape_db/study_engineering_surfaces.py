"""Sequential, resumable BCH/RM/random Q1 engineering surfaces.

Random means describe a uniform rate-half subspace sampled once and reused.
Conditional 60-bit simultaneous spectrum caps are reported separately.
All numerical outputs are local diagnostics, never certificates.
"""
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import activation_q1_refresh as refresh
import activation_refresh_native as native
import bch_growth_model as onset
import random_spectrum_variance as variance
import read_grid_receipts as receipts
import run_occupation_grid as source

HERE=Path(__file__).resolve().parent
DIRECTORY=HERE/'engineering_surface_v1'
LN2=math.log(2)


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def random_mean_counts(block):
    probability=((1 << (block//2))-1)/((1 << block)-1)
    return {w:math.comb(block,w)*probability for w in range(1,block+1)}


def summarize_coefficients(block,rows,counts,best,witnesses):
    weights=sorted(counts)
    terms=np.array([math.log(rows)+math.log(counts[w])+best[w] for w in weights])
    total=float(np.logaddexp.reduce(terms));dominant=weights[int(np.argmax(terms))]
    return dict(margin_bits=-total/LN2,intercept_bits=-total/LN2+math.log2(rows),
                dominant_weight=dominant,dominant_log_tilt=float(witnesses[dominant]),
                dominant_share=math.exp(float(max(terms))-total))


def evaluate(kernel,block,t,s,exponent,ac,spectra):
    rows=(1<<exponent)//(block//2);cutoff=block*rows//10
    grid=np.unique(np.r_[np.arange(-160,1)/10,np.arange(-120,141)/20-math.log(rows)])

    def values(logs):
        lam=np.exp(logs)
        moments=kernel.coefficients(*refresh.epoch_logs(t,s,ac,lam),rows//t,block)
        return np.minimum(0.,moments+cutoff*lam[:,None])

    initial=values(grid);indices=np.argmin(initial,axis=0)
    best=initial[indices,np.arange(block+1)];witnesses=grid[indices]
    relevant=set()
    for counts in spectra.values():
        terms={w:math.log(n)+best[w] for w,n in counts.items()}
        relevant.update(w for w,v in terms.items() if v>=max(terms.values())-24*LN2)
    extra=np.unique(np.concatenate([witnesses[w]+np.arange(-10,11)/200 for w in relevant]))
    refined=values(extra);indices=np.argmin(refined,axis=0)
    candidate=refined[indices,np.arange(block+1)];improve=candidate<best
    best[improve]=candidate[improve];witnesses[improve]=extra[indices[improve]]
    lo=min(min(grid),min(extra));hi=max(max(grid),max(extra))
    result=[]
    for family,counts in spectra.items():
        row=dict(family=family,block_bits=block,dimension=block//2,step_bits=t,state_bits=s,
                 message_exponent=exponent,outer_rows=rows,epochs_per_region=rows//t,
                 **summarize_coefficients(block,rows,counts,best,witnesses))
        row['witness_at_edge']=int(row['dominant_log_tilt'] in (lo,hi))
        row['setup_charged_margin_bits']=(-float(np.logaddexp(-row['margin_bits']*LN2,-60*LN2))/LN2
                                          if family=='random_caps60' else '')
        result.append(row)
    return result


def main():
    _,observations,_=receipts.snapshot()
    spectra={};maps={};dependencies={}
    for name in source.grid.BASES:
        dependencies.update(source.grid.verify(HERE/name)['source_sha256'])
    for row in observations.values():
        family='bch' if 'BCH' in row['series'] else 'rm' if row['series'].startswith('RM(') else None
        if family is None: continue
        b,t,s=int(row['block_bits']),int(row['step_bits']),int(row['state_bits'])
        if family not in spectra.setdefault(b,{}): spectra[b][family]=source.exact_counts(row)
        if (t,s) not in maps:
            path=source.grid.pilot.ROOT/row['map_source']
            maps[t,s]=json.loads(path.read_text());dependencies[row['map_source']]=sha(path)
    for b in (8,32):
        if spectra[b]['bch']!=spectra[b]['rm']: raise ValueError('shared exact anchor differs')
    for b in (8,32,64,128,256,512):
        spectra.setdefault(b,{})['random_mean']=random_mean_counts(b)
        spectra[b]['random_caps60']=variance.caps(b,b//2,60)
    for path in (Path(__file__),Path(refresh.__file__),Path(native.__file__),native.LIBRARY,
                 HERE/'activation_refresh_kernel.cpp',Path(onset.__file__),Path(variance.__file__),
                 Path(receipts.__file__),Path(source.__file__)):
        dependencies[path.relative_to(source.grid.pilot.ROOT).as_posix()]=sha(path)
    fingerprint=hashlib.sha256(json.dumps(dependencies,sort_keys=True).encode()).hexdigest()
    keys={(b,64,20,e) for b in spectra for e in range(12,27)}
    keys|={(b,t,s,20) for b in spectra for t,s in maps}
    keys|={(b,64,s,e) for b in spectra for s in (10,12,16,20) for e in (16,18,20,22,24)}
    keys=sorted((b,t,s,e) for b,t,s,e in keys if (1<<e)//(b//2)>=t)
    DIRECTORY.mkdir(exist_ok=True);kernel=native.RefreshKernel();results=[]
    for index,(b,t,s,e) in enumerate(keys,1):
        path=DIRECTORY/f'b{b}_t{t}_s{s}_e{e}.json'
        if path.exists():
            cached=json.loads(path.read_text())
            if cached['input_fingerprint']!=fingerprint: raise ValueError(f'stale checkpoint {path}')
            rows=cached['rows']
        else:
            ac={w:n for w,n in enumerate(maps[t,s]['a_counts']) if w and n}
            rows=evaluate(kernel,b,t,s,e,ac,spectra[b])
            path.write_text(json.dumps(dict(input_fingerprint=fingerprint,rows=rows),indent=2)+'\n')
        results.extend(rows)
        print(f'{index}/{len(keys)} B{b} t{t} s{s} e{e}: '+', '.join(f'{r["family"]}={r["margin_bits"]:.4f}' for r in rows),flush=True)
    csv_path=HERE/'engineering_surfaces.csv'
    with csv_path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
    models={}
    for b,entries in spectra.items():
        for family,counts in entries.items():
            value=onset.spectrum_model(b,counts)
            for shell in value['shells']:
                if not math.isfinite(shell['log_upper']): shell['log_upper']=None
            models[f'{family}_{b}']=value
    payload=dict(status='Q1_ENGINEERING_SURFACES_DIAGNOSTIC',rows=len(results),geometries=len(keys),
                 csv_sha256=sha(csv_path),source_sha256=dependencies,models=models,
                 limitations=['Q1 only; full occupation bounds unchanged.',
                              'Random mean averages over one reused random outer and encoder setup.',
                              'Random caps60 conditional margins exclude the separately reported setup charge.',
                              'Onset idealization omits cancellations and finite-epoch fluctuations.'])
    (HERE/'engineering_surfaces.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')


if __name__=='__main__': main()
