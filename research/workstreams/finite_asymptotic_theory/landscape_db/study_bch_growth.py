"""Sequential BCH growth study; all generated results remain local and ignored."""
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import activation_q1 as old
import activation_q1_native as native
import activation_q1_refresh as refresh
import bch_growth_model as model
import read_grid_receipts as receipts
import run_occupation_grid as source

HERE=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    _,observations,_=receipts.snapshot()
    chosen=[r for r in observations.values() if 'BCH' in r['series']]
    spectra={};maps={};templates={};dependencies={}
    for name in source.grid.BASES:
        dependencies.update(source.grid.verify(HERE/name)['source_sha256'])
    for row in chosen:
        b,t,s=int(row['block_bits']),int(row['step_bits']),int(row['state_bits'])
        if b not in spectra:spectra[b]=source.exact_counts(row)
        if (t,s) not in maps:
            path=source.grid.pilot.ROOT/row['map_source'];maps[t,s]=json.loads(path.read_text())
            dependencies[row['map_source']]=sha(path)
        templates[b,t,s,int(row['message_exponent'])]=row
    keys={(b,64,20,e) for b in spectra for e in range(12,27)}
    keys|={(b,t,s,20) for b in spectra for t,s in maps}
    results=[];old.coefficient_logs=native.Q1Kernel().coefficients
    for b,t,s,e in sorted(keys):
        length=(1<<e)//(b//2)
        if length%t or length<t:continue
        config=maps[t,s];ac={w:n for w,n in enumerate(config['a_counts']) if w and n}
        tilts=np.unique(np.r_[np.arange(-160,1)/10,np.arange(-120,141)/20-math.log(length)])
        after=refresh.screen(t,s,ac,b,length,spectra[b],tilts,refine=True)
        tilts=after.pop('log_tilts_evaluated')
        before=old.screen(t,s,ac,b,b//2,e,{w:math.log(n) for w,n in spectra[b].items()},tilts)
        if after['margin_bits']<before['margin_bits']-1e-8:
            raise ArithmeticError('uniform-refresh refinement weakened the old bound')
        published=templates.get((b,t,s,e),{})
        results.append(dict(block_bits=b,dimension=b//2,step_bits=t,state_bits=s,message_exponent=e,
            outer_rows=length,epochs_per_region=length//t,margin_bits=after['margin_bits'],
            old_same_witness_margin_bits=before['margin_bits'],published_margin_bits=published.get('margin_bits',''),
            intercept_bits=after['margin_bits']+math.log2(length),
            old_intercept_bits=before['margin_bits']+math.log2(length),
            dominant_weight=after['dominant_weight'],dominant_log_surprisal=after['dominant_log_surprisal'],
            scaled_tilt=math.exp(after['dominant_log_surprisal'])*length,
            witness_at_edge=int(after['dominant_log_surprisal'] in (min(tilts),max(tilts))),
            density_penalty_ceiling_bits=(b*length/t)*math.log1p(1/((1<<s)-2))/math.log(2),
            evidence='Q1 binary64 diagnostic; four-state refresh; no full-distance claim'))
        print(f'B{b} t{t} s{s} k2^{e}: {after["margin_bits"]:.6f} bits; improvement {after["margin_bits"]-before["margin_bits"]:.6f}',flush=True)
    csv_path=HERE/'bch_growth.csv'
    with csv_path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
    models={str(b):model.spectrum_model(b,counts) for b,counts in spectra.items()}
    # JSON uses null for the zero-probability idealization's infinite margin.
    for value in models.values():
        for shell in value['shells']:
            if not math.isfinite(shell['log_upper']):shell['log_upper']=None
    fixed=[r for r in results if r['step_bits']==64 and r['state_bits']==20 and 16<=r['message_exponent']<=24]
    summary={}
    for b in spectra:
        rows=[r for r in fixed if r['block_bits']==b]
        x=np.array([r['message_exponent'] for r in rows]);y=np.array([r['margin_bits'] for r in rows])
        summary[str(b)]=dict(slope_bits_per_k_doubling=float(np.polyfit(x,y,1)[0]),
            intercept_min=min(r['intercept_bits'] for r in rows),intercept_max=max(r['intercept_bits'] for r in rows))
    for path in (Path(__file__),Path(refresh.__file__),Path(model.__file__),Path(old.__file__),
                 Path(source.__file__),Path(receipts.__file__),native.LIBRARY):
        dependencies[path.relative_to(source.grid.pilot.ROOT).as_posix()]=sha(path)
    payload=dict(status='BCH_Q1_GROWTH_DIAGNOSTIC',rows=len(results),models=models,
        fixed_t64_s20_trends=summary,source_sha256=dependencies,csv_sha256=sha(csv_path),
        limitations=['Q1 only; higher occupations require separate refinement.',
            'Persistent-state models neglect cancellations and finite-epoch output fluctuations.',
            'No fit to BCH-256 or claim of an exact larger BCH spectrum.'])
    (HERE/'bch_growth.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
