"""Bounded factored-macroblock Q1 screen; exact maps, binary64 bounds only.

The expansion repeats the certified 128-position A two/four times. B has
distinct weight-three columns and is independent of A. No existing producer
or certificate is modified. Single-row screens are not full certificates.
"""
import hashlib
import json
import math
from pathlib import Path
import random
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
ROOT=PARENT.parents[1]
sys.path.insert(0,str(PARENT/'asymmetric'))
import search_feedback as search
wm=search.wm


def spectrum(columns,s):
    """Enumerate all states in 64-position slices, without uint8 overflow."""
    weights=np.zeros(1<<s,dtype=np.uint16)
    for start in range(0,len(columns),64):
        rows=[sum(((c>>j)&1)<<p for p,c in enumerate(columns[start:start+64])) for j in range(s)]
        words=np.zeros(1,dtype=np.uint64)
        for row in rows:words=np.concatenate((words,words^np.uint64(row)))
        weights+=np.bitwise_count(words)
    return {int(w):int(n) for w,n in zip(*np.unique(weights,return_counts=True))}


def record(a,b,s):
    assert len(a)==len(b) and search.rank(b)==s and all(b) and len(set(b))==len(b)
    spec=spectrum(a,s);del spec[0]
    weights={q:sum((c&q).bit_count()&1 for c in a) for q in b}
    cancellation=[[weights[q],weights[q]+1-2*((a[p]&q).bit_count()&1)] for p,q in enumerate(b)]
    return dict(t=len(a),s=s,spectrum=spec,levels=sorted(spec),cancellation=cancellation)


def screen(data):
    tilts=np.arange(-1200,-699,dtype=float)/100
    counts={w:n for w,n in enumerate(search.g.tv.smaller_outer.spectrum()) if w and n}
    epochs=32768//data['t'];assert epochs*data['t']==32768
    output=[]
    for rounds in (1,2,None):
        rz,ra=wm.regions(*wm.transfers(data,np.exp(tilts),rounds),epochs)
        moments=wm.coefficients(rz,ra-math.log(epochs),128)
        for delta in (search.g.F(33,200),search.g.F(19,100)):
            cutoff=(1<<22)*delta.numerator//delta.denominator
            costs=np.minimum(0,moments+cutoff*np.exp(tilts)[:,None])
            indices=np.argmin(costs,axis=0);terms=[math.log(32768*n)+costs[indices[w],w] for w,n in counts.items()]
            dominant=list(counts)[int(np.argmax(terms))]
            output.append(dict(rounds=rounds,distance_target=str(delta),
                q1_margin_bits=-float(np.logaddexp.reduce(terms))/math.log(2),
                dominant_outer_weight=dominant,dominant_log_tilt=float(tilts[indices[dominant]]),
                grid_edge=bool(indices[dominant] in (0,len(tilts)-1))))
    return output


def main():
    path=PARENT/'NO_CONSTANT_MAP.json';base=json.loads(path.read_text());a0=base['columns']
    bank=json.loads((PARENT/'asymmetric/FEEDBACK_SCREEN.json').read_text())
    baseline=next(r['b_columns'] for r in bank['candidates'] if r['name']=='greedy3_2')
    candidates=[('control128',1,baseline)]
    for repeat in (2,4):
        for seed in (1,2):
            candidates.append((f'repeat{repeat}_random3_seed{seed}',repeat,random.Random(seed).sample(search.pool(19,3),128*repeat)))
    results=[]
    for name,repeat,b in candidates:
        a=a0*repeat;data=record(a,b,19);dual=spectrum(b,19)
        row=dict(name=name,t=len(a),s=19,repeats=repeat,a_columns=a,b_columns=b,
            a_spectrum=data['spectrum'],b_dual_spectrum=dual,low_kernel=search.low_kernel(b),
            direct_B_transpose_xors=sum(c.bit_count()-1 for c in b),
            extra_chunk_combine_xors=128*(repeat-1),q1=screen(data))
        results.append(row);print(json.dumps({k:row[k] for k in ('name','low_kernel','q1')}),flush=True)
    sources=[Path(__file__),path,PARENT/'asymmetric/FEEDBACK_SCREEN.json',Path(search.__file__),Path(wm.__file__),Path(search.g.tv.smaller_outer.__file__)]
    payload=dict(status='EXACT_MACRO_MAPS_AND_BINARY64_Q1_NOT_CERTIFICATE',message_bits=1<<20,
        outer=[128,32,32],results=results,
        source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (HERE/'MACRO_Q1.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
