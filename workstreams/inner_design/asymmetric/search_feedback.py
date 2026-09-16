"""Bounded independent-B search; exact audits and binary64 Q1 diagnostics.

No certified sources are edited. Sparse column costs are arithmetic screens,
not timing estimates. Matrices below use explicit independent A and B columns.
"""
from collections import Counter,defaultdict
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import general_occupancies as g
import weight_memory as wm


def image(a_columns,q):
    return sum(((c&q).bit_count()&1)<<p for p,c in enumerate(a_columns))


def inject(b_columns,x):
    q=0
    while x:
        bit=x&-x;x^=bit;q^=b_columns[bit.bit_length()-1]
    return q


def rank(words):
    pivots={}
    for word in words:
        while word:
            p=word.bit_length()-1
            if p not in pivots:pivots[p]=word;break
            word^=pivots[p]
    return len(pivots)


def low_kernel(columns):
    assert all(columns) and len(set(columns))==len(columns)
    pairs=Counter(a^b for a,b in itertools.combinations(columns,2))
    triples=sum(pairs[c] for c in columns)
    quads=sum(n*(n-1)//2 for n in pairs.values())
    assert triples%3==quads%3==0
    return dict(weight1=0,weight2=0,weight3=triples//3,weight4=quads//3,
                maximum_pair_fiber=max(pairs.values()))


def low_cancellation(a_columns,b_columns,s,maximum=2):
    weights=wm.all_weights(a_columns,s)
    images=[image(a_columns,b) for b in b_columns];result={}
    for j in range(1,maximum+1):
        by_syndrome=defaultdict(Counter);by_weight=defaultdict(Counter)
        for support in itertools.combinations(range(len(b_columns)),j):
            syndrome=0;emitted=0
            for p in support:syndrome^=b_columns[p];emitted^=images[p]^(1<<p)
            if not syndrome:continue
            w=emitted.bit_count();by_syndrome[syndrome][w]+=1
            by_weight[int(weights[syndrome])][w]+=1
        result[j]=dict(patterns=sorted({tuple(sorted(v.items())) for v in by_syndrome.values()}),
                      by_weight=dict(by_weight),
                      nonzero_input_count=sum(sum(v.values()) for v in by_syndrome.values()))
    return result


def pool(s,weight):return [sum(1<<j for j in support) for support in itertools.combinations(range(s),weight)]


def greedy_triples(s,t,seed):
    rng=random.Random(seed);remaining=pool(s,3);chosen=[];pairs=Counter();degrees=[0]*s
    for _ in range(t):
        trial=rng.sample(remaining,min(128,len(remaining)))
        def score(c):
            return (sum(pairs[c^b] for b in chosen),sum(degrees[j] for j in range(s) if c>>j&1))
        c=min(trial,key=score)
        for b in chosen:pairs[c^b]+=1
        chosen.append(c);remaining.remove(c)
        for j in range(s):degrees[j]+=(c>>j)&1
    rng.shuffle(chosen)
    return chosen


def q1(a_columns,b_columns,s):
    assert all(b_columns) and len(set(b_columns))==len(b_columns)
    weights=wm.all_weights(a_columns,s)
    spectrum={int(w):int(n) for w,n in zip(*np.unique(weights,return_counts=True)) if w}
    record=dict(t=len(a_columns),s=s,spectrum=spectrum,levels=sorted(spectrum),
                cancellation=[[int(weights[b]),(image(a_columns,b)^(1<<p)).bit_count()] for p,b in enumerate(b_columns)])
    # Fixed common tilts from the balanced certificate; discovery only.
    tilts=np.array([-9.39,-9.66]);lam=np.exp(tilts);length=32768;epochs=length//128
    rz,ra=wm.regions(*wm.transfers(record,lam,1),epochs)
    moments=wm.coefficients(rz,ra-math.log(epochs),128)
    outer={w:n for w,n in enumerate(g.tv.smaller_outer.spectrum()) if w and n}
    result=[]
    for i,delta in enumerate((g.F(33,200),g.F(19,100))):
        cutoff=128*length*delta.numerator//delta.denominator
        bound=np.logaddexp.reduce([math.log(length*n)+moments[i,w]+cutoff*lam[i] for w,n in outer.items()])
        result.append(dict(distance_target=str(delta),log_surprisal=float(tilts[i]),margin_bits=-float(bound)/math.log(2)))
    return result


def audit(a_columns,b_columns,s,name):
    t=len(a_columns);assert len(b_columns)==t and rank(b_columns)==s
    b_weights=wm.all_weights(b_columns,s)
    spectrum={int(w):int(n) for w,n in zip(*np.unique(b_weights,return_counts=True))}
    kernel=g.tv.fixed.maps.dual_spectrum(spectrum,t,s)
    low=low_kernel(b_columns)
    for j in range(1,5):assert kernel.get(j,0)==low[f'weight{j}']
    ba=[inject(b_columns,image(a_columns,1<<j)) for j in range(s)]
    return dict(name=name,t=t,s=s,b_columns=b_columns,rank=s,
                column_weight_histogram=dict(sorted(Counter(c.bit_count() for c in b_columns).items())),
                direct_transpose_state_sum_xors=sum(c.bit_count()-1 for c in b_columns),
                ba_zero=not any(ba),ba_rank=rank(ba),low_kernel=low,
                dual_minimum_weight=min(w for w in spectrum if w),dual_spectrum=spectrum,
                kernel_minimum_weight=min(w for w,n in kernel.items() if w and n),kernel_spectrum=kernel,
                q1=q1(a_columns,b_columns,s))


def main():
    path=HERE.parent/'NO_CONSTANT_MAP.json';a_columns=json.loads(path.read_text())['columns'];s=19;t=128
    candidates=[('symmetric_control',a_columns)]
    rng=random.Random(1)
    candidates.append(('weight2_plus_one',rng.sample(pool(s,2),127)+[1]))
    candidates.append(('mixed2_3',rng.sample(pool(s,2),64)+rng.sample(pool(s,3),64)))
    for seed in (1,2):
        candidates.append((f'random3_{seed}',random.Random(seed).sample(pool(s,3),t)))
        candidates.append((f'greedy3_{seed}',greedy_triples(s,t,seed)))
    results=[]
    for name,columns in candidates:
        row=audit(a_columns,columns,s,name);results.append(row)
        print(json.dumps({k:row[k] for k in ('name','direct_transpose_state_sum_xors','low_kernel','dual_minimum_weight','q1')}),flush=True)
    sources=[Path(__file__),path,Path(g.__file__),Path(wm.__file__),Path(g.tv.__file__),Path(g.tv.fixed.maps.__file__)]
    payload=dict(status='EXACT_MAP_AUDITS_AND_BINARY64_Q1_ONLY',fixed_a='NO_CONSTANT_MAP.json',
                 candidates=results,source_sha256={p.relative_to(g.tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (HERE/'FEEDBACK_SCREEN.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
