"""Search expansion subspaces with weights in [24,40], keeping feedback fixed."""
from collections import Counter
from functools import lru_cache
from pathlib import Path
import argparse
import json
import random
import math
from fractions import Fraction as F
import model
import refine
from flint import arb,ctx


def nullspace(rows, width):
    pivots={}
    for value in rows:
        while value:
            j=value.bit_length()-1
            if j in pivots: value^=pivots[j]
            else:
                pivots[j]=value
                break
    result=[]
    for j in range(width):
        if j in pivots: continue
        vector=1 << j
        for k,row in sorted(pivots.items()):
            if (vector&row).bit_count()&1: vector^=1 << k
        assert all(not ((vector&row).bit_count()&1) for row in rows)
        result.append(vector)
    return result


def record_from(a,b,s,**metadata):
    maps=model.base.independent.g.tv.fixed.maps
    sa=maps.spectrum(maps.generators(a,s))
    assert sa[0]==1 and sum(sa.values())==1 << s and 64 not in sa
    assert len(a)==len(b)==len(set(a))==len(set(b))==64 and all(a) and all(b)
    assert model.base.independent.search.rank(a)==model.base.independent.search.rank(b)==s
    images=[sum(((c&q).bit_count()&1)<<j for j,c in enumerate(a)) for q in b]
    spectrum={w:n for w,n in sa.items() if w}
    return dict(t=64,s=s,spectrum=spectrum,levels=sorted(spectrum),
        expansion_columns=a,feedback_columns=b,transvection_rounds=1,
        cancellation=[[v.bit_count(),(v^(1<<j)).bit_count()] for j,v in enumerate(images)],
        **metadata)


def search(s,parent_s=15,limit=20000):
    assert 6<s<parent_s<=20
    old=model.study.core.inner(64,parent_s)
    a=old['expansion_columns']
    weights=model.base.independent.search.wm.all_weights(a,parent_s)
    bad=[q for q in range(1,1 << parent_s) if not 24<=weights[q]<=40]
    # Preserve all six linear generators. Constraints affect only quadratics.
    bad_quadratic=sorted(set(q>>6 for q in bad))
    assert 0 not in bad_quadratic
    rng=random.Random(20260920+100*s+parent_s)
    codimension=parent_s-s
    for trial in range(limit):
        constraints=[rng.randrange(1,1 << (parent_s-6)) for _ in range(codimension)]
        if any(all(not ((q&c).bit_count()&1) for c in constraints) for q in bad_quadratic):
            continue
        basis=nullspace(constraints,parent_s-6)
        if len(basis)!=s-6: continue
        basis=[1 << j for j in range(6)]+[v<<6 for v in basis]
        new_a=[sum(((column&v).bit_count()&1)<<j for j,v in enumerate(basis)) for column in a]
        if not all(new_a): continue
        record=record_from(new_a,model.study.core.inner(64,s)['feedback_columns'],s,
            parent_s=parent_s,parent_basis=basis,constraints=constraints,trial=trial,
            parent_bad_states=len(bad),parent_bad_quadratics=len(bad_quadratic))
        assert min(record['spectrum'])>=24 and max(record['spectrum'])<=40
        return record
    return None


class Engine(refine.Engine):
    def __init__(self, record):
        # Set up the explicit fixed maps, with no mutation of old module globals.
        self.record=record
        self.t,self.s=record['t'],record['s']
        assert self.t==64
        self.a_columns,self.columns=record['expansion_columns'],record['feedback_columns']
        exact=model.base.independent.g.tv.fixed.maps
        sa=exact.spectrum(exact.generators(self.a_columns,self.s))
        sb=exact.spectrum(exact.generators(self.columns,self.s))
        kernel=exact.dual_spectrum(sb,self.t,self.s)
        self.spectrum={w:n for w,n in sa.items() if w}
        assert self.spectrum=={int(w):n for w,n in record['spectrum'].items()}
        assert sa[0]==sb[0]==1 and 64 not in sa
        assert sum(sa.values())==sum(sb.values())==1 << self.s
        self.b_spectrum={w:n for w,n in sb.items() if w}
        self.kernel=[kernel.get(j,0) for j in range(self.t+1)]
        assert self.kernel[:3]==[1,0,0] and sum(self.kernel)==1 << (self.t-self.s)
        self.caps=model.base.independent.g.fiber_caps(self.t,self.s,self.b_spectrum,self.kernel)
        self.low=model.base.independent.search.low_cancellation(self.a_columns,self.columns,self.s)
        self.levels=sorted(self.spectrum)
        self.n,self.m=len(self.levels)+2,(1 << self.s)-1
        self.exponent,self.length=16,512
        self.outer_length,self.output_bits=256,131072
        self.cutoff=13107


def run(a):
    if a.output.exists(): raise FileExistsError(a.output)
    ctx.prec=256
    rows=[]
    for s in a.s:
        record=search(s,a.parent_s,a.limit)
        if record is None:
            print('no subspace found',s,flush=True)
            rows.append(dict(s=s,found=False));continue
        q1=model.study.core.evaluate(record,256,16,sharp=True)
        row=dict(s=s,found=True,inner=record,q1=q1)
        print('subspace',s,'trial',record['trial'],'spectrum',record['spectrum'],
              'diagnostic Q1',q1['q1_margin_bits'],flush=True)
        rows.append(row)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    model.base.base.write_new(a.output,dict(status='K16_EXACT_SUBSPACE_SEARCH_WITH_DIAGNOSTIC_Q1',
        rows=rows,full_distance_proved=False,source_sha256=model.base.sources()))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--s',type=int,nargs='+',default=[12,13,14])
    p.add_argument('--parent-s',type=int,default=15)
    p.add_argument('--limit',type=int,default=20000)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())
