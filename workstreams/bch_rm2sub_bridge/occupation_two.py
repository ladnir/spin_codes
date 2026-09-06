"""Activation-aware occupation-two bound for the fixed BCH/RM2Sub model.

Exact region support counts, Arb epoch/region bounds, upward binary64 pair
counts, and exact rational spectrum aggregation. No M22 probability reuse.
"""
import argparse
import math
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from flint import arb,ctx
import bridge as base
import activation_bridge as q1

sys.path.insert(0,str(base.BCH/'code'))
from certify_random_inner_threshold_outward import add_up,mul_up,upper_float
from certify_bch_m25_closure import deterministic_caps


def epochs(t,s,spectrum,z,number):
    zero,one=q1.matrices(t,s,spectrum,z,number)
    den=(1<<s)-1; kappa=number(den)/(den-1); distance=min(spectrum)
    m2=number(0)
    for w,count in spectrum.items():
        for overlap in range(3):
            if overlap<=w and 2-overlap<=t-w:
                m2+=number(count*math.comb(w,overlap)*math.comb(t-w,2-overlap))*z**(w+2-2*overlap)/(den*math.comb(t,2))
    two=(number(0),z*z,number(0),z**(distance-2)/den,number(0),z**(distance-2),
         kappa*m2/den,number(0),kappa*m2)
    return zero,one,two


def regions(t,s,spectrum,z,number,length=8192):
    assert length%t==0
    kernels=epochs(t,s,spectrum,z,number)
    weighted=[tuple(number(math.comb(t,j))*x for x in kernels[j]) for j in range(3)]
    current=[tuple(number(int(i==j)) for i in range(3) for j in range(3)),(number(0),)*9,(number(0),)*9]
    for _ in range(length//t):
        updated=[]
        for degree in range(3):
            terms=[q1.positive_mul(current[degree-j],weighted[j]) for j in range(degree+1)]
            updated.append(tuple(sum((term[k] for term in terms),number(0)) for k in range(9)))
        current=updated
    return [tuple(v/math.comb(length,j) for v in current[j]) for j in range(3)]


def row_product(rows,matrix):
    result=np.empty_like(rows)
    # Fixed-width inner kernel; each primitive nonnegative operation rounds up.
    result[...,0]=add_up(add_up(mul_up(rows[...,0],matrix[0]),mul_up(rows[...,1],matrix[3])),mul_up(rows[...,2],matrix[6]))
    result[...,1]=add_up(add_up(mul_up(rows[...,0],matrix[1]),mul_up(rows[...,1],matrix[4])),mul_up(rows[...,2],matrix[7]))
    result[...,2]=add_up(add_up(mul_up(rows[...,0],matrix[2]),mul_up(rows[...,1],matrix[5])),mul_up(rows[...,2],matrix[8]))
    return result


def pair_counts(region,positions=256):
    current=np.zeros((positions+1,positions+1,3)); scratch=np.zeros_like(current)
    current[0,0,0]=1.
    for n in range(positions):
        size=n+1; old=current[:size,:size]; updated=scratch[:size+1,:size+1]
        updated.fill(0)
        updated[:size,:size]=row_product(old,region[0])
        mixed=row_product(old,region[1])
        updated[1:,:size]=add_up(updated[1:,:size],mixed)
        updated[:size,1:]=add_up(updated[:size,1:],mixed)
        updated[1:,1:]=add_up(updated[1:,1:],row_product(old,region[2]))
        current,scratch=scratch,current
    assert np.isfinite(current).all() and (current>=0).all()
    inverse=np.array([upper_float(arb(1)/math.comb(positions,w)) for w in range(positions+1)])
    total=add_up(add_up(current[...,0],current[...,1]),current[...,2])
    return mul_up(mul_up(total,inverse[:,None]),inverse[None,:])


def aggregate(conditional):
    caps=deterministic_caps(); terms=[]; total=F(0)
    placements=math.comb(base.ROWS,2)
    for a in base.WEIGHTS:
        for b in base.WEIGHTS:
            term=placements*caps[a]*caps[b]*F.from_float(float(conditional[a,b]))
            total+=term; terms.append((term,a,b))
    return total,sorted(terms,reverse=True)[:10]


def source_paths():
    local=[Path(__file__),Path(base.__file__),Path(q1.__file__),base.HERE/'inputs/manifest.json']
    outer=[base.BCH/p for p in ('code/certify_random_inner_threshold_outward.py',
           'code/certify_bch_m25_closure.py','generated/bch256_closure_deterministic_envelope.json',
           'generated/johnson_n256_w52_d38.json','generated/johnson_n256_w54_d38.json')]
    return local,outer


def run(name,tenth):
    t,s,spectrum=base.load_map(name)
    # load_map independently checks all B columns are nonzero and distinct.
    ctx.prec=256; lam=(arb(tenth)/10).exp(); z=(-lam).exp()
    region=regions(t,s,spectrum,z,arb)
    rounded=np.array([[upper_float(v) for v in matrix] for matrix in region])
    values=pair_counts(rounded)
    conditional=np.minimum(1.,mul_up(values,upper_float((base.CUTOFF*lam).exp())))
    total,dominant=aggregate(conditional)
    stem=base.HERE/'generated'/f'{name}_q2_j{tenth}'
    output=stem.with_suffix('.json'); cache=stem.with_suffix('.npz')
    assert not output.exists() and not cache.exists(),'Write-once output already exists'
    with cache.open('xb') as stream:
        np.savez_compressed(stream,conditional_pair_upper=conditional,region_upper=rounded)
    local,outer=source_paths()
    payload=dict(status='ACTIVATION_AWARE_OUTWARD_Q2_BOUND',configuration=name,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,
                        distance_cutoff=209716,occupation=2,step_bits=t,state_bits=s),
        witness_tenth=tenth,precision_bits=ctx.prec,Q2_upper=base.encode(total),
        margin_bits_diagnostic=math.log2(total.denominator)-math.log2(total.numerator),
        Q2_below_2_to_minus_60=total<F(1,1<<60),all_occupations_certified=False,
        dominant_pairs=[dict(weights=[a,b],upper=base.encode(v)) for v,a,b in dominant],
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local+[cache]},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer},
        model='Fresh nonzero multiplier each epoch; output precedes update; initial state zero; independent row and region permutations',
        support_law='Two independent outer support sets; a region with two bits is a uniform two-subset of 8192 positions, not two draws with replacement')
    base.write_new(output,payload)
    print(name,'Q2 margin',payload['margin_bits_diagnostic'],'below 2^-60',payload['Q2_below_2_to_minus_60'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--configuration',choices=base.CONFIGS,default='t128_s15')
    parser.add_argument('--tenth',type=int,default=-75)
    args=parser.parse_args()
    assert -120<=args.tenth<=0
    run(args.configuration,args.tenth)
