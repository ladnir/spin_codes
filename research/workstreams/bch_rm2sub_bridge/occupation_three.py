"""Compact Q3 bound: exact conditioning from auxiliary Bernoulli supports.

The construction still uses fixed-weight uniform supports. Bernoulli rows
are only a coefficient upper-bound device; every conditioning cost is paid.
No monotonicity assumption, high-weight deletion, or cubic tensor cache.
"""
import argparse
import itertools
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import occupation_two as q2

NAME='t128_s15'
BANDS=((38,40,42),(44,46,48),(50,52,54),tuple(range(56,66,2)),
       tuple(range(66,80,2)),tuple(range(80,102,2)),tuple(range(102,130,2)),
       tuple(range(130,156,2)),tuple(range(156,178,2)),tuple(range(178,200,2)),
       tuple(range(200,212,2)),tuple(range(212,220,2)),(256,))


def region_matrices(t,s,spectrum,z,number,length=8192):
    epoch=list(q2.epochs(t,s,spectrum,z,number))
    den=(1<<s)-1; kappa=number(den)/(den-1); d=min(spectrum)
    m3=number(0)
    for w,c in spectrum.items():
        for v in range(4):
            if v<=w and 3-v<=t-w:
                m3+=number(c*math.comb(w,v)*math.comb(t-w,3-v))*z**(w+3-2*v)/(den*math.comb(t,3))
    epoch.append((number(0),z**3,number(0),z**(d-3)/den,number(0),z**(d-3),
                  kappa*m3/den,number(0),kappa*m3))
    weighted=[tuple(v*math.comb(t,j) for v in matrix) for j,matrix in enumerate(epoch)]
    current=[tuple(number(int(i==j)) for i in range(3) for j in range(3))]+[(number(0),)*9 for _ in range(3)]
    for _ in range(length//t):
        updated=[]
        for degree in range(4):
            terms=[q1.positive_mul(current[j],weighted[degree-j]) for j in range(degree+1)]
            updated.append(tuple(sum((v[k] for v in terms),number(0)) for k in range(9)))
        current=updated
    return [tuple(v/math.comb(length,j) for v in current[j]) for j in range(4)]


def probabilities(ps,number):
    values=[number(1)]
    for p in ps:
        p=number(p.numerator)/p.denominator
        updated=[number(0)]*(len(values)+1)
        for j,v in enumerate(values):
            updated[j]+=v*(1-p); updated[j+1]+=v*p
        values=updated
    return values


def condition_cost(band,p,caps):
    return sum((F(caps[w],math.comb(256,w))/(p**w*(1-p)**(256-w)) for w in band),F(0))


def bands():
    assert sorted(w for band in BANDS for w in band)==list(base.WEIGHTS)
    caps=q2.deterministic_caps(); result=[]
    for band in BANDS:
        p=F(band[0]+band[-1],512)
        result.append((p,condition_cost(band,p,caps)))
    return result


def coefficient(region,ps,number,positions=256):
    probability=probabilities(ps,number)
    matrix=tuple(sum((p*r[k] for p,r in zip(probability,region)),number(0)) for k in range(9))
    result=tuple(number(int(i==j)) for i in range(3) for j in range(3))
    exponent=positions
    while exponent:
        if exponent&1: result=q1.positive_mul(result,matrix)
        exponent>>=1
        if exponent: matrix=q1.positive_mul(matrix,matrix)
    return sum(result[:3],number(0))


def multiplicity(indices):
    return 1 if indices[0]==indices[2] else (3 if len(set(indices))==2 else 6)


def log_coefficient(region,ps):
    probability=probabilities(ps,float)
    matrix=tuple(sum(p*r[k] for p,r in zip(probability,region)) for k in range(9))
    scale=0.
    # Exactly 256 regions: eight squarings, scaled to prevent underflow.
    for _ in range(8):
        matrix=q1.positive_mul(matrix,matrix)
        largest=max(matrix)
        assert largest>0 and math.isfinite(largest)
        matrix=tuple(v/largest for v in matrix)
        scale=2*scale+math.log(largest)
    return math.log(sum(matrix[:3]))+scale


def screen():
    t,s,spectrum=base.load_map(NAME)
    kernel=base.read(base.HERE/'inputs'/f'{NAME}_b_kernel_spectrum.json')['by_total_weight']
    assert all(row['kernel_words']==0 for row in kernel if 1<=row['total_weight']<=3)
    groups=bands(); boxes=list(itertools.combinations_with_replacement(range(len(groups)),3))
    best={box:(math.inf,None) for box in boxes}
    for tenth in range(-85,-54,5):
        lam=math.exp(tenth/10); region=region_matrices(t,s,spectrum,math.exp(-lam),float)
        for box in boxes:
            log=log_coefficient(region,[groups[i][0] for i in box])+209716*lam+sum(math.log(groups[i][1].numerator)-math.log(groups[i][1].denominator) for i in box)
            if log<best[box][0]: best[box]=(log,tenth)
    logsum=base.lse([v[0]+math.log(multiplicity(box)) for box,v in best.items()])+math.log(math.comb(8192,3))
    payload=dict(status='Q3_BERNOULLI_CONDITIONING_SCREEN_ONLY',configuration=NAME,bands=BANDS,
        margin_bits_diagnostic=-logsum/math.log(2),
        boxes=[dict(bands=box,witness_tenth=best[box][1]) for box in boxes],
        source_sha256={'occupation_three.py':base.sha(Path(__file__))})
    base.write_new(base.HERE/'generated'/'t128_s15_q3_conditioning_screen.json',payload)
    print('Q3 conditioning screen margin',payload['margin_bits_diagnostic'],'boxes',len(boxes),flush=True)


def certify():
    import sys
    sys.path.insert(0,str(base.BCH/'code'))
    from audit_bch_q1_full_arb import rational
    path=base.HERE/'generated'/'t128_s15_q3_conditioning_screen.json'; discovery=base.read(path)
    assert discovery['source_sha256']['occupation_three.py']==base.sha(Path(__file__))
    t,s,spectrum=base.load_map(NAME); groups=bands(); ctx.prec=256
    kernel=base.read(base.HERE/'inputs'/f'{NAME}_b_kernel_spectrum.json')['by_total_weight']
    assert all(row['kernel_words']==0 for row in kernel if 1<=row['total_weight']<=3)
    witnesses={tuple(row['bands']):row['witness_tenth'] for row in discovery['boxes']}
    assert set(witnesses)==set(itertools.combinations_with_replacement(range(len(BANDS)),3))
    tables={}; rows=[]; total=F(0)
    for box,tenth in sorted(witnesses.items()):
        assert isinstance(tenth,int) and -120<=tenth<=0
        if tenth not in tables:
            lam=(arb(tenth)/10).exp()
            tables[tenth]=(region_matrices(t,s,spectrum,(-lam).exp(),arb),(209716*lam).exp())
        region,correction=tables[tenth]
        moment=coefficient(region,[groups[i][0] for i in box],arb)*correction
        upper=rational(moment.upper())
        contribution=math.comb(8192,3)*multiplicity(box)*math.prod(groups[i][1] for i in box)*upper
        assert contribution>0
        total+=contribution
        rows.append(dict(bands=box,witness_tenth=tenth,weighted_moment_upper=base.encode(upper),contribution_upper=base.encode(contribution)))
    local,outer=q2.source_paths()
    local += [Path(__file__),path]
    payload=dict(status='OUTWARD_Q3_EXACT_BERNOULLI_CONDITIONING_BOUND',configuration=NAME,
        parameters=dict(occupation=3,message_bits=1<<20,output_bits=1<<21,outer_rows=8192,
                        outer_length=256,distance_cutoff=209716,step_bits=t,state_bits=s),
        bands=[dict(weights=band,p=base.encode(p),conditioning_cost=base.encode(cost)) for band,(p,cost) in zip(BANDS,groups)],
        boxes=rows,Q3_upper=base.encode(total),margin_bits_diagnostic=math.log2(total.denominator)-math.log2(total.numerator),
        Q3_below_2_to_minus_60=total<F(1,1<<60),all_occupations_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer},
        probability_scope='Original fixed-weight row and region permutations; Bernoulli rows are auxiliary and exact binomial conditioning costs are paid',
        monotonicity_assumed=False,weight_truncation_used=False)
    base.write_new(base.HERE/'generated'/'t128_s15_q3_conditioning_outward.json',payload)
    print('OUTWARD Q3 margin',payload['margin_bits_diagnostic'],'below 2^-60',payload['Q3_below_2_to_minus_60'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('mode',choices=('screen','certify'))
    args=parser.parse_args(); (screen if args.mode=='screen' else certify)()
