"""Outward adaptive-group range bound; no composition enumeration."""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path

from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import general_occupancy as general
from general_batch_certificate import density_cost

SCREEN=base.HERE/'generated'/'adaptive_q17_q64_refined_screen.json'
OUTPUT=base.HERE/'generated'/'adaptive_q17_q64_outward.json'


def adaptive_kernel(region,ps,roots):
    current=region
    left=[r*(1-p) for p,r in zip(ps,roots)]
    right=[r*p for p,r in zip(ps,roots)]
    for _ in range(len(region)-1):
        updated=[]
        for a,b in zip(current,current[1:]):
            # upper() is an exact singleton, so max has no interval-order ambiguity.
            updated.append(tuple(max((x*a[k]+y*b[k]).upper() for x,y in zip(left,right)) for k in range(9)))
        current=updated
    return current[0]


def build(verify=False):
    if not verify:assert not OUTPUT.exists()
    old=base.read(OUTPUT) if verify else None
    if old:
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        for name,digest in old['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
    screen=base.read(SCREEN)
    for name,digest in screen['source_sha256'].items():assert base.sha(base.HERE/name)==digest
    assert [r['occupation'] for r in screen['rows']]==list(range(17,65))
    maximum={}
    for row in screen['rows']:
        j=row['witness_tenth'];assert isinstance(j,int) and -120<=j<=0
        maximum[j]=max(maximum.get(j,0),row['occupation'])
    t,s,spectrum=base.load_map(general.NAME);kernel=general.kernel_spectrum();caps=general.q2.deterministic_caps()
    ctx.prec=512 if verify else 256;cache={};total=F(0);rows=[]
    from audit_bch_q1_full_arb import rational
    for row in screen['rows']:
        q=row['occupation'];j=row['witness_tenth']
        if j not in cache:
            lam=(arb(j)/10).exp()
            cache[j]=(general.regions(t,s,spectrum,kernel,(-lam).exp(),arb,maximum[j]),(209716*lam).exp())
        region,correction=cache[j]
        ps=[base.decode(p) for p in row['p']];assert len(ps)==5 and all(0<p<1 for p in ps)
        costs=[density_cost(band,p,caps) for band,p in zip(general.BANDS,ps)]
        roots=[((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
        probabilities=[arb(p.numerator)/p.denominator for p in ps]
        matrix=adaptive_kernel(region[:q+1],probabilities,roots)
        for _ in range(8):matrix=q1.positive_mul(matrix,matrix)
        bound=math.comb(8192,q)*5**q*rational((sum(matrix[:3],arb(0))*correction).upper())
        assert 0<bound<F(1,1<<142)
        if old:assert bound<=base.decode(old['rows'][q-17]['upper'])
        total+=bound
        rows.append(dict(occupation=q,upper=base.encode(bound),margin_bits_diagnostic=math.log2(bound.denominator)-math.log2(bound.numerator),
                         witness_tenth=j,p=row['p']))
        if q%8==0 or q==17:print('Adaptive Q',q,'verified' if verify else 'certified','margin',rows[-1]['margin_bits_diagnostic'],flush=True)
    previous=base.HERE/'generated'/'general_q4_q16_outward.json'
    partial=total+base.decode(base.read(previous)['Q1_through_Q16_upper'])
    assert total<F(1,1<<138) and partial<F(1,1<<49)
    if old:
        assert total<=base.decode(old['Q17_through_Q64_upper']) and partial<=base.decode(old['Q1_through_Q64_upper'])
        print('All adaptive occupancies 17..64 replayed at 512 bits; total 1..64 <2^-49',flush=True);return
    local,outer=general.q2.source_paths()
    local += [Path(__file__),Path(general.__file__),base.HERE/'general_batch_certificate.py',
              base.HERE/'certify_q3_compact.py',SCREEN,previous]
    local += [base.HERE/name for name in screen['source_sha256']]
    payload=dict(status='OUTWARD_ADAPTIVE_GROUP_BOUND_Q17_THROUGH_Q64',configuration=general.NAME,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,
                        step_bits=t,state_bits=s,distance_cutoff=209716),bands=general.BANDS,rows=rows,
        Q17_through_Q64_upper=base.encode(total),Q1_through_Q64_upper=base.encode(partial),
        Q1_through_Q64_below_2_to_minus_49=True,all_occupations_certified=False,
        method='Entrywise adaptive choice upper-bounds every fixed ordered group assignment; multiply by 5^Q',
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer})
    base.write_new(OUTPUT,payload)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();build(args.verify)
