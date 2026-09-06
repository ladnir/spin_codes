"""Small batch receipt: all composition witnesses live in hashed screens."""
import argparse
import math
import sys
from fractions import Fraction as F
from pathlib import Path

from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import general_occupancy as model
from certify_q3_compact import dyadic_upper

sys.path.insert(0,str(base.BCH/'code'))
from audit_bch_q1_full_arb import rational
OUTPUT=base.HERE/'generated'/'general_q4_q16_outward.json'
SCREENS=[base.HERE/'generated'/f'general_q{a}_q{b}_screen.json' for a,b in ((4,16),(5,15))]


def density_cost(band,p,caps):
    p=arb(p.numerator)/p.denominator
    values=[arb(caps[w])/(math.comb(256,w)*p**w*(1-p)**(256-w)) for w in band]
    return dyadic_upper(max(rational(v.upper()) for v in values))


def build(verify=False):
    if not verify:assert not OUTPUT.exists()
    saved=base.read(OUTPUT) if verify else None
    if saved:
        for filename,digest in saved['local_sha256'].items():assert base.sha(base.HERE/filename)==digest
        for filename,digest in saved['outer_sha256'].items():assert base.sha(base.BCH/filename)==digest
    discovery={}
    for path in SCREENS:
        data=base.read(path)
        assert data['source_sha256']['general_occupancy.py']==base.sha(Path(model.__file__))
        for row in data['rows']:
            assert row['occupation'] not in discovery
            discovery[row['occupation']]=row
    assert sorted(discovery)==list(range(4,17))
    t,s,spectrum=base.load_map(model.NAME);kernel=model.kernel_spectrum();caps=model.q2.deterministic_caps()
    ctx.prec=512 if verify else 256;tables={};results=[];combined=F(0)
    for occupancy,row in sorted(discovery.items()):
        assert [tuple(r['counts']) for r in row['compositions']]==list(model.compositions(occupancy))
        total=F(0);worst=F(0);worst_counts=None
        for witness in row['compositions']:
            counts=witness['counts'];j=witness['witness_tenth'];ps=[base.decode(v) for v in witness['p']]
            assert len(counts)==len(ps)==5 and sum(counts)==occupancy
            assert all(0<p<1 for p in ps) and isinstance(j,int) and -120<=j<=0
            if j not in tables:
                lam=(arb(j)/10).exp()
                tables[j]=(model.regions(t,s,spectrum,kernel,(-lam).exp(),arb,16),(209716*lam).exp())
            region,correction=tables[j]
            dist=model.distribution(counts,ps,arb)
            matrix=tuple(sum((p*r[k] for p,r in zip(dist,region)),arb(0)) for k in range(9))
            for _ in range(8):matrix=q1.positive_mul(matrix,matrix)
            moment=rational((sum(matrix[:3],arb(0))*correction).upper())
            costs=[density_cost(band,p,caps) if count else F(1) for band,p,count in zip(model.BANDS,ps,counts)]
            contribution=math.comb(8192,occupancy)*model.multiplicity(counts)*moment*math.prod(c**n for c,n in zip(costs,counts))
            assert contribution>0;total+=contribution
            if contribution>worst:worst=contribution;worst_counts=counts
        margin=math.log2(total.denominator)-math.log2(total.numerator)
        if verify:
            old=next(r for r in saved['rows'] if r['occupation']==occupancy)
            assert total<=base.decode(old['upper'])
        assert total<F(1,1<<150)
        results.append(dict(occupation=occupancy,composition_count=len(row['compositions']),upper=base.encode(total),
            margin_bits_diagnostic=margin,worst_counts=worst_counts,worst_contribution=base.encode(worst)))
        combined+=total
        print('Q',occupancy,'verified' if verify else 'certified','margin',margin,flush=True)
    old_paths=[base.HERE/'generated'/f for f in ('t128_s15_activation_q1_outward.json','t128_s15_q2_j-75.json','t128_s15_q3_compact_outward.json')]
    partial=combined+sum((base.decode(base.read(p)[key]) for p,key in zip(old_paths,('Q1_upper','Q2_upper','Q3_upper'))),F(0))
    assert partial<F(1,1<<49)
    if verify:
        assert combined<=base.decode(saved['Q4_through_Q16_upper'])
        assert partial<=base.decode(saved['Q1_through_Q16_upper'])
        print('All 4..16 compositions replayed at 512 bits; partial 1..16 bound <2^-49',flush=True);return
    local,outer=model.q2.source_paths()
    local += [Path(__file__),Path(model.__file__),base.HERE/'certify_q3_compact.py']+SCREENS+old_paths
    payload=dict(status='GENERAL_KERNEL_OUTWARD_OCCUPANCIES_4_THROUGH_16',configuration=model.NAME,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,
                        distance_cutoff=209716,step_bits=t,state_bits=s),
        bands=model.BANDS,rows=results,Q4_through_Q16_upper=base.encode(combined),
        Q1_through_Q16_upper=base.encode(partial),Q1_through_Q16_below_2_to_minus_49=True,
        general_formula_occupancy_range=[1,8192],numerically_certified_occupancy_range=[1,16],
        all_occupations_certified=False,local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer})
    base.write_new(OUTPUT,payload)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();build(args.verify)
