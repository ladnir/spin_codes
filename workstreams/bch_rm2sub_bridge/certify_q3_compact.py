"""Compact dyadic Q3 aggregation; rational conditioning rounded only upward."""
import argparse
import itertools
import math
import sys
from fractions import Fraction as F
from pathlib import Path

from flint import arb,ctx
import bridge as base
import occupation_three as q3
import optimized_conditioning as search

sys.path.insert(0,str(base.BCH/'code'))
from audit_bch_q1_full_arb import rational
OUTPUT=base.HERE/'generated'/'t128_s15_q3_compact_outward.json'


def dyadic_upper(value,bits=192):
    assert value>0
    exponent=value.numerator.bit_length()-value.denominator.bit_length()
    shift=bits-exponent
    scale=F(1<<shift) if shift>=0 else F(1,1<<-shift)
    scaled=value*scale
    return F(-(-scaled.numerator//scaled.denominator))/scale


def build():
    assert not OUTPUT.exists()
    path=base.HERE/'generated'/'t128_s15_q3_optimized_screen.json'; screen=base.read(path)
    for filename,digest in screen['source_sha256'].items():assert base.sha(base.HERE/filename)==digest
    t,s,spectrum=base.load_map(q3.NAME); caps=q3.q2.deterministic_caps(); ctx.prec=256
    kernel=base.read(base.HERE/'inputs'/f'{q3.NAME}_b_kernel_spectrum.json')['by_total_weight']
    assert all(row['kernel_words']==0 for row in kernel if 1<=row['total_weight']<=3)
    assert sorted(w for band in q3.BANDS for w in band)==list(base.WEIGHTS)
    expected=list(itertools.combinations_with_replacement(range(len(q3.BANDS)),3))
    assert [tuple(r['bands']) for r in screen['boxes']]==expected
    tables={};rows=[];total=F(0)
    for row in screen['boxes']:
        box=tuple(row['bands']); tenth=row['witness_tenth']; ps=[base.decode(v) for v in row['p']]
        assert len(ps)==3 and -120<=tenth<=0
        assert all(0<p<=1 and (p<1 or q3.BANDS[i]==(256,)) for i,p in zip(box,ps))
        if tenth not in tables:
            lam=(arb(tenth)/10).exp()
            tables[tenth]=(q3.region_matrices(t,s,spectrum,(-lam).exp(),arb),(209716*lam).exp())
        region,correction=tables[tenth]
        moment=rational((q3.coefficient(region,ps,arb)*correction).upper())
        costs=[dyadic_upper(q3.condition_cost(q3.BANDS[i],p,caps)) for i,p in zip(box,ps)]
        contribution=math.comb(8192,3)*q3.multiplicity(box)*math.prod(costs)*moment
        assert contribution>0;total+=contribution
        rows.append(dict(bands=box,witness_tenth=tenth,p=row['p'],moment_upper=base.encode(moment),
                         conditioning_upper=[base.encode(v) for v in costs],contribution_upper=base.encode(contribution)))
    first=base.read(base.HERE/'generated'/'t128_s15_activation_q1_outward.json')
    second=base.read(base.HERE/'generated'/'t128_s15_q2_j-75.json')
    partial=base.decode(first['Q1_upper'])+base.decode(second['Q2_upper'])+total
    local,outer=q3.q2.source_paths()
    local += [Path(__file__),Path(q3.__file__),Path(search.__file__),path,
              base.HERE/'generated'/'t128_s15_activation_q1_outward.json',base.HERE/'generated'/'t128_s15_q2_j-75.json']
    payload=dict(status='OUTWARD_Q3_COMPACT_CONDITIONING_BOUND',configuration=q3.NAME,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,
                        step_bits=t,state_bits=s,occupation=3,distance_cutoff=209716),bands=q3.BANDS,boxes=rows,
        Q3_upper=base.encode(total),Q123_upper=base.encode(partial),
        margin_bits_diagnostic=math.log2(total.denominator)-math.log2(total.numerator),
        Q3_below_2_to_minus_126=total<F(1,1<<126),Q123_below_2_to_minus_49=partial<F(1,1<<49),
        all_occupations_certified=False,conditioning_rounding='192-bit relative dyadic upper bounds',
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer})
    base.write_new(OUTPUT,payload)
    print('OUTWARD Q3 margin',payload['margin_bits_diagnostic'],'Q123 < 2^-49',payload['Q123_below_2_to_minus_49'],flush=True)


if __name__=='__main__':
    build()
