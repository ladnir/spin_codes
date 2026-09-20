"""Independent degeneracy and combinatorial cross-checks for Q2/Q3 receipts."""
from __future__ import annotations

import itertools
import json
import math
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode=True
import numpy as np
from flint import arb,ctx
from audit_bch_q1_full_arb import regions,support_moments,rational,decode,encode
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
ctx.prec=256


def q1_conditionals(u):
    exact=Fraction.from_float(math.exp(u))
    s=arb(exact.numerator)/exact.denominator
    b=(1+(-s).exp())/2
    q=arb(1)/(1<<22)
    zero=(arb(1),arb(0),q*b,(1-q)*b)
    active=(q*b,(1-q)*b,q*b,(1-q)*b)
    rz,ra=regions(zero,active,8192,arb)
    moment=support_moments(rz,ra,256,arb)
    correction=(209716*s).exp()
    return [min(Fraction(1),rational((x*correction).upper())) for x in moment]


def checked(path):
    receipt=json.loads(path.read_text())
    for name,expected in receipt['source_sha256'].items():
        assert sha(ROOT/name)==expected
    return receipt


def main():
    output=GEN/'bch256_low_occupation_crosschecks.json'
    assert not output.exists()
    q2=checked(GEN/'bch256_q2_positive_outward.json')
    q3=checked(GEN/'bch256_q3_capped_outward.json')
    with np.load(GEN/'bch256_q2_positive_outward.npz') as data:
        pair=data['conditional_pair_upper'].copy()
    with np.load(GEN/'bch256_q3_capped_outward.npz') as data:
        triple=data['conditional_upper'].copy()
    single=q1_conditionals(-7.5)
    pair_axis_checks=0
    for w in range(257):
        for value in (pair[w,0],pair[0,w]):
            assert Fraction.from_float(value)>=single[w]
            pair_axis_checks+=1
    first,second=q1_conditionals(-7.25),q1_conditionals(-7.0)
    triple_axis_checks=0
    for w in range(71):
        target=min(first[w],second[w])
        for value in (triple[w,0,0],triple[0,w,0],triple[0,0,w]):
            assert Fraction.from_float(value)>=target
            triple_axis_checks+=1
    assert pair[0,0]==triple[0,0,0]==1.0
    pair_caps={int(w):cap for w,cap in q2['spectrum_envelope'].items() if int(w) and cap}
    pair_sum=Fraction(0)
    for a,b in itertools.combinations_with_replacement(sorted(pair_caps),2):
        value=Fraction.from_float(max(pair[a,b],pair[b,a]))
        multiplicity=1 if a==b else 2
        pair_sum+=multiplicity*math.comb(8192,2)*pair_caps[a]*pair_caps[b]*value
    assert decode(q2['Q2_upper'])<=pair_sum<decode(q2['Q2_upper'])*Fraction(1000001,1000000)
    triple_caps={int(w):cap for w,cap in q3['weight_class_caps'].items()}
    triple_sum=Fraction(0)
    for weights in itertools.combinations_with_replacement(sorted(triple_caps),3):
        permutations=set(itertools.permutations(weights))
        multiplicity=math.factorial(3)//math.prod(math.factorial(n) for n in Counter(weights).values())
        assert multiplicity==len(permutations)
        value=Fraction.from_float(max(triple[p] for p in permutations))
        triple_sum+=multiplicity*math.comb(8192,3)*math.prod(triple_caps[w] for w in weights)*value
    assert decode(q3['Q3_upper'])<=triple_sum<decode(q3['Q3_upper'])*Fraction(1000001,1000000)
    assert pair_sum+triple_sum<Fraction(1,1<<60)
    result=dict(classification='Independent degeneracy checks and symmetrized combinatorial reaggregation',
        Arb_precision_bits=ctx.prec,Q2_zero_weight_axis_checks=pair_axis_checks,
        Q3_two_zero_weight_axis_checks=triple_axis_checks,
        all_zero_input_probability_one=True,
        unordered_pair_aggregate_upper=encode(pair_sum),unordered_triple_aggregate_upper=encode(triple_sum),
        both_symmetrized_aggregates_below_2_to_minus_60=True,
        symmetry_rounding_inflation_below_one_part_per_million=True,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),
            ROOT/'code/audit_bch_q1_full_arb.py',GEN/'bch256_q2_positive_outward.json',
            GEN/'bch256_q3_capped_outward.json',GEN/'bch256_q2_positive_outward.npz',
            GEN/'bch256_q3_capped_outward.npz')})
    write_new(output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
