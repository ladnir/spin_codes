"""Outward Q3 bound with exact low weights and a weight-70 upper class.

Deleting input ones increases the RandomStepConv PGF in expectation, by its
two-state Markov coupling. This is not pointwise monotonicity of a fixed code.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode=True
import numpy as np
from flint import arb,ctx
from audit_bch_q1_full_arb import identity,matrix_mul,matrix_add,row_mul,encode
from certify_bch_q2_positive import row_product_up
from certify_random_inner_threshold_outward import mul_up,add_up,upper_float
from evaluate_bch_q2_envelopes import envelopes
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
OUTPUT=GEN/'bch256_q3_capped_outward.json'
ctx.prec=192


def region_matrices(zero,active,length,scalar):
    current=[identity(scalar)]+[(scalar(0),)*4 for _ in range(3)]
    for _ in range(length):
        current=[matrix_mul(current[0],zero)]+[
            matrix_add(matrix_mul(current[j],zero),matrix_mul(current[j-1],active)) for j in range(1,4)]
    return [tuple(x/math.comb(length,j) for x in current[j]) for j in range(4)]


def triple_upper(regions,positions,maximum_weight):
    side=maximum_weight+1
    current=np.zeros((side,side,side,2))
    scratch=np.zeros_like(current)
    current[0,0,0,0]=1
    patterns=[p for p in itertools.product((0,1),repeat=3)]
    for completed in range(positions):
        size=min(completed+1,side)
        outsize=min(completed+2,side)
        old=current[:size,:size,:size]
        updated=scratch[:outsize,:outsize,:outsize]
        updated.fill(0)
        for selected in range(4):
            term=row_product_up(old,regions[selected])
            for bits in patterns:
                if sum(bits)!=selected:
                    continue
                lengths=[min(size,outsize-b) for b in bits]
                src=tuple(slice(0,n) for n in lengths)
                dst=tuple(slice(b,b+n) for b,n in zip(bits,lengths))
                updated[dst]=add_up(updated[dst],term[src])
        current,scratch=scratch,current
    assert np.isfinite(current).all() and (current>=0).all()
    inverse=np.array([upper_float(arb(1)/math.comb(positions,w)) for w in range(side)])
    totals=add_up(current[...,0],current[...,1])
    totals=mul_up(totals,inverse[:,None,None])
    totals=mul_up(totals,inverse[None,:,None])
    return mul_up(totals,inverse[None,None,:])


def toy_check():
    q,b=Fraction(1,16),Fraction(9,10)
    zero=(Fraction(1),Fraction(0),q*b,(1-q)*b)
    active=(q*b,(1-q)*b,q*b,(1-q)*b)
    regions=region_matrices(zero,active,4,Fraction)
    for count in range(4):
        total=(Fraction(0),)*4
        for selected in itertools.combinations(range(4),count):
            product=identity(Fraction)
            for j in range(4):
                product=matrix_mul(product,active if j in selected else zero)
            total=matrix_add(total,product)
        assert regions[count]==tuple(x/math.comb(4,count) for x in total)
    rounded=np.array([[upper_float(arb(x.numerator)/x.denominator) for x in row] for row in regions])
    full=triple_upper(rounded,3,3)
    truncated=triple_upper(rounded,3,2)
    assert np.array_equal(full[:3,:3,:3],truncated)
    for weights in itertools.product(range(4),repeat=3):
        total=Fraction(0)
        families=[list(itertools.combinations(range(3),w)) for w in weights]
        for supports in itertools.product(*families):
            row=(Fraction(1),Fraction(0))
            for j in range(3):
                row=row_mul(row,regions[sum(j in s for s in supports)])
            total+=sum(row)
        exact=total/math.prod(math.comb(3,w) for w in weights)
        assert Fraction.from_float(full[weights])>=exact
    return dict(all_64_weight_triples_checked_against_exact_rationals=True,
                truncation_matches_full_recurrence=True)


def main():
    assert not OUTPUT.exists()
    checks=toy_check()
    spectrum=envelopes()['LP_only']
    caps={w:spectrum[w] for w in range(38,70,2)}
    caps[70]=(1<<128)-1
    weights=sorted(caps)
    best=np.ones((71,71,71))
    witnesses=np.zeros_like(best)
    for u in (-7.25,-7.0):
        s=Fraction.from_float(math.exp(u))
        tilt=arb(s.numerator)/s.denominator
        bit=(1+(-tilt).exp())/2
        q=arb(1)/(1<<22)
        zero=(arb(1),arb(0),q*bit,(1-q)*bit)
        active=(q*bit,(1-q)*bit,q*bit,(1-q)*bit)
        regions=region_matrices(zero,active,8192,arb)
        upper=np.array([[upper_float(x) for x in row] for row in regions])
        moment=triple_upper(upper,256,70)
        tail=np.minimum(1.0,mul_up(moment,upper_float((209716*tilt).exp())))
        better=tail<best
        witnesses[better]=u
        best=np.minimum(best,tail)
        print('Q3 outward tilt complete',u,flush=True)
    total=Fraction(0)
    terms=[]
    placements=math.comb(8192,3)
    for a,b,c in itertools.product(weights,repeat=3):
        term=placements*caps[a]*caps[b]*caps[c]*Fraction.from_float(best[a,b,c])
        total+=term
        terms.append((term,a,b,c))
    cache=OUTPUT.with_suffix('.npz')
    assert not cache.exists()
    np.savez_compressed(cache,conditional_upper=best,log_tilt_witness=witnesses)
    receipt=dict(classification='Outward Q3 bound via positive support-count recurrence and monotone weight truncation',
        parameters=dict(message_bits=1<<20,outer_rows=8192,output_bits=1<<21,
            distance_cutoff=209716,memory_bits=22,occupation=3,outer_length=256),
        toy_checks=checks,weight_class_caps=caps,
        upper_class='The class at 70 represents every actual weight >=70, including 256, with total mass bounded by 2^128-1',
        monotonicity_scope='Expected RandomStepConv PGF only. The two-state chain admits an order-preserving coin coupling when input ones are deleted. Uniform row supports couple by a shared random ordering; region permutations preserve input-set inclusion.',
        low_cap_sources='Existing deterministic LP, Johnson, packing and total-mass caps; no statistical shell acceptance needed',
        counting_law='Choose three row positions, then sum over ordered triples of weight classes. Independent row-coordinate permutations suffice even if local messages coincide.',
        Q3_upper=encode(total),Q3_margin_bits_diagnostic=-encode(total)['log2_diagnostic'],
        exact_below_2_to_minus_60=total<Fraction(1,1<<60),
        dominant_triples=[dict(weights=[a,b,c],upper=encode(term)) for term,a,b,c in sorted(terms,reverse=True)[:10]],
        arithmetic_argument='Unnormalized support counts are bounded in exact arithmetic by 8^256=2^768. The upward binary64 recurrence remains finite and repairs positive underflow. Exact binomial reciprocals and Chernoff factors are enclosed by Arb. Final aggregation uses Fractions.',
        limitations=['Only Q=3, not Q>=4','RandomStepConv-M22, not RM2Sub'],
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'code/certify_bch_q2_positive.py',
            ROOT/'code/audit_bch_q1_full_arb.py',ROOT/'code/certify_random_inner_threshold_outward.py',
            ROOT/'code/evaluate_bch_q2_envelopes.py',GEN/'exact_bch_lp_caps.json',
            GEN/'johnson_n256_w52_d38.json',GEN/'johnson_n256_w54_d38.json',cache)})
    write_new(OUTPUT,receipt)
    print(json.dumps({k:receipt[k] for k in ('classification','Q3_margin_bits_diagnostic','exact_below_2_to_minus_60','dominant_triples')},indent=2),flush=True)


if __name__=='__main__':
    main()
