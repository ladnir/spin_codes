"""Outward Q2 certificate using exact-count positive recurrences.

Region entries use Arb; the support-pair recurrence uses upward binary64
arithmetic with positive-underflow repair; final summation uses Fractions.
The existing LP-only caps suffice if the resulting exact comparison passes.
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
from certify_random_inner_threshold_outward import mul_up,add_up,upper_float
from evaluate_bch_q2_envelopes import envelopes
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
OUTPUT=GEN/'bch256_q2_positive_outward.json'
ctx.prec=192


def region_matrices(zero,active,length,scalar):
    current=[identity(scalar),(scalar(0),)*4,(scalar(0),)*4]
    for _ in range(length):
        current=[matrix_mul(current[0],zero),
                 matrix_add(matrix_mul(current[1],zero),matrix_mul(current[0],active)),
                 matrix_add(matrix_mul(current[2],zero),matrix_mul(current[1],active))]
    return [tuple(x/math.comb(length,j) for x in current[j]) for j in range(3)]


def row_product_up(rows,matrix):
    result=np.empty_like(rows)
    result[...,0]=add_up(mul_up(rows[...,0],matrix[0]),mul_up(rows[...,1],matrix[2]))
    result[...,1]=add_up(mul_up(rows[...,0],matrix[1]),mul_up(rows[...,1],matrix[3]))
    return result


def pair_counts_upper(regions,positions):
    # No conditioning weights enter this recurrence. Each of the four choices
    # appends one outer coordinate to an ordered pair of support sets.
    current=np.zeros((positions+1,positions+1,2),dtype=np.float64)
    scratch=np.zeros_like(current)
    current[0,0,0]=1
    for completed in range(positions):
        size=completed+1
        old=current[:size,:size]
        updated=scratch[:size+1,:size+1]
        updated.fill(0)
        updated[:size,:size]=row_product_up(old,regions[0])
        mixed=row_product_up(old,regions[1])
        updated[1:,:size]=add_up(updated[1:,:size],mixed)
        updated[:size,1:]=add_up(updated[:size,1:],mixed)
        both=row_product_up(old,regions[2])
        updated[1:,1:]=add_up(updated[1:,1:],both)
        current,scratch=scratch,current
    assert np.isfinite(current).all() and (current>=0).all()
    reciprocal=np.array([upper_float(arb(1)/math.comb(positions,w)) for w in range(positions+1)])
    totals=add_up(current[...,0],current[...,1])
    return mul_up(mul_up(totals,reciprocal[:,None]),reciprocal[None,:])


def toy_check():
    q,b=Fraction(1,16),Fraction(9,10)
    zero=(Fraction(1),Fraction(0),q*b,(1-q)*b)
    active=(q*b,(1-q)*b,q*b,(1-q)*b)
    regions=region_matrices(zero,active,5,Fraction)
    for count in range(3):
        total=(Fraction(0),)*4
        for selected in itertools.combinations(range(5),count):
            product=identity(Fraction)
            for j in range(5):
                product=matrix_mul(product,active if j in selected else zero)
            total=matrix_add(total,product)
        assert regions[count]==tuple(x/math.comb(5,count) for x in total)
    rounded=np.array([[upper_float(arb(x.numerator)/x.denominator) for x in row] for row in regions])
    upper=pair_counts_upper(rounded,4)
    maximum_excess=Fraction(0)
    for a in range(5):
        for b in range(5):
            total=Fraction(0)
            for first in itertools.combinations(range(4),a):
                for second in itertools.combinations(range(4),b):
                    row=(Fraction(1),Fraction(0))
                    for j in range(4):
                        row=row_mul(row,regions[int(j in first)+int(j in second)])
                    total+=sum(row)
            exact=total/(math.comb(4,a)*math.comb(4,b))
            excess=Fraction.from_float(upper[a,b])-exact
            assert excess>=0
            maximum_excess=max(maximum_excess,excess)
    # Explicitly exercise the path that must not round a positive product to zero.
    tiny=np.array([np.nextafter(0.0,np.inf)])
    assert mul_up(tiny,np.array([0.5]))[0]>0
    return dict(all_25_pairs_checked_against_exact_rational_enumeration=True,
                maximum_absolute_outward_excess=float(maximum_excess),positive_underflow_check_passed=True)


def main():
    assert not OUTPUT.exists(), 'Refusing to overwrite a certificate'
    checks=toy_check()
    s=Fraction.from_float(math.exp(-7.5))
    tilt=arb(s.numerator)/s.denominator
    bit=(1+(-tilt).exp())/2
    q=arb(1)/(1<<22)
    zero=(arb(1),arb(0),q*bit,(1-q)*bit)
    active=(q*bit,(1-q)*bit,q*bit,(1-q)*bit)
    regions=region_matrices(zero,active,8192,arb)
    region_upper=np.array([[upper_float(x) for x in row] for row in regions])
    moments=pair_counts_upper(region_upper,256)
    correction=upper_float((209716*tilt).exp())
    conditional=np.minimum(1.0,mul_up(moments,correction))
    spectrum=envelopes()['LP_only']
    weights=sorted(w for w,cap in spectrum.items() if w and cap)
    total=Fraction(0)
    pair_terms=[]
    placements=math.comb(8192,2)
    for a in weights:
        for b in weights:
            term=placements*spectrum[a]*spectrum[b]*Fraction.from_float(conditional[a,b])
            total+=term
            pair_terms.append((term,a,b))
    assert total>0
    cache=OUTPUT.with_suffix('.npz')
    assert not cache.exists()
    np.savez_compressed(cache,conditional_pair_upper=conditional,region_matrices_upper=region_upper)
    receipt=dict(classification='Outward Q2 bound: Arb region enclosures, positive binary64 count recurrence, exact rational aggregation',
        parameters=dict(message_bits=1<<20,outer_rows=8192,output_bits=1<<21,
                        distance_cutoff=209716,memory_bits=22,outer_length=256,occupation=2),
        fixed_surprisal=encode(s),arithmetic_precision_bits=ctx.prec,toy_checks=checks,
        probability_space='Fixed repeated C; independent row-coordinate permutations, independent region permutations, and RandomStepConv-M22 setup',
        counting_law='Choose two row positions; for each ordered weight pair (a,b), multiply A_a*A_b by the conditional support-pair bound. Local messages may coincide.',
        envelope_source='Existing deterministic LP, Johnson and packing caps only; no statistical shell caps needed',
        spectrum_envelope=spectrum,Q2_upper=encode(total),Q2_margin_bits_diagnostic=-encode(total)['log2_diagnostic'],
        exact_below_2_to_minus_60=total<Fraction(1,1<<60),
        dominant_pairs=[dict(weights=[a,b],upper=encode(term)) for term,a,b in sorted(pair_terms,reverse=True)[:10]],
        arithmetic_argument='All recurrences sum nonnegative products, each rounded upward. Products that positively underflow are replaced by the least positive subnormal. Unnormalized support counts are at most 4^256 in exact arithmetic; the observed upward computation remains finite. Exact binomial reciprocals are enclosed by Arb. The final sum is rational.',
        limitation='Occupation two only; Q>=3 and the RM2Sub transfer remain separate',
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),
            ROOT/'code/audit_bch_q1_full_arb.py',ROOT/'code/certify_random_inner_threshold_outward.py',
            ROOT/'code/evaluate_bch_q2_envelopes.py',GEN/'exact_bch_lp_caps.json',
            GEN/'johnson_n256_w52_d38.json',GEN/'johnson_n256_w54_d38.json',cache)})
    write_new(OUTPUT,receipt)
    print(json.dumps({k:receipt[k] for k in ('classification','Q2_margin_bits_diagnostic','exact_below_2_to_minus_60','toy_checks')},indent=2),flush=True)


if __name__=='__main__':
    main()
