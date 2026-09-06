"""Independent Arb transfer recurrence and evidence-aware full Q1 audit.

Uses unnormalized support-count recurrences and exact binomial denominators,
including weight 256. Does not certify Q>=2 or replace RandomStepConv by RM2Sub.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True
from flint import arb, ctx
from run_higher_endpoint_preflight import sha, write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
TRANSFER=GEN/'bch256_q1_full_arb_transfer.json'
AUDIT=GEN/'bch256_q1_full_arb_evidence_audit.json'
ctx.prec=192


def identity(scalar):
    return (scalar(1),scalar(0),scalar(0),scalar(1))


def matrix_mul(a,b):
    return (a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3])


def matrix_add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def row_mul(a,b):
    return (a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3])


def regions(zero,active,length,scalar):
    # Unnormalized sum of all products with exactly one active position.
    z=identity(scalar)
    a=(scalar(0),)*4
    for _ in range(length):
        a=matrix_add(matrix_mul(a,zero),matrix_mul(z,active))
        z=matrix_mul(z,zero)
    return z,tuple(x/length for x in a)


def support_moments(zero,active,positions,scalar):
    # Count every support before dividing by its exact binomial multiplicity.
    rows=[(scalar(1),scalar(0))]
    for n in range(positions):
        updated=[]
        for w in range(n+2):
            a=row_mul(rows[w],zero) if w<=n else (scalar(0),scalar(0))
            b=row_mul(rows[w-1],active) if w else (scalar(0),scalar(0))
            updated.append((a[0]+b[0],a[1]+b[1]))
        rows=updated
    return [(a+b)/math.comb(positions,w) for w,(a,b) in enumerate(rows)]


def rational_toy_check():
    q,b=Fraction(1,16),Fraction(9,10)
    zero=(Fraction(1),Fraction(0),q*b,(1-q)*b)
    active=(q*b,(1-q)*b,q*b,(1-q)*b)
    z,a=regions(zero,active,5,Fraction)
    exact=(Fraction(0),)*4
    for selected in range(5):
        product=identity(Fraction)
        for j in range(5):
            product=matrix_mul(product,active if j==selected else zero)
        exact=matrix_add(exact,product)
    assert a==tuple(x/5 for x in exact)
    actual=support_moments(z,a,6,Fraction)
    for weight in range(7):
        total=Fraction(0)
        for support in itertools.combinations(range(6),weight):
            row=(Fraction(1),Fraction(0))
            for j in range(6):
                row=row_mul(row,a if j in support else z)
            total+=sum(row)
        assert actual[weight]==total/math.comb(6,weight)
    return dict(exact_rational_exhaustive_check_passed=True,region_length=5,
                outer_positions=6,all_supports_checked=64)


def rational(value):
    mantissa,exponent=value.man_exp()
    mantissa,exponent=int(mantissa),int(exponent)
    return Fraction(mantissa<<exponent) if exponent>=0 else Fraction(mantissa,1<<-exponent)


def encode(value):
    return dict(numerator=str(value.numerator),denominator=str(value.denominator),
                log2_diagnostic=math.log2(value.numerator)-math.log2(value.denominator) if value>0 else None)


def decode(value):
    return Fraction(int(value['numerator']),int(value['denominator']))


def transfer_inputs():
    names=['random_inner_oa15_majorant_diagnostic.json','random_inner_threshold_outward.json']
    return {str((GEN/name).relative_to(ROOT)):sha(GEN/name) for name in names} | {
        str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}


def build_transfer():
    assert not TRANSFER.exists(), 'Refusing to overwrite transfer receipt'
    checks=rational_toy_check()
    threshold=json.loads((GEN/'random_inner_threshold_outward.json').read_text())
    diagnostic=json.loads((GEN/'random_inner_oa15_majorant_diagnostic.json').read_text())
    params=threshold['parameters']
    assert (params['outer_rows'],params['output_bits'],params['distance'],params['memory_bits'])==(8192,2097152,209716,22)
    assert threshold['source_sha256']['diagnostic']==sha(GEN/'random_inner_oa15_majorant_diagnostic.json')
    groups={}
    for w in range(38,220,2):
        log_tilt=float(diagnostic['randomstepconv']['best_tilt_by_weight'][w])
        s=Fraction.from_float(math.exp(log_tilt))
        groups.setdefault(s,[]).append(w)
    # A separate explicit rational witness suffices for the unique all-one word.
    groups[Fraction(1,128)]=[256]
    rows={}
    for s,weights in sorted(groups.items()):
        t=arb(s.numerator)/s.denominator
        z=(-t).exp()
        bit=(1+z)/2
        collision=arb(1)/(1<<22)
        terminate=collision*bit
        survive=(1-collision)*bit
        zero=(arb(1),arb(0),terminate,survive)
        active=(terminate,survive,terminate,survive)
        rz,ra=regions(zero,active,8192,arb)
        moments=support_moments(rz,ra,256,arb)
        chernoff=(209716*t).exp()
        for w in weights:
            value=8192*moments[w]*chernoff
            upper=min(Fraction(8192),rational(value.upper()))
            assert upper>0
            row=dict(weight=w,surprisal=encode(s),coefficient_upper=encode(upper))
            if w!=256:
                # Recover the earlier coefficient using Arb, without binary64 underflow.
                old_log=float(threshold['coefficient_log2_upper'][str(w)])
                old_ratio=Fraction.from_float(old_log)
                old=arb(old_ratio.numerator)/old_ratio.denominator
                old_lower=rational((old*arb(2).log()).exp().lower())
                row['below_previous_directed_log_coefficient']=upper<=old_lower
                assert row['below_previous_directed_log_coefficient'], f'Unexpected directed discrepancy at weight {w}'
            rows[str(w)]=row
        print('ARB transfer',str(s),weights,flush=True)
    assert decode(rows['256']['coefficient_upper'])<Fraction(1,1<<128)
    payload=dict(classification='Full Arb interval transfer computation; exact support counts and binomial division',
        precision_bits=ctx.prec,parameters=params,toy_check=checks,
        initial_state='inactive, row vector (1,0)',
        shell_coverage='Every possible nonzero shell: even weights 38 through 218 and the all-one shell 256',
        coefficient_meaning='8192 times a Chernoff upper bound for one fixed codeword after independent coordinate/region permutations and RandomStepConv setup',
        coefficient_rows=rows,source_sha256=transfer_inputs(),
        limitation='This certifies numerical transfer bounds, not multiplicity caps or Q>=2')
    write_new(TRANSFER,payload)
    print('TRANSFER COMPLETE',TRANSFER,flush=True)


def audit():
    assert not AUDIT.exists(), 'Refusing to overwrite evidence audit'
    from run_endpoint_certificate import verify as verify38
    from run_higher_endpoint_certificate import verify as verify_higher
    old_dir=GEN/'endpoint_certificate_20260904'
    new_dir=GEN/'higher_endpoint_certificate_20260904'
    old=verify38(old_dir)
    new=verify_higher(new_dir)
    for directory,value in ((old_dir,old),(new_dir,new)):
        assert value==json.loads((directory/'certificate.json').read_text())
    transfer=json.loads(TRANSFER.read_text())
    assert transfer['source_sha256']==transfer_inputs()
    lp=json.loads((GEN/'exact_bch_lp_caps.json').read_text())
    caps={row['weight']:row['exact_Aw_C_cap'] for row in lp['shells']}
    sources={w:'exact BCH-sandwich LP cap' for w in caps}
    accepted={}
    for weight,result in [(38,old)]+[(int(w),r) for w,r in new['results'].items()]:
        if result['status']=='conditional_statistical_acceptance':
            cap=result['a38_cap'] if weight==38 else result['target_cap']
            caps[weight]=min(caps[weight],cap)
            sources[weight]='conditionally accepted fixed-budget endpoint test'
            accepted[str(weight)]=cap
    for w in range(52,130,2):
        caps[w]=min(1<<128,math.comb(256,w-18)//math.comb(w,w-18))
        sources[w]='constant-weight packing and total mass'
        if w in (52,54):
            johnson=json.loads((GEN/f'johnson_n256_w{w}_d38.json').read_text())
            caps[w]=min(caps[w],johnson['integer_shell_upper'])
            sources[w]='Johnson integer shell cap'
    coefficients={int(w):decode(row['coefficient_upper']) for w,row in transfer['coefficient_rows'].items()}
    total=coefficients[256]
    terms=[dict(weights=[256],cap=1,source='at most one all-one word',upper=encode(total))]
    for w,cap in sorted(caps.items()):
        pair=coefficients[w]+(coefficients[256-w] if w!=128 else 0)
        contribution=cap*pair
        total+=contribution
        terms.append(dict(weights=[w,256-w] if w!=128 else [128],cap=cap,source=sources[w],upper=encode(contribution)))
    assert {x for term in terms for x in term['weights']}==set(coefficients)
    target=Fraction(1,1<<40)
    residual=target-total
    q2=json.loads((GEN/'bch256_q2_envelope_diagnostic.json').read_text())
    receipt=dict(classification='Exact rational full-Q1 aggregation of independently recomputed Arb bounds; shell evidence is conditional statistical acceptance',
        accepted_shell_caps=accepted,full_Q1_upper=encode(total),
        Q1_margin_bits_diagnostic=-encode(total)['log2_diagnostic'],
        exact_Q1_comparison_passes=total<target,
        residual_for_Q_ge_2=encode(residual),residual_fraction_of_target=float(residual/target),
        all_one_included=True,all_one_upper=encode(coefficients[256]),terms=terms,
        statistical_familywise_error_under_ideal_IID='at most 2^-39 for the three predeclared shell tests',
        randomness_qualification=new['randomness_qualification'],
        probability_space='Fixed repeated outer C; independent row-coordinate permutations and region permutations; RandomStepConv-M22 setup',
        parameters=transfer['parameters'],
        Q2_status=dict(classification='binary64 diagnostic only; no outward certificate',
            existing_LP_envelope_margin_bits=q2['results']['LP_only']['q2_margin_bits']),
        remaining_obligations=['Certify and sum all Q>=2 contributions within the exact residual',
            'Q2 numerical pass is not yet outward-certified; Q>=3 remains uncontrolled here',
            'RM2Sub requires its own transfer analysis; this audit uses RandomStepConv-M22',
            'A deterministic BCH spectrum theorem is not established by the sampling tests'],
        full_SPIN_40_bit_theorem_established=False,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (TRANSFER,GEN/'exact_bch_lp_caps.json',GEN/'johnson_n256_w52_d38.json',
             GEN/'johnson_n256_w54_d38.json',old_dir/'certificate.json',new_dir/'certificate.json',
             GEN/'bch256_q2_envelope_diagnostic.json',Path(__file__))})
    write_new(AUDIT,receipt)
    print(json.dumps({k:receipt[k] for k in ('accepted_shell_caps','Q1_margin_bits_diagnostic',
        'exact_Q1_comparison_passes','residual_fraction_of_target','all_one_upper',
        'full_SPIN_40_bit_theorem_established')},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('transfer','audit'))
    args=parser.parse_args()
    build_transfer() if args.operation=='transfer' else audit()
