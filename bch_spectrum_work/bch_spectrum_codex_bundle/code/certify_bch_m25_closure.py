"""Deterministic-spectrum M25 variant; does not close the original M22 target.

Primary: unnormalized full-Arb support-count recurrence.
Independent check: normalized full-Arb recurrence and binary polynomial powers.
Higher occupations use the existing M22 bounds by marginal monotonicity, with
an explicit allowance for removing all empirical orbit-search lower bounds.
"""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
import numpy as np
from flint import arb,ctx
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
from audit_bch_q1_full_arb import regions,support_moments,rational,encode,decode
from screen_bch_random_inner_closure import region,base,uniform_coefficients
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
OUTPUT=GEN/'bch256_m25_unconditional_distance.json'


def deterministic_caps():
    receipt=json.loads((GEN/'bch256_closure_deterministic_envelope.json').read_text())
    result={0:0,256:1}
    low={row['weight']:row['cap'] for row in receipt['shells']}
    for w in range(38,130,2):
        value=min(1<<128,math.comb(256,w-18)//math.comb(w,w-18))
        if w in low:
            value=min(value,low[w])
        if w in (52,54):
            value=min(value,json.loads((GEN/f'johnson_n256_w{w}_d38.json').read_text())['integer_shell_upper'])
        result[w]=result[256-w]=value
    return result


def source_paths():
    return [Path(__file__),ROOT/'code/audit_bch_q1_full_arb.py',
        ROOT/'code/screen_bch_random_inner_closure.py',
        GEN/'bch256_random_inner_closure_screen.json',
        GEN/'bch256_closure_deterministic_envelope.json',
        GEN/'johnson_n256_w52_d38.json',GEN/'johnson_n256_w54_d38.json',
        GEN/'bch256_q2_positive_outward.json',GEN/'bch256_q3_capped_outward.json',
        GEN/'bch256_occupation_tail_outward.json',GEN/'bch256_low_occupation_crosschecks.json',
        GEN/'bch256_tail_transfer_crosschecks.json',GEN/'bch256_tail_envelope_checks.json']


def old_higher_occupations():
    q2=json.loads((GEN/'bch256_q2_positive_outward.json').read_text())
    q3=json.loads((GEN/'bch256_q3_capped_outward.json').read_text())
    tail=json.loads((GEN/'bch256_occupation_tail_outward.json').read_text())
    # These receipts used no statistical upper caps. Their original shell
    # envelopes are comparison measures, enlarged by r^q in the new proof.
    def upper(value):
        for key in ('Q2_upper','Q3_upper','q2_upper','q3_upper','total_upper'):
            if key in value:
                return decode(value[key])
        raise KeyError(list(value))
    return upper(q2),upper(q3),decode(tail['tail_upper'])


def build():
    assert not OUTPUT.exists()
    caps=deterministic_caps()
    saved=json.loads((GEN/'bch256_random_inner_closure_screen.json').read_text())
    candidates=sorted({r['s'] for c in saved['cases'] for r in c['coefficient_rows']})
    best=np.full(257,np.inf)
    witnesses=np.zeros(257)
    for s in candidates:
        z,o=base.step_matrices(math.exp(-s),25)
        rz,ra=region(z,o,8192)
        a=uniform_coefficients(rz,ra,256,256)
        values=np.logaddexp(a[:,0,0],a[:,0,1])+209716*s
        improve=values<best
        best=np.minimum(best,values)
        witnesses[improve]=s
    groups={}
    for w,c in caps.items():
        if c:
            groups.setdefault(Fraction.from_float(float(witnesses[w])),[]).append(w)
    ctx.prec=192
    coefficients={}
    for s,weights in sorted(groups.items()):
        t=arb(s.numerator)/s.denominator
        b=(1+(-t).exp())/2
        epsilon=arb(1)/(1<<25)
        d,f=epsilon*b,(1-epsilon)*b
        z,o=(arb(1),arb(0),d,f),(d,f,d,f)
        rz,ra=regions(z,o,8192,arb)
        moments=support_moments(rz,ra,256,arb)
        correction=(209716*t).exp()*8192
        for w in weights:
            v=correction*moments[w]
            upper=min(Fraction(8192),rational(v.upper())*Fraction((1<<64)+1,1<<64))
            coefficients[str(w)]=dict(weight=w,s=encode(s),coefficient_upper=encode(upper),cap=caps[w])
        print('CERTIFIED M25 SHELLS',sorted(weights),flush=True)
    q1=sum((decode(row['coefficient_upper'])*row['cap'] for row in coefficients.values()),Fraction(0))
    q2,q3,tail=old_higher_occupations()
    # For r=1+2^-20 and q<=8192:
    # r^q <= sum_{j>=0}(8192/2^20)^j = 128/127.
    # This is an analytic bound, not a floating evaluation of r^8192.
    inflation=Fraction(128,127)
    total=q1+inflation*(q2+q3+tail)
    assert total<Fraction(1,1<<40)
    result=dict(classification='All-occupation ideal RandomStepConv-M25 distance certificate using deterministic BCH bounds',
        parameters=dict(message_bits=1<<20,outer_rows=8192,outer_length=256,
                        output_bits=1<<21,distance_cutoff=209716,memory_bits=25),
        coefficient_rows=coefficients,Q1_upper=encode(q1),
        old_M22_Q2_upper=encode(q2),old_M22_Q3_upper=encode(q3),old_M22_tail_upper=encode(tail),
        higher_occupation_inflation=encode(inflation),
        full_first_moment_upper=encode(total),
        full_margin_bits_diagnostic=-encode(total)['log2_diagnostic'],
        statistical_shell_assumptions_used=False,empirical_orbit_lower_bounds_used=False,
        exact_full_bound_below_2_to_minus_40=True,all_occupations_covered=True,
        original_M22_target_closed=False,paper_random_convolution_transfer_established=False,
        RM2Sub_transfer_established=False,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in source_paths()})
    write_new(OUTPUT,result)
    print(json.dumps({k:result[k] for k in ('full_margin_bits_diagnostic','exact_full_bound_below_2_to_minus_40','original_M22_target_closed')},indent=2))


def normalized_check(s):
    # Independent matrix representation and multiplication; no primary helpers.
    t=arb(s.numerator)/s.denominator
    bit=(1+(-t).exp())/2
    d=bit/(1<<25)
    f=bit-d
    zero=[[arb(1),arb(0)],[d,f]]
    one=[[d,f],[d,f]]
    def mul(a,b):
        return [[sum((a[i][k]*b[k][j] for k in range(2)),arb(0)) for j in range(2)] for i in range(2)]
    def add(a,b):
        return [[a[i][j]+b[i][j] for j in range(2)] for i in range(2)]
    ident=[[arb(1),arb(0)],[arb(0),arb(1)]]
    empty=[[arb(0),arb(0)],[arb(0),arb(0)]]
    def polymul(a,b):
        return [mul(a[0],b[0]),add(mul(a[0],b[1]),mul(a[1],b[0]))]
    power=[zero,one]
    result=[ident,empty]
    exponent=8192
    while exponent:
        if exponent&1:
            result=polymul(result,power)
        exponent>>=1
        if exponent:
            power=polymul(power,power)
    rz=result[0]
    ra=[[v/8192 for v in row] for row in result[1]]
    rows=[[arb(1),arb(0)]]
    for length in range(1,257):
        updated=[]
        for w in range(length+1):
            item=[arb(0),arb(0)]
            if w<length:
                item=[sum((rows[w][k]*rz[k][j] for k in range(2)),arb(0))*(length-w)/length for j in range(2)]
            if w:
                item=[item[j]+sum((rows[w-1][k]*ra[k][j] for k in range(2)),arb(0))*w/length for j in range(2)]
            updated.append(item)
        rows=updated
    return [sum(row)*(209716*t).exp()*8192 for row in rows]


def verify():
    value=json.loads(OUTPUT.read_text())
    assert value['source_sha256']=={str(p.relative_to(ROOT)):sha(p) for p in source_paths()}
    caps=deterministic_caps()
    groups={}
    for name,row in value['coefficient_rows'].items():
        assert int(name)==row['weight'] and row['cap']==caps[int(name)]
        groups.setdefault(decode(row['s']),[]).append(row)
    assert set(map(int,value['coefficient_rows']))=={w for w,c in caps.items() if c}
    ctx.prec=256
    checked=0
    for s,rows in groups.items():
        moments=normalized_check(s)
        for row in rows:
            assert rational(moments[row['weight']].upper())<=decode(row['coefficient_upper'])
            checked+=1
    q1=sum((decode(row['coefficient_upper'])*row['cap'] for row in value['coefficient_rows'].values()),Fraction(0))
    q2,q3,tail=old_higher_occupations()
    assert q1==decode(value['Q1_upper'])
    assert (q2,q3,tail)==tuple(decode(value[k]) for k in ('old_M22_Q2_upper','old_M22_Q3_upper','old_M22_tail_upper'))
    assert decode(value['higher_occupation_inflation'])==Fraction(128,127)
    total=q1+Fraction(128,127)*(q2+q3+tail)
    assert total==decode(value['full_first_moment_upper'])<Fraction(1,1<<40)
    output=GEN/'bch256_m25_independent_transfer_check.json'
    result=dict(classification='Independent 256-bit Arb normalized recurrence, region polynomial powers, and exact aggregate',
        all_92_shells_independently_checked=checked==92,shells_checked=checked,
        exact_full_bound_below_2_to_minus_40=True,
        source_sha256={str(OUTPUT.relative_to(ROOT)):sha(OUTPUT),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))})
    if output.exists():
        assert json.loads(output.read_text())==result
    else:
        write_new(output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    verify() if '--verify' in sys.argv else build()
