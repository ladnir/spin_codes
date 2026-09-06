"""Outward certification of every occupation q=4,...,8192.

Low occupations use positive matrix recurrences. Higher occupations use
Arb-enclosed Perron-vector and coefficient bounds at saved dyadic witnesses.
"""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import numpy as np
from flint import arb,ctx
from bch_tail_reference import ROOT,GEN,L,B,D,contiguous
from audit_bch_q1_full_arb import matrix_mul,identity,rational,encode,decode
from certify_random_inner_threshold_outward import (upper_float,uniform_coefficients_upper,arb_of_float)
from run_higher_endpoint_preflight import sha,write_new
ctx.prec=192


def as_arb(x):
    if isinstance(x,Fraction):
        return arb(x.numerator)/x.denominator
    return arb(x)


def power(matrix,n):
    result=identity(arb)
    while n:
        if n&1:
            result=matrix_mul(result,matrix)
        n//=2
        if n:
            matrix=matrix_mul(matrix,matrix)
    return result


def log_record(bound):
    upper=rational(bound.upper())
    exponent=-(-upper.numerator//upper.denominator)
    return dict(log2_upper=encode(upper),dyadic_upper_exponent=exponent)


def perron_bound(q,factor,rho,s,p):
    bit=(1+(-s).exp())/2
    epsilon=arb(1)/(1<<22)
    eta=rho*p
    a=1-eta+eta*epsilon*bit
    b=(1-epsilon)*bit
    c=eta*b
    d=epsilon*bit
    eigen=(a+b+((a-b)**2+4*c*d).sqrt())/2
    gap=eigen-b
    assert gap.lower()>0
    prefactor=gap/d
    # The positive eigenvector is (1,d/(lambda-b)). Taking max(1,1/v_1)
    # bounds the terminal all-ones vector even at a rounding boundary.
    prefactor=as_arb(max(Fraction(1),rational(prefactor.upper())))
    choose=arb(math.comb(L,q)).log()
    point=choose+q*p.log()
    if q<L:
        point+=(L-q)*(1-p).log()
    full=(choose+q*arb(factor).log()+prefactor.log()+L*B*eigen.log()+D*s-B*point)/arb(2).log()
    return full


def toy_check():
    # Exact Perron witness inequality, with an independent finite power check.
    checks=0
    for rho in (Fraction(1,3),Fraction(1,2)):
        for p in (Fraction(1,5),Fraction(2,3),Fraction(1)):
            bit=arb(9)/10
            epsilon=arb(1)/16
            eta=as_arb(rho*p)
            a=1-eta+eta*epsilon*bit
            b=(1-epsilon)*bit
            c=eta*b
            d=epsilon*bit
            eigen=(a+b+((a-b)**2+4*c*d).sqrt())/2
            prefactor=max(Fraction(1),rational(((eigen-b)/d).upper()))
            actual=power((a,c,d,b),12)
            upper=as_arb(prefactor)*eigen**12
            assert (upper-(actual[0]+actual[1])).lower()>0
            checks+=1
    return dict(six_independent_finite_power_Perron_checks=checks)


def main():
    output=GEN/'bch256_occupation_tail_outward.json'
    assert not output.exists()
    moment_path=GEN/'bch256_tail_moment_envelopes.json'
    witness_path=GEN/'bch256_tail_reference_refined.json'
    moments=json.loads(moment_path.read_text())
    diagnostic=json.loads(witness_path.read_text())
    for document in (moments,diagnostic):
        for name,expected in document['source_sha256'].items():
            assert sha(ROOT/name)==expected
    factors={r['rho']:int(r['factor']) for r in moments['references']}
    checks=toy_check()
    maximum=86
    best={}
    for rho_text,tilts in [('1/2',[-7.25+j*.25 for j in range(12)]),('9/20',[-4.75,-4.5])]:
        rho=as_arb(Fraction(rho_text))
        factor=factors[rho_text]
        for u in tilts:
            s=arb_of_float(math.exp(u))
            bit=(1+(-s).exp())/2
            epsilon=arb(1)/(1<<22)
            end=epsilon*bit
            stay=(1-epsilon)*bit
            zero=np.array([[1.,0.],[upper_float(end),upper_float(stay)]])
            candidate=np.array([[upper_float(1-rho+rho*end),upper_float(rho*stay)],
                                [upper_float(end),upper_float(stay)]])
            mantissas,exponents=uniform_coefficients_upper(zero,candidate,L,maximum)
            for q in range(4,maximum+1):
                matrix=tuple(arb_of_float(float(x))*arb(2)**int(exponents[q]) for x in mantissas[q].flat)
                complete=power(matrix,B)
                moment=complete[0]+complete[1]
                log_bound=(arb(math.comb(L,q)).log()+q*arb(factor).log()+moment.log()+D*s)/arb(2).log()
                row=log_record(log_bound)
                row.update(q=q,method='outward_positive_region',rho=rho_text,
                           s_dyadic=encode(Fraction.from_float(math.exp(u))))
                if q not in best or decode(row['log2_upper'])<decode(best[q]['log2_upper']):
                    best[q]=row
            print('OUTWARD sparse',rho_text,u,flush=True)
    for diagnostic_row in diagnostic['rows']:
        q=diagnostic_row['q']
        if q<=maximum:
            continue
        w=diagnostic_row['witness']
        assert w['method']=='perron_conditioned'
        rho_text=w['rho']
        s,p=Fraction.from_float(w['s']),Fraction.from_float(w['p'])
        assert s>0 and 0<p<=1 and (p<1 or q==L)
        log_bound=perron_bound(q,factors[rho_text],as_arb(Fraction(rho_text)),as_arb(s),as_arb(p))
        row=log_record(log_bound)
        row.update(q=q,method='Arb_Perron_conditioned',rho=rho_text,s_dyadic=encode(s),p_dyadic=encode(p))
        best[q]=row
        if q%1024==0:
            print('OUTWARD dense',q,flush=True)
    assert set(best)==set(range(4,L+1))
    rows=[best[q] for q in range(4,L+1)]
    total=sum((Fraction(1<<r['dyadic_upper_exponent']) if r['dyadic_upper_exponent']>=0
               else Fraction(1,1<<-r['dyadic_upper_exponent']) for r in rows),Fraction(0))
    passing=total<=Fraction(1,1<<45)
    payload=dict(classification='Outward-certified occupation tail using exact moment envelopes and saved numerical witnesses',
        parameters=dict(outer_rows=L,outer_length=B,message_bits=1<<20,output_bits=L*B,distance_cutoff=D,memory_bits=22),
        rows=rows,tail_upper=encode(total),tail_margin_bits_diagnostic=-encode(total)['log2_diagnostic'],
        exact_tail_at_most_2_to_minus_45=passing,
        all_q_at_most_2_to_minus_58=all(r['dyadic_upper_exponent']<=-58 for r in rows),
        unresolved_occupations=contiguous(r['q'] for r in rows if r['dyadic_upper_exponent']>-58),
        toy_checks=checks,Arb_precision_bits=ctx.prec,
        scope='Fixed repeated BCH-derived C, independent coordinate/region permutations, RandomStepConv-M22; no RM2Sub transfer. Tail uses deterministic caps and OA15 only, not shell sampling.',
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),moment_path,witness_path,
            ROOT/'code/bch_tail_reference.py',ROOT/'code/audit_bch_q1_full_arb.py',
            ROOT/'code/certify_random_inner_threshold_outward.py')})
    write_new(output,payload)
    print(json.dumps({k:payload[k] for k in ('tail_margin_bits_diagnostic','exact_tail_at_most_2_to_minus_45',
        'all_q_at_most_2_to_minus_58','unresolved_occupations')},indent=2))


if __name__=='__main__':
    main()
