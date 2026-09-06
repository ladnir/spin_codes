"""Diagnostic probe: cumulative shell envelopes versus biased IID supports.

Exact rational prefix comparisons justify a scalar envelope for decreasing
functions of row weight. Transfer values below use nearest binary64 only.
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode=True
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
sys.path.insert(0,str(GEN/'bch256_q2_envelope_diagnostic_sources'))
import evaluate_ebch128_randomstepconv_g1 as base
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients
from evaluate_bch_q2_envelopes import envelopes
from run_higher_endpoint_preflight import sha,write_new


def prefix_envelope(spectrum,rho):
    a,b=rho.numerator,rho.denominator
    cdf=0
    mass=0
    candidates=[]
    for w in range(257):
        cdf+=math.comb(256,w)*a**w*(b-a)**(256-w)
        mass=min((1<<128)-1,mass+spectrum.get(w,0))
        candidates.append((Fraction(mass*b**256,cdf),w))
    ratio,weight=max(candidates)
    factor=-(-ratio.numerator//ratio.denominator)
    return factor,weight


def main():
    output=GEN/'cumulative_reference_q3_diagnostic.json'
    assert not output.exists()
    spectrum=envelopes()['LP_only']
    results=[]
    q_values=(3,4,8,16)
    for rho in (Fraction(3,10),Fraction(7,20),Fraction(2,5),Fraction(9,20),Fraction(1,2)):
        factor,w=prefix_envelope(spectrum,rho)
        best={q:math.inf for q in q_values}
        witnesses={}
        for u in (-9,-8,-7.5,-7,-6.5,-6,-5.5,-5,-4,-3,-2,-1,0,0.75):
            s=math.exp(u)
            zero,active=base.step_matrices(math.exp(-s),22)
            candidate=(1-float(rho))*zero+float(rho)*active
            region=uniform_coefficients(base.log_entries(zero),base.log_entries(candidate),8192,16)
            powered=base.log_power(region[list(q_values)],256)
            moments=np.logaddexp(powered[:,0,0],powered[:,0,1])
            for q,moment in zip(q_values,moments):
                value=math.log(math.comb(8192,q))+q*math.log(factor)+min(0,float(moment)+209716*s)
                if value<best[q]:
                    best[q]=value
                    witnesses[q]=u
        row=dict(rho=str(rho),integer_prefix_envelope=str(factor),
                 envelope_log2=math.log2(factor),maximizing_prefix_weight=w,
                 margin_bits={str(q):-v/math.log(2) for q,v in best.items()},best_log_tilts=witnesses)
        results.append(row)
        print(json.dumps(row),flush=True)
    write_new(output,dict(classification='Exact cumulative-envelope comparisons; binary64 transfer diagnostic only',
        results=results,scope='Only Q in {3,4,8,16}; no all-Q certificate',
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'code/evaluate_bch_q2_envelopes.py',
            GEN/'bch256_q2_envelope_diagnostic_sources/evaluate_ebch128_randomstepconv_g1.py',
            GEN/'bch256_q2_envelope_diagnostic_sources/evaluate_ebch128_randomstepconv_q1_exact.py')}))


if __name__=='__main__':
    main()
