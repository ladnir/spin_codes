"""Exact cumulative caps from squared Krawtchouk polynomials of degree <=7.

The degree-15 orthogonal-array identities and complement symmetry are the
same retained BCH-sandwich assumptions used by the earlier exact LP model.
"""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
from bch_tail_reference import ROOT,GEN,envelopes
from run_higher_endpoint_preflight import sha,write_new


def kraw(n,j,w):
    return sum((-1)**i*math.comb(w,i)*math.comb(n-w,j-i)
               for i in range(max(0,j-(n-w)),min(j,w)+1))


def main():
    output=GEN/'bch256_tail_moment_envelopes.json'
    assert not output.exists()
    n=256
    mass=1<<128
    ks=[[kraw(n,j,w) for w in range(n+1)] for j in range(8)]
    norms=[math.comb(n,j) for j in range(8)]
    # Verify all orthogonality identities used for the polynomial expectation.
    for i in range(8):
        for j in range(8):
            assert sum(math.comb(n,w)*ks[i][w]*ks[j][w] for w in range(n+1))==((1<<n)*norms[i] if i==j else 0)
    spectrum=envelopes()['LP_only']
    cumulative=[]
    previous=0
    witnesses=[]
    for w in range(n+1):
        previous=min(mass-1,previous+spectrum.get(w,0))
        cap=previous
        witness=dict(method='cumulative existing shell caps')
        if 38<=w<128:
            cap=min(cap,(mass-2)//2)
            for degree in range(1,8):
                # Single orthogonal square and same-parity reproducing kernel.
                for indices in ([degree],list(range(degree%2,degree+1,2))):
                    scale=math.lcm(*(norms[j] for j in indices))
                    coefficients={j:ks[j][w]*scale//norms[j] for j in indices}
                    values=[sum(coefficients[j]*ks[j][v] for j in indices) for v in range(n+1)]
                    assert all(values[n-v]==(-1)**degree*values[v] for v in range(n+1))
                    minimum=min(values[v]**2 for v in range(38,w+1,2))
                    if not minimum:
                        continue
                    expected_scaled=sum(coefficients[j]**2*norms[j] for j in indices)
                    numerator=mass*expected_scaled-2*values[0]**2
                    assert numerator>=0
                    bound=numerator//(2*minimum)
                    if bound<cap:
                        cap=bound
                        witness=dict(method='symmetric polynomial square',degree=degree,
                            coefficients={str(j):str(c) for j,c in coefficients.items()},
                            expectation_scaled=str(expected_scaled),zero_value_squared=str(values[0]**2),
                            minimum_squared_on_prefix=str(minimum))
        cumulative.append(cap)
        witnesses.append(witness)
    # A cap at a later prefix also bounds every earlier prefix.
    for w in range(255,-1,-1):
        if cumulative[w+1]<cumulative[w]:
            cumulative[w]=cumulative[w+1]
            witnesses[w]=dict(method='later prefix',prefix=w+1)
    references=[]
    for a in range(30,51):
        rho=Fraction(a,100)
        b=rho.denominator
        a0=rho.numerator
        cdf=0
        candidates=[]
        for w in range(n+1):
            cdf+=math.comb(n,w)*a0**w*(b-a0)**(n-w)
            candidates.append((Fraction(cumulative[w]*b**n,cdf),w))
        ratio,w=max(candidates)
        factor=-(-ratio.numerator//ratio.denominator)
        references.append(dict(rho=str(rho),factor=str(factor),factor_log2=math.log2(factor),maximizing_prefix=w))
    payload=dict(classification='Exact integer polynomial-square cumulative caps and rational prefix envelopes',
        assumptions=['Fixed dimension-128 even complement-symmetric C with minimum distance at least 38',
                     'The retained OA15 identities: expectations of all degree<=15 polynomials match Bin(256,1/2)',
                     'Existing deterministic shell caps; no statistical caps used'],
        cumulative_nonzero_caps=cumulative,witnesses=witnesses,references=references,
        all_64_Krawtchouk_orthogonality_checks_passed=True,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'code/evaluate_bch_q2_envelopes.py',
            GEN/'coupled_lp_exact.json',GEN/'exact_bch_lp_caps.json',GEN/'johnson_n256_w52_d38.json',GEN/'johnson_n256_w54_d38.json')})
    write_new(output,payload)
    print(json.dumps(references,indent=2))


if __name__=='__main__':
    main()
