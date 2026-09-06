"""Independent exact checks of the published OA anchor and prefix witnesses."""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
from pypdf import PdfReader
from anchors import L71_HALF
from spectrum_exact import full_from_symmetric_half,macwilliams,dual_min_distance,krawtchouk
from bch_quotient import generator_polynomial,binary_poly_degree,binary_poly_divmod
from evaluate_bch_q2_envelopes import envelopes
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'


def main():
    output=GEN/'bch256_tail_envelope_checks.json'
    assert not output.exists()
    paper=ROOT/'paper/e104-a_9_1321.pdf'
    text=PdfReader(paper).pages[5].extract_text().split('Table 7',1)[1].split('Table 8',1)[0]
    extracted={}
    for line in text.splitlines():
        parts=line.split()
        if len(parts)>1 and parts[0].isdigit() and parts[-1].isdigit():
            extracted[int(parts[0])]=int(parts[-1])
    assert extracted==L71_HALF
    low=full_from_symmetric_half(extracted)
    dual=macwilliams(low,71)
    assert min(dual)>=0 and sum(dual)==1<<185 and dual_min_distance(dual)==16
    g_low,g_q=generator_polynomial(59),generator_polynomial(39)
    assert 255-binary_poly_degree(g_low)==71
    assert 255-binary_poly_degree(g_q)==123
    assert binary_poly_divmod(g_low,g_q)[1]==0
    assert binary_poly_divmod((1<<255)-1,g_q)[1]==0
    for degree in range(16):
        assert sum(count*math.comb(w,degree) for w,count in enumerate(low))==(1<<(71-degree))*math.comb(256,degree)
    path=GEN/'bch256_tail_moment_envelopes.json'
    receipt=json.loads(path.read_text())
    for name,expected in receipt['source_sha256'].items():
        assert sha(ROOT/name)==expected
    caps=receipt['cumulative_nonzero_caps']
    witnesses=receipt['witnesses']
    assert len(caps)==len(witnesses)==257
    spectrum=envelopes()['LP_only']
    original=[]
    value=0
    for w in range(257):
        value=min((1<<128)-1,value+spectrum.get(w,0))
        original.append(min(value,(1<<127)-1) if 38<=w<128 else value)
    squares=0
    for w in range(256,-1,-1):
        witness=witnesses[w]
        if witness['method']=='later prefix':
            assert witness['prefix']==w+1
            assert caps[w]==caps[w+1]<=original[w]
        elif witness['method']=='cumulative existing shell caps':
            assert caps[w]==original[w]
        else:
            assert witness['method']=='symmetric polynomial square' and 38<=w<128
            coefficients={int(j):int(c) for j,c in witness['coefficients'].items()}
            assert all(0<=j<=7 for j in coefficients)
            values=[sum(c*krawtchouk(256,j,v) for j,c in coefficients.items()) for v in range(257)]
            square=[v*v for v in values]
            assert square==square[::-1]
            minimum=min(square[38:w+1:2])
            assert minimum>0 and minimum==int(witness['minimum_squared_on_prefix'])
            expectation=Fraction(sum(math.comb(256,v)*square[v] for v in range(257)),1<<256)
            assert expectation.denominator==1 and expectation==int(witness['expectation_scaled'])
            assert sum(a*p for a,p in zip(low,square))==(1<<71)*expectation
            assert square[0]==int(witness['zero_value_squared'])
            bound=((1<<128)*expectation-2*square[0])//(2*minimum)
            assert caps[w]==bound<original[w]
            squares+=1
    assert all(0<=caps[w]<=caps[w+1] for w in range(256))
    prefix_checks=0
    for reference in receipt['references']:
        rho=Fraction(reference['rho'])
        factor=int(reference['factor'])
        probabilities=[Fraction(math.comb(256,w))*rho**w*(1-rho)**(256-w) for w in range(257)]
        cdf=Fraction(0)
        ratios=[]
        for cap,probability in zip(caps,probabilities):
            cdf+=probability
            assert cap<=factor*cdf
            ratios.append(Fraction(cap,cdf))
            prefix_checks+=1
        maximum=max(ratios)
        assert factor==-(-maximum.numerator//maximum.denominator)
        assert ratios[reference['maximizing_prefix']]==maximum
    paths=[Path(__file__),paper,path,ROOT/'code/anchors.py',ROOT/'code/spectrum_exact.py',
           ROOT/'code/bch_quotient.py',ROOT/'code/affine_wambach.py',ROOT/'code/evaluate_bch_q2_envelopes.py']
    result=dict(classification='Independent exact anchor, polynomial-square and prefix-envelope checks',
        published_Table7_k71_matches_every_anchor_entry=True,
        anchor_dual_minimum_distance=16,
        generator_containment_L71_in_Q123_verified=True,
        Q_contains_all_one_verified=True,
        OA15_moments_checked=16,polynomial_square_witnesses_checked=squares,
        exact_prefix_comparisons_checked=prefix_checks,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
        proof_dependency='The published Table 7 is the external mathematical anchor, not a newly enumerated 2^71-word spectrum. L subset Q subset C implies C dual subset L dual and hence OA15.')
    write_new(output,result)
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':
    main()
