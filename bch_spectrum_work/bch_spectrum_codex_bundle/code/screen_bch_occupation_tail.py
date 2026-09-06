"""Screen every remaining occupation with exact-prefix reference envelopes."""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from scipy.special import logsumexp
from bch_tail_reference import (ROOT,GEN,L,B,D,base,envelopes,prefix_envelope,
    normalized_regions,log_choose_array,dense_moment_log,contiguous,toy_check)
from run_higher_endpoint_preflight import sha,write_new


def main():
    output=GEN/'bch256_tail_reference_screen.json'
    assert not output.exists()
    checks=toy_check()
    qs=np.arange(4,L+1)
    log_choose=log_choose_array(L,qs)
    best=np.full(len(qs),math.inf)
    witness=[None]*len(qs)
    references=[]
    spectrum=envelopes()['LP_only']
    for rho in (Fraction(3,10),Fraction(7,20),Fraction(2,5),Fraction(9,20),Fraction(1,2)):
        factor,w=prefix_envelope(spectrum,rho)
        references.append(dict(rho=str(rho),factor=str(factor),maximizing_prefix=w))
        outer=log_choose+qs*math.log(factor)
        for u in np.arange(-8.,1.001,.25):
            s=math.exp(float(u))
            for scale in (1.,1.5,2.):
                p=np.minimum(.999999,scale*qs/L)
                p[-1]=1.
                moments=dense_moment_log(s,float(rho),qs,p,log_choose)
                bound=outer+np.minimum(0,moments+D*s)
                improved=np.flatnonzero(bound<best)
                for i in improved:
                    witness[i]=dict(method='conditioned',rho=str(rho),u=float(u),p=float(p[i]))
                best=np.minimum(best,bound)
        print('DENSE',str(rho),'unresolved',contiguous(qs[best>-58*math.log(2)]),flush=True)
    # Refine exact small occupations; no large-q inference from these values.
    maximum=128
    for rho in (Fraction(9,20),Fraction(1,2)):
        factor,_=prefix_envelope(spectrum,rho)
        small_q=np.arange(4,maximum+1)
        outer=log_choose[:len(small_q)]+small_q*math.log(factor)
        for u in np.arange(-8.,-4.499,.05):
            s=math.exp(float(u))
            zero,active=base.step_matrices(math.exp(-s),22)
            matrices=normalized_regions(zero,(1-float(rho))*zero+float(rho)*active,L,maximum)[4:]
            powered=base.log_power(base.log_entries(matrices),B)
            moments=np.logaddexp(powered[:,0,0],powered[:,0,1])
            values=outer+np.minimum(0,moments+D*s)
            improved=np.flatnonzero(values<best[:len(values)])
            for i in improved:
                witness[i]=dict(method='exact_region',rho=str(rho),u=float(u))
            best[:len(values)]=np.minimum(best[:len(values)],values)
            print('EXACT',str(rho),round(float(u),2),flush=True)
    rows=[dict(q=int(q),log2_upper=float(v/math.log(2)),witness=w) for q,v,w in zip(qs,best,witness)]
    payload=dict(classification='Binary64 diagnostic of mathematically justified bounds; not outward-certified',
        parameters=dict(outer_rows=L,outer_length=B,distance_cutoff=D,memory_bits=22),
        references=references,toy_checks=checks,rows=rows,
        tail_margin_bits_diagnostic=-float(logsumexp(best))/math.log(2),
        unresolved_at_per_q_2_to_minus_58=contiguous(qs[best>-58*math.log(2)]),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'code/bch_tail_reference.py',
            ROOT/'code/probe_cumulative_reference_q3.py',ROOT/'code/evaluate_bch_q2_envelopes.py',
            GEN/'exact_bch_lp_caps.json',GEN/'johnson_n256_w52_d38.json',GEN/'johnson_n256_w54_d38.json')})
    write_new(output,payload)
    print(json.dumps({k:payload[k] for k in ('tail_margin_bits_diagnostic','unresolved_at_per_q_2_to_minus_58')},indent=2))


if __name__=='__main__':
    main()
