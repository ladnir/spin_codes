"""Diagnostic refinement using moment envelopes and a Perron-vector bound."""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from scipy.special import logsumexp
from bch_tail_reference import ROOT,GEN,L,B,D,contiguous,log_choose_array
from run_higher_endpoint_preflight import sha,write_new


def perron_log_bound(q,choose,s,rho,p):
    bit=(1+np.exp(-s))/2
    epsilon=2.**-22
    eta=rho*p
    a=1-eta+eta*epsilon*bit
    b=(1-epsilon)*bit
    c=eta*b
    d=epsilon*bit
    eigen=(a+b+np.hypot(a-b,2*np.sqrt(c*d)))/2
    factor=np.maximum(1,(eigen-b)/d)
    log_point=choose+q*np.log(p)
    interior=q<L
    log_point[interior]+=(L-q[interior])*np.log1p(-p[interior])
    return np.log(factor)+L*B*np.log(eigen)+D*s-B*log_point


def golden(function,lower,upper,steps=42):
    alpha=(math.sqrt(5)-1)/2
    a=np.broadcast_to(lower,upper.shape).copy()
    b=upper.copy()
    x=b-alpha*(b-a)
    y=a+alpha*(b-a)
    fx,fy=function(x),function(y)
    for _ in range(steps):
        left=fx<=fy
        b=np.where(left,y,b)
        a=np.where(left,a,x)
        x=b-alpha*(b-a)
        y=a+alpha*(b-a)
        fx,fy=function(x),function(y)
    return np.where(fx<fy,x,y)


def main():
    output=GEN/'bch256_tail_reference_refined.json'
    assert not output.exists()
    source=json.loads((GEN/'bch256_tail_reference_screen.json').read_text())
    moments=json.loads((GEN/'bch256_tail_moment_envelopes.json').read_text())
    refs={r['rho']:r for r in moments['references']}
    old_refs={r['rho']:r for r in source['references']}
    qs=np.arange(4,L+1)
    choose=log_choose_array(L,qs)
    best=np.full(len(qs),math.inf)
    witness=[None]*len(qs)
    for i,row in enumerate(source['rows']):
        w=row['witness']
        if w['method']=='exact_region':
            rho=w['rho']
            best[i]=row['log2_upper']*math.log(2)+int(qs[i])*(math.log(int(refs[rho]['factor']))-math.log(int(old_refs[rho]['factor'])))
            witness[i]=w.copy()
    for reference in moments['references']:
        rho=float(Fraction(reference['rho']))
        factor=int(reference['factor'])
        r=qs/L
        # Zero-collision limiting saddle, used only to initialize diagnostics.
        lo=r.copy()
        hi=np.minimum(1.,(4/9)/rho)*np.ones_like(r)
        for _ in range(48):
            p=(lo+hi)/2
            derivative=-rho/(1-rho*p)+.2*rho/(1-2*rho*p)
            expected=p+p*(1-p)*derivative
            lo=np.where(expected<r,p,lo)
            hi=np.where(expected>=r,p,hi)
        p=np.clip((lo+hi)/2,1e-9,1-1e-9)
        p[-1]=1.
        eta=np.minimum(rho*p,4/9)
        u=np.log(-np.log1p(-2*eta))
        for _ in range(10):
            u=golden(lambda v:perron_log_bound(qs,choose,np.exp(v),rho,p),u-.7,u+.7)
            ell=np.log(np.minimum(p,1-1e-12)/(1-np.minimum(p,1-1e-12)))
            def with_logits(x):
                selected=1/(1+np.exp(-x))
                selected[-1]=1.
                return perron_log_bound(qs,choose,np.exp(u),rho,selected)
            ell=golden(with_logits,ell-1.5,ell+1.5)
            p=1/(1+np.exp(-ell))
            p[-1]=1.
        values=choose+qs*math.log(factor)+np.minimum(0,perron_log_bound(qs,choose,np.exp(u),rho,p))
        improved=np.flatnonzero(values<best)
        for i in improved:
            witness[i]=dict(method='perron_conditioned',rho=reference['rho'],s=float(math.exp(u[i])),p=float(p[i]))
        best=np.minimum(best,values)
        print('REFINED',reference['rho'],'unresolved',contiguous(qs[best>-58*math.log(2)]),flush=True)
    payload=dict(classification='Binary64 diagnostic with exact moment-prefix envelopes; not outward-certified',
        rows=[dict(q=int(q),log2_upper=float(v/math.log(2)),witness=w) for q,v,w in zip(qs,best,witness)],
        unresolved_at_per_q_2_to_minus_58=contiguous(qs[best>-58*math.log(2)]),
        tail_margin_bits_diagnostic=-float(logsumexp(best))/math.log(2),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),ROOT/'code/bch_tail_reference.py',
            GEN/'bch256_tail_reference_screen.json',GEN/'bch256_tail_moment_envelopes.json')})
    write_new(output,payload)
    print(json.dumps({k:payload[k] for k in ('unresolved_at_per_q_2_to_minus_58','tail_margin_bits_diagnostic')},indent=2))


if __name__=='__main__':
    main()
