"""Same adaptive inequality and Arb coefficients, with log-domain recurrence.

This is a binary64 numerical diagnostic, not an outward replay receipt.
"""
import argparse
import math
from pathlib import Path
import time
import numpy as np
from scipy.special import logsumexp
from flint import arb,ctx

import endpoint_band_diagnostic as bands

core=bands.core;backend=bands.backend;base=core.base


def adaptive_matrix(logs,left,right):
    with np.errstate(divide='ignore'):la,lb=np.log(left),np.log(right)
    current=logs.copy()
    for _ in range(1,len(current)):
        a,b=current[:-1],current[1:]
        updated=np.logaddexp(la[0]+a,lb[0]+b)
        for x,y in zip(la[1:],lb[1:]):np.maximum(updated,np.logaddexp(x+a,y+b),out=updated)
        current=updated
    return current[0]


def run(certificate,output):
    data=base.read(certificate);core.authenticate(data,base.ROOT)
    spec=data['task']['instance'];core.require(spec==core.instance(spec['configuration'],spec['message_exponent']),'Wrong instance')
    rows=spec['rows'];core.require(data['task']['interval']==[rows,rows],'Endpoint required')
    t,s,spectrum,kernel=core.inputs.load(spec['configuration']);ctx.prec=256;started=time.monotonic()
    witness=data['result']['witness'];ps=[base.decode(v) for v in witness['p']]
    costs=backend.sparse.costs_for(backend.sparse.BANDS,ps,core.inputs.caps_module.caps())
    # Use exactly the frozen outward coefficient pairs, not merely nearby roots.
    roots=[((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probs=[arb(p.numerator)/p.denominator for p in ps]
    left=np.nextafter(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probs)]),np.inf)
    right=np.nextafter(np.array([float((r*p).upper()) for r,p in zip(roots,probs)]),np.inf)
    keep=backend.sparse.hull.indices(left,right);lam=(arb(witness['tilt'])/10).exp()
    region=backend.sparse.poly.regions(t,s,spectrum,kernel,(-lam).exp(),rows,length=rows)
    logs=np.array([[float(v.log()) if v>0 else -np.inf for v in row] for row in region]).reshape(-1,3,3)
    print('log-domain region ready',time.monotonic()-started,flush=True)
    matrix=adaptive_matrix(logs,left[keep],right[keep]);final=bands.logpower(matrix,256)[0]
    moment=float(logsumexp(final));correction=rows*math.log(len(ps))+spec['cutoff']*float(lam)
    margin=-(moment+correction)/math.log(2)
    result=dict(status='ADAPTIVE_LOG_DOMAIN_DIAGNOSTIC_ONLY',instance=spec,full_coverage=False,
        certificate_sha256=base.sha(certificate),original_scaled_margin_bits=data['result']['diagnostic_margins'][0],
        log_domain_margin_bits=margin,final_state_fractions=np.exp(final-moment).tolist(),
        final_region_log_matrix=matrix.tolist(),hull_band_indices=keep,seconds=time.monotonic()-started,
        source_sha256={**backend.sources(),**{Path(p).relative_to(base.ROOT).as_posix():base.sha(Path(p))
            for p in (__file__,bands.__file__)}})
    base.write_new(output,result);print({k:v for k,v in result.items() if k not in ('source_sha256','instance')},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificate',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.certificate,a.output)
