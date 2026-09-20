"""Log-domain zero-path diagnostic; supersedes underflow-prone scalar screens."""
import argparse
import math
from pathlib import Path
import time
import numpy as np
from flint import arb,ctx
import endpoint_zero_exact_diagnostic as exact

core=exact.core;backend=exact.backend;base=exact.base


def adaptive_log(values,left,right):
    current=np.array([float(v.log()) if v>0 else -np.inf for v in values])
    with np.errstate(divide='ignore'):la,lb=np.log(left),np.log(right)
    for _ in range(1,len(values)):
        a,b=current[:-1],current[1:]
        updated=np.logaddexp(la[0]+a,lb[0]+b)
        for x,y in zip(la[1:],lb[1:]):np.maximum(updated,np.logaddexp(x+a,y+b),out=updated)
        current=updated
    return float(current[0])


def run(certificate,output):
    data=base.read(certificate);core.authenticate(data,base.ROOT)
    spec=data['task']['instance'];core.require(spec==core.instance(spec['configuration'],spec['message_exponent']),'Wrong instance')
    rows=spec['rows'];core.require(data['task']['interval']==[rows,rows],'Endpoint required')
    t,s,spectrum,kernel=core.inputs.load(spec['configuration']);ctx.prec=256;started=time.monotonic()
    witness=data['result']['witness'];ps=[base.decode(v) for v in witness['p']]
    costs=backend.sparse.costs_for(backend.sparse.BANDS,ps,core.inputs.caps_module.caps())
    roots=np.array([math.exp((math.log(c.numerator)-math.log(c.denominator))/256) for c in costs])
    probs=np.array(list(map(float,ps)));left,right=roots*(1-probs),roots*probs
    keep=backend.sparse.hull.indices(left,right);lam=(arb(witness['tilt'])/10).exp()
    value=adaptive_log(exact.zero_region(t,kernel,(-lam).exp(),rows),left[keep],right[keep])
    log_bound=256*value+rows*math.log(len(ps))+spec['cutoff']*float(lam)
    full_margin=data['result']['diagnostic_margins'][0];zero_margin=-log_bound/math.log(2)
    result=dict(status='LOG_DOMAIN_ZERO_PATH_DIAGNOSTIC_ONLY',instance=spec,
        certificate_sha256=base.sha(certificate),full_adaptive_margin_bits=full_margin,
        never_activated_margin_bits=zero_margin,zero_to_full_log_ratio_bits=full_margin-zero_margin,
        zero_path_cannot_be_bad=38*rows>spec['cutoff'],seconds=time.monotonic()-started,
        source_sha256={**backend.sources(),**{Path(p).relative_to(base.ROOT).as_posix():base.sha(Path(p))
            for p in (__file__,exact.__file__)}})
    base.write_new(output,result)
    print({k:v for k,v in result.items() if k not in ('source_sha256','instance')},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificate',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.certificate,a.output)
