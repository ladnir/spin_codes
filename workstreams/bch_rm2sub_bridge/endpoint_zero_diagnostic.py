"""Scalar adaptive envelope of paths that never activate, for diagnosis only."""
import argparse
import math
from pathlib import Path
import time
import numpy as np
from flint import arb,arb_poly,ctx

import certificate_search_core as core
import certificate_search_backend as backend

base=core.base


def zero_region(t,kernel,z,rows):
    power=arb_poly([kernel.get(j,0)*z**j for j in range(t+1)])
    current=arb_poly([1]);n=rows//t
    while n:
        if n&1:current=current*power
        n>>=1
        if n:power=power*power
    return [max(arb(0),(current[j]/math.comb(rows,j)).upper()) for j in range(rows+1)]


def adaptive_scalar(values,left,right):
    mantissas=[];powers=[]
    for value in values:
        value=value.upper()
        m,e=value.man_exp();power=int(e)+int(m).bit_length() if m else 0
        mantissas.append(float(value*arb(2)**(-power)));powers.append(power)
    current=np.array(mantissas);powers=np.array(powers,dtype=np.int64)
    for _ in range(1,len(values)):
        common=np.maximum(powers[:-1],powers[1:])
        a=np.ldexp(current[:-1],powers[:-1]-common)
        b=np.ldexp(current[1:],powers[1:]-common)
        updated=left[0]*a+right[0]*b
        for x,y in zip(left[1:],right[1:]):np.maximum(updated,x*a+y*b,out=updated)
        current,shift=np.frexp(updated);powers=common+shift
    return math.log(float(current[0]))+int(powers[0])*math.log(2)


def run(certificate,output):
    data=base.read(certificate);core.authenticate(data,base.ROOT)
    spec=data['task']['instance'];core.require(spec==core.instance(spec['configuration'],spec['message_exponent']),'Wrong instance')
    rows=spec['rows'];core.require(data['task']['interval']==[rows,rows],'Endpoint required')
    t,s,spectrum,kernel=core.inputs.load(spec['configuration']);ctx.prec=256;started=time.monotonic()
    witness=data['result']['witness'];ps=[base.decode(v) for v in witness['p']]
    costs=backend.sparse.costs_for(backend.sparse.BANDS,ps,core.inputs.caps_module.caps())
    roots=np.array([math.exp((math.log(c.numerator)-math.log(c.denominator))/256) for c in costs])
    probs=np.array(list(map(float,ps)));left,right=roots*(1-probs),roots*probs
    keep=backend.sparse.hull.indices(left,right)
    lam=(arb(witness['tilt'])/10).exp()
    value=adaptive_scalar(zero_region(t,kernel,(-lam).exp(),rows),left[keep],right[keep])
    log_bound=256*value+rows*math.log(len(ps))+spec['cutoff']*float(lam)
    full_margin=data['result']['diagnostic_margins'][0];zero_margin=-log_bound/math.log(2)
    result=dict(status='ADAPTIVE_NEVER_ACTIVATED_DIAGNOSTIC_ONLY',instance=spec,
        certificate_sha256=base.sha(certificate),full_adaptive_margin_bits=full_margin,
        never_activated_margin_bits=zero_margin,zero_to_full_log_ratio_bits=full_margin-zero_margin,
        zero_path_cannot_be_bad=38*rows>spec['cutoff'],seconds=time.monotonic()-started,
        source_sha256={**backend.sources(),Path(__file__).relative_to(base.ROOT).as_posix():base.sha(Path(__file__))})
    base.write_new(output,result)
    print({k:v for k,v in result.items() if k not in ('source_sha256','instance')},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificate',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.certificate,a.output)
