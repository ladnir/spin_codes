"""Never-active envelope with exact integer kernel coefficients before tilting."""
import argparse
import math
from pathlib import Path
import time
import numpy as np
from flint import arb,fmpz_poly,ctx
import endpoint_zero_diagnostic as original

core=original.core;backend=original.backend;base=core.base


def zero_region(t,kernel,z,rows):
    # Factoring z^j out avoids absolute Arb polynomial errors in tiny coefficients.
    polynomial=fmpz_poly([kernel.get(j,0) for j in range(t+1)])**(rows//t)
    return [(arb(polynomial[j])*z**j/math.comb(rows,j)).upper() for j in range(rows+1)]


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
    value=original.adaptive_scalar(zero_region(t,kernel,(-lam).exp(),rows),left[keep],right[keep])
    log_bound=256*value+rows*math.log(len(ps))+spec['cutoff']*float(lam)
    full_margin=data['result']['diagnostic_margins'][0];zero_margin=-log_bound/math.log(2)
    result=dict(status='EXACT_KERNEL_COEFFICIENT_ZERO_PATH_DIAGNOSTIC_ONLY',instance=spec,
        certificate_sha256=base.sha(certificate),full_adaptive_margin_bits=full_margin,
        never_activated_margin_bits=zero_margin,zero_to_full_log_ratio_bits=full_margin-zero_margin,
        zero_path_cannot_be_bad=38*rows>spec['cutoff'],seconds=time.monotonic()-started,
        source_sha256={**backend.sources(),**{Path(p).relative_to(base.ROOT).as_posix():base.sha(Path(p))
            for p in (__file__,original.__file__)}})
    base.write_new(output,result)
    print({k:v for k,v in result.items() if k not in ('source_sha256','instance')},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificate',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.certificate,a.output)
