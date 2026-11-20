"""Diagnostic decomposition of a dense witness into pure-group contributions."""
import argparse
import math
import numpy as np
from scipy.stats import binom
from flint import arb,ctx
import bridge as base
import polynomial_regions as poly
import christoffel_caps as outer
import general_occupancy as general
import scaled_adaptive as scaled

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--screen',required=True);a=p.parse_args()
    screen=base.read(base.HERE/'generated'/a.screen);ctx.prec=192;caps=outer.deterministic_caps()
    t,s,spectrum=base.load_map('t128_s15')
    for row in screen['rows']:
        q=row['occupation'];lam=(arb(row['witness_tenth'])/10).exp()
        region=poly.regions(t,s,spectrum,poly.tight.kernel_spectrum(),(-lam).exp(),q)
        mantissas,exponents=scaled.initial(region,outward=False)
        with np.errstate(divide='ignore'):logregion=np.log(mantissas)+exponents[:,None,None]*math.log(2)
        for band,encoded in zip(screen['bands'],row['p']):
            p=float(base.decode(encoded))
            gamma=max(math.log(caps[w])-math.log(math.comb(256,w))-w*math.log(p)-(256-w)*math.log1p(-p) for w in band)
            terms=binom.logpmf(np.arange(q+1),q,p)[:,None,None]+logregion
            shift=float(terms.max());matrix=np.exp(terms-shift).sum(axis=0)
            value=general.log_power(matrix)+256*shift+q*gamma+209716*float(lam)+math.log(math.comb(8192,q))
            print('Q',q,'band',band[0],band[-1],'p',round(p,5),'pure margin',round(-value/math.log(2),3),flush=True)
