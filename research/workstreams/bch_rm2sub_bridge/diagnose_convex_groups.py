"""Pure-group tilt sweep and cost of convexifying region coefficients."""
import argparse
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import binom
from flint import arb,ctx
import bridge as base
import general_occupancy as general
import tightened_occupancy as tight
import polynomial_regions as poly
import occupation_three as groups
import christoffel_caps as caps
import scaled_adaptive as scaled
import convex_region_majorant as convex


def logs(region):
    mantissa,exponent=scaled.initial(region,outward=False)
    with np.errstate(divide='ignore'):return np.log(mantissa)+exponent[:,None,None]*math.log(2)


def log_moment(logregion,q,p):
    terms=binom.logpmf(np.arange(q+1),q,p)[:,None,None]+logregion[:q+1]
    shift=float(terms.max());matrix=np.exp(terms-shift).sum(axis=0)
    return general.log_power(matrix)+256*shift


def run(q,tilts,tag):
    t,s,spectrum=base.load_map(tight.NAME);kernel=tight.kernel_spectrum();upper=caps.deterministic_caps()
    ctx.prec=256;rows=[]
    for tenth in tilts:
        lam=(arb(tenth)/10).exp();region=poly.regions(t,s,spectrum,kernel,(-lam).exp(),q)
        majorant=convex.convexify(region);convex.assert_convex(majorant)
        original=logs(region);envelope=logs(majorant);results=[]
        for band in groups.BANDS:
            def objective(theta,source=envelope):
                p=1/(1+math.exp(-theta))
                gamma=max(math.log(upper[w])-math.log(math.comb(256,w))-w*math.log(p)-(256-w)*math.log1p(-p) for w in band)
                return log_moment(source,q,p)+q*gamma+209716*float(lam)+math.log(math.comb(8192,q))
            opt=minimize_scalar(objective,bounds=(-4,12),method='bounded',options={'xatol':1e-5})
            loss=(objective(float(opt.x))-objective(float(opt.x),original))/math.log(2)
            results.append(dict(band=list(band),p=1/(1+math.exp(-float(opt.x))),margin_bits=-float(opt.fun)/math.log(2),convex_loss_bits=loss))
        rows.append(dict(tilt=tenth,groups=results))
        print('Convex Q',q,'tilt',tenth,'margins',[round(r['margin_bits']) for r in results],
              'maximum convex loss',max(r['convex_loss_bits'] for r in results),flush=True)
    base.write_new(base.HERE/'generated'/f'convex_groups_{tag}.json',dict(status='PURE_GROUP_DIAGNOSTIC_ONLY',occupation=q,rows=rows,
        local_sha256={p.name:base.sha(p) for p in (Path(__file__),Path(convex.__file__),Path(poly.__file__),Path(caps.__file__))}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--q',type=int,required=True);p.add_argument('--tilts',nargs='+',type=int,required=True)
    p.add_argument('--tag',required=True);a=p.parse_args();run(a.q,a.tilts,a.tag)
