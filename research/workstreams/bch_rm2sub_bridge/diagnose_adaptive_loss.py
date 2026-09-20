"""Compare pure-group moments with the adaptive bound (diagnostic only)."""
import argparse
import math
import numpy as np
from scipy.stats import binom
from scipy.optimize import minimize_scalar
import bridge as base
import tightened_occupancy as tight
import general_occupancy as general
import adaptive_range as adaptive
import occupation_three as q3

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--q',type=int,required=True)
    parser.add_argument('--tilt',type=int,required=True);args=parser.parse_args();q=args.q
    t,s,spectrum=base.load_map(tight.NAME);caps=general.q2.deterministic_caps();lam=math.exp(args.tilt/10)
    epoch=np.array(tight.epoch_matrices(t,s,spectrum,tight.kernel_spectrum(),math.exp(-lam),float,q)).reshape(-1,3,3)
    region=adaptive.normalized_regions(epoch,q)
    for band in q3.BANDS:
        def objective(theta):
            p=1/(1+math.exp(-theta))
            gamma=max(math.log(caps[w])-math.log(math.comb(256,w))-w*math.log(p)-(256-w)*math.log1p(-p) for w in band)
            matrix=np.sum(binom.pmf(np.arange(q+1),q,p)[:,None,None]*region,axis=0)
            scale=float(matrix.max())
            if scale==0:return math.inf
            return general.log_power(matrix/scale)+256*math.log(scale)+q*gamma+209716*lam+math.log(math.comb(8192,q))
        opt=minimize_scalar(objective,bounds=(-4,3),method='bounded')
        print('Pure',band[0],band[-1],'p',1/(1+math.exp(-opt.x)),'margin',-opt.fun/math.log(2),flush=True)
