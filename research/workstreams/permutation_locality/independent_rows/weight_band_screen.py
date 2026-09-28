"""Binary64 aggregate-input-weight band proposals, not certificates.

Every original row weight lies in [64,192]. Therefore the sum of the 4q
row weights lies in [256q,768q]. Each disjoint band can use a different
Bernoulli witness. The spectrum and moment arithmetic here is diagnostic;
an outward certificate must replay both the scalar and inner calculation.
"""
import argparse
from math import comb, log

import numpy as np
from scipy.special import logsumexp
from row_verify import float_inner
from row_counts import row_gamma_function
from shape_inner import prepare_shared, build
from averaged_windows import AveragedHighInner
from bch_joint_support import authenticated_caps
from aggregate_weight import AggregateWeights


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=80)
    parser.add_argument('--tilts',nargs='+',default=['.052','.06'])
    parser.add_argument('--probabilities',type=float,nargs='+',default=[.36,.4,.44,.48,.5,.52,.56,.6,.64,.7,.8])
    parser.add_argument('--bands',type=int,default=32)
    parser.add_argument('--cut',type=int,default=8)
    parser.add_argument('--full-feedback',type=int,default=6)
    parser.add_argument('--window-histogram',type=int,default=8)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--column-density',action='store_true')
    parser.add_argument('--feedback-density',type=int,choices=range(0,7),default=0)
    parser.set_defaults(class_tail=True,joint_cancellation=True,penalty='1')
    args = parser.parse_args()
    if not 1<=args.groups<=2048 or args.bands<1 or any(not 0<p<1 for p in args.probabilities):
        parser.error('invalid group count, band count, or witness probability')
    caps = authenticated_caps()
    gamma = row_gamma_function(caps,64,192,exclude_zero=True)
    shared = prepare_shared(args)
    q,rows = args.groups,4*args.groups
    scalar_model = AggregateWeights(caps,64,192,rows)
    lower,upper = 64*rows,192*rows
    width = max(1,(upper-lower+1+args.bands-1)//args.bands)
    bands = [(lo,min(upper,lo+width-1)) for lo in range(lower,upper+1,width)]
    best = [(np.inf,None) for _ in bands]
    witnesses = sorted(set(args.probabilities+[.5]))
    location = log(comb(2048,q))
    for tilt in args.tilts:
        args.tilt = tilt
        model = AveragedHighInner(build(args,shared))
        for p in args.probabilities:
            moment = float_inner(model,q,tilt,p,4,clamp=True)
            homogeneous = (moment+location+rows*gamma(p))/log(2)
            print('DIAGNOSTIC reference tilt/p/log2moment/log2homogeneous:',
                  tilt,p,(moment-float(tilt)*209715)/log(2),homogeneous,flush=True)
            for i,(lo,hi) in enumerate(bands):
                scalar,s = min((scalar_model.majorant_log(p,lo,hi,s),s) for s in witnesses)
                candidate = ((moment+location+scalar)/log(2),(tilt,p,s))
                if candidate[0] < best[i][0]:
                    best[i] = candidate
        score = logsumexp([x[0]*log(2) for x in best])/log(2)
        print('DIAGNOSTIC complete aggregate-band sum log2 upper:',score,flush=True)
    for band,(score,witness) in zip(bands,best):
        print('DIAGNOSTIC band',band,'log2 upper',score,'tilt/p/majorant',witness,flush=True)
    print('No outward arithmetic or full occupancy cover is asserted.',flush=True)


if __name__ == '__main__':
    main()
