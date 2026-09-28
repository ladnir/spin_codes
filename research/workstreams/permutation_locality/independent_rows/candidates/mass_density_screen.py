"""Compare whole-column mass-density choices with the refined baseline.

This is a selected-point binary64 screen, never a certificate. Local
matrices are regenerated from authenticated maps and exact censuses.
No snapshot is read, and no production verifier or cache is modified.
"""
import argparse
from fractions import Fraction as Q
from math import comb, log
from pathlib import Path
import sys

import numpy as np
from flint import ctx
from scipy.optimize import minimize_scalar

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import verify as baseline
import cancellation_joint
import feedback_density
import full_feedback_refinement
import universal_density
import window_histogram
from mass_density import blend, conditioned_coefficients
from occupancy_sensitivity import float_placement
from occupancy_memory import Z, C


def epoch_grid(tilts, penalty, precision=192, *, rounds=2):
    """Complete-shape baselines sharing the expensive exact censuses."""
    if type(rounds) is not int or not 1<=rounds<=32:
        raise ValueError('integer update count in 1..32 required')
    tilts=tuple(dict.fromkeys(map(str,tilts)))
    if not tilts or any(Q(tilt)<=0 for tilt in tilts):
        raise ValueError('at least one positive output tilt required')
    ctx.prec=precision
    prepared=baseline.prepare(baseline.local_data(4))
    fresh=baseline.fresh_census(prepared)
    inputs=baseline.prepare_inputs()
    tails=baseline.tail_census(inputs,prepared)
    pairs=baseline.pair_census(inputs,prepared)
    zeros=baseline.zero_census()
    feedback=full_feedback_refinement.census(6)
    full_feedback_refinement.check_pairs(feedback)
    histograms=window_histogram.census(prepared)
    cancellations=cancellation_joint.census()
    cancellation_joint.check_fresh(cancellations,prepared)
    cancellation_joint.check_feedback(cancellations,feedback)
    density=feedback_density.build(6,tilts,precision,rounds)
    ctx.prec=precision
    result={}
    for tilt in tilts:
        windows=baseline.averages(inputs,tilt)
        moments=window_histogram.moments(histograms,tilt,8)
        ops=baseline.fresh_refine(prepared,fresh,tilt,32,penalty,'1')
        ops=baseline.window_refine(prepared,windows,tilt,32,penalty,ops,'1')
        ops=baseline.multi_refine(ops,fresh,tilt,penalty,'1')
        ops=baseline.collision_refine(ops,prepared,tilt,penalty,'1')
        ops=baseline.zero_refine(ops,zeros,tilt,penalty,'1')
        ops=baseline.lift(ops,inputs[3],tails,tilt,penalty,64,pairs,input_penalty='1')
        ops=baseline.transform(ops,inputs[3],tilt,rounds)
        ops=window_histogram.refine(ops,histograms,moments,penalty,rounds,'1')
        ops=full_feedback_refinement.refine(ops,feedback,inputs[3],tilt,penalty,rounds,'1')
        ops=cancellation_joint.refine(ops,cancellations,tilt,penalty,rounds,'1')
        ops=universal_density.refine(ops,inputs[3],tilt,penalty,'1',
                                   window_averages=windows,feedback=density,rounds=rounds)
        baseline.tail_test(inputs,ops,tilt,penalty,rounds,'1')
        result[tilt]=ops
    return result,feedback


def epochs(tilt, penalty, precision=192, *, rounds=2):
    """The current complete-shape baseline; retain its integer census."""
    result,feedback=epoch_grid([tilt],penalty,precision,rounds=rounds)
    return result[str(tilt)],feedback


def as_array(ops):
    return np.array([[[float(t[i,j]) for j in range(11)] for i in range(11)] for t in ops])


def score(region, q, support, count, tilt):
    def objective(z):
        p=1/(1+np.exp(-z))
        return (baseline.log_power_moment(baseline.matrix_for_probabilities(region,[p]*q),256,
                                          baseline.TAIL_TERMINAL)
                +float(tilt)*209715-q*baseline.log_binomial_mass(256,support,p))
    fit=minimize_scalar(objective,bounds=(-8.,16.),method='bounded')
    outer=q*(log(count.numerator)-log(count.denominator))+log(comb(2048,q))
    return (fit.fun+outer)/log(2),1/(1+np.exp(-fit.x))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilt',default='.056')
    parser.add_argument('--penalty',default='.9')
    parser.add_argument('--groups',type=int,nargs='+',default=[64,80])
    parser.add_argument('--supports',type=int,nargs='+',default=[192,200,208])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--condition-through',type=int,default=6,choices=range(6,33),help='Extend the six-window census by conditioning on exposed packets')
    parser.add_argument('--zero-returns',action='store_true',help='Apply the alternative to mature lazy zero returns, preserving refresh coefficients')
    args=parser.parse_args()
    if (Q(args.tilt)<=0 or not 0<Q(args.penalty)<=1 or args.precision<128
            or any(not 1<=q<=2048 for q in args.groups)
            or any(not 38<=u<256 for u in args.supports)):
        parser.error('positive tilt, penalty in (0,1], precision >=128, q in 1..2048, support in 38..255 required')
    print('BINARY64 SELECTED-POINT DIAGNOSTIC ONLY; no complete coverage or certificate',flush=True)
    base,data=epochs(args.tilt,args.penalty,args.precision)
    candidate=conditioned_coefficients(data,6,args.condition_through,args.tilt,args.penalty)
    spectrum=baseline.authenticated_caps()
    cdf=baseline.integer_cdf(baseline.weighted_cdf_upper(spectrum,1<<128,full_weight=1/Q(args.penalty)))
    shells=baseline.weighted_union_shells(spectrum,full_weight=1/Q(args.penalty))
    counts={u:min(Q(cdf[u]),shells[u]) for u in args.supports}
    choices=[('baseline',{})]
    choices += [(f'all-{fraction}',{j:fraction for j in range(1,7)}) for fraction in (Q(1,4),Q(1,2),Q(1))]
    choices += [(f'window-{j}',{j:Q(1)}) for j in range(1,7)]
    if args.condition_through>6:
        choices += [(f'window-{j}',{j:Q(1)}) for j in range(7,args.condition_through+1)]
        choices += [(f'from-{start}',{j:Q(1) for j in range(start,args.condition_through+1)})
                    for start in (5,6,7)]
    for name,fractions in choices:
        changed=blend(base,candidate,fractions,target=Z if args.zero_returns else C)
        region=float_placement(as_array(changed),max(args.groups))
        for q in args.groups:
            for u in args.supports:
                value,p=score(region,q,u,counts[u],args.tilt)
                print('DIAGNOSTIC NOT CERTIFICATE:',name,'target','zero' if args.zero_returns else 'density',
                      'q/u',q,u,'tilt/rho',args.tilt,args.penalty,
                      'log2 score',value,'p',p,flush=True)


if __name__=='__main__':
    main()
