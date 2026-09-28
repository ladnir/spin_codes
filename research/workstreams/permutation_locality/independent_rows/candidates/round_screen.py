"""Four-bit route: compare actual local envelopes for different update counts.

Selected homogeneous support vectors and binary64 proposals only. This
does not reuse a two-update numerical certificate for another ensemble.
"""
import argparse
from fractions import Fraction as Q

from mass_density_screen import epoch_grid,as_array,float_placement,score,baseline


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rounds',type=int,nargs='+',default=[2,3])
    parser.add_argument('--tilts',nargs='+',default=['.048','.056','.064'])
    parser.add_argument('--penalty',default='.9')
    parser.add_argument('--groups',type=int,nargs='+',default=[64,80])
    parser.add_argument('--supports',type=int,nargs='+',default=[192,200,208])
    args=parser.parse_args()
    if (any(not 1<=r<=32 for r in args.rounds) or any(Q(t)<=0 for t in args.tilts)
            or not 0<Q(args.penalty)<=1 or any(not 1<=q<=2048 for q in args.groups)
            or any(not 38<=u<256 for u in args.supports)):
        parser.error('valid update counts, tilts, penalty, groups, and supports required')
    print('Selected four-bit support events, cutoff 209715; binary64 screen, NOT a full certificate.',flush=True)
    caps=baseline.authenticated_caps()
    cdf=baseline.integer_cdf(baseline.weighted_cdf_upper(caps,1<<128,full_weight=1/Q(args.penalty)))
    shells=baseline.weighted_union_shells(caps,full_weight=1/Q(args.penalty))
    for rounds in args.rounds:
        print('Building actual operators for update count',rounds,flush=True)
        operators,_=epoch_grid(args.tilts,args.penalty,rounds=rounds)
        for tilt,ops in operators.items():
            region=float_placement(as_array(ops),max(args.groups))
            for q in args.groups:
                for u in args.supports:
                    value,p=score(region,q,u,min(Q(cdf[u]),shells[u]),tilt)
                    print('FOUR-BIT updates/q/u',rounds,q,u,'tilt/rho',tilt,args.penalty,
                          'log2 proposal',value,'support tilt',p,flush=True)


if __name__=='__main__':main()
