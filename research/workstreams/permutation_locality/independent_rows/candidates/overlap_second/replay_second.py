"""Selected q/support replay of the quadratic cancellation refinement."""
import argparse
from fractions import Fraction as Q
from flint import ctx

from second import build
from quadratic import zero_refine
from spectral_feedback import build as feedback_build
from replay import refined_family,optimize
from mass_density_screen import baseline,as_array


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--epoch-cache',required=True)
    parser.add_argument('--groups',type=int,default=80)
    parser.add_argument('--support',type=int,default=200)
    parser.add_argument('--tilt',default='.072')
    parser.add_argument('--penalty',default='.9')
    parser.add_argument('--iterations',type=int,default=30)
    parser.add_argument('--density-chord',action='store_true')
    parser.add_argument('--quadratic-density',action='store_true')
    args=parser.parse_args();ctx.prec=192
    if not 1<=args.groups<=2048 or not 38<=args.support<256 or args.iterations<1:
        parser.error('valid selected occupancy and support required')
    base,coefficients=refined_family(args.epoch_cache,args.tilt,args.penalty,args.density_chord)
    second=build(10)
    if args.quadratic_density:
        from density_second import build as density_build,refine as density_refine
        from group_moment import maps
        records,feedback=density_build(10,[args.tilt],second)
    else:
        feedback=feedback_build(10)
    changed=zero_refine(base,second,feedback,args.tilt,args.penalty)
    if args.quadratic_density:
        changed=density_refine(changed,records[args.tilt],maps()[2],args.penalty,minimum=5,maximum=10)
    caps=baseline.authenticated_caps()
    cdf=baseline.integer_cdf(baseline.weighted_cdf_upper(caps,1<<128,full_weight=1/Q(args.penalty)))
    shells=baseline.weighted_union_shells(caps,full_weight=1/Q(args.penalty));shells[0]-=1
    count=min(Q(cdf[args.support]),shells[args.support])
    selected=optimize(changed,coefficients,args.groups,args.support,count,args.tilt,args.iterations)
    exact=baseline.placement(selected,rounding=baseline.rounded,maximum_groups=args.groups)
    label=args.penalty+':quadratic-zero:optimized-mass'
    if args.density_chord:label+=':density-chord'
    if args.quadratic_density:label+=':quadratic-density'
    args.probe_supports=[args.support];args.target_bits=52
    baseline.shell_points(args,{(args.tilt,label):(exact,as_array(exact))},{label:cdf},{label:shells})
    print('Selected support vector only; no full occupancy or full-code certificate.',flush=True)


if __name__=='__main__':main()
