"""Outward replay for one equal-support class, not a full occupancy proof."""
import argparse
from fractions import Fraction
from math import comb
import numpy as np
from flint import arb,ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from occupancy_allones import full_column_caps,self_test
from occupancy_weighted import allocate
from occupancy_memory import prepare,epoch_operators,rounded,TERMINAL
from occupancy_model import local_data,placement
from occupancy_screen import optimize
from occupancy_memory_verify import replay,DENOMINATOR
from group_rank_one_verify import up


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=16,choices=range(4,33))
    parser.add_argument('--support',type=int,default=80)
    parser.add_argument('--tilt',default='.0064')
    parser.add_argument('--penalty',default='.5')
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--fresh',action='store_true')
    parser.add_argument('--window-average',action='store_true')
    args=parser.parse_args()
    rho=Fraction(args.penalty)
    # The reused replay may cap its tilted bound at one. This is valid
    # here because rho^J<=1; do not silently allow rho>1 in this verifier.
    assert 0<rho<=1 and args.precision>=128 and 38<=args.support<=256
    self_test()
    spectrum=authenticated_caps()
    dimensions=dimension_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    total=sum(row[args.support] for row in caps)
    shells=full_column_caps(spectrum,caps,dimensions,args.support,total)
    weighted_count=sum((c*rho**(-j) for j,c in allocate(shells,total,rho)),Fraction(0))
    prepared=prepare(local_data(4))
    if args.fresh:
        from occupancy_fresh_moment import fresh_census,refine
        fresh=fresh_census(prepared)
    if args.window_average:
        import occupancy_window_average as window_average
        windows=window_average.prepare_inputs()
    ctx.prec=args.precision
    if args.fresh:
        ops=refine(prepared,fresh,args.tilt,args.groups,args.penalty)
    else:
        ops=epoch_operators(prepared,args.tilt,maximum_groups=args.groups,pair_conditioned=True,full_penalty=args.penalty)
    if args.window_average:
        values=window_average.averages(windows,args.tilt)
        window_average.self_test(windows,values,args.tilt)
        ops=window_average.refine(prepared,values,args.tilt,args.groups,args.penalty,ops)
    regions=placement(ops,rounding=rounded)
    arrays=[np.array([[float(t[i,j]) for j in range(9)] for i in range(9)]) for t in regions]
    _,ps=optimize(arrays,(args.support,)*args.groups,float(args.tilt),TERMINAL)
    nums=[DENOMINATOR if args.support==256 else max(1,min(DENOMINATOR-1,round(p*DENOMINATOR))) for p in ps]
    moment=replay(regions,args.tilt,nums,[(args.support,args.support,0)]*args.groups)
    count=arb(weighted_count.numerator)/weighted_count.denominator
    upper=up(moment*comb(2048,args.groups)*count**args.groups)
    print('OUTWARD',args.precision,'groups',args.groups,'each support',args.support,
          'tilt',args.tilt,'all-one penalty',args.penalty,'fresh',args.fresh,'window-average',args.window_average,flush=True)
    print('Witness probability numerators over',DENOMINATOR,nums,flush=True)
    print('Contribution upper',upper,'margin',-upper.log()/arb(2).log(),flush=True)
    assert 0<upper<arb(2)**-40
    print('VERIFIED: this single support class, all ranks and group locations, <2^-40.')
    print('Other support vectors at this occupancy remain uncovered; NOT a full-code certificate.')


if __name__=='__main__':
    main()
