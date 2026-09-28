"""Independent-row SPIN: stronger inner bounds and separate outward covers.

K=2^20, BCH[256,128], four-row packets, IMT(128,19), two updates.
Imports universal inner/cover helpers without changing the shared-route
entry point. Only a returned outward cover establishes an occupancy bound.
"""
import argparse
from fractions import Fraction as Q
from math import comb, log
from pathlib import Path
import sys

import numpy as np
from flint import arb, arb_mat, ctx
from scipy.optimize import minimize_scalar
from support import weighted_cdf_upper, weighted_union_shells, integer_cdf, tilted_support_caps

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from bch_joint_support import authenticated_caps
from occupancy_model import local_data, placement
from occupancy_memory import prepare, rounded
from occupancy_fresh_moment import fresh_census, refine as fresh_refine
from occupancy_window_average import averages, refine as window_refine
from occupancy_multi_average import refine as multi_refine
from fresh_collision import refine as collision_refine
from zero_moment import census as zero_census, refine as zero_refine
from mature_tail import prepare_inputs, census as tail_census, lift, TAIL_TERMINAL, self_test as tail_test
from pair_tail import census as pair_census
from mixing_rounds import transform, self_test as mixing_test
from occupancy_cdf_cover import (cover, self_test as folding_test, geometry_test,
                                retained_test, placement_prefix_test)
from group_rank_one_verify import up
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment, log_binomial_mass


def shell_points(args,operators,count_sets,shell_sets):
    """Outward replay for selected homogeneous supports, not a cover."""
    q=args.groups
    for u in args.probe_supports:
        best=None
        for choice,(_,region) in operators.items():
            tilt,penalty=choice
            count=min(Q(count_sets[penalty][u]),shell_sets[penalty][u])
            if count==0:
                continue
            def objective(p):
                value=log_power_moment(matrix_for_probabilities(region,[p]*q),256,TAIL_TERMINAL)
                return value+float(tilt)*209715-q*log_binomial_mass(256,u,p)
            if u==256:
                numerator,denominator=1,1
            else:
                fit=minimize_scalar(lambda z:objective(1/(1+np.exp(-z))),bounds=(-12,18),method='bounded')
                denominator=10**9
                numerator=max(1,min(denominator-1,round(denominator/(1+np.exp(-fit.x)))))
            value=objective(numerator/denominator)+q*(log(count.numerator)-log(count.denominator))+log(comb(2048,q))
            if best is None or value<best[0]:
                best=value,choice,numerator,denominator,count
        if best is None:
            print('Independent-row shell point excluded by exact count:',q,u,flush=True)
            continue
        _,choice,numerator,denominator,count=best
        exact=operators[choice][0]
        p=arb(numerator)/denominator
        size=len(TAIL_TERMINAL)
        matrix=sum((arb(comb(q,j))*p**j*(1-p)**(q-j)*exact[j]
                    for j in range(q+1)),arb_mat(size,size))**256
        term=sum((matrix[0,j] for j in range(size) if TAIL_TERMINAL[j]),arb(0))
        mass=arb(comb(256,u))*p**u*(1-p)**(256-u)
        assert mass>0
        term=up(term*(arb(choice[0])*209715).exp()*comb(2048,q)
                *((arb(count.numerator)/count.denominator)/mass)**q)
        assert term>0
        print('OUTWARD INDEPENDENT-ROW SHELL POINT q/u',q,u,'choice',choice,
              'p',f'{numerator}/{denominator}','log2 upper',term.log()/arb(2).log(),
              'meets point budget',bool(term<arb(2)**(-args.target_bits)),flush=True)
    print('Only stated homogeneous support events are bounded; no support cover or full-code certificate.',flush=True)


def build_operators(args):
    ctx.prec=args.precision
    import operator_cache
    cache_directory=getattr(args,'operator_cache',None)
    sources=operator_cache.source_digest() if cache_directory else None
    size=len(TAIL_TERMINAL)
    result={}
    def arrays(exact):
        return [np.array([[float(t[i,j]) for j in range(size)] for i in range(size)]) for t in exact]
    if cache_directory:
        for tilt in args.tilts:
            for penalty in args.penalties:
                key=operator_cache.parameters(args,tilt,penalty,sources)
                exact=operator_cache.load(cache_directory,key)
                if exact is not None:
                    result[tilt,penalty]=exact,arrays(exact)
                    print('Reused exact dyadic operator memo:',tilt,penalty,'degree',args.groups,flush=True)
        if len(result)==len(set(args.tilts))*len(set(args.penalties)):
            return result
    prepared=prepare(local_data(4))
    fresh=fresh_census(prepared)
    inputs=prepare_inputs()
    tails=tail_census(inputs,prepared)
    pairs=pair_census(inputs,prepared)
    zeros=zero_census()
    if args.full_feedback:
        import full_feedback_refinement as feedback_refine
        feedback=feedback_refine.census(args.full_feedback)
        feedback_refine.check_pairs(feedback)
    if args.window_histogram:
        import window_histogram
        window_histogram.self_test()
        histograms=window_histogram.census(prepared)
    if args.joint_cancellation:
        import cancellation_joint
        cancellations=cancellation_joint.census()
        cancellation_joint.check_fresh(cancellations,prepared)
        if args.full_feedback>=3:
            cancellation_joint.check_feedback(cancellations,feedback)
    density_data=None
    if args.feedback_density:
        import feedback_density
        density_data=feedback_density.build(args.feedback_density,args.tilts,args.precision,2)
    ctx.prec=args.precision
    local=max(4,min(args.groups,32))
    for tilt in args.tilts:
        if all((tilt,penalty) in result for penalty in args.penalties):
            continue
        windows=averages(inputs,tilt)
        if args.window_histogram:
            moments=window_histogram.moments(histograms,tilt,args.window_histogram)
        for penalty in args.penalties:
            if (tilt,penalty) in result:
                continue
            a=args.weight_tilt
            ops=fresh_refine(prepared,fresh,tilt,local,penalty,a)
            ops=window_refine(prepared,windows,tilt,local,penalty,ops,a)
            ops=multi_refine(ops,fresh,tilt,penalty,a)
            ops=collision_refine(ops,prepared,tilt,penalty,a)
            ops=zero_refine(ops,zeros,tilt,penalty,a)
            ops=lift(ops,inputs[3],tails,tilt,penalty,64,pairs,input_penalty=a)
            ops=transform(ops,inputs[3],tilt,2)
            if args.window_histogram:
                ops=window_histogram.refine(ops,histograms,moments,penalty,2,a)
            if args.full_feedback:
                ops=feedback_refine.refine(ops,feedback,inputs[3],tilt,penalty,2,a)
            if args.joint_cancellation:
                ops=cancellation_joint.refine(ops,cancellations,tilt,penalty,2,a)
            if args.column_density or args.feedback_density:
                import universal_density
                ops=universal_density.refine(ops,inputs[3],tilt,penalty,a,
                    window_averages=windows if args.column_density else None,
                    feedback=density_data,rounds=2)
            tail_test(inputs,ops,tilt,penalty,2,a)
            exact=placement(ops,rounding=rounded,maximum_groups=args.groups)
            result[tilt,penalty]=exact,arrays(exact)
            if cache_directory:
                key=operator_cache.parameters(args,tilt,penalty,sources)
                operator_cache.save(cache_directory,key,exact)
            print('Independent-row inner operators built:',tilt,penalty,'degree',args.groups,flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=64,choices=range(1,2049),metavar='Q')
    parser.add_argument('--occupancies',type=int,nargs='+',choices=range(1,2049),metavar='Q')
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--operator-cache',help='Optional local exact-dyadic memo directory; omit for independent regeneration')
    parser.add_argument('--tilts',nargs='+',default=['.032','.048'])
    parser.add_argument('--penalties',nargs='+',default=['.95'])
    parser.add_argument('--weight-tilt',default='1',help='Coupled total-input-weight witness a>0 (not a construction change)')
    parser.add_argument('--full-feedback',type=int,default=0,choices=[0,*range(2,11)])
    parser.add_argument('--window-histogram',type=int,default=0,choices=range(13))
    parser.add_argument('--joint-cancellation',action='store_true')
    parser.add_argument('--column-density',action='store_true',help='Complete shape maxima of conditioned one-window density bounds')
    parser.add_argument('--feedback-density',type=int,default=0,choices=range(7),help='Complete exact feedback-density census through this local occupancy')
    parser.add_argument('--probe-supports',type=int,nargs='+',default=[])
    parser.add_argument('--probe-vector',type=int,nargs='+',action='append',default=[])
    parser.add_argument('--shell-points',action='store_true',help='Use exact shell caps and outward replay for homogeneous probe supports only')
    parser.add_argument('--shell-cover',action='store_true',help='Use both exact shell caps and CDF caps in every support interval')
    parser.add_argument('--prefix-rank',action='store_true',help='Combine nested prefix and shell constraints; implies --shell-cover')
    parser.add_argument('--max-splits',type=int,default=200)
    parser.add_argument('--target-bits',type=int,default=40)
    parser.add_argument('--screen-only',action='store_true')
    args=parser.parse_args()
    args.shell_cover=args.shell_cover or args.prefix_rank
    if args.precision<128 or any(Q(t)<=0 for t in args.tilts) or any(not 0<Q(p)<=1 for p in args.penalties) or Q(args.weight_tilt)<=0:
        parser.error('precision >=128, positive tilts, penalties in (0,1], and positive weight tilt required')
    if any(not 38<=u<=256 for u in args.probe_supports):
        parser.error('probe supports must be in [38,256]')
    if args.shell_points and (not args.probe_supports or args.probe_vector):
        parser.error('shell points require homogeneous probe supports and exclude probe vectors')
    occupancies=sorted(set(args.occupancies or [args.groups]))
    args.groups=max(occupancies)
    # These affect only the imported general support-cover routine.
    args.retain_parents=True
    args.joint_witness=True
    args.joint_top=2
    print('ENSEMBLE: INDEPENDENT ROW SHUFFLES; two inner updates. Old certificates are not inputs.',flush=True)
    print('Total-input-weight coefficient witness:',args.weight_tilt,flush=True)
    mixing_test();folding_test();geometry_test();retained_test();placement_prefix_test()
    spectrum=authenticated_caps()
    count_sets={}
    shell_sets={}
    for p in args.penalties:
        if Q(args.weight_tilt)!=1:
            shell_sets[p],exact_counts=tilted_support_caps(spectrum,1<<128,
                                      full_weight=1/Q(p),input_weight=Q(args.weight_tilt))
        else:
            exact_counts=weighted_cdf_upper(spectrum,1<<128,full_weight=1/Q(p))
        count_sets[p]=integer_cdf(exact_counts)
        if (args.shell_points or args.shell_cover) and p not in shell_sets:
            shell_sets[p]=weighted_union_shells(spectrum,full_weight=1/Q(p))
            shell_sets[p][0]-=1
        print('Exact averaged outer CDF, rounded upward to integers:',p,flush=True)
    operators=build_operators(args)
    totals=[]
    for q in occupancies:
        args.groups=q
        print('INDEPENDENT-ROW occupancy',q,flush=True)
        if args.shell_points:
            result=shell_points(args,operators,count_sets,shell_sets)
        elif args.shell_cover:
            from shell_cover import cover as combined_cover
            result=combined_cover(args,operators,count_sets,TAIL_TERMINAL,shell_sets,
                                  prefix_rank=args.prefix_rank)
        else:
            result=cover(args,operators,count_sets,TAIL_TERMINAL)
        totals.append(result)
    if all(value is not None for value in totals):
        total=arb(0)
        for value in totals:
            total=up(total+value)
        print('VERIFIED INDEPENDENT-ROW OCCUPANCIES',occupancies,'precision',args.precision,
              'upper',total,'margin',-total.log()/arb(2).log(),flush=True)
        if occupancies!=list(range(1,2049)):
            print('Other occupancies remain; NOT a full-code certificate.',flush=True)
    else:
        print('Independent-row run has no aggregate outward certificate.',flush=True)


if __name__=='__main__':
    main()
