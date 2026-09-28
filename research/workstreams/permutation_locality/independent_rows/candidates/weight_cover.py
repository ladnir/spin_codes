"""Complete support covers combining different total-input-weight witnesses.

Each operator choice retains the outer measure with its matching reciprocal
weight. This controller uses the existing unrestricted verifier and cover;
it does not import the experimental mass-density refinements. Only the
outward cover result certifies the requested occupancies.
"""
import argparse
from fractions import Fraction as Q
from pathlib import Path
import sys
from types import SimpleNamespace

from flint import arb

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from verify import (build_operators, authenticated_caps, weighted_cdf_upper,
                    weighted_union_shells, integer_cdf, tilted_support_caps,
                    TAIL_TERMINAL, up)
from shell_cover import cover


def measure_key(penalty, weight):
    """An opaque cover label binding both reciprocal outer factors."""
    rho,a=Q(penalty),Q(weight)
    if not 0<rho<=1 or a<=0:
        raise ValueError('require 0<rho<=1 and positive input-weight witness')
    return f'rho={rho};a={a}'


def append_grid(destination, count_sets, shell_sets, caps, args):
    """Build and label one grid from the same argument snapshot."""
    args=SimpleNamespace(**vars(args))
    weight=Q(args.weight_tilt)
    built=build_operators(args)
    for (tilt,penalty),value in built.items():
        label=measure_key(penalty,weight)
        if (tilt,label) in destination:
            raise ValueError('duplicate operator witness')
        if label not in count_sets:
            rho=Q(penalty)
            if weight==1:
                counts=weighted_cdf_upper(caps,1<<128,full_weight=1/rho)
                shells=weighted_union_shells(caps,full_weight=1/rho)
                shells[0]-=1
            else:
                shells,counts=tilted_support_caps(caps,1<<128,full_weight=1/rho,input_weight=weight)
            count_sets[label]=integer_cdf(counts)
            shell_sets[label]=shells
        destination[tilt,label]=value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occupancies',type=int,nargs='+',default=[64])
    parser.add_argument('--weight-tilts',nargs='+',default=['.98'])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--max-splits',type=int,default=400)
    parser.add_argument('--target-bits',type=int,default=52)
    parser.add_argument('--operator-cache')
    parser.add_argument('--screen-only',action='store_true')
    parser.add_argument('--no-prefix-rank',action='store_true')
    args=parser.parse_args()
    if (any(not 1<=q<=2048 for q in args.occupancies) or args.precision<128
            or args.max_splits<0 or args.target_bits<1
            or any(Q(a)<=0 or Q(a)==1 for a in args.weight_tilts)
            or len({Q(a) for a in args.weight_tilts})!=len(args.weight_tilts)):
        parser.error('require valid occupancies, precision >=128, positive target, and distinct positive extra weight tilts other than one')
    args.probe_supports=[]; args.probe_vector=[]
    args.retain_parents=True; args.joint_witness=True; args.joint_top=2
    args.full_feedback=6; args.window_histogram=8; args.joint_cancellation=True
    args.column_density=True; args.weight_tilt='1'
    caps=authenticated_caps()
    operators,counts,shells={},{},{}
    args.groups=max(64,max(args.occupancies))
    args.tilts=['.016','.024']; args.penalties=['.75','1']; args.feedback_density=0
    append_grid(operators,counts,shells,caps,args)
    args.tilts=['.032','.04','.048']; args.penalties=['.75','.9','1']; args.feedback_density=6
    append_grid(operators,counts,shells,caps,args)
    args.tilts=['.052','.056','.06']; args.penalties=['.9']
    append_grid(operators,counts,shells,caps,args)
    args.groups=max(80,max(args.occupancies))
    args.tilts=['.052','.056','.064','.072']
    for weight in args.weight_tilts:
        args.weight_tilt=weight
        append_grid(operators,counts,shells,caps,args)
    totals=[]
    for q in sorted(set(args.occupancies)):
        args.groups=q
        print('COMPLETE MIXED-WEIGHT COVER:',q,'groups;',len(operators),'operator choices;',
              'outer keys',sorted(counts),flush=True)
        result=cover(args,operators,counts,TAIL_TERMINAL,shells,prefix_rank=not args.no_prefix_rank)
        totals.append(result)
    if all(value is not None for value in totals):
        total=up(sum(totals,arb(0)))
        print('VERIFIED REQUESTED OCCUPANCIES',sorted(set(args.occupancies)),
              'upper',total,'margin',-total.log()/arb(2).log(),flush=True)
    else:
        print('No aggregate certificate from this run.',flush=True)
    print('This controller does not assemble the full q=1..2048 certificate.',flush=True)


if __name__=='__main__':
    main()
