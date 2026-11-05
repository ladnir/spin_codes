"""Fixed h6/a0/c6 no-full-reset branch; no construction changes.

Reuse verified local factors, but retain five episodes and arbitrary alignment.
All arithmetic bounds are integer/rational; no generated files.
"""
import argparse
import itertools
import time
from fractions import Fraction as F
from functools import lru_cache

from diagnose_pair_activation_route import M, NEG, choose
from diagnose_pair_aligned_frontier import local_envelope, single_table
from diagnose_pair_activation_route import episode_table
from diagnose_pair_region_rank import count_cap, early_joint, joint_cell
from diagnose_pair_defects import profile
from diagnose_pair_shared_regions import sharing_price

ENDPOINT_BITS = 156  # six intervals: (4097)^12 < 2^156


def legacy_branch(report=True,initial=0,episodes=6):
    began=time.monotonic()
    residual=[]; closed=NEG; worst=(NEG,0)
    assert initial in (0,1) and 1<=episodes<=6+initial
    for live in range(episodes,4096):
        base=profile(6,initial,episodes,live,460,True)
        if base is None:continue
        value=base['log2_upper']
        if value+186>-80:
            local=NEG
            for early,used,regions in early_joint(6,initial,episodes):
                cell=joint_cell(6,initial,episodes,live,early,used,regions,275)
                if cell is not None:
                    local=max(local,value if 'unresolved' in cell else min(value,cell['upper']))
            value=local
        if value+186>-80:residual.append(live)
        else:closed=max(closed,value+186)
        worst=max(worst,(value+186,live))
        if time.monotonic()-began>40:raise TimeoutError('legacy 40s limit')
    result=dict(branch=(6,initial,episodes),residual_live_counts=residual,closed_good_route_upper=closed,
                worst=worst,seconds=time.monotonic()-began)
    if report:print(result,flush=True)
    return result


def test():
    checks=0
    for active in (2,3,4):
        for crossings in range(6):
            upper=combined(active,2,crossings,3)
            for lengths in itertools.product((1,2),repeat=5):
                dp=[0]+[NEG]*(2*active)
                for index,length in enumerate(lengths):
                    table=(episode_table if index<crossings else single_table)(active,length)
                    nxt=[NEG]*len(dp)
                    for used,previous in enumerate(dp):
                        if previous==NEG:continue
                        for extra,value in enumerate(table[:len(dp)-used]):
                            if value!=NEG:nxt[used+extra]=max(nxt[used+extra],previous+value)
                    dp=nxt
                exact=tuple(itertools.accumulate(dp,max))
                assert all(a>=b for a,b in zip(upper,exact))
                checks+=1
    print(dict(five_episode_length_envelope_checks=checks),flush=True)


@lru_cache(maxsize=2048)
def combined(active,maximum,crossings,virtual_dimension):
    budget=(virtual_dimension-1)*active
    dp=[0]+[NEG]*budget
    for index in range(5):
        table=local_envelope(active,maximum,index<crossings)
        nxt=[NEG]*(budget+1)
        for used,previous in enumerate(dp):
            if previous==NEG:continue
            for extra,value in enumerate(table[:budget-used+1]):
                if value!=NEG:nxt[used+extra]=max(nxt[used+extra],previous+value)
        dp=nxt
    return tuple(itertools.accumulate(dp,max))


@lru_cache(maxsize=1024)
def eligible_cap(early,crossings,tail):
    fragments=5+crossings
    need=129-tail
    if fragments<need:return 0
    assert need>0
    mean=F(early+5,16*fragments)
    beta=min([F(1)]+[F(choose(fragments,j),choose(need,j))*mean**j
                     for j in range(1,need+1)])
    return count_cap(beta,ENDPOINT_BITS+100)


@lru_cache(maxsize=8192)
def source_term(active,maximum,crossings,virtual_dimension,sharing):
    actual_dimension=virtual_dimension-sharing
    envelope=combined(active,maximum,crossings,virtual_dimension)
    values=[]
    for zeros in range((actual_dimension-1)*active+1):
        route=envelope[zeros+sharing*active]
        if route!=NEG:
            values.append((choose((128+actual_dimension)*active,zeros)-1).bit_length()
                          +32*(actual_dimension*active-zeros)+route)
    if not values:return NEG
    return ((choose(M,active)-1).bit_length()+max(values)
            +(len(values)-1).bit_length()+sharing_price(active,sharing))


def run(live,screen=False,report=True):
    began=time.monotonic()
    worst=NEG; witness=None; cells=excluded=0; largest_cap=largest_dimension=0
    for maximum in range(1,9):
        for early in range(maximum+4,5*maximum+1):
            tail=(live-early+16)//16
            for crossings in range(6):
                virtual_dimension=5+crossings+tail-128
                if virtual_dimension<=0:continue
                largest_dimension=max(largest_dimension,virtual_dimension)
                cap=eligible_cap(early,crossings,tail)
                largest_cap=max(largest_cap,cap)
                for sharing in range(virtual_dimension):
                    cells+=1
                    if cap<1<<maximum:
                        excluded+=1;continue
                    if not screen:
                        terms=[source_term(s,maximum,crossings,virtual_dimension,sharing)
                               for s in range(1<<maximum,cap+1)]
                        value=max(terms)+1+(len(terms)-1).bit_length()+ENDPOINT_BITS
                        if value>worst:
                            worst=value
                            witness=dict(maximum=maximum,early=early,crossings=crossings,
                                         tail_regions=tail,sharing=sharing,cap=cap,
                                         virtual_dimension=virtual_dimension)
                    if time.monotonic()-began>40:raise TimeoutError('new count 40s limit')
    result=dict(profile=(6,0,6,live),cells=cells,excluded=excluded,
                largest_cap=largest_cap,largest_dimension=largest_dimension,
                good_route_upper=None if screen else worst,worst=witness,
                whole_code_certificate=False,seconds=time.monotonic()-began)
    if report:print(result,flush=True)
    return result


def integrate(clean_early_branch=False):
    began=time.monotonic()
    old=legacy_branch(False,initial=int(clean_early_branch),episodes=6)
    worst=NEG; witness=None; cells=0
    for live in old['residual_live_counts']:
        result=run(live,report=False)
        cells+=result['cells']
        if result['good_route_upper']>worst:worst,witness=result['good_route_upper'],live
        if time.monotonic()-began>40:raise TimeoutError('integration 40s limit')
    total=max(old['closed_good_route_upper'],-80,worst,-100)+2
    residual=old['residual_live_counts']
    kind=('h6_a1_c6_no_full_reset_clean_early_subclass' if clean_early_branch
          else 'integrated_h6_a0_c6_no_full_reset')
    print(dict(kind=kind,formerly_residual=len(residual),
               residual_range=(min(residual),max(residual)) if residual else None,
               new_count_cells=cells,new_good_route_upper=worst,worst_live=witness,
               old_remaining_good_route_upper=old['closed_good_route_upper'],
               old_route_exception=-80,new_route_exception=-100,covered_family_upper=total,
               complete_branch=not clean_early_branch,
               closed=total<=-40,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--legacy-branch',action='store_true')
    p.add_argument('--integrate',action='store_true')
    p.add_argument('--clean-early-branch',action='store_true')
    p.add_argument('--screen',action='store_true')
    p.add_argument('--test',action='store_true')
    p.add_argument('--initial',type=int,default=0)
    p.add_argument('--episodes',type=int,default=6)
    p.add_argument('--live',type=int,default=1938)
    args=p.parse_args()
    if args.test:test()
    elif args.legacy_branch:legacy_branch(initial=args.initial,episodes=args.episodes)
    elif args.clean_early_branch:integrate(True)
    elif args.integrate:integrate()
    else:run(args.live,args.screen)
