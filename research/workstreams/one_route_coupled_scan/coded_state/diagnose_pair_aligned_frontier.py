"""Actual h6/a1/c7 clean-episode frontier, including noncrossing intervals.

Uses upper envelopes over lengths, exact RS-zero costs, and sharing penalties.
Restricted global profile only; old proof helpers remain untouched.
"""
import argparse
import itertools
import math
import time
from fractions import Fraction as F
from functools import lru_cache

from diagnose_pair_activation_route import M,B,Q,NEG,choose,logceil_ratio,episode_table
from diagnose_pair_region_rank import count_cap
from diagnose_pair_shared_regions import sharing_price


@lru_cache(maxsize=512)
def single_table(active,length):
    out=[NEG]*(active+1)
    for count in range(1,active+1):
        numerator=sum(choose(B//2,w//2)*choose(length*w,count-w)
                      *length**w*Q**(active-w)
                      for w in range(1<<length,min(count,B)+1,2))
        if numerator:
            out[active-count]=logceil_ratio(numerator,choose(M,count)*Q**active)
    return tuple(out)


@lru_cache(maxsize=512)
def local_envelope(active,maximum,cross):
    tables=[episode_table(active,length) if cross else single_table(active,length)
            for length in range(1,maximum+1)]
    return tuple(max(x) for x in zip(*tables))


@lru_cache(maxsize=2048)
def combined(active,maximum,crossings):
    budget=3*active
    dp=[0]+[NEG]*budget
    for index in range(6):
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
    fragments=6+crossings
    need=129-tail
    if fragments<need:return 0
    mean=F(early+6,16*fragments)
    beta=min([F(1)]+[F(choose(fragments,j),choose(need,j))*mean**j
                     for j in range(1,need+1)])
    return count_cap(beta,282)  # 182 endpoint bits, 100-bit route exception


@lru_cache(maxsize=8192)
def source_term(active,maximum,crossings,virtual_dimension,sharing):
    actual_dimension=virtual_dimension-sharing
    envelope=combined(active,maximum,crossings)
    values=[]
    for zeros in range((actual_dimension-1)*active+1):
        route=envelope[zeros+sharing*active]
        if route!=NEG:
            values.append(((choose((128+actual_dimension)*active,zeros)-1).bit_length()
                           +32*(actual_dimension*active-zeros)+route,zeros))
    if not values:return NEG
    return ((choose(M,active)-1).bit_length()+max(values)[0]
            +(len(values)-1).bit_length()+sharing_price(active,sharing))


def run(live=1938,screen=False,maximum_limit=8,report=True):
    began=time.monotonic()
    worst=NEG;witness=None;cells=excluded=0;largest_cap=0
    for maximum in range(1,maximum_limit+1):
        for early in range(maximum+5,6*maximum+1):
            tail=(live-early+16)//16
            for crossings in range(6):  # first episode starts at batch zero
                virtual_dimension=6+crossings+tail-128
                if virtual_dimension<=0:continue
                assert virtual_dimension<=4
                cap=eligible_cap(early,crossings,tail)
                largest_cap=max(largest_cap,cap)
                for sharing in range(virtual_dimension):
                    cells+=1
                    if cap < 1<<maximum:
                        excluded+=1;continue
                    if not screen:
                        terms=[source_term(s,maximum,crossings,virtual_dimension,sharing)
                               for s in range(1<<maximum,cap+1)]
                        value=max(terms)+1+(len(terms)-1).bit_length()+182
                        if value>worst:
                            worst=value
                            witness=dict(maximum=maximum,early=early,crossings=crossings,
                                         tail_regions=tail,sharing=sharing,cap=cap,
                                         virtual_dimension=virtual_dimension)
                    if time.monotonic()-began>40:
                        print(dict(incomplete=True,cells=cells,worst=worst,witness=witness),flush=True)
                        raise TimeoutError('40-second aligned-frontier cap')
    result=dict(profile=(6,1,7,live),maximum_limit=maximum_limit,cells=cells,
               excluded=excluded,largest_eligible_cap=largest_cap,
               good_route_upper=None if screen else worst,route_exception_upper=-100,
               family_upper=None if screen else max(worst,-100)+1,worst=witness,
               whole_code_certificate=False,seconds=time.monotonic()-began)
    if report:print(result,flush=True)
    return result


def legacy_branch(report=True):
    from diagnose_pair_defects import profile
    from diagnose_pair_region_rank import early_joint,joint_cell
    began=time.monotonic()
    residual=[];worst=NEG;witness=None;closed_worst=NEG
    for live in range(7,4096):
        base=profile(6,1,7,live,460,True)
        if base is None:continue
        value=base['log2_upper']
        if value+186>-80:
            local=NEG
            for early,used,regions in early_joint(6,1,7):
                cell=joint_cell(6,1,7,live,early,used,regions,275)
                if cell is not None:
                    local=max(local,value if 'unresolved' in cell else min(value,cell['upper']))
            value=local
        if value+186>-80:residual.append(live)
        else:closed_worst=max(closed_worst,value+186)
        if value>worst:worst,witness=value,live
        if time.monotonic()-began>40:raise TimeoutError('40-second legacy branch cap')
    result=dict(kind='legacy_h6_a1_c7',residual_live_counts=residual,
               worst_live=witness,good_route_upper=worst+186,
               closed_good_route_upper=closed_worst,
               whole_code_certificate=False,seconds=time.monotonic()-began)
    if report:print(result,flush=True)
    return result


def integrate():
    began=time.monotonic()
    old=legacy_branch(False)
    worst=NEG;witness=None;cells=0
    for live in old['residual_live_counts']:
        result=run(live=live,report=False)
        cells+=result['cells']
        if result['good_route_upper']>worst:
            worst,witness=result['good_route_upper'],live
        if time.monotonic()-began>40:raise TimeoutError('40-second integrated branch cap')
    # The endpoint union already covers all lengths and count classes.
    # Each realized endpoint tuple fixes T; take a maximum, not another T union.
    total=max(old['closed_good_route_upper'],-80,worst,-100)+2
    print(dict(kind='integrated_h6_a1_c7_no_full_reset',
               formerly_residual=len(old['residual_live_counts']),
               residual_range=(min(old['residual_live_counts']),max(old['residual_live_counts'])),
               new_count_cells=cells,new_good_route_upper=worst,worst_live=witness,
               old_remaining_good_route_upper=old['closed_good_route_upper'],
               old_route_exception=-80,new_route_exception=-100,
               complete_branch_upper=total,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)
    assert total<=-40


def test():
    checks=multi=allocations=0
    for width in (2,4):
        for length in (1,2,3):
            for count in range(1,min(8,width*(length+1))+1):
                by_weight={}
                for subset in itertools.combinations(range(width*(length+1)),count):
                    initial={j for j in subset if j<width}
                    if not initial or any((j^1) not in initial for j in initial):continue
                    if any(j%width not in initial for j in subset):continue
                    by_weight[len(initial)]=by_weight.get(len(initial),0)+1
                for weight in range(2,width+1,2):
                    assert by_weight.get(weight,0)==choose(width//2,weight//2)*choose(length*weight,count-weight)
                    checks+=1
    for population in range(3,10):
        for active in range(1,population):
            for pieces in (2,3,4):
                correction=2**((pieces-1)*active)*F(population,population-active+1)**((pieces-1)*active)
                for counts in itertools.product(range(active+1),repeat=pieces):
                    total=sum(counts)
                    if total>active:continue
                    ratio=F(math.prod(choose(population,n) for n in counts),choose(population,total))
                    assert ratio<=correction
                    multi+=1
                assert choose(active+pieces-1,pieces-1)<=(active+1)**(pieces-1)
                allocations+=1
    print(dict(single_region_image_checks=checks,multifragment_ratio_checks=multi,
               allocation_count_checks=allocations),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--test',action='store_true')
    p.add_argument('--screen',action='store_true')
    p.add_argument('--maximum-limit',type=int,default=8)
    p.add_argument('--legacy-branch',action='store_true')
    p.add_argument('--integrate',action='store_true')
    p.add_argument('--live',type=int,default=1938)
    args=p.parse_args()
    assert 1<=args.maximum_limit<=8
    if args.test:test()
    elif args.legacy_branch:legacy_branch()
    elif args.integrate:integrate()
    else:run(live=args.live,screen=args.screen,maximum_limit=args.maximum_limit)
