"""Share the internal-defect budget between early episodes and full growth.

Fixed three-shear model. Old localized and episode verifiers are preserved.
"""
import argparse
import itertools
import math
import time
from fractions import Fraction as F
from functools import lru_cache

from diagnose_pair_defects import profile
from diagnose_pair_region_rank import early_joint,joint_cell,count_cap
from diagnose_pair_full_reset_frontier import localized_cap

NEG=-(1<<50)


@lru_cache(maxsize=2048)
def restricted_tail(mean,partial,threshold,minimum_full,maximum_regions):
    worst=F(0); witness=None
    for j in range(partial+1):
        lo=max(minimum_full,math.ceil(mean-j))
        hi=min(maximum_regions-j,math.floor(mean))
        for complete in range(lo,hi+1):
            remainder=mean-complete
            if j==0:
                if remainder:continue
                value=F(int(complete>=threshold)); order=0
            else:
                if not F(j,16)<=remainder<=F(15*j,16):continue
                need=threshold-complete
                if need<=0:value,order=F(1),0
                elif need>j:value,order=F(0),0
                else:
                    value,order=min([(F(1),0)]+[
                        (F(math.comb(j,k),math.comb(need,k))*(remainder/j)**k,k)
                        for k in range(1,need+1)])
            if value>worst:worst,witness=value,(j,complete,order)
    return worst,witness


@lru_cache(maxsize=2048)
def terminal_cap(h,initial,episodes,live,early,early_regions):
    bits=26*(h+2)+(2*(h+2)-1).bit_length()+12+80
    batches=live-early+1
    full,rem=divmod(batches,16)
    regions=early_regions+full+int(bool(rem))
    mean=F(live+episodes,16)
    partial=2*episodes-initial-1
    beta=F(0);witness=None
    # A terminal batch is in a completely included region, or the initial
    # terminal region containing at least rem/16 of the total support set.
    for slots in range(rem if rem else 16,17):
        rest,case=restricted_tail(mean-F(slots,16),partial-int(slots<16),128,
                                  max(0,full-int(slots==16)),regions-1)
        if rest/16>beta:beta,witness=rest/16,(slots,case)
    return dict(cap=count_cap(beta,bits),witness=witness)


@lru_cache(maxsize=2048)
def split_terminal_caps(h,initial,episodes,live,early,early_regions,coarse_endpoints=False):
    # One extra bit pays for the two cap events at each endpoint tuple/batch.
    if coarse_endpoints:
        bits=26*(h+2)+(2*(h+2)-1).bit_length()+12+80+1
    else:
        # Terminal end is fixed, as is the first start when initial=1.
        # The endpoint tuple determines lengths and their minimum defect cost.
        # Union c,a classes, chosen batches, and the two cap events.
        bits=13*(2*episodes-initial-1)+(2*(h+1)-1).bit_length()+12+80+1
    full,rem=divmod(live-early+1,16)
    regions=early_regions+full+int(bool(rem))
    mean=F(live+episodes,16)
    partial=2*episodes-initial-1
    full_beta,_=restricted_tail(mean-1,partial,128,max(0,full-1),regions-1)
    first_beta=F(0)
    if rem:
        for slots in range(rem,17):
            # The initial terminal region is distinct from all `full`
            # guaranteed complete terminal regions, even if earlier inputs
            # fill it to occupancy one in I.
            rest,_=restricted_tail(mean-F(slots,16),partial-int(slots<16),128,
                                   full,regions-1)
            first_beta=max(first_beta,rest/16)
    return dict(first=count_cap(first_beta,bits) if rem else 0,
                full=count_cap(full_beta/16,bits),remainder=rem)


def ordered_growth(extra,first,full,remainder):
    if remainder==0:return ((1<<(extra+1))-1)*full
    # Activation uses one of the `remainder` initial-region batches.
    # All initial-region defects precede all complete-region defects.
    return max((1<<(extra-k))*((1<<(k+1))-1)*first
               +((1<<(extra-k))-1)*full
               for k in range(min(extra,remainder-1)+1))


def growth_required(cap,width=512):
    if cap==0:return width
    needed=0
    while ((1<<(needed+1))-1)*cap<width:needed+=1
    return needed


def selected(h=8,initial=1,episodes=7,live=1955,report=True,terminal=False,split=False,coarse_endpoints=False):
    began=time.monotonic()
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    old_bits=placement+9+80
    local_bits=26*(h+2)+(2*(h+2)-1).bit_length()+12+80
    base=profile(h,initial,episodes,live,460,True)
    assert base is not None and live>base['partial_cap']
    cap=localized_cap(F(live+episodes,16),2*episodes-initial-1,local_bits)['cap']
    extra=h-episodes+initial
    required=growth_required(cap)
    old=new=NEG; old_witness=new_witness=None; discarded=0
    for early,used,regions in early_joint(h,initial,episodes):
        if time.monotonic()-began>30:raise TimeoutError('selected growth 30s cap')
        cell=joint_cell(h,initial,episodes,live,early,used,regions,old_bits)
        if cell is None:continue
        value=base['log2_upper'] if 'unresolved' in cell else min(base['log2_upper'],cell['upper'])
        if value>old:old,old_witness=value,(early,used,regions)
        local_cap=terminal_cap(h,initial,episodes,live,early,regions)['cap'] if terminal else cap
        needed=growth_required(local_cap)
        feasible=used+needed<=extra
        split_caps=None
        if split:
            split_caps=split_terminal_caps(h,initial,episodes,live,early,regions,coarse_endpoints)
            feasible=ordered_growth(extra-used,split_caps['first'],split_caps['full'],
                                    split_caps['remainder'])>=512
        if not feasible:
            discarded+=1;continue
        if value>new:new,new_witness=value,(early,used,regions,local_cap,needed,split_caps)
        if time.monotonic()-began>30:raise TimeoutError('selected growth 30s cap')
    result=dict(profile=(h,initial,episodes,live),batch_cap=cap,internal_budget=extra,
                growth_required=required,old_upper=old+placement,new_upper=new+placement,
                terminal_specific=terminal,split_caps=split,
                old_witness=old_witness,new_witness=new_witness,discarded_cells=discarded,
                whole_code_certificate=False,seconds=time.monotonic()-began)
    if report:print(result,flush=True)
    return result


def frontier(terminal=False,split=False,coarse_endpoints=False):
    """Recheck all 45 retained h8 residual T values, without global claims."""
    began=time.monotonic(); results=[]
    for initial,episodes,lo,hi in ((0,5,1950,1951),(1,5,1940,1960),
                                  (1,6,1943,1959),(1,7,1955,1959)):
        for live in range(lo,hi+1):
            result=selected(8,initial,episodes,live,False,terminal,split,coarse_endpoints)
            results.append(result)
            if time.monotonic()-began>45:raise TimeoutError('frontier 45s cap')
    residual=[r for r in results if r['new_upper']>-80]
    print(dict(kind='retained_h8_frontier',checked=len(results),residual=len(residual),
               worst=max(results,key=lambda r:r['new_upper']),
               residual_examples=residual[:4],whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


def cover(h=8):
    began=time.monotonic(); checked=old_residual=local_rescued=new_rescued=0
    worst=NEG; witness=None; residual=[]
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    old_bits=placement+9+80
    local_bits=26*(h+2)+(2*(h+2)-1).bit_length()+12+80
    for initial in (0,1):
        for episodes in range(1,h+initial+1):
            for live in range(episodes,4096):
                base=profile(h,initial,episodes,live,460,True)
                checked+=1
                if base is None:continue
                value=base['log2_upper']
                if value+placement>-80:
                    local=NEG
                    for early,used,regions in early_joint(h,initial,episodes):
                        cell=joint_cell(h,initial,episodes,live,early,used,regions,old_bits)
                        if cell is not None:
                            local=max(local,value if 'unresolved' in cell else min(value,cell['upper']))
                        if time.monotonic()-began>60:raise TimeoutError('complete cover 60s limit')
                    value=local
                if value+placement>-80:
                    old_residual+=1
                    assert live>base['partial_cap']
                    extra=h-episodes+initial
                    cap=localized_cap(F(live+episodes,16),2*episodes-initial-1,local_bits)['cap']
                    if ((1<<(extra+1))-1)*cap<512:
                        local_rescued+=1;value=NEG
                    else:
                        result=selected(h,initial,episodes,live,False,split=True)
                        value=result['new_upper']-placement
                        if value+placement<=-80:new_rescued+=1
                        else:residual.append(result)
                if value>worst:worst,witness=value,(initial,episodes,live)
                if time.monotonic()-began>60:raise TimeoutError('complete cover 60s limit')
    upper=max(worst+placement,-80)+2
    print(dict(kind='coupled_growth_no_full_reset',h=h,checked=checked,
               old_residual=old_residual,old_local_rescued=local_rescued,
               new_rescued=new_rescued,residual=len(residual),examples=residual[:3],
               good_route_upper=worst+placement,worst=witness,
               route_exceptions=(-80,-80,-80),family_upper=upper,closed=upper<=-40,
               whole_code_certificate=False,seconds=time.monotonic()-began),flush=True)


def test():
    began=time.monotonic(); geometry=growth=moment=0
    batches=12; region=4
    for mask in range(1,1<<(batches-1)):
        if not (mask>>(batches-2)&1):continue
        active=[i for i in range(1,batches) if mask>>(i-1)&1]
        intervals=[]
        for i in active:
            if not intervals or i>intervals[-1][1]+1:intervals.append([i,i])
            else:intervals[-1][1]=i
        inputs=set().union(*(set(range(lo-1,hi+1)) for lo,hi in intervals))
        terminal=set(range(intervals[-1][0]-1,batches))
        full,rem=divmod(len(terminal),region)
        counts=[len(inputs.intersection(range(lo,lo+region))) for lo in range(0,batches,region)]
        early_regions=sum((hi//region)-((lo-1)//region)+1 for lo,hi in intervals[:-1])
        initial=int(active[0]==1); episodes=len(intervals)
        assert len(inputs)==len(active)+episodes
        assert sum(n==region for n in counts)>=full
        assert sum(n>0 for n in counts)<=early_regions+full+bool(rem)
        assert sum(0<n<region for n in counts)<=2*episodes-initial-1
        for batch in terminal:
            r=batch//region; slots=counts[r]
            is_first=rem and batch<(min(terminal)+rem)
            rest=counts[:r]+counts[r+1:]
            assert sum(n==region for n in rest)>=full-int(not is_first)
            assert slots>=rem if is_first else slots==region
            geometry+=1
    for extra in range(5):
        for first,full in itertools.product(range(5),repeat=2):
            for rem in range(1,6):
                upper=ordered_growth(extra,first,full,rem)
                for k in range(min(extra,rem-1)+1):
                    weight=first
                    for _ in range(k):weight=2*weight+first
                    for _ in range(extra-k):weight=2*weight+full
                    assert weight<=upper
                    growth+=1
    for sizes in itertools.product(range(5),repeat=3):
        ps=[F(n,4) for n in sizes]
        complete=sum(n==4 for n in sizes)
        partial=sum(0<n<4 for n in sizes)
        touched=sum(n>0 for n in sizes)
        for threshold in (1,2,3):
            beta,_=restricted_tail(sum(ps,F(0)),partial,threshold,complete,touched)
            # Exact one-row tail under independent regional membership.
            exact=sum((math.prod(ps[i] if bits[i] else 1-ps[i] for i in range(3))
                       for bits in itertools.product((0,1),repeat=3) if sum(bits)>=threshold),F(0))
            assert exact<=beta
            moment+=1
    print(dict(terminal_geometry_checks=geometry,ordered_growth_checks=growth,
               restricted_moment_checks=moment,seconds=time.monotonic()-began),flush=True)
    assert time.monotonic()-began<20


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--frontier',action='store_true')
    p.add_argument('--cover',action='store_true')
    p.add_argument('--test',action='store_true')
    p.add_argument('--defects',type=int,default=8)
    p.add_argument('--terminal',action='store_true')
    p.add_argument('--split',action='store_true')
    p.add_argument('--coarse-endpoints',action='store_true')
    p.add_argument('--initial',type=int,default=1)
    p.add_argument('--episodes',type=int,default=7)
    p.add_argument('--live',type=int,default=1955)
    args=p.parse_args()
    assert 0<=args.defects<=16
    if args.test:test()
    elif args.cover:cover(args.defects)
    elif args.frontier:
        assert args.defects==8
        frontier(args.terminal,args.split,args.coarse_endpoints)
    else:selected(h=args.defects,initial=args.initial,episodes=args.episodes,live=args.live,
                  terminal=args.terminal,split=args.split,coarse_endpoints=args.coarse_endpoints)
