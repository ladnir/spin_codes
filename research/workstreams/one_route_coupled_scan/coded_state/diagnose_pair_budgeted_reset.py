"""Input-capped proper states and full runs share the reset defect budget.

Fixed SL2 construction; bounded integer/rational diagnostics, no output files.
"""
import argparse
import time
from fractions import Fraction as F

from diagnose_pair_defects import profile,volume_bound
from diagnose_pair_episodes import volume_bits
from diagnose_pair_region_rank import dimension_cap,rank_cap
from diagnose_pair_full_reset_frontier import localized_cap,legacy_map
from diagnose_pair_growth_budget import growth_required
from diagnose_pair_defect_envelope import partner

NEG=-(1<<50)


def growth_precursors(cap,width=512):
    """Pair-image counts at preceding defects, latest first, before first full.

    Clean steps cannot grow support. To attain R coordinates at a defective
    update, at least R-cap coordinates must come from live pairs; their count
    is even. The incoming state has at least ceil((R-cap)/2) coordinates.
    Repeat while this exceeds the activation cap. No distribution assumption.
    """
    assert 0<cap<width
    required=width; result=[]
    while required>cap:
        prior=(required-cap+1)//2
        result.append(2*prior)
        required=prior
    return tuple(result)


def selected(h=8,initial=1,episodes=8,live=1984,report=True,
             exact=True,precursors=True,prior_upper=None):
    began=time.monotonic()
    base=profile(h,initial,episodes,live,460)
    assert base is not None
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    # State-interval endpoints, first start fixed iff initial=1. The free
    # terminal input adds no variable endpoint. Two types of new route event.
    endpoints=13*(2*episodes-initial)+(2*(h+1)-1).bit_length()
    local_bits=endpoints+12+1+80
    dimension_bits=endpoints+7+1+80
    rank_bits=26*(h+2)+(2*(h+2)-1).bit_length()+9+80
    mean=F(live+episodes+1,16)
    partial=2*episodes+1-initial
    cap=localized_cap(mean,partial,local_bits)['cap']
    dimension=dimension_cap(mean,partial,dimension_bits)[0]
    extra=h-episodes+initial
    required=growth_required(cap)
    result=dict(h=h,a=initial,c=episodes,T=live,batch_cap=cap,dimension=dimension,
                internal_budget=extra,growth_required=required)
    if required>extra:
        result.update(excluded=True,upper=NEG)
    elif cap<512:
        clean_length=cap.bit_length()-1
        proper=min(live,episodes*clean_length+9*extra)
        full=live-proper
        assert extra>=1
        regions=min(256,(full+30*extra)//16)
        rank=(full*512 if regions<=128 else
              rank_cap(full,extra-1,rank_bits)['rank_lower'])
        weight=460*(1<<21)//1000
        inverse=1+volume_bound(base['slots'],weight)-32*rank
        decay=9*(extra-required+1)
        reset=1+32*dimension+12+(decay**512-1).bit_length()-32*512
        coarse_upper=min(base['log2_upper'],inverse,reset)+placement
        growth=growth_precursors(cap)
        assert len(growth)==required
        # Before the first full state these g incoming batches are proper,
        # hence disjoint from the selected full batches. Each adds at most
        # one touched region. On <=128 regions RS parity columns are independent.
        bonus=sum(growth) if precursors and regions+len(growth)<=128 else 0
        if bonus:
            rank=full*512+bonus
            inverse=1+volume_bound(base['slots'],weight)-32*rank
        used_exact=False
        if (exact and (prior_upper is None or prior_upper>-80)
                and min(base['log2_upper'],inverse,reset)+placement>-80):
            inverse=1+volume_bits(base['slots'],weight)-32*rank
            used_exact=True
        result.update(excluded=False,proper=proper,full=full,full_run_cap=extra,
                      rank=rank,full_region_cap=regions,decay_window=decay,
                      coarse_upper=coarse_upper,growth_precursors=growth,
                      precursor_rank_bonus=bonus,exact_volume=used_exact,
                      inverse=inverse+placement,reset=reset+placement,
                      upper=min(base['log2_upper'],inverse,reset)+placement)
    else:
        # No stronger state-count assertion is justified without a proper
        # activation. Keep the previously certified direct bound.
        result.update(excluded=False,upper=base['log2_upper']+placement)
    result.update(whole_code_certificate=False,seconds=time.monotonic()-began)
    if report:print(result,flush=True)
    return result


def test():
    began=time.monotonic(); checks=0
    width=4; full=15; d=2
    for cap in (1,2,3):
        # Values retain the largest number of proper incoming states for
        # each (current support, live episodes, internal defects, full runs).
        states={(mask,1,0,0):0 for mask in range(1,full) if mask.bit_count()<=cap}
        for t in range(7):
            nxt_states={}
            for (current,episodes,extra,runs),proper in states.items():
                live=current|partner(current,t%d,width)
                for nxt in range(full+1):
                    if (nxt&~live).bit_count()>cap:continue
                    if current==0 and nxt.bit_count()>cap:continue
                    new_episodes=episodes+int(current==0 and nxt!=0)
                    defective=int(current!=0 and (live!=current or bool(nxt&~live)))
                    new_extra=extra+defective
                    new_runs=runs+int(current!=full and nxt==full)
                    if new_extra>3 or new_episodes>3:continue
                    new_proper=proper+int(0<current<full)
                    assert new_runs<=new_extra
                    assert new_proper<=new_episodes*(cap.bit_length()-1)+d*new_extra
                    key=(nxt,new_episodes,new_extra,new_runs)
                    nxt_states[key]=max(nxt_states.get(key,-1),new_proper)
                    checks+=1
            states=nxt_states
            if time.monotonic()-began>20:raise TimeoutError('state-count test 20s cap')
    print(dict(input_capped_state_transition_checks=checks,seconds=time.monotonic()-began),flush=True)
    checks=0
    for cap in (1,2,3):
        needed=growth_precursors(cap,width)
        assert len(needed)==growth_required(cap,width)
        states={(mask,()) for mask in range(1,full) if mask.bit_count()<=cap}
        for t in range(7):
            nxt_states=set()
            for current,recent in states:
                live=current|partner(current,t%d,width)
                for nxt in range(1,full+1):
                    if (nxt&~live).bit_count()>cap:continue
                    defective=live!=current or bool(nxt&~live)
                    tail=(recent+(live.bit_count(),))[-len(needed):] if defective else recent
                    if nxt==full:
                        assert len(tail)==len(needed)
                        assert all(x>=y for x,y in zip(reversed(tail),needed))
                        checks+=1
                    else:nxt_states.add((nxt,tail))
            states=nxt_states
            if time.monotonic()-began>20:raise TimeoutError('growth precursor test 20s cap')
    print(dict(first_full_precursor_checks=checks,seconds=time.monotonic()-began),flush=True)


def integrate(h=8,exact=True,precursors=True):
    began=time.monotonic()
    old=legacy_map(h,False)
    worst=old['old_closed_good_route']; witness=None; residual=[]; improved=0
    exact_cells=0; old_window_closed=0
    for cell in old['old_residual']:
        if cell['new']<=-80:
            improved+=1; old_window_closed+=1
            if cell['new']>worst:worst,witness=cell['new'],cell
            continue
        new=selected(h,cell['a'],cell['c'],cell['T'],False,exact,precursors,cell['new'])
        exact_cells+=int(new.get('exact_volume',False))
        value=min(cell['new'],new['upper'])
        if value<=-80:improved+=1
        else:residual.append(new)
        if value>worst:worst,witness=value,new
        if time.monotonic()-began>60:raise TimeoutError('budgeted reset integration 60s cap')
    upper=max(worst,-80)+2
    print(dict(kind='budgeted_full_reset',h=h,old_residual=len(old['old_residual']),
               covered=improved,residual=len(residual),examples=residual[:4],
               old_window_closed=old_window_closed,exact_cells=exact_cells,
               good_route_upper=worst,worst=witness,
               route_exceptions=(-80,-80,-80),family_upper=upper,closed=upper<=-40,
               whole_code_certificate=False,seconds=time.monotonic()-began),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--test',action='store_true')
    p.add_argument('--integrate',action='store_true')
    p.add_argument('--coarse',action='store_true')
    p.add_argument('--no-precursors',action='store_true')
    p.add_argument('--defects',type=int,default=8)
    p.add_argument('--initial',type=int,default=1)
    p.add_argument('--episodes',type=int,default=8)
    p.add_argument('--live',type=int,default=1984)
    args=p.parse_args()
    if args.test:test()
    elif args.integrate:integrate(args.defects,not args.coarse,not args.no_precursors)
    else:selected(args.defects,args.initial,args.episodes,args.live,
                  exact=not args.coarse,precursors=not args.no_precursors)
