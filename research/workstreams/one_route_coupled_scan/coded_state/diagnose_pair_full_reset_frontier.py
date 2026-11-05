"""Map/refine the fixed SL2 full-reset frontier; preserve old verifiers.

All bounds use integer/rational arithmetic. Bounded serial jobs, no files.
"""
import argparse
import itertools
import math
import time
from fractions import Fraction as F
from diagnose_pair_defects import profile,volume_bound,ceildiv
from diagnose_pair_region_rank import rank_cap,input_cap,dimension_cap,row_tail,count_cap
from diagnose_pair_defect_envelope import partner

B=512
N=1<<21
NEG=-(1<<50)


def localized_cap(mean,partial,bits):
    """Eligible rows whose routed coordinate lies in one specified batch.

    Condition only by multiplying the joint batch membership moment; routing
    in all other regions stays independent. Maximize the containing region's
    allowed fraction before selecting the row cap.
    """
    beta=F(0); witness=None
    for slots in range(1,17):
        rest_partial=partial-int(slots<16)
        if rest_partial<0 or mean<F(slots,16):continue
        rest,case=row_tail(mean-F(slots,16),rest_partial,128)
        candidate=rest/16
        if candidate>beta:beta,witness=candidate,(slots,case)
    return dict(cap=count_cap(beta,bits),beta=beta,witness=witness)


def local_selected(h=7,initial=1,episodes=7,live=1968):
    began=time.monotonic()
    bits=26*(h+2)+(2*(h+2)-1).bit_length()+12+80
    result=localized_cap(F(live+episodes+1,16),2*episodes+1-initial,bits)
    extra=h-episodes+initial
    state_cap=((1<<(extra+1))-1)*result['cap']
    print(dict(h=h,a=initial,c=episodes,T=live,eligible_per_batch=result['cap'],
               witness=result['witness'],extra=extra,state_cap=state_cap,
               excludes_full_state=state_cap<B,route_exception=-80,
               whole_code_certificate=False,seconds=time.monotonic()-began),flush=True)


def legacy_map(h=6,report=True):
    began=time.monotonic(); residual=[]; closed=NEG; checked=0
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    bits=26*(h+2)+(2*(h+2)-1).bit_length()+9+80
    for initial in (0,1):
        for episodes in range(1,h+initial+1):
            for live in range(episodes,4096):
                base=profile(h,initial,episodes,live,460)
                if base is None:continue
                value=base['log2_upper']; checked+=1
                if value+placement>-80:
                    full=live-base['partial_cap']
                    rank=rank_cap(full//8*8,h,bits)['rank_lower']
                    padded=min(4096,ceildiv(live+episodes+1,8)*8)
                    dimension=min(base['input_dimension'],input_cap(padded,episodes+1,bits))
                    inverse=1+volume_bound(base['slots'],460*N//1000)-32*rank
                    reset=1+32*dimension+B*12-32*B
                    value=min(value,inverse,reset)
                    if value+placement>-80:
                        extra=h-episodes+initial
                        window=9*(extra+1)
                        new_reset=1+32*dimension+12+(window**B-1).bit_length()-32*B
                        residual.append(dict(a=initial,c=episodes,T=live,extra=extra,
                                             input_dimension=dimension,old=value+placement,
                                             decay_window=window,
                                             new=min(value,new_reset)+placement))
                    else:closed=max(closed,value+placement)
                else:closed=max(closed,value+placement)
                if time.monotonic()-began>45:raise TimeoutError('full-reset map 45s limit')
    result=dict(h=h,checked=checked,old_residual=residual,old_closed_good_route=closed,
                decay_survivors=sum(cell['new']>-80 for cell in residual),
                whole_code_certificate=False,seconds=time.monotonic()-began)
    if report:print(result,flush=True)
    return result


def integrate(h=6,detail=False):
    began=time.monotonic()
    old=legacy_map(h,False)
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    bits=26*(h+2)+(2*(h+2)-1).bit_length()+9+80
    local_bits=26*(h+2)+(2*(h+2)-1).bit_length()+12+80
    refined=[]; good=old['old_closed_good_route']
    for cell in old['old_residual']:
        bound=cell['new']; kind='decay_window'
        if bound>-80:
            assert cell['T']+cell['c']+1<4096, 'exact padding requires an available gap'
            mass=F(cell['T']+cell['c']+1,16)
            partial=2*cell['c']+1-cell['a']
            dimension,layers=dimension_cap(mass,partial,bits)
            eligible=layers[0][1] if layers else 0
            cell=dict(cell,aligned_dimension=dimension,eligible_rows=eligible)
            if cell['extra']==0 and eligible<B:
                bound=NEG; kind='full_activation_excluded'
            else:
                bound=min(bound,1+32*dimension+12
                          +(cell['decay_window']**B-1).bit_length()-32*B+placement)
                kind='aligned_dimension_and_decay'
            if bound>-80:
                local=localized_cap(mass,partial,local_bits)
                state_cap=((1<<(cell['extra']+1))-1)*local['cap']
                cell=dict(cell,eligible_per_batch=local['cap'],state_cap=state_cap)
                if state_cap<B:
                    bound=NEG;kind='localized_full_state_excluded'
        good=max(good,bound)
        refined.append(dict(a=cell['a'],c=cell['c'],T=cell['T'],kind=kind,upper=bound,
                            aligned_dimension=cell.get('aligned_dimension'),
                            eligible_rows=cell.get('eligible_rows'),
                            eligible_per_batch=cell.get('eligible_per_batch'),
                            state_cap=cell.get('state_cap')))
        if time.monotonic()-began>45:raise TimeoutError('integration 45s limit')
    upper=max(good,-80)+2  # good-route, old route, exact-alignment route, local route
    survivors=[cell for cell in refined if cell['upper']>-80]
    groups={}
    for cell in survivors:
        key=(cell['a'],cell['c'])
        group=groups.setdefault(key,dict(count=0,first_T=cell['T'],last_T=cell['T'],worst=NEG))
        group['count']+=1; group['last_T']=cell['T'];group['worst']=max(group['worst'],cell['upper'])
    print(dict(kind='full_reset_refinement',h=h,
               old_residual_count=len(refined),residual_count=len(survivors),
               residual_groups=groups,refined=refined if detail else None,
               good_route_upper=good,route_exceptions=(-80,-80,-80),
               complete_family_upper=upper,closed=upper<=-40,
               combined_with_no_full_reset_upper=max(upper,-69)+1 if h==6 else None,
               whole_code_certificate=False,seconds=time.monotonic()-began),flush=True)


def test():
    began=time.monotonic(); geometry=decays=localized=0
    # Exhaust every incoming support sequence on four coordinates, then use
    # the last full state preceding the first subsequent zero.
    width=4; full=(1<<width)-1; phases=2
    for length in range(1,6):
        for middle in itertools.product(range(1,full+1),repeat=length-1):
            history=(full,)+middle+(0,)
            last=max(i for i,mask in enumerate(history[:-1]) if mask==full)
            extra=0
            for t in range(last,length):
                current,nxt=history[t:t+2]
                live=current|partner(current,t%phases,width)
                extra+=int(live!=current or bool(nxt&~live))
            assert length-last<=phases*(extra+1)
            decays+=1
        if time.monotonic()-began>20:raise TimeoutError('tests 20s limit')
    # Exact input set for each binary incoming-state pattern, including the
    # unconstrained terminal input. The first/last anchoring saves boundaries.
    batches=12; region=4
    for mask in range(1,1<<(batches-1)):
        states=[bool(mask>>(i-1)&1) for i in range(1,batches)]
        initial=int(states[0]); episodes=sum(x and (i==0 or not states[i-1]) for i,x in enumerate(states))
        live=sum(states)
        inputs={batches-1}
        for i,active in enumerate(states,1):
            if active:inputs.update((i-1,i))
        partial=sum(0<len(inputs.intersection(range(start,start+region)))<region
                    for start in range(0,batches,region))
        assert len(inputs)<=live+episodes+1
        assert partial<=2*episodes+1-initial
        target=live+episodes+1
        if len(inputs)<target<batches:
            assert target-len(inputs)==1
            candidates=[inputs|{j} for j in range(batches) if j not in inputs]
            assert any(sum(0<len(padded.intersection(range(start,start+region)))<region
                           for start in range(0,batches,region))<=2*episodes+1-initial
                       for padded in candidates)
        geometry+=1
    population=4
    for batch_size in (1,2):
        batch_masks=[mask for mask in range(1<<population) if mask.bit_count()==batch_size]
        for sizes in itertools.product(range(population+1),repeat=2):
            ps=[F(size,population) for size in sizes]
            masks=[[mask for mask in range(1<<population) if mask.bit_count()==size] for size in sizes]
            outcomes=list(itertools.product(batch_masks,*masks))
            for need in (1,2):
                beta=min([F(1)]+[sum((math.prod(ps[i] for i in subset)
                                     for subset in itertools.combinations(range(2),k)),F(0))
                                  /math.comb(need,k) for k in range(1,need+1)])
                beta*=F(batch_size,population)
                for rows in range(1,population+1):
                    actual=F(sum(all((batch>>r)&1 and sum((mask>>r)&1 for mask in regions)>=need
                                     for r in range(rows)) for batch,*regions in outcomes),len(outcomes))
                    assert actual<=beta**rows
                    localized+=1
    print(dict(last_full_decay_checks=decays,anchored_input_geometry_checks=geometry,
               exact_localized_row_moment_checks=localized,
               seconds=time.monotonic()-began),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--defects',type=int,default=6)
    p.add_argument('--integrate',action='store_true')
    p.add_argument('--test',action='store_true')
    p.add_argument('--detail',action='store_true')
    p.add_argument('--local-cap',action='store_true')
    p.add_argument('--initial',type=int,default=1)
    p.add_argument('--episodes',type=int,default=7)
    p.add_argument('--live',type=int,default=1968)
    args=p.parse_args()
    assert 0<=args.defects<=16
    if args.test:test()
    elif args.local_cap:local_selected(args.defects,args.initial,args.episodes,args.live)
    elif args.integrate:integrate(args.defects,args.detail)
    else:legacy_map(args.defects)
