"""Single defective batch inside a closed episode, fixed SL2 construction.

Bounded combinatorial checks and actual a1/c6 frontier geometry. No files.
"""
import argparse
import itertools
import time
from fractions import Fraction as F
from functools import lru_cache

from diagnose_pair_activation_route import B,M,Q,NEG,choose,logceil_ratio
from diagnose_pair_defect_envelope import partner,statistics,initial_masks
from diagnose_pair_region_rank import count_cap
from diagnose_pair_aligned_frontier import local_envelope
from diagnose_pair_shared_regions import sharing_price

BANDS=((1,1),(2,2),(3,4),(5,8),(9,16),(17,32),
       (33,64),(65,128),(129,256),(257,512))


def one_batch_histories(width,length,weight,position,extra):
    remaining=length-position-1
    return min(1 << (width*(length-1)),
               length**weight*sum(choose(width,b)*remaining**b
                                  for b in range(min(extra,width)+1)))


def test():
    began=time.monotonic(); checked=counts=denominators=vandermonde=0
    assert BANDS[0][0]==1 and BANDS[-1][1]==B
    assert all(high+1==low for (_,high),(low,_) in zip(BANDS,BANDS[1:]))
    for width,length in ((2,1),(4,1),(4,2),(4,3),(8,1),(8,2),(8,3),(8,4),(8,5)):
        phases=(width-1).bit_length()
        full=(1<<width)-1
        totals={}
        for position in range(length):
            allowed=[tuple(mask for mask in range(1,full)
                           if t==position or partner(mask,t%phases,width)==mask)
                     for t in range(length)]+[(0,)]

            def histories(prefix):
                t=len(prefix)-1
                current=prefix[-1]
                live=current|partner(current,t%phases,width)
                for nxt in allowed[t+1]:
                    bad=(live!=current or bool(nxt&~live))
                    if bad!=(t==position):continue
                    extended=prefix+(nxt,)
                    if t+1==length:yield extended
                    else:yield from histories(extended)

            for history in itertools.chain.from_iterable(histories((initial,)) for initial in allowed[0]):
                initial=history[0]
                remaining=length-position-1
                item=statistics(history,width)
                assert position<phases and remaining<phases
                assert initial.bit_count()>=1<<position
                if remaining:
                    assert history[position+1].bit_count()>=1<<remaining
                assert item['initial']+item['extra']>=1<<remaining
                key=(initial,position,item['extra'])
                totals[key]=totals.get(key,0)+1
                checked+=1
                if checked%1000==0 and time.monotonic()-began>20:
                    raise TimeoutError('test 20s limit')
        for initial,position in {(u,j) for u,j,e in totals}:
            cumulative=0
            for extra in range(1,width+1):
                cumulative+=totals.get((initial,position,extra),0)
                assert cumulative<=one_batch_histories(width,length,initial.bit_count(),position,extra)
                counts+=1
        if time.monotonic()-began>20:raise TimeoutError('test 20s limit')
    for population in range(3,11):
        for active in range(1,min(5,population-1)+1):
            for total in range(1,2*active+1):
                for weight in range(1,min(active,total)+1):
                    lo=max(weight,total-active); hi=min(active,total)
                    lower=min(choose(population,n)*choose(population,total-n) for n in (lo,hi))
                    for nleft in range(lo,hi+1):
                        assert lower<=choose(population,nleft)*choose(population,total-nleft)
                        denominators+=1
                        for length in range(1,5):
                            for left in range(1,length+1):
                                for size in range(weight,weight+3):
                                    actual=choose((left-1)*size,nleft-weight)*choose((length+1-left)*size,total-nleft)
                                    assert actual<=choose(length*size,total-weight)
                                    vandermonde+=1
    print(dict(single_defect_support_checks=checked,history_count_checks=counts,
               denominator_endpoint_checks=denominators,
               vandermonde_allocation_checks=vandermonde,
               seconds=time.monotonic()-began),flush=True)


@lru_cache(maxsize=8192)
def cap_for(early,fragments,tail):
    need=129-tail
    if fragments<need:return 0
    assert need>0
    mean=F(early+5,16*fragments)
    beta=min([F(1)]+[F(choose(fragments,j),choose(need,j))*mean**j
                     for j in range(1,need+1)])
    # Endpoints identify all fragments, not the internally defective episode.
    # Pay for all five designations here as a conservative route union.
    return count_cap(beta,259)


def cells(live):
    for maximum in range(1,9):
        for clean_length in range(maximum+3,4*maximum+1):
            for length in range(1,18):
                early=clean_length+length
                tail=(live-early+16)//16
                for first in (False,True):
                    fs=((length+16)//16,) if first else range(1,4)
                    for fragments in fs:
                        if not first and not (fragments<=length+1<=16*fragments
                                              and (fragments<3 or length==17)):
                            continue
                        for crossings in range(5 if first else 4):
                            m=4+crossings+fragments
                            dimension=m+tail-128
                            if dimension<=0:continue
                            cap=cap_for(early,m,tail)
                            minimum=max(1<<maximum,1<<max(0,length-9))
                            for sharing in range(dimension):
                                yield dict(maximum=maximum,clean_length=clean_length,
                                           length=length,first=first,fragments=fragments,
                                           crossings=crossings,virtual_dimension=dimension,
                                           sharing=sharing,cap=cap,minimum=minimum,
                                           tail=tail)


def screen(live):
    began=time.monotonic(); total=excluded=0; largest=0; survivors={}
    for cell in cells(live):
        total+=1; largest=max(largest,cell['cap'])
        if cell['cap']<cell['minimum']:excluded+=1;continue
        key=cell['length']
        value=(cell['virtual_dimension'],cell['cap'],cell['maximum'])
        survivors[key]=max(survivors.get(key,(0,0,0)),value)
        if time.monotonic()-began>30:raise TimeoutError('screen 30s limit')
    print(dict(live=live,cells=total,excluded=excluded,largest_cap=largest,
               survivors_by_length=survivors,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


@lru_cache(maxsize=2048)
def birth_descriptions(length,position,high):
    remaining=length-position-1
    return sum(choose(B,b)*remaining**b for b in range(min(high,B)+1))


@lru_cache(maxsize=8192)
def history_coefficient(length,weight,high):
    result=0
    for position in range(max(0,length-9),min(8,length-1)+1):
        remaining=length-position-1
        if weight<1<<position or weight+high<1<<remaining:continue
        masks=(initial_masks(B,weight,high) if position==0 else
               choose(B//2,weight//2) if weight%2==0 else 0)
        histories=min(1 << (B*(length-1)),
                      length**weight*birth_descriptions(length,position,high))
        result+=masks*histories
    return result


@lru_cache(maxsize=1024)
def defective_table(active,length,fragments):
    """Vandermonde envelope; sum all charges and exceptional positions.

    For two fragments, relax the allocation after the activation batch.
    The denominator minimum is attained at an endpoint by log-concavity.
    """
    assert fragments in (1,2)
    out=[NEG]*(fragments*active+1)
    for total in range(1,fragments*active+1):
        bands=[]
        for low,high in BANDS:
            terms=[]
            for weight in range(1,min(active,total)+1):
                coefficient=history_coefficient(length,weight,high)
                if not coefficient:continue
                numerator=coefficient*choose(length*min(B,weight+high),total-weight)
                if not numerator:continue
                if fragments==1:
                    denominator=choose(M,total)
                else:
                    left_min=max(weight,total-active)
                    left_max=min(active,total)
                    denominator=min(choose(M,left)*choose(M,total-left)
                                    for left in (left_min,left_max))
                terms.append(logceil_ratio(numerator,denominator*Q**(weight+low)))
            if terms:bands.append(max(terms)+(len(terms)-1).bit_length())
        if bands:out[fragments*active-total]=max(bands)+(len(bands)-1).bit_length()
    return tuple(out)


@lru_cache(maxsize=2048)
def combined(active,maximum,crossings,length,fragments,dimension):
    budget=(dimension-1)*active
    dp=[0]+[NEG]*budget
    tables=[local_envelope(active,maximum,index<crossings) for index in range(4)]
    tables.append(defective_table(active,length,fragments))
    for table in tables:
        nxt=[NEG]*(budget+1)
        for used,previous in enumerate(dp):
            if previous==NEG:continue
            for extra,value in enumerate(table[:budget-used+1]):
                if value!=NEG:nxt[used+extra]=max(nxt[used+extra],previous+value)
        dp=nxt
    return tuple(itertools.accumulate(dp,max))


@lru_cache(maxsize=8192)
def source_term(active,maximum,crossings,length,fragments,dimension,sharing):
    actual=dimension-sharing
    envelope=combined(active,maximum,crossings,length,fragments,dimension)
    values=[]
    for zeros in range((actual-1)*active+1):
        route=envelope[zeros+sharing*active]
        if route!=NEG:
            values.append((choose((128+actual)*active,zeros)-1).bit_length()
                          +32*(actual*active-zeros)+route)
    if not values:return NEG
    return ((choose(M,active)-1).bit_length()+max(values)
            +(len(values)-1).bit_length()+sharing_price(active,sharing))


def bound(live,length_limit=17,report=True):
    began=time.monotonic(); worst=NEG; witness=None; total=excluded=0
    for cell in cells(live):
        if cell['length']>length_limit:continue
        total+=1
        if cell['cap']<cell['minimum']:excluded+=1;continue
        assert cell['fragments']<=2, 'three-region survivor requires its own factor'
        terms=[source_term(s,cell['maximum'],cell['crossings'],cell['length'],
                           cell['fragments'],cell['virtual_dimension'],cell['sharing'])
               for s in range(cell['minimum'],cell['cap']+1)]
        # Endpoints, five choices of defective episode, global density,
        # and active-source counts. Positions and charges are in the factor.
        value=max(terms)+156+3+1+(len(terms)-1).bit_length()
        if value>worst:worst,witness=value,cell
        if time.monotonic()-began>40:
            print(dict(incomplete=True,cells=total,worst=worst,witness=witness),flush=True)
            raise TimeoutError('bound 40s limit')
    result=dict(live=live,length_limit=length_limit,cells=total,excluded=excluded,
                good_route_upper=worst,route_exception=-100,
                family_upper=max(worst,-100)+1,witness=witness,
                whole_code_certificate=False,seconds=time.monotonic()-began)
    if report:print(result,flush=True)
    return result


def integrate():
    from diagnose_pair_unaligned_frontier import legacy_branch,run as clean_bound
    began=time.monotonic()
    old=legacy_branch(False,initial=1,episodes=6)
    bad_worst=clean_worst=NEG; bad_witness=clean_witness=None; total=0
    for live in old['residual_live_counts']:
        bad=bound(live,report=False)
        clean=clean_bound(live,report=False)
        total+=bad['cells']
        if bad['good_route_upper']>bad_worst:bad_worst,bad_witness=bad['good_route_upper'],live
        if clean['good_route_upper']>clean_worst:clean_worst,clean_witness=clean['good_route_upper'],live
        if time.monotonic()-began>60:
            print(dict(incomplete=True,last_live=live,bad_worst=bad_worst,
                       clean_worst=clean_worst,cells=total),flush=True)
            raise TimeoutError('integration 60s limit')
    # Old good-route part, old route exception, two new good-route parts,
    # and their two route exceptions. Six contributions cost <=3 bits.
    upper=max(old['closed_good_route_upper'],-80,bad_worst,clean_worst,-100)+3
    print(dict(kind='complete_h6_a1_c6_no_full_reset',
               residual_count=len(old['residual_live_counts']),defective_count_cells=total,
               defective_good_route_upper=bad_worst,defective_worst_live=bad_witness,
               clean_good_route_upper=clean_worst,clean_worst_live=clean_witness,
               old_remaining_good_route_upper=old['closed_good_route_upper'],
               old_route_exception=-80,new_route_exceptions=(-100,-100),
               complete_branch_upper=upper,closed=upper<=-40,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


def audit_six():
    from diagnose_pair_unaligned_frontier import legacy_branch
    began=time.monotonic(); results=[]
    for initial in (0,1):
        for episodes in range(1,6):
            result=legacy_branch(False,initial,episodes)
            results.append(dict(branch=(6,initial,episodes),
                                residual_count=len(result['residual_live_counts']),
                                remaining_good_route_upper=result['closed_good_route_upper']))
            if result['residual_live_counts']:
                print(dict(unexpected_residual=result),flush=True)
                raise AssertionError('uncovered six-defect branch')
            if time.monotonic()-began>60:raise TimeoutError('remaining branch audit 60s limit')
    print(dict(kind='other_six_defect_no_full_reset_branches',branches=results,
               worst_good_route_upper=max(r['remaining_good_route_upper'] for r in results),
               common_old_route_exception=-80,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--test',action='store_true')
    p.add_argument('--live',type=int,default=1941)
    p.add_argument('--bound',action='store_true')
    p.add_argument('--integrate',action='store_true')
    p.add_argument('--audit-six',action='store_true')
    p.add_argument('--length-limit',type=int,default=17)
    args=p.parse_args()
    if args.test:test()
    elif args.audit_six:audit_six()
    elif args.integrate:integrate()
    elif args.bound:bound(args.live,args.length_limit)
    else:screen(args.live)
