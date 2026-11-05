"""Closed-episode support/zero identities for the fixed SL2 model.

Combinatorial exhaustive checks, not a whole-code certificate. No output files.
"""
import argparse
import itertools
import math
import time
from functools import lru_cache

from diagnose_pair_activation_route import M, B, Q, NEG, choose, logceil_ratio, episode_table
from diagnose_pair_eligible_rows import eligible_majorant
from diagnose_pair_region_rank import count_cap


def partner(mask,phase,width):
    return sum(1 << (j ^ (1 << phase)) for j in range(width) if mask>>j&1)


def statistics(history,width):
    assert history[0] and history[-1]==0
    phases=(width-1).bit_length()
    initial=history[0]
    empty_births=singleton=births=deaths=forced=0
    envelope=initial
    input_union=initial
    for t,(current,nxt) in enumerate(zip(history,history[1:])):
        partners=partner(current,t%phases,width)
        live=current|partners
        ghost=partners&~current
        born=nxt&~current
        empty=born&~live
        singleton+=ghost.bit_count()
        empty_births+=empty.bit_count()
        births+=born.bit_count()
        deaths+=(current&~nxt).bit_count()
        forced+=(live&~nxt).bit_count()
        envelope|=empty|ghost
        input_union|=live|nxt
    charge=initial.bit_count()+empty_births+singleton
    return dict(charge=charge,forced=forced,initial=initial.bit_count(),
                births=births,deaths=deaths,extra=empty_births+singleton,
                envelope=envelope,input_union=input_union)


def initial_masks(width,weight,extra):
    # Number of supports having at most extra singleton pairs in one matching.
    pairs=width//2
    return sum(math.comb(pairs,singles)*2**singles
               * math.comb(pairs-singles,(weight-singles)//2)
               for singles in range(min(extra,weight,pairs)+1)
               if (weight-singles)%2==0
               and 0 <= (weight-singles)//2 <= pairs-singles)


def history_cap(width,length,weight,extra):
    # B births, each identified by coordinate+time; each initial/born
    # active component dies at one of <=length times. No free amplitude
    # variables are counted in this forward support-history bound.
    return length**weight * sum(math.comb(width*length,b)*length**b
                               for b in range(min(extra,width*length)+1))


def test():
    began=time.monotonic()
    identities=caps=masks=0
    for width,length in ((2,5),(4,4),(8,2)):
        full=(1<<width)-1
        counts={}
        for initial in range(1,full+1):
            for middle in itertools.product(range(1,full+1),repeat=length-1):
                history=(initial,)+middle+(0,)
                item=statistics(history,width)
                assert item['forced']==item['charge']
                assert item['deaths']==item['initial']+item['births']
                assert item['births']<=item['extra']
                assert item['input_union']&~item['envelope']==0
                assert item['envelope'].bit_count()<=item['charge']
                if item['extra']<=2:
                    key=(initial,item['extra'])
                    counts[key]=counts.get(key,0)+1
                identities+=1
        for (initial,extra),count in counts.items():
            assert count<=history_cap(width,length,initial.bit_count(),extra)
            caps+=1
        for weight in range(width+1):
            for extra in range(width//2+1):
                exact=sum(mask.bit_count()==weight and
                          (partner(mask,0,width)&~mask).bit_count()<=extra
                          for mask in range(full+1))
                assert exact==initial_masks(width,weight,extra)
                masks+=1
        if time.monotonic()-began>20:
            raise TimeoutError('20-second support-identity test cap')
    print(dict(support_history_identities=identities,history_count_bounds=caps,
               initial_mask_enumerators=masks,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


def costs():
    # Upper bit costs for one and two extra charges before the route count.
    # For fixed activation U, the clean factor is (length/q)^weight.
    for length in (4,8,17):
        for extra in (1,2):
            multiplicity=sum(math.comb(512*length,b)*length**b
                             for b in range(extra+1))
            loss=(multiplicity-1).bit_length()
            print(dict(length=length,extra=extra,history_bits=loss,
                       net_probability_bits=32*extra-loss,
                       excludes_route_and_activation_choices=True),flush=True)


@lru_cache(maxsize=128)
def defective_table(active,length,extra):
    assert 1<=extra<=2 and 1<=length<=8+9*extra
    best=[NEG]*(2*active+1)
    for left in range(1,length+1):
        right=length+1-left
        for nleft in range(1,active+1):
            for nright in range(active+1):
                numerator=0
                for weight in range(1,min(nleft,B)+1):
                    size=min(B,weight+extra)
                    numerator += (initial_masks(B,weight,extra)
                                  * choose((left-1)*size,nleft-weight)
                                  * choose(right*size,nright)
                                  * history_cap(B,length,weight,extra)
                                  * Q**(active-weight))
                if numerator:
                    denominator=choose(M,nleft)*choose(M,nright)*Q**(active+extra)
                    deficit=2*active-nleft-nright
                    best[deficit]=max(best[deficit],logceil_ratio(numerator,denominator))
    return tuple(best)


def profile(extra=1):
    assert 0<=extra<=2
    began=time.monotonic()
    lengths=(3,3,3,3,3)
    dimension=3
    placements=156
    cap=count_cap(eligible_majorant(lengths,dimension),placements+100)
    terms=[]
    for active in range(8,cap+1):
        budget=(dimension-1)*active
        dp=[0]+[NEG]*budget
        for index,length in enumerate(lengths):
            table=(defective_table(active,length,extra) if index==0 and extra
                   else episode_table(active,length))
            nxt=[NEG]*(budget+1)
            for used,previous in enumerate(dp):
                if previous==NEG:continue
                for added,value in enumerate(table[:budget-used+1]):
                    if value!=NEG:nxt[used+added]=max(nxt[used+added],previous+value)
            dp=nxt
        envelope=tuple(itertools.accumulate(dp,max))
        zeros=[]
        for count in range(budget+1):
            if envelope[count]!=NEG:
                zeros.append(((choose(131*active,count)-1).bit_length()
                              +32*(dimension*active-count)+envelope[count],count))
        if zeros:
            best=max(zeros)
            terms.append(((choose(M,active)-1).bit_length()+best[0]
                          +(len(zeros)-1).bit_length(),active,best[1]))
        if time.monotonic()-began>35:
            raise TimeoutError('35-second defective profile cap; incomplete')
    worst=max(terms)
    # Five possible exceptional episodes. The actual interval endpoints
    # determine order but not which interval carries the extra charge.
    upper=worst[0]+1+(len(terms)-1).bit_length()+placements+(3 if extra else 0)
    print(dict(lengths=lengths,regions=131,defective_episodes=int(bool(extra)),extra_charge=extra,
               eligible_cap=cap,source_counts=len(terms),worst=worst,
               good_route_upper=upper,route_exception_upper=-100,
               family_upper=max(upper,-100)+1,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--test',action='store_true')
    p.add_argument('--profile',action='store_true')
    p.add_argument('--extra',type=int,default=1)
    args=p.parse_args()
    if args.test:test()
    elif args.profile:profile(args.extra)
    else:costs()
