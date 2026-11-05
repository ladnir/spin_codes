"""Joint activation containment, pair closure, route, and reset diagnostic.

Fixed construction. Restricted disjoint crossing episodes, no output claim.
Exact integer arithmetic; capped serial run, no numerical output files.
"""
import argparse
import itertools
import math
import time
from functools import lru_cache

from diagnose_pair_eligible_rows import eligible_majorant
from diagnose_pair_region_rank import count_cap

M, B, Q = 8192, 512, 1 << 32
NEG = -(1 << 50)


@lru_cache(maxsize=65536)
def choose(n,k):
    return math.comb(n,k) if 0 <= k <= n else 0


def logceil_ratio(n,d):
    assert n > 0 and d > 0
    return n.bit_length()-d.bit_length()+1


def episode_numerator(width, left, right, nleft, nright, minimum, length, q, power):
    """Weight summed over activation subsets U; denominator is external.

    There are width/2 matched pairs. The left region has left batches,
    beginning with activation; the right region has right batches.
    All later inputs must use coordinates in U. q^power is a common
    denominator for (length/q)^|U|.
    """
    total = 0
    for weight in range(minimum+(minimum&1), min(nleft,width)+1, 2):
        if weight > power:
            break
        total += (choose(width//2,weight//2)
                  * choose((left-1)*weight,nleft-weight)
                  * choose(right*weight,nright)
                  * length**weight * q**(power-weight))
    return total


@lru_cache(maxsize=256)
def episode_table(active, length):
    # Entry d is an upper logarithm for 2*active-d nonzero coordinates
    # across the episode's two regions. Maximize over splits and their
    # distribution; these are determined by a fixed word, not union choices.
    best = [NEG]*(2*active+1)
    minimum = 1 << length
    for left in range(1,length+1):
        right = length+1-left
        for nleft in range(minimum,active+1):
            for nright in range(active+1):
                numerator = episode_numerator(B,left,right,nleft,nright,
                                             minimum,length,Q,active)
                if not numerator:
                    continue
                denominator = choose(M,nleft)*choose(M,nright)*Q**active
                deficit = 2*active-nleft-nright
                best[deficit] = max(best[deficit],logceil_ratio(numerator,denominator))
    return tuple(best)


def joint_reset_route(active,lengths,dimension):
    budget = (dimension-1)*active
    dp = [0]+[NEG]*budget
    for length in lengths:
        table = episode_table(active,length)
        nxt = [NEG]*(budget+1)
        for used,previous in enumerate(dp):
            if previous == NEG:
                continue
            for extra,value in enumerate(table[:budget-used+1]):
                if value != NEG:
                    nxt[used+extra] = max(nxt[used+extra],previous+value)
        dp = nxt
    return max(dp)


def run(lengths=(3,3,3,3,3),dimension=3,report=True):
    began=time.monotonic()
    placements=26*(len(lengths)+1)
    cap=count_cap(eligible_majorant(lengths,dimension),placements+100)
    minimum=1<<max(lengths)
    terms=[]
    for active in range(minimum,cap+1):
        joint=joint_reset_route(active,lengths,dimension)
        if joint != NEG:
            terms.append(((choose(M,active)-1).bit_length()+32*dimension*active+joint,active))
        if time.monotonic()-began>40:
            raise TimeoutError('40-second activation-route cap; incomplete')
    if terms:
        worst=max(terms)
        upper=worst[0]+1+(len(terms)-1).bit_length()+placements
    else:
        worst=None
        upper=NEG
    result=dict(lengths=lengths,regions=128+dimension,eligible_cap=cap,
               source_counts=len(terms),worst=worst,good_route_upper=upper,
               route_exception_upper=-100,family_upper=max(upper,-100)+1,
               whole_code_certificate=False,seconds=time.monotonic()-began)
    if report:print(result,flush=True)
    return result


def cover():
    began=time.monotonic()
    worst=NEG
    witness=None
    cases=excluded=0
    for lengths in itertools.combinations_with_replacement(range(1,9),5):
        result=run(lengths,3,False)
        cases+=1
        excluded+=result['source_counts']==0
        if result['good_route_upper']>worst:
            worst,witness=result['good_route_upper'],result
        if time.monotonic()-began>40:
            raise TimeoutError('40-second family cap; incomplete')
    print(dict(regions=131,length_multisets=cases,excluded=excluded,worst=witness,
               good_route_upper=worst,route_exception_upper=-100,
               family_upper=max(worst,-100)+1,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


def test():
    began=time.monotonic()
    checks=0
    # Enumerate both regional images for a width-four toy, including
    # activation pair closure and containment of all later inputs in U.
    width=4
    for left,right in ((1,1),(1,2),(2,1),(2,2)):
        for nleft in range(1,5):
            for nright in range(5):
                counts={}
                for a in itertools.combinations(range(left*width),nleft):
                    initial={x for x in a if x<width}
                    if not initial or any((x^1) not in initial for x in initial):
                        continue
                    if any(x%width not in initial for x in a):
                        continue
                    for b in itertools.combinations(range(right*width),nright):
                        if all(x%width in initial for x in b):
                            w=len(initial)
                            counts[w]=counts.get(w,0)+1
                direct=sum(count*3**w*16**(4-w) for w,count in counts.items())
                formula=episode_numerator(width,left,right,nleft,nright,2,3,16,4)
                assert direct==formula
                checks+=1
    print(dict(exact_joint_route_checks=checks,seconds=time.monotonic()-began),flush=True)
    assert time.monotonic()-began<10


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--test',action='store_true')
    parser.add_argument('--cover',action='store_true')
    args=parser.parse_args()
    if args.test:test()
    elif args.cover:cover()
    else:run()
