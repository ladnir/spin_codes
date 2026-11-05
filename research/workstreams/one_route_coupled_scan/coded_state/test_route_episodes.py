"""Independent combinatorial checks for the routed episode certificate."""
from bisect import bisect_left, insort
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import comb
from route_episode_certificate import cover, feasible_count, N


def intervals_test(n):
    counts = Counter()
    def extend(end, volume, number):
        for start in range(end,n):
            for stop in range(start+1,n+1):
                terminal = stop==n
                length = stop-start
                if length<2 and not terminal:
                    continue
                key = terminal,volume+length,number+1
                counts[key] += 1
                if not terminal:
                    extend(stop,volume+length,number+1)
    extend(0,0,0)
    checked = 0
    for terminal in (False,True):
        delta = int(terminal)
        for c in range(1,n+1):
            for v in range(1,n+1):
                if not c+(c-delta)<=v<=n-(1-delta):
                    expected = 0
                else:
                    expected = comb(v-(c-delta)-1,c-1)*comb(n-v+c-1,c-delta)
                assert counts[terminal,v,c] == expected
                checked += 1
    return checked


def route_test():
    # Four regions and four labels. Exhaust every within-region selected subset.
    subsets = tuple(range(16))
    checks = 0
    for tau in (Fraction(3,2),Fraction(2)):
        for p in (1,2):
            numerators = defaultdict(Fraction)
            for masks in product(subsets,repeat=4):
                sizes = tuple(mask.bit_count() for mask in masks)
                value = Fraction(1)
                for label in range(p):
                    occupied = sum(bool(mask>>label & 1) for mask in masks)
                    value *= tau**(occupied-2) if occupied>=3 else 0
                numerators[sizes] += value
            for sizes,total in numerators.items():
                denominator = 1
                distribution = [Fraction(1)]
                for size in sizes:
                    denominator *= comb(4,size)
                    probability = Fraction(size,4)
                    result = [Fraction(0)]*(len(distribution)+1)
                    for t,value in enumerate(distribution):
                        result[t] += value*(1-probability)
                        result[t+1] += value*probability
                    distribution = result
                positive = sum(distribution[t]*tau**(t-2) for t in range(3,5))
                ordinary = sum(value*tau**max(0,t-2) for t,value in enumerate(distribution))
                factor = min(ordinary,(ordinary-1)*tau/(tau-1))
                assert total/denominator <= positive**p <= factor**p
                checks += 1
    return checks


def cover_test(boxes, expected):
    # Rectangles must be pairwise disjoint on each source-count slice.
    slices = 0
    for terminal in (False,True):
        for p in range(32,65):
            events = []
            for index,(kind,box,_,_) in enumerate(boxes):
                pl,ph,vl,vh,cl,ch = box
                if kind==terminal and pl<=p<=ph:
                    interval = cl,ch,index
                    events.extend(((vl,1,interval),(vh+1,-1,interval)))
            active = []
            for _,kind,interval in sorted(events):
                position = bisect_left(active,interval)
                if kind<0:
                    assert active[position] == interval
                    active.pop(position)
                else:
                    if position:
                        assert active[position-1][1] < interval[0]
                    if position<len(active):
                        assert interval[1] < active[position][0]
                    insort(active,interval)
            assert not active
            slices += 1
    # Independent triangular-domain count, rather than the box counting helper.
    n = N//128
    direct = 33*(sum(n-2*c for c in range(1,n//2+1))
                 +sum(n-2*c+2 for c in range(1,(n+1)//2+1)))
    actual = sum(feasible_count(box,kind) for kind,box,_,_ in boxes)
    assert actual == expected == direct == 4429185024
    return slices


if __name__=='__main__':
    print(f'Interval/composition checks: {intervals_test(8)} PASS',flush=True)
    print(f'Exact positive route-moment checks: {route_test()} PASS',flush=True)
    boxes,expected = cover(491*N//1000)
    print(f'Disjoint cover slices: {cover_test(boxes,expected)} PASS; complete count={expected}',flush=True)
