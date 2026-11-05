"""Exact-integer selected-class gate using a uniform matching-graph envelope.

This is not full coverage of state histories. No files or raw masks are saved.
"""
from fractions import Fraction
from math import comb, isqrt
import heapq
import random
import time


def graph_envelope(M, regions):
    choose_bits = []
    coefficient = 1
    for s in range(M//2+1):
        if s:
            coefficient = coefficient*(M-s+1)//s
        choose_bits.append(coefficient.bit_length())

    def ceil_sqrt_ratio(num,den):
        root = isqrt(num//den)
        return root+(den*root*root < num)

    bounds = []
    for s in range(M+1):
        small = min(s, M-s)
        n = regions*(small//2)
        pnum = small-1 if small % 2 == 0 else small
        if not n:
            upper = 0
        else:
            mean_up = (n*pnum+M-2)//(M-1)
            # ln(2)<7/10; union subsets at fixed cardinality, <=8192 sizes,
            # and <=256 suffixes. Compare two proved binomial-mgf tails.
            Lbits = choose_bits[small]+77
            hoeffding = ceil_sqrt_ratio(7*n*Lbits,20)
            bernstein = (ceil_sqrt_ratio(7*n*pnum*Lbits,5*(M-1))
                         +(7*Lbits+14)//15)
            slack = min(hoeffding,bernstein)
            upper = min(n, mean_up+slack)
        if s > M//2:
            upper += regions*(s-M//2)
        bounds.append(upper)
    return bounds


def free_dimension(M, regions, k, D, envelope):
    assert D % 2 == 0
    return min(D//2, max(0, max(min(D, D//2+envelope[s], regions*s)-k*s
                               for s in range(M+1))))


def dimension_sweep(M, regions, k, limit, envelope, seconds=15):
    """Exact maximum of the three affine pieces, without storing all masks."""
    started = time.monotonic()
    events = []
    for s, E in enumerate(envelope):
        assert 2*E <= regions*s
        for start,end,slope,intercept in (
                (0,E,2,-k*s),
                (E,regions*s-E,1,E-k*s),
                (regions*s-E,limit,0,(regions-k)*s)):
            if start <= min(end,limit):
                events.append((start,min(end,limit),slope,intercept))
    events.sort()
    heaps = [[],[],[]]
    cursor = 0
    for m in range(limit+1):
        if m % 4096 == 0 and time.monotonic()-started > seconds:
            raise RuntimeError('Exact sweep time cap; no partial coverage claim')
        while cursor < len(events) and events[cursor][0] <= m:
            start,end,slope,intercept = events[cursor]
            heapq.heappush(heaps[slope],(-intercept,end))
            cursor += 1
        values = [0]
        for slope,heap in enumerate(heaps):
            while heap and heap[0][1] < m:
                heapq.heappop(heap)
            if heap:
                values.append(slope*m-heap[0][0])
        yield min(m,max(values))


def suffix_distance_cap(regions):
    M, B, bits, batch = 8192, 256, 32, 512
    N, u = M*B, M*regions
    P = (u-batch)//2
    envelope = graph_envelope(M,regions)
    target = -80
    lo = max(0,(target-P-1)//bits-max(0,u-N//2)+P+1)
    assert lo <= P
    gap,m = min_dimension_gap(regions,B//2,lo,P,envelope)
    cap = P+gap+(target-u-P-1)//bits
    e = P-m
    dim = free_dimension(M,regions,B//2,2*m,envelope)
    forward = P+bits*(max(0,u-N//2)-e)+1
    assert forward > target and m-dim == gap
    assert u+P+bits*(cap-P-m+dim)+1 <= target
    return dict(kind='analytic_complete_dense_suffix_cover', regions=regions,
                exception_counts=P+1, weight_threshold=cap,
                distance_percent=100*cap/N, conditional_log2_per_class=target,
                route_failure_below_power_of_two=-56, bottleneck=(e,m,dim,forward),
                whole_code_certificate=False)


def min_dimension_gap(regions,k,lo,hi,envelope):
    # m-dim(m)=max(0,min_s[max(-m,-E_s,m-R*s)+k*s]).
    # For fixed s, the bracket is minimized on [E_s,R*s-E_s].
    best = None
    for s,E in enumerate(envelope):
        m = max(lo,min(hi,E))
        value = max(0,max(-m,-E,m-regions*s)+k*s)
        if best is None or value < best[0]:
            best = value,m
    return best


def selected_gate(regions, exceptions):
    M, B, b, bits = 8192, 256, 512, 32
    N = M*B
    W = 45*N//100
    u = regions*M
    live = u-b
    pairs = live//2
    D = live-2*exceptions
    assert 0 <= exceptions <= pairs
    envelope = graph_envelope(M, regions)
    dim = free_dimension(M, regions, B//2, D, envelope)
    # One prescribed activation batch; all subsequent incoming coordinates
    # are nonzero. At most 2^u output supports and 2^pairs branch masks.
    inverse = u+pairs+bits*(W-D-exceptions+dim)+1
    forward = pairs+bits*(max(0, u-N//2)-exceptions)+1
    return dict(regions=regions, exceptions=exceptions, free_coordinates=D,
                uniform_shortened_dimension=dim,
                conditional_inverse_log2_strict_upper=inverse,
                conditional_forward_log2_strict_upper=forward,
                route_failure_below_power_of_two=-56,
                whole_code_certificate=False)


def tests():
    factorial_cases = 0
    for M in range(2, 15, 2):
        for s in range(M+1):
            n = s//2
            p = Fraction((s-1 if s % 2 == 0 else s), M-1)
            actual = Fraction(1)
            for j in range(n+1):
                if j:
                    i = j-1
                    actual *= Fraction((s-2*i)*(s-2*i-1), 2*j*(M-2*i-1))
                assert actual <= comb(n,j)*p**j
                factorial_cases += 1
    rng = random.Random(51)
    graph_cases = 0
    for _ in range(12):
        M, regions, k = 6, 3, 2
        edges = []
        for region in range(regions):
            order = list(range(M))
            rng.shuffle(order)
            edges += list(zip(order[::2],order[1::2]))
        envelope = [0]*(M+1)
        for mask in range(1 << M):
            internal = sum(bool(mask & (1 << a)) and bool(mask & (1 << b)) for a,b in edges)
            size = mask.bit_count()
            envelope[size] = max(envelope[size], internal)
        caps = [free_dimension(M,regions,k,2*m,envelope) for m in range(len(edges)+1)]
        assert caps == list(dimension_sweep(M,regions,k,len(edges),envelope))
        for lo in range(len(edges)+1):
            gap,m = min_dimension_gap(regions,k,lo,len(edges),envelope)
            assert gap == min(j-caps[j] for j in range(lo,len(edges)+1))
            assert gap == m-caps[m]
        for mask in range(1 << len(edges)):
            degree = [0]*M
            for j,(a,b) in enumerate(edges):
                if mask & (1 << j):
                    degree[a] += 1
                    degree[b] += 1
            dim = sum(max(0,d-k) for d in degree)
            assert dim <= caps[mask.bit_count()]
            graph_cases += 1
    return dict(factorial_moment_checks=factorial_cases, selected_edge_graph_checks=graph_cases)


def main():
    start = time.monotonic()
    print(tests())
    N = 2**21
    for R,e in ((192,3*N//20),(160,N//10),(224,N//5),(256,N//4)):
        print(selected_gate(R,e))
    for R in (160,192,224,256):
        print(suffix_distance_cap(R),flush=True)
    worst = None
    covered = 0
    for R in range(129,257):
        result = suffix_distance_cap(R)
        covered += result['exception_counts']
        if worst is None or result['weight_threshold'] < worst['weight_threshold']:
            worst = result
        if time.monotonic()-start > 50:
            raise RuntimeError('Full-suffix cap; do not claim full coverage')
    assert covered < 2**27
    print(dict(all_region_aligned_suffixes=128,exception_count_classes=covered,
               worst=worst,conditional_union_log2_upper=-53,
               combined_failure_below_power_of_two=-52))
    holes = 6144
    threshold = worst['weight_threshold']-2*holes
    assert threshold == 842604 and 100*threshold > 40*N
    print(dict(nonfull_pair_budget=holes,weight_threshold=threshold,
               distance_percent=100*threshold/N,
               arbitrary_first_batch_within_first_region=True,
               combined_failure_below_power_of_two=-52,
               whole_code_certificate=False))
    assert time.monotonic()-start < 65


if __name__ == '__main__':
    main()
