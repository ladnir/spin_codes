"""Actual-regional rank bounds for full-batch intervals, fixed SL2 model.

Standard library, bounded serial work, no files. Integer/rational proofs;
not a whole-code certificate. Prior episode and defect helpers stay frozen.
"""
import argparse
import itertools
import math
import time
from fractions import Fraction as F
from functools import lru_cache

from diagnose_pair_defects import profile, volume_bound, ceildiv

M, B, N, K = 8192, 512, 1 << 21, 1 << 20
SCALE = 256


def upper_log_fraction(value):
    assert value > 0
    num, den = value.numerator ** SCALE, value.denominator ** SCALE
    return num.bit_length() - den.bit_length() + 1


@lru_cache(maxsize=32768)
def row_tail(mean, partial, threshold):
    """Worst joint-row positive-polynomial majorant, not independence of rows.

    Mean is the total number of included coordinates per RS row. There are
    <=partial nonconstant regions; every nonconstant occupancy is 1..15/16.
    Full regions shift the threshold. For fixed j,C, optimize the factorial
    moment order before taking the worst occupancy profile.
    """
    worst = F(0)
    case = None
    for j in range(partial + 1):
        for complete in range(max(0, math.ceil(mean - j)), min(256-j, math.floor(mean)) + 1):
            rem = mean - complete
            if j == 0:
                if rem != 0:
                    continue
                candidate = F(int(complete >= threshold))
                order = 0
            else:
                if not F(j, 16) <= rem <= F(15*j, 16):
                    continue
                need = threshold - complete
                if need <= 0:
                    candidate, order = F(1), 0
                elif need > j:
                    candidate, order = F(0), 0
                else:
                    p = rem / j
                    options = [(F(math.comb(j, k), math.comb(need, k)) * p ** k, k)
                               for k in range(1, need + 1)]
                    candidate, order = min([(F(1), 0)] + options)
            if candidate > worst:
                worst, case = candidate, (j, complete, order)
    return worst, case


def count_cap(beta, bits):
    """Except with probability 2^-bits, <=cap rows satisfy the event."""
    if beta == 0:
        return 0
    if beta >= 1:
        return M
    price = upper_log_fraction(beta)
    if price >= 0:
        return M

    def passes(s):
        return SCALE * (math.comb(M, s) - 1).bit_length() + s * price <= -SCALE * bits

    if not passes(M):
        return M
    lo, hi = 0, M
    # The binomial moment bound is log-concave in s and starts at 1. Its
    # negative sublevel is a suffix; the conservative ceil logs can shift
    # a crossing by <1 bit. Binary search need only return a verified cap.
    while hi - lo > 1:
        mid = (hi + lo) // 2
        if passes(mid):
            hi = mid
        else:
            lo = mid
    assert passes(hi)
    return hi - 1


def dimension_cap(mean, partial, bits):
    total = 0
    layers = []
    for threshold in range(129, 257):
        beta, case = row_tail(mean, partial, threshold)
        if beta == 0:
            break
        cap = count_cap(beta, bits)
        total += cap
        layers.append((threshold, cap, case))
    return total, layers


@lru_cache(maxsize=1024)
def rank_cap(full_batches, defects, bits=240):
    mean = F(full_batches, 16)
    partial = min(256, 2 * (defects + 1))
    dim, layers = dimension_cap(mean, partial, bits)
    codim, complement = dimension_cap(256 - mean, partial, bits)
    rank = max(0, full_batches * B - dim, K - codim)
    return dict(full_batches=full_batches, partial_regions=partial,
                dimension_cap=dim, complement_dimension_cap=codim,
                rank_lower=rank, dimension_layers=layers, complement_layers=complement)


@lru_cache(maxsize=1024)
def input_cap(batches, intervals, bits):
    return dimension_cap(F(batches, 16), min(256, 2*intervals), bits)[0]


@lru_cache(maxsize=128)
def early_joint(h, a, c, d=9, region=16):
    """(early length, necessary extra defects) -> upper early region cover.

    Each early interval has no full state and length <=8+9*internal_defects.
    Count the minimum defects its length requires; remaining defects can
    enlarge the proper tail of the final, terminal-live interval.
    """
    extra = h - c + a
    states = {(0, 0): 0}
    for i in range(c-1):
        nxt = {}
        for (length, cost), r in states.items():
            for ell in range(1, d + d*(extra-cost)):
                used = ceildiv(max(0, ell-(d-1)), d)
                regions = ceildiv(ell+1, region) if i == 0 and a else (ell+2*region-1)//region
                key = length+ell, cost+used
                nxt[key] = max(nxt.get(key, -1), r+regions)
        states = nxt
    return tuple((s, e, r) for (s, e), r in sorted(states.items()))


def joint_cell(h, a, c, live, early, used, early_regions, bits):
    extra = h-c+a
    if early >= live:
        return None
    partial = min(live, early + 8 + 9*(extra-used))
    if partial == live:
        # This could have no full state. The old count handles these small
        # lengths directly; callers must not silently discard this branch.
        return dict(unresolved='no guaranteed full state')
    regions = min(256, early_regions + ceildiv(live-early+1,16))
    if regions <= 128:
        return None
    slots = min(N, B*(live+c))
    random_min = B*(live-partial)+2*partial
    full = live-partial
    # Truncate full batches downward, pad allowed inputs upward. Both preserve
    # the number of intervals, so the good event for all interval sets applies.
    grid_full = full//8*8
    grid_input = min(4096, ceildiv(live+c,8)*8)
    rank = rank_cap(grid_full,h,bits)['rank_lower']
    dim = min(slots*(regions-128)//regions, input_cap(grid_input,c,bits))
    w = 460*N//1000
    inverse = 1+volume_bound(slots,w)-32*rank
    forward = 1+volume_bound(random_min,w)+32*(dim-random_min)
    return dict(h=h,a=a,c=c,live=live,early=early,used=used,partial=partial,
                regions=regions,rank=rank,input_dimension=dim,
                inverse=inverse,forward=forward,upper=min(inverse,forward))


def joint_cover(h=5, seconds=45):
    began=time.monotonic()
    # Covers both dimensions and complements on all <=h+1 interval sets,
    # every length already included. Three uses, <=128 tail levels each.
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    bits=placement+9+80
    worst=-(1<<60); witness=None; base_cells=refined=excluded=0
    residual=[]
    for a in (0,1):
        for c in range(1,h+a+1):
            for live in range(c,4096):
                base=profile(h,a,c,live,460,True)
                base_cells+=1
                if base is None:
                    excluded+=1;continue
                bound=base['log2_upper']
                if bound+placement>-80:
                    local=-(1<<60);local_witness=None
                    for early,used,r in early_joint(h,a,c):
                        item=joint_cell(h,a,c,live,early,used,r,bits)
                        refined+=1
                        if item is None:continue
                        if 'unresolved' in item:
                            value=bound
                        else:
                            value=min(bound,item['upper'])
                        if value>local:
                            local,local_witness=value,item
                        if time.monotonic()-began>seconds:
                            raise TimeoutError('joint-cover time cap; incomplete, no certificate')
                    bound=local
                    if bound+placement>-80:
                        residual.append(local_witness)
                if bound>worst:
                    worst=bound;witness=(a,c,live)
    result=dict(kind='no_full_reset_joint_cover',defects=h,base_cells=base_cells,
                excluded=excluded,refined_cells=refined,worst=witness,
                good_route_family_upper=worst+placement,route_failure_upper=-80,
                family_upper=max(worst+placement,-80)+1,
                residual_count=len(residual),residual_first=residual[:1],
                whole_code_certificate=False,seconds=time.monotonic()-began)
    print(result,flush=True)


def full_reset_cover(h=5, seconds=45):
    began=time.monotonic()
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    # I can include an isolated final batch, hence <=h+2 intervals.
    route_placements=26*(h+2)+(2*(h+2)-1).bit_length()
    bits=route_placements+9+80
    worst=-(1<<60); witness=None; cells=refined=excluded=0; residual=[]
    for a in (0,1):
        for c in range(1,h+a+1):
            for live in range(c,4096):
                base=profile(h,a,c,live,460)
                cells+=1
                if base is None:
                    excluded+=1;continue
                bound=base['log2_upper']
                if bound+placement>-80:
                    refined+=1
                    full=live-base['partial_cap']
                    rank=rank_cap(full//8*8,h,bits)['rank_lower']
                    ninput=min(4096,ceildiv(live+c+1,8)*8)
                    dim=min(base['input_dimension'],input_cap(ninput,c+1,bits))
                    inverse=1+volume_bound(base['slots'],460*N//1000)-32*rank
                    reset=1+32*dim+512*12-32*512
                    bound=min(bound,inverse,reset)
                    if bound+placement>-80:
                        residual.append((a,c,live,rank,dim,inverse,reset))
                if bound>worst:
                    worst=bound;witness=(a,c,live)
                if time.monotonic()-began>seconds:
                    raise TimeoutError('full-reset cover cap; incomplete, no certificate')
    print(dict(kind='full_reset_regional_cover',defects=h,cells=cells,
               excluded=excluded,refined_cells=refined,worst=witness,
               good_route_family_upper=worst+placement,route_failure_upper=-80,
               family_upper=max(worst+placement,-80)+1,
               residual_count=len(residual),residual_first=residual[:1],
               whole_code_certificate=False,seconds=time.monotonic()-began),flush=True)


def sparse_reset_profile():
    """First six-defect residual, using all five short resets jointly.

    Restricted no-full-reset profile: h<=6,a=0,c=6,T=1893. All six defects
    are births from zero, so every incoming state has only empty/full pairs.
    Each early interval has length<=8. RS support forces them all to length1.
    """
    lengths_checked=0
    for lengths in itertools.product(range(1,9),repeat=5):
        early=sum(lengths)
        regions=10+ceildiv(1893-early+1,16)
        if regions>=129:
            assert lengths==(1,1,1,1,1)
        lengths_checked+=1
    # With 129 regions and <=11 partial ones, the exact input mass forces
    # j=11 and every partial region to have occupancy 1/16.
    mean=F(1893+6,16)
    possible=[j for j in range(12) if F(129-j)+F(j,16)<=mean]
    assert possible==[11] and mean==118+F(11,16)
    q=1<<32
    # One row: <=q message choices, route containment 2^-44, and five
    # length-one resets of its s activation coordinates cost q^(-5s).
    # A single active row would activate a one-sided pair, impossible here.
    x=F(1,1<<44)*F(1,q**4)
    y=M*x
    assert y==F(1,1<<159)
    # Factor two for pair density; s!>=2 for s>=2 gives y^2/(1-y).
    per_placement=y*y/(1-y)
    assert per_placement < F(1,1<<317)
    print(dict(kind='six_defect_sparse_reset_profile',h=6,a=0,c=6,live=1893,
               early_length_vectors=lengths_checked,eligible_partial_regions=11,
               minimum_active_rows=2,closed_short_resets=5,
               per_placement_upper=-317,placement_union_bits=186,
               profile_upper=-131,whole_code_certificate=False),flush=True)


def selected(defects=5, live=2058, episodes=6):
    began = time.monotonic()
    base = profile(defects, 1, episodes, live, 460, True)
    assert base is not None
    full = live - base['partial_cap']
    rank = rank_cap(full, defects)
    good_bound = 1 + volume_bound(base['slots'], 460*N//1000) - 32 * rank['rank_lower']
    # All full-state unions of <=h+1 intervals; lengths already included.
    # Both rank sides and <=128 threshold levels use a further 8 bits.
    route_union = 26 * (defects + 1) + (defects + 1).bit_length() + 8
    print(dict(base=base, regional_rank=rank, good_route_class_upper=good_bound,
               route_failure_upper=-240+route_union,
               whole_code_certificate=False, seconds=time.monotonic()-began), flush=True)
    assert time.monotonic()-began < 30


def test():
    began = time.monotonic()
    checks = 0
    # Exact joint factorial domination on a small regional permutation.
    for sizes in itertools.product(range(5), repeat=3):
        masks = [[m for m in range(16) if m.bit_count() == size] for size in sizes]
        outcomes = list(itertools.product(*masks))
        p = tuple(F(size, 4) for size in sizes)
        for need in range(1, 4):
            probabilities = []
            for order in range(1, need + 1):
                e = sum((math.prod(p[j] for j in subset)
                         for subset in itertools.combinations(range(3), order)), F(0))
                probabilities.append(e / math.comb(need, order))
            beta = min([F(1)] + probabilities)
            for rows in range(1, 5):
                actual = F(sum(all(sum((m >> row) & 1 for m in regions) >= need
                                   for row in range(rows)) for regions in outcomes), len(outcomes))
                assert actual <= beta ** rows
                checks += 1
    caps=0
    for beta in (F(0),F(1),F(1,2),F(1,100),F(1,8192),F(1,10**6)):
        for bits in (40,240,300):
            cap=count_cap(beta,bits)
            if cap<M:
                s=cap+1
                assert math.comb(M,s)*beta**s <= F(1,1<<bits)
            caps+=1
    geometry=0
    for h in range(5):
        for a in (0,1):
            for c in range(1,min(h+a,3)+1):
                extra=h-c+a
                expected={}
                for lengths in itertools.product(range(1,2+2*extra),repeat=c-1):
                    cost=sum(ceildiv(max(0,ell-1),2) for ell in lengths)
                    if cost>extra:continue
                    regions=sum(ceildiv(ell+1,4) if i==0 and a else (ell+7)//4
                                for i,ell in enumerate(lengths))
                    key=sum(lengths),cost
                    expected[key]=max(expected.get(key,-1),regions)
                assert early_joint(h,a,c,2,4)==tuple((s,e,r) for (s,e),r in sorted(expected.items()))
                geometry+=1
    print(dict(exact_joint_factorial_checks=checks,exact_tail_cap_checks=caps,
               exhaustive_joint_geometry_tables=geometry,
               seconds=time.monotonic()-began), flush=True)
    assert time.monotonic()-began < 20


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--test', action='store_true')
    p.add_argument('--defects', type=int, default=5)
    p.add_argument('--live', type=int, default=2058)
    p.add_argument('--episodes', type=int, default=6)
    p.add_argument('--joint-cover',action='store_true')
    p.add_argument('--full-reset-cover',action='store_true')
    p.add_argument('--sparse-reset',action='store_true')
    args = p.parse_args()
    assert 0 <= args.defects <= 16
    if args.test:
        test()
    elif args.joint_cover:
        joint_cover(args.defects)
    elif args.full_reset_cover:
        full_reset_cover(args.defects)
    elif args.sparse_reset:
        sparse_reset_profile()
    else:
        selected(args.defects, args.live, args.episodes)


if __name__ == '__main__':
    main()
