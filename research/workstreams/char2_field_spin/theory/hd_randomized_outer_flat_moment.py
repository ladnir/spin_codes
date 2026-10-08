"""Positive zero-mask bound for a uniform nonzero incoming state.

The DP retains state constraints from occupied tiles. Its coefficients are
positive; removing the zero state never subtracts an upper bound from another.
No whole-code distance certificate or native implementation is provided.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, log2


def state_rank_join(previous, local, state=4, *, distinct_odd=True):
    """Safe span lower bound for constraints contributed by distinct tiles."""
    if distinct_odd and previous == local == 1 and state >= 2:
        return 2
    if distinct_odd and previous == local == 3 and state == 4:
        return 4
    return max(previous, local)


def profile_rank_lower(inputs, zeros, *, width=8, state=4, distinct_odd=True):
    """Counts-only lower bound on the extra rank supplied by the common state."""
    if len(inputs) != len(zeros) or any(not 0 <= x <= width for x in inputs+zeros):
        raise ValueError('Invalid tile count profile')
    rank = untouched = gap = 0
    for a, v in zip(inputs, zeros):
        rank = state_rank_join(rank, min(state, max(v-a, 0)), state,
                               distinct_odd=distinct_odd)
        untouched += v if a == 0 else 0
        gap += width-max(a, v)
    length, active = width*len(inputs), sum(inputs)
    global_rank = max(0, min(state, length-active)-gap)
    return max(rank, min(state, untouched), global_rank)


@lru_cache(None)
def local_terms(mark, field=65536, width=8, state=4):
    """Exact positive local coefficients, scaled by field**a."""
    mark = F(mark)
    if not 0 < mark < 1 or not 1 <= state <= width:
        raise ValueError('Need 0<z<1 and a valid state size')
    grouped = {}
    for a in range(width+1):
        for v in range(width+1):
            local_rank = min(state, max(v-a, 0))
            if local_rank == state:
                continue  # No nonzero state survives this tile's constraints.
            key = a, local_rank, min(state, v if a == 0 else 0), min(state, width-max(a, v))
            coefficient = (comb(width, a)*comb(width, v)*mark**(width-v)
                           *(1-mark)**v*field**max(a-v, 0))
            grouped[key] = grouped.get(key, F(0))+coefficient
    return tuple((*key, coefficient) for key, coefficient in sorted(grouped.items()))


def flat_moments(mark, field=65536, *, width=8, tiles=8, state=4,
                 directed=False, distinct_odd=True, global_rank=True,
                 return_stats=False):
    """Bounds for every j; Arb enclosures or FLOAT base-two logarithms.

    The actual target uses width8, tiles8, state4, with beta64's disjoint
    input/output cosets. The width2/state1 toy is independently tested.
    """
    if directed:
        from flint import arb
        zero, one = arb(0), arb(1)
        convert = lambda value: arb(value.numerator)/value.denominator
        algebraic_mark = isinstance(mark, arb)
    else:
        zero, one = 0.0, 1.0
        convert = float
        algebraic_mark = False
    if not algebraic_mark:
        mark = F(mark)
    length = width*tiles
    if mark == 1:
        result = (one if directed else 0.0,)*(length+1)
        return (result, dict(peak_states=1)) if return_stats else result
    if algebraic_mark:
        if not 0 < mark < 1:
            raise ValueError('The mark enclosure must lie strictly between zero and one')
        local = tuple((a, min(state, max(v-a, 0)), min(state, v if a == 0 else 0),
                       min(state, width-max(a, v)),
                       comb(width, a)*comb(width, v)*mark**(width-v)
                       *(1-mark)**v*field**max(a-v, 0))
                      for a in range(width+1) for v in range(width+1)
                      if max(v-a, 0) < state)
    else:
        local = tuple((a, r, b, g, convert(coefficient))
                      for a, r, b, g, coefficient in local_terms(mark, field, width, state))
    # r is ONLY the span bound from whole tile constraint spaces. Do not fold
    # b into r: the 3+3 rule requires two genuine three-dimensional tile spaces.
    dp = {(0, 0, 0, 0):one}
    peak = 1
    for unused in range(tiles):
        updated = {}
        for (j, rank, untouched, gap), value in dp.items():
            for a, local_rank, local_untouched, local_gap, coefficient in local:
                successor = state_rank_join(rank, local_rank, state,
                                            distinct_odd=distinct_odd)
                new_untouched = min(state, untouched+local_untouched)
                if successor == state or new_untouched == state:
                    continue
                key = j+a, successor, new_untouched, min(state, gap+local_gap)
                updated[key] = updated.get(key, zero)+value*coefficient
        dp = updated
        peak = max(peak, len(dp))
    weighted = [zero for unused in range(length+1)]
    denominator_state = field**state-1
    for (j, rank, untouched, gap), value in dp.items():
        global_lower = (max(0, min(state, length-j)-gap) if global_rank else 0)
        lower = max(rank, untouched, global_lower)
        if lower < state:
            numerator = field**(state-lower)-1
            factor = (arb(numerator)/denominator_state if directed
                      else numerator/denominator_state)
            weighted[j] += value*factor
    if directed:
        result = tuple(value/(field**j*comb(length, j)) for j, value in enumerate(weighted))
    else:
        result = tuple(log2(value)-j*log2(field)-log2(comb(length, j))
                       if value else float('-inf') for j, value in enumerate(weighted))
    return (result, dict(peak_states=peak)) if return_stats else result


def compare_selected(mark=F(1, 2), positions=(0, 1, 8, 16, 24, 32, 40, 48, 56, 60, 64)):
    """Bounded local diagnostic only; the comparison has no distance claim."""
    from hd_randomized_outer_state_screen import NEG, total, mds_moment
    # Reconstruct the old V_j bound directly; the cone matrix combines it with
    # another column term, so reading that matrix would not isolate V_j.
    log_z = log2(mark)
    from hd_randomized_outer_state_profile import punctured_state_moment
    from hd_randomized_outer_state_screen import add
    dp = {(0, 0):0.0}
    for unused in range(8):
        following = {}
        for (j, u), value in dp.items():
            for a in range(9):
                key = j+a, u+(a != 0)
                term = value+log2(comb(8, a))+mds_moment(a, log_z)
                following[key] = add(following.get(key, NEG), term)
        dp = following
    previous = tuple(total(value+punctured_state_moment(64-8*u, log_z)
                           for (jj, u), value in dp.items() if jj == j)-log2(comb(64, j))
                     for j in range(65))
    new, stats = flat_moments(mark, return_stats=True)
    return dict(mark=str(mark), directed=False, peak_states=stats['peak_states'],
                positions=[dict(j=j, old_log2=previous[j], new_log2=new[j],
                                gain_bits=previous[j]-new[j]) for j in positions])


if __name__ == '__main__':
    import json
    for mark in (F(1, 2), F(1, 16), F(1, 256)):
        print(json.dumps(compare_selected(mark), indent=2), flush=True)
