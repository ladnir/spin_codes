"""Active-step Bellman bound retaining old-state and augmented-source ranks.

Supports one or two active coordinates at widths2/4. Independent equations
cost1/tau; dependent roots cost1/q. One common global zero-form space persists.
This nonlinear operator must precede fixed-gap linear envelopes, not be
averaged over random route positions. No calculations or files at import.
"""

from functools import lru_cache

from coupled_rank_transfer import rotate_state_mask


def active_step(values, width, input_mask, z, tau, coefficient_count=256,
                cycle=None, fresh=False):
    if width not in (2, 4) or input_mask <= 0 or input_mask >= 1 << width:
        raise ValueError("need width2/4 and a nonempty input mask")
    k = input_mask.bit_count()
    if k > 2 or not 0 < z <= 1 or tau < 1 or coefficient_count < 2:
        raise ValueError("need at most two input forms and valid markers")
    size = 1 << width
    if len(values) != width + 1 or any(len(row) != size for row in values):
        raise ValueError("terminal must be indexed by outgoing rank and support")
    if any(value < 0 for row in values for value in row):
        raise ValueError("terminal values must be nonnegative")
    zero, one = z * 0, z * 0 + 1
    price, p = one / tau, one / coefficient_count
    result = [[zero for _ in range(size)] for _ in range(width + 1)]
    result[0][0] = z**k * (values[k][input_mask] if fresh
                           else max(values[s][input_mask] for s in range(1, k + 1)))

    for state in range(1, size):
        rotated = rotate_state_mask(width, state, cycle)

        @lru_cache(maxsize=None)
        def branch(index, a, b, outside, following):
            # a=rank(old Q modulo W), b=rank(old Q and all X modulo W).
            if index == width:
                lower = max(b - k, b - a, 0)
                if not following:
                    return values[0][0] if lower == 0 else zero
                lower = max(1, lower)
                upper = min(width, b, following.bit_count())
                return max((values[s][following] for s in range(lower, upper + 1)), default=zero)
            bit = 1 << index
            active = bool(input_mask & bit)
            current, prior = bool(state & bit), bool(rotated & bit)

            def child(output_nonzero, next_nonzero, aa=a, bb=b, pp=outside):
                value = branch(index + 1, aa, bb, pp, following | (bit if next_nonzero else 0))
                return z * value if output_nonzero else value

            # An equation outside the old span consumes one input-image rank.
            # With <=2 inputs, their nonzero image set determines its update.
            outside_after = 0 if b - a == 1 else outside & ~bit

            if not active:
                if not current:
                    return child(False, prior)
                if not prior:
                    return child(True, True)
                nonzero, cancelled = child(True, True), child(True, False)
                best = max(nonzero, (one - p) * nonzero + p * cancelled)
                if a >= 2:
                    best = max(best, nonzero + price * child(True, False, a - 1, b - 1))
                return best

            is_outside = bool(outside & bit)
            if not current:
                if not prior:
                    return child(True, True)
                nonzero = child(True, True)
                if is_outside:
                    return nonzero + price * child(True, False, a, b - 1, outside_after)
                # x=RQ, unequal proportional forms, or independent old forms.
                best = max(nonzero, child(True, False))
                if a >= 2:
                    best = max(best, nonzero + price * child(True, False, a - 1, b - 1))
                return best

            if not prior:
                nonzero, cancelled = child(True, True), child(False, False)
                if is_outside:
                    return nonzero + price * child(False, False, a, b - 1, outside_after)
                best = max(nonzero, (one - p) * nonzero + p * cancelled)
                if a >= 2:
                    best = max(best, nonzero + price * child(False, False, a - 1, b - 1))
                return best

            # Both old state coordinates are nonzero, so y=Q'=0 is impossible.
            both = child(True, True)
            if is_outside:
                return (both + price * child(False, True, a, b - 1, outside_after)
                        + price * child(True, False, a, b - 1, outside_after))
            output_zero, next_zero = child(False, True), child(True, False)
            best = max(both, (one - p) * both + p * output_zero,
                       (one - p) * both + p * next_zero,
                       (one - 2 * p) * both + p * (output_zero + next_zero))
            if a >= 2:
                independent_y = price * child(False, True, a - 1, b - 1)
                independent_p = price * child(True, False, a - 1, b - 1)
                best = max(best,
                           both + independent_y,
                           both + independent_p,
                           p * output_zero + (one - p) * (both + independent_p),
                           p * next_zero + (one - p) * (both + independent_y),
                           both + independent_y + independent_p)
            return best

        for rank in range(1, state.bit_count() + 1):
            if fresh:
                patterns = ((k, input_mask),)
            else:
                patterns = [(0, 0)]
                patterns.extend((1, mask) for mask in range(1, size) if mask & ~input_mask == 0)
                if k == 2:
                    patterns.append((2, input_mask))
            result[rank][state] = max(branch(0, rank, rank + ell, outside, 0)
                                      for ell, outside in patterns)
    return result
