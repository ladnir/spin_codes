"""Small Bellman upper bound for a genuinely empty input interval.

Entries are indexed by EXACT residual state-span rank and current support.
The nonlinearity is a maximum over paired coefficient-dependence outcomes,
not a transition probability. Do not average this operator over random
route positions without a separate proof permitting that interchange.

Works with fractions for exact local tests or floats for diagnostic screens.
No calculations at import, file output, or width larger than eight.
"""

from functools import lru_cache

from coupled_rank_transfer import rotate_state_mask


def terminal_ones(width, one=1):
    size = 1 << width
    values = [[one * 0 for _ in range(size)] for _ in range(width + 1)]
    values[0][0] = one
    for rank in range(1, width + 1):
        for mask in range(1, size):
            if rank <= mask.bit_count():
                values[rank][mask] = one
    return values


def empty_step(values, width, z, tau, coefficient_count=256, cycle=None):
    """One empty-step continuation; values[r][H] are nonnegative bounds.

    A zero equation independent of the current coordinate costs one global
    rank and lowers local rank. A dependent cancellation costs at most
    1/coefficient_count probability without lowering rank. Nonzero and zero
    children must remain paired when maximizing over those cases.
    """
    # The rank price and root probability are distinct. Unlike the older
    # fixed-mask bound, this paired Bellman rule permits tau > q.
    if width not in (2, 4, 8) or not 0 < z <= 1 or tau < 1 or coefficient_count < 2:
        raise ValueError("invalid width or rank/output marker")
    size = 1 << width
    if len(values) != width + 1 or any(len(row) != size for row in values):
        raise ValueError("expected (width+1) by 2**width table")
    if any(value < 0 for row in values for value in row):
        raise ValueError("continuation values must be nonnegative")
    zero = z * 0
    one = zero + 1
    inverse_rank, root_probability = one / tau, one / coefficient_count
    result = [[zero for _ in range(size)] for _ in range(width + 1)]
    result[0][0] = values[0][0]
    for mask in range(1, size):
        rotated = rotate_state_mask(width, mask, cycle)
        mandatory = mask ^ rotated
        optional = tuple(j for j in range(width) if ((mask & rotated) >> j) & 1)

        @lru_cache(maxsize=None)
        def branch(index, rank, following):
            if index == len(optional):
                if not following or rank > following.bit_count():
                    return zero
                return values[rank][following]
            bit = 1 << optional[index]
            nonzero = branch(index + 1, rank, following | bit)
            dependent_zero = branch(index + 1, rank, following)
            # Outside the coefficient set there is no exceptional root.
            best = max(nonzero, (1 - root_probability) * nonzero + root_probability * dependent_zero)
            if rank >= 2:
                independent_zero = branch(index + 1, rank - 1, following)
                best = max(best, nonzero + inverse_rank * independent_zero)
            return best

        for rank in range(1, mask.bit_count() + 1):
            result[rank][mask] = z**mask.bit_count() * branch(0, rank, mandatory)
    return result


def diagnostic():
    """Small float-only comparison; not a route-averaged certificate."""
    import argparse
    import math

    from fractions import Fraction
    from coupled_rank_transfer import local_matrix, pair_antipodal_cycle

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", type=str, default="0.994")
    parser.add_argument("--steps", type=int, default=2048)
    args = parser.parse_args()
    if not 1 <= args.steps <= 8192:
        parser.error("steps must be between1 and8192")
    z = Fraction(args.z)
    if not 0 < z < 1:
        parser.error("need 0<z<1")
    width, tau = 4, 256
    cycle = pair_antipodal_cycle(width)
    matrix = [[float(v) for v in row] for row in local_matrix(width, 0, z, tau, cycle)]
    blind = [1.0] * (1 << width)
    ranked = terminal_ones(width, 1.0)
    print(f"FLOATING EMPTY-GAP DIAGNOSTIC, z={z}, tau={tau}; no route averaging or whole-code claim")
    for step in range(1, args.steps + 1):
        blind = [sum(weight * v for weight, v in zip(row, blind)) for row in matrix]
        ranked = empty_step(ranked, width, float(z), float(tau), cycle=cycle)
        if step in (1, 4, 16, 64, 256, args.steps):
            logs = [math.log2(max(ranked[r])) for r in range(1, width + 1)]
            print(f"L={step} blind_log2={math.log2(max(blind[1:])):.6f} rank1..4_log2="
                  + ",".join(f"{value:.6f}" for value in logs))


if __name__ == "__main__":
    diagnostic()
