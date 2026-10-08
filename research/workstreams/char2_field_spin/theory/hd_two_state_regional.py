"""Exact small regional rank-envelope contraction; not a distance certificate.

The outer is coordinate-major, both pre-prefixes carry globally, and the
post-inner has two additive states with globally carried state. Regional
occupancies are retained separately. No global binomial normalization is used.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import comb, log2


def sharp_factor(tau, q):
    """Fresh affine equation, with its pivot form required to remain live."""
    tau = F(tau)
    if tau < 1 or q < 4:
        raise ValueError("Need tau>=1 and q>=4")
    return max(1 / tau, F(1, q - 1))


def zero_matrix():
    return [[F(0) for _ in range(4)] for _ in range(4)]


def identity():
    return [[F(i == j) for j in range(4)] for i in range(4)]


def matmul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(4)), F(0))
             for j in range(4)] for i in range(4)]


def vecmul(v, a):
    return [sum((v[i] * a[i][j] for i in range(4)), F(0)) for j in range(4)]


def pre_transfers(g):
    def options(old, active):
        if not active:
            return ((old, F(1)),)
        if not old:
            return ((1, F(1)),)
        return ((0, g), (1, F(1)))
    result = []
    for active in (0, 1):
        matrix = zero_matrix()
        for first, second in product((0, 1), repeat=2):
            for new_first, a in options(first, active):
                for new_second, b in options(second, new_first):
                    matrix[2 * first + second][2 * new_first + new_second] += a * b
        result.append(matrix)
    return tuple(result)


def post_transfers(y, g):
    y, g = F(y), F(g)
    result = []
    for active in (0, 1):
        matrix = zero_matrix()
        for old, new in product(range(4), repeat=2):
            if not active and new != old:
                continue
            if active and any(not (old & bit) and not (new & bit) for bit in (1, 2)):
                continue
            returns = (old & ~new).bit_count() if active else 0
            for out in (0, 1):
                if not out:
                    if (not active and old in (1, 2)) or (active and not old):
                        continue
                    if (active or old) and not new:
                        continue
                elif not active and not old:
                    continue
                pivots = returns + int(not out and bool(active or old))
                matrix[old][new] += g ** pivots * y ** out
        result.append(matrix)
    return tuple(result)


def polynomial_power(t0, t1, rows):
    if not 0 <= rows <= 32:
        raise ValueError("Exact reference is bounded to at most32 rows")
    coefficients = [identity()]
    for _ in range(rows):
        following = [zero_matrix() for _ in range(len(coefficients) + 1)]
        for h, matrix in enumerate(coefficients):
            for active, factor in enumerate((t0, t1)):
                term = matmul(matrix, factor)
                for i, j in product(range(4), repeat=2):
                    following[h + active][i][j] += term[i][j]
        coefficients = following
    return coefficients


def regional_post_kernels(rows, y, g):
    coefficients = polynomial_power(*post_transfers(y, g), rows)
    return [[[entry / comb(rows, h) for entry in row] for row in matrix]
            for h, matrix in enumerate(coefficients)]


def post_profile(profile, kernels):
    vector = [F(1), F(0), F(0), F(0)]
    for h in profile:
        vector = vecmul(vector, kernels[h])
    return sum(vector)


def source_rows(b, k):
    if not 1 <= k < b:
        raise ValueError("Invalid outer dimension")
    return [(mask, max(0, mask.bit_count() - b + k)) for mask in range(1 << b)
            if mask == 0 or mask.bit_count() > b - k]


def precoder_profiles(b, k, rows, tau, g):
    """Exact coefficients of the four-flag rank envelope, not exact ranks."""
    choices = source_rows(b, k)
    if b * rows > 16 or len(choices) ** rows > 100000:
        raise ValueError("Use only bounded small exact proof objects")
    transitions = pre_transfers(g)
    result = defaultdict(F)
    for sources in product(choices, repeat=rows):
        if not any(mask for mask, _ in sources):
            continue  # Exclude the all-zero message before the final1/tau.
        dimension = sum(d for _, d in sources)
        states = {(0, (0,) * b): F(tau) ** dimension}
        # This is the actual chronological order, including boundary carry.
        for region in range(b):
            for row in range(rows):
                active = bool(sources[row][0] & (1 << region))
                following = defaultdict(F)
                for (old, profile), mass in states.items():
                    for new, weight in enumerate(transitions[active][old]):
                        if weight:
                            updated = list(profile)
                            updated[region] += new & 1
                            following[new, tuple(updated)] += mass * weight
                states = following
        for (_, profile), mass in states.items():
            result[profile] += mass
    return dict(result)


def exact_bound(b=4, k=2, rows=2, q=4, tau=3, y=F(1, 2), cutoff=3):
    if not 0 < y <= 1 or not 0 <= cutoff < b * rows:
        raise ValueError("Invalid output mark or cutoff")
    g = sharp_factor(tau, q)
    profiles = precoder_profiles(b, k, rows, tau, g)
    kernels = regional_post_kernels(rows, y, g)
    moment = sum((mass * post_profile(profile, kernels)
                  for profile, mass in profiles.items()), F(0))
    z = F(q * q - 3 * q + 3, (q - 1) ** 2)
    upper = moment / F(tau) / F(y) ** cutoff / z ** (b * rows)
    return dict(outer=[b, k], rows=rows, q=q, tau=F(tau), y=F(y), cutoff=cutoff,
                g=g, profile_count=len(profiles), rank_moment=moment,
                probability_upper_unclipped=upper, probability_upper=min(F(1), upper),
                margin_bits=-log2(float(upper)), exact_arithmetic=True,
                nontrivial_upper_bound=upper < 1, whole_code_certified=False,
                route="independent uniform regional")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer", type=int, nargs=2, default=(4, 2))
    parser.add_argument("--rows", type=int, default=2)
    parser.add_argument("--q", type=int, default=4)
    parser.add_argument("--tau", type=F, default=F(3))
    parser.add_argument("--y", type=F, default=F(1, 2))
    parser.add_argument("--cutoff", type=int, default=3)
    args = parser.parse_args()
    result = exact_bound(*args.outer, args.rows, args.q, args.tau, args.y, args.cutoff)
    for key, value in result.items():
        print(key, str(value))


if __name__ == "__main__":
    main()
