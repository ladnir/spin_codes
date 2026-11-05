"""Expose a uniform pair route online, retaining the whole local rank table.

The suffix bound is uniform over remaining labeled forms and common quotients.
No lookahead equations are allowed. See ONLINE_ROUTE.md; this is not an
interchange of the older fixed-route nonlinear bound and route averaging.
No jobs, file output, or generated data at import.
"""

from check_two_groups import fresh_terms
from coupled_rank_transfer import local_matrix, pair_antipodal_cycle
from empty_span_rank import empty_step


VALID = [(0, 0)] + [(r, h) for h in range(1, 16) for r in range(1, h.bit_count() + 1)]


def zeros(z):
    return [[z * 0 for _ in range(16)] for _ in range(5)]


def active_terms(mask, z, tau, *, fresh=False, outer_randomized=False):
    """Positive support terms; retain rank gains only for actual fresh inputs."""
    from fractions import Fraction as F
    from decimal import Decimal
    zz, tt = F(z), F(tau)
    def convert(x):
        if isinstance(z, Decimal):
            return Decimal(x.numerator) / Decimal(x.denominator)
        return x if isinstance(z, F) else float(x)
    if fresh or outer_randomized:
        price = tt if fresh else min(tt, 256)
        return {h: [(j, convert(v), gain) for j, v, gain in fresh_terms(mask, h, zz, price)]
                for h in range(1, 16)}
    matrix = local_matrix(4, mask, zz, min(tt, 256), pair_antipodal_cycle(4))
    return {h: [(j, convert(v), 0) for j, v in enumerate(matrix[h]) if v] for h in range(1, 16)}


def active(values, mask, z, tau, *, fresh=False, outer_randomized=False, terms=None, rank_floor=False):
    if rank_floor and not (fresh or outer_randomized):
        raise ValueError("rank-floor branch needs unaggregated input charge counts")
    terms = terms if terms is not None else active_terms(mask, z, tau, fresh=fresh,
                                                        outer_randomized=outer_randomized)
    out = zeros(z)
    k = mask.bit_count()
    out[0][0] = z**k * (values[k][mask] if fresh else max(values[s][mask] for s in range(1, k + 1)))
    for r, h in VALID[1:]:
        for following, weight, gained in terms[h]:
            floor = 0
            if rank_floor:
                old_equations = (h & ~mask & ~following).bit_count()
                floor = max(k - gained, r - gained - old_equations) if fresh else r - k - gained - old_equations
                floor = max(0, floor)
            if not following:
                if not floor:
                    out[r][h] += weight * values[0][0]
                continue
            cap = min(r + k - (gained if fresh else 0), following.bit_count())
            floor = max(1, floor)
            if floor <= cap:
                out[r][h] += weight * max(values[s][following] for s in range(floor, cap + 1))
    return out


def combine(z, terms):
    """Positive sum of (weight, table) pairs."""
    out = zeros(z)
    for weight, table in terms:
        for r, h in VALID:
            out[r][h] += weight * table[r][h]
    return out


def counts(m, active_pairs):
    """Exact next-batch counts for a singleton in each of active_pairs pairs."""
    if m < 1 or active_pairs not in (0, 1, 2) or active_pairs > 2 * m:
        raise ValueError("invalid remaining route counts")
    from math import comb
    denominator = comb(2 * m, 2)
    return [(c, comb(active_pairs, c) * comb(2 * m - active_pairs, 2 - c), denominator)
            for c in range(3) if c <= active_pairs and 2 - c <= 2 * m - active_pairs]


def region(values, z, tau, steps=2048, *, fresh=False, outer_randomized=False, rank_floor=False):
    if not 1 <= steps <= 2048 or (fresh and outer_randomized):
        raise ValueError("unsupported route or simultaneous input assumptions")
    masks = ((3, 12), (1, 2, 4, 8), (5, 6, 9, 10))
    compiled = {mask: active_terms(mask, z, tau, fresh=fresh, outer_randomized=outer_randomized)
                for group in masks for mask in group}
    one = z * 0 + 1

    def mean_active(table, group):
        return combine(z, [(one / len(group), active(table, mask, z, tau, fresh=fresh,
                            outer_randomized=outer_randomized, terms=compiled[mask], rank_floor=rank_floor)) for mask in group])

    empty = values
    pair, single, split = zeros(z), zeros(z), zeros(z)
    for m in range(1, steps + 1):
        new_empty = empty_step(empty, 4, z, tau, cycle=pair_antipodal_cycle(4))
        new_pair = combine(z, [(one / m, mean_active(empty, masks[0]))])
        new_single = combine(z, [(one / m, mean_active(empty, masks[1]))])
        if m > 1:
            new_pair = combine(z, [(one, new_pair), ((one * (m - 1)) / m,
                                empty_step(pair, 4, z, tau, cycle=pair_antipodal_cycle(4)))])
            new_single = combine(z, [(one, new_single), ((one * (m - 1)) / m,
                                  empty_step(single, 4, z, tau, cycle=pair_antipodal_cycle(4)))])
        split_terms = []
        for consumed, numerator, denominator in counts(m, 2):
            if not numerator:
                continue
            if consumed == 0:
                table = empty_step(split, 4, z, tau, cycle=pair_antipodal_cycle(4))
            elif consumed == 1:
                table = mean_active(single, masks[1])
            else:
                table = mean_active(empty, masks[2])
            split_terms.append((one * numerator / denominator, table))
        empty, pair, single, split = new_empty, new_pair, new_single, combine(z, split_terms)
    return pair, split


def envelope(vector, z, tau, steps=2048, *, fresh=False, outer_randomized=False):
    """Collapse only after an entire region, separately for its two terminals."""
    reset = zeros(z)
    reset[0][0] = z * 0 + 1
    live = [[value for value in row] for row in vector]
    live[0][0] = z * 0
    r_pair, r_split = region(reset, z, tau, steps, fresh=fresh, outer_randomized=outer_randomized)
    l_pair, l_split = region(live, z, tau, steps, fresh=fresh, outer_randomized=outer_randomized)
    def collapse(reset_table, live_table):
        return [[reset_table[0][0], live_table[0][0]],
                [max(reset_table[r][h] / vector[r][h] for r, h in VALID[1:]),
                 max(live_table[r][h] / vector[r][h] for r, h in VALID[1:])]]
    return collapse(r_pair, l_pair), collapse(r_split, l_split)


def main():
    import argparse
    import math
    from decimal import Decimal as D, ROUND_CEILING, localcontext
    from fractions import Fraction as F
    from check_empty_span_potential import check
    from check_sparse_full_gate import multiply, upward
    from check_two_groups import evaluate_fresh_columns
    from lazy_outer_step import pivot_envelope

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, default=0, choices=(0, 80, 120))
    args = parser.parse_args()
    v, z, tau, lam = check(F(499, 500), F(993, 1000), tau=F(1024))
    with localcontext() as ctx:
        ctx.prec = args.precision or 80
        ctx.rounding = ROUND_CEILING
        convert = upward if args.precision else float
        vector = [[convert(x) for x in row] for row in v]
        zz, tt = convert(z), convert(tau)
        late = envelope(vector, zz, tt, outer_randomized=True)
        print("online late regional envelopes", late, flush=True)
        fresh = envelope(vector, zz, tt, fresh=True)
        print("online fresh regional envelopes", fresh, flush=True)
        # Float screens are diagnostic only; converting them to Decimal does
        # not make them rigorous. Directed runs construct the tables in Decimal.
        matrix = lambda a: [[x if isinstance(x, D) else D(str(x)) for x in row] for row in a]
        p = [[upward(x) for x in row] for row in pivot_envelope(v, tau)]
        bounds = evaluate_fresh_columns(*(matrix(a) for a in late), *(matrix(a) for a in fresh),
                                       min(v[r][h] for r, h in VALID[1:]), z, tau, lam,
                                       late_prefix=multiply(p, p))
        factor = upward(F(256, 255)**64)
        bounds = {h: x * factor for h, x in bounds.items()}
        print("UPWARD" if args.precision else "FLOAT DIAGNOSTIC", "precision", args.precision,
              "log2 union", math.log2(sum(bounds.values(), D(0))))
        print("selected log2 bounds", {h: math.log2(bounds[h]) for h in (33, 48, 64)})
        print("Partial two-group class only; no whole-code certificate.")


if __name__ == "__main__":
    main()
