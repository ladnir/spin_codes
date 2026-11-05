"""Exact bounded quotient checks for ACTIVE_RANK_FLOOR.md; no distance claim.

Reuses the independent GF16/nibble-form Gaussian oracle and its 12 fixed
form fixtures. No files, random search, or benchmark; run serially with -B.
"""

from collections import Counter
from itertools import product

from test_paired_rank_step import (COEFFICIENTS, fixtures, image_data, reduced,
                                   require, rotated, scale, span, support)


def main():
    coverage = Counter()
    for fixture in fixtures():
        name, variables, state, inputs, history, initial_zeros, fresh = fixture
        width = len(state)
        cycle = (0, 1) if width == 2 else (0, 2, 1, 3)
        initial = span(initial_zeros, variables)
        state, inputs, history = (reduced(words, initial, variables)
                                  for words in (state, inputs, history))
        h, mask = support(state), support(inputs)
        r, k = len(span(state, variables)), support(inputs).bit_count()
        old_required = tuple(x for x in state + inputs if x) + history
        prior = rotated(state, cycle)
        rotated_h = support(prior)
        for coefficients in product(COEFFICIENTS, repeat=width):
            output = tuple(x ^ scale(a, q, variables)
                           for x, a, q in zip(inputs, coefficients, state))
            following = tuple(p ^ y for p, y in zip(prior, output))
            for j, hp in product(range(1 << width), repeat=2):
                if j & ~support(output) or hp & ~support(following):
                    continue
                zero_forms = tuple(x for index, x in enumerate(output) if not j & (1 << index))
                zero_forms += tuple(x for index, x in enumerate(following) if not hp & (1 << index))
                required = old_required
                required += tuple(x for index, x in enumerate(output) if j & (1 << index))
                required += tuple(x for index, x in enumerate(following) if hp & (1 << index))
                basis = span(zero_forms, variables)
                if any(not x for x in reduced(required, basis, variables)):
                    continue
                gain = len(basis)
                s = len(span(reduced(following, basis, variables), variables))
                c = (mask & ~(j & hp)).bit_count()
                d = (h & ~mask & ~hp).bit_count()
                require(gain <= c + d, (name, "zero-equation count", gain, c, d))
                floor = max(0, r - k - c - d)
                require(floor <= s <= min(hp.bit_count(), r + k),
                        (name, "ordinary floor", r, k, c, d, s, hp))
                c_max = (mask & (~hp | (h & rotated_h))).bit_count()
                require(c <= c_max, (name, "aggregated charge cap", c, c_max))
                require(max(0, r - k - c_max - d) <= s, (name, "aggregate floor"))
                if floor:
                    coverage["positive ordinary floor"] += 1
                if fresh:
                    a, b, outside = image_data(state, inputs, basis, variables)
                    require(gain >= c, (name, "fresh pivot rank", gain, c))
                    require(a == r - (gain - c) and b - a == k - c,
                            (name, "exact fresh projection rank", a, b, r, gain, c, k))
                    require((outside == 0) == (k == c), (name, "fresh outside support"))
                    fresh_floor = max(0, k - c, r - c - d)
                    require(fresh_floor <= s <= min(hp.bit_count(), r + k - c),
                            (name, "fresh floor", r, k, c, d, s, hp))
                    if fresh_floor > floor:
                        coverage["strict fresh improvement"] += 1
                    if k == 2 and c == 0:
                        require(s >= 2, (name, "fresh pair cannot collapse to rank1"))
                        coverage["fresh pair excludes rank1"] += 1
                    if not r:
                        require(s == k and j == mask and hp == mask,
                                (name, "fresh quiet activation"))
                        coverage["fresh quiet activation"] += 1
                    if not hp:
                        require(fresh_floor == 0, (name, "fresh reset floor"))
                if not hp:
                    require(floor == 0, (name, "ordinary reset floor"))
                    coverage["reset branch"] += 1
                coverage["retained branches"] += 1
        coverage["form fixtures"] += 1
    require(coverage["form fixtures"] == 12, "fixture count")
    for key in ("positive ordinary floor", "strict fresh improvement", "fresh pair excludes rank1",
                "fresh quiet activation", "reset branch"):
        require(coverage[key] > 0, ("missing coverage", key))
    print("PASS: exact active-rank floors, projection ranks, and aggregate charge caps:", dict(sorted(coverage.items())))
    print("Local fixture checks only; no whole-code or numerical margin claim.")


if __name__ == "__main__":
    main()
