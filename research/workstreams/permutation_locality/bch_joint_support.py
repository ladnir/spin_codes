"""Apply authenticated BCH shell caps to grouped-support/prefix bounds.

This checks one active 16-row group, not full SPIN distance. No files are
written. Run from any directory with python-flint installed.
"""

from fractions import Fraction as F
from itertools import accumulate, combinations, product
from math import comb, log2, prod
from pathlib import Path
import random
import re
import sys

from flint import fmpz_poly
from joint_support import gf2_rank, joint_bound, self_test, span
from shortened_bound import dimension_caps

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
REPO = RESEARCH.parent


def authenticated_caps():
    sys.path.insert(0, str(RESEARCH / 'workstreams/inner_design/asymmetric/bch256'))
    import verify_progress as verifier
    dependencies = verifier.outer_caps()
    counts = verifier.model.outer.inputs.caps_module.caps()
    spectrum = [counts.get(w, 0) for w in range(257)]
    spectrum[0] = 1
    assert next(w for w in range(1, 257) if spectrum[w]) == 38
    print(f'BCH caps authenticated: {len(dependencies)} dependency files', flush=True)
    return spectrum


def zero_feedback_probability():
    source = (REPO / 'spin/src/kernels/generated/SelectedMaps.h').read_text()
    body = source.split('struct Map128S19 {', 1)[1]
    raw = re.search(r'feedbackColumns\{([^}]+)\}', body).group(1)
    columns = [int(value.strip(), 0) for value in raw.split(',')]
    # Current production uses weight5_seed0, not the older greedy3_2 map.
    # Reconstruct that candidate and enumerate the actual production columns.
    pool = [sum(1 << j for j in support) for support in combinations(range(19), 5)]
    expected = random.Random(0).sample(pool, 128)
    assert columns == expected and len(columns) == 128
    zero = [0] * 17
    for start in range(0, 128, 16):
        syndromes = [0] * (1 << 16)
        for mask in range(1, 1 << 16):
            bit = (mask & -mask).bit_length() - 1
            syndromes[mask] = syndromes[mask & (mask - 1)] ^ columns[start + bit]
            zero[mask.bit_count()] += syndromes[mask] == 0
    assert zero == [0] * 8 + [2] + [0] * 8
    q = max(F(zero[a], 8 * comb(16, a)) for a in range(1, 17))
    assert q == F(1, 51480)
    print(f'Production feedback census: worst nonempty zero probability {q}', flush=True)
    return q


def rank_total(k, g, h):
    """Number of g-by-k binary matrices of rank h, with exact division."""
    numerator = prod(((1 << k) - (1 << j)) * ((1 << g) - (1 << j)) for j in range(h))
    denominator = prod((1 << h) - (1 << j) for j in range(h))
    quotient, remainder = divmod(numerator, denominator)
    assert remainder == 0
    return quotient


def support_caps(spectrum, g=16, dimensions=None):
    polynomial = fmpz_poly([0] + spectrum[1:])
    power = fmpz_poly([1])
    result = []
    for h in range(1, g + 1):
        power *= polynomial
        cumulative = list(accumulate(int(v) for v in power.coeffs()))
        total = rank_total(128, g, h)
        bounds = [min(total, joint_bound(spectrum, g, h, u, cumulative=cumulative)[0])
                  for u in range(257)]
        if dimensions is not None:
            for u in range(257):
                # Every tuple with support <=u lies in a shortened code on
                # some u-subset. Overcount those subsets, using a uniform
                # upper dimension on each shortened code.
                cap = (comb(256,u) * rank_total(dimensions[u],g,h)
                       if dimensions[u] >= h else 0)
                bounds[u] = min(bounds[u], cap)
        # Monotonicity of the unknown true CDF permits caps at larger u to
        # improve earlier ones. Never interpret differences as true shell caps.
        for u in range(255, -1, -1):
            bounds[u] = min(bounds[u], bounds[u + 1])
        assert bounds[0] == 0 and bounds[-1] == total
        result.append(bounds)
        selected = ', '.join(f'u{u}:{log2(bounds[u]):.2f}' if bounds[u] else f'u{u}:zero'
                             for u in (38, 57, 80, 96, 128, 192, 256))
        minimum = next(u for u, v in enumerate(bounds) if v)
        print(f'rank {h:2}: support lower bound {minimum:2}; log2 CDF caps {selected}', flush=True)
    return result


def prefix_probabilities(n, length, q):
    """E[q^J], J=intersection of a uniform u-subset with a fixed prefix."""
    if not 0 <= length <= n or not 0 <= q <= 1:
        raise ValueError('invalid prefix or probability')
    # Common denominator avoids a rational operation for every summand.
    a, b = q.numerator, q.denominator
    weighted = [comb(length, j) * a**j * b**(length-j) for j in range(length + 1)]
    free = [comb(n-length, j) for j in range(n-length + 1)]
    polynomial = fmpz_poly(weighted) * fmpz_poly(free)
    numerators = [polynomial[u] for u in range(n+1)]
    result = [F(int(v), b**length * comb(n, u)) for u, v in enumerate(numerators)]
    assert len(result) == n + 1 and result[0] == 1
    assert all(a >= b for a, b in zip(result, result[1:]))
    return result


def weighted_cdf_bound(caps, decreasing_weights):
    """Summation by parts: only nonnegative multiples of upper CDFs."""
    assert len(caps) == len(decreasing_weights)
    assert all(a <= b for a, b in zip(caps, caps[1:]))
    assert all(a >= b >= 0 for a, b in zip(decreasing_weights, decreasing_weights[1:]))
    return (caps[-1] * decreasing_weights[-1]
            + sum((caps[u] * (decreasing_weights[u] - decreasing_weights[u + 1])
                   for u in range(len(caps) - 1)), F(0)))


def prefix_tests():
    for n in range(1, 9):
        for length in range(n + 1):
            q = F(1, 3)
            probabilities = prefix_probabilities(n, length, q)
            for u in range(n + 1):
                actual = sum((q**sum(j < length for j in support)
                              for support in combinations(range(n), u)), F(0)) / comb(n, u)
                assert probabilities[u] == actual
            counts = list(range(n + 1))
            caps = list(accumulate(counts))
            direct = sum((a * b for a, b in zip(counts, probabilities)), F(0))
            assert weighted_cdf_bound(caps, probabilities) == direct
            assert weighted_cdf_bound([a + 1 for a in caps], probabilities) >= direct
            for q in (F(0), F(1)):
                edge = prefix_probabilities(n,length,q)
                for u in range(n+1):
                    expected = (F(comb(n-length,u),comb(n,u)) if u <= n-length else F(0)) if not q else F(1)
                    assert edge[u] == expected
    print('Prefix and CDF checks passed for all subset sizes at n=1..8', flush=True)


def shortened_count_tests():
    for rows, n, d in (([0b1111000,0b1100110,0b1010101],7,4),
                       ([0b1001,0b1010,0b1100],4,2)):
        words = span(rows)
        assert all(w.bit_count() % 2 == 0 for w in words)
        caps = dimension_caps(n,len(rows),d,last_lp=n)
        for u in range(n+1):
            for support in combinations(range(n),u):
                mask = sum(1 << j for j in support)
                size = sum(w & ~mask == 0 for w in words)
                assert size <= 1 << caps[u]
        actual = [[0] * (n+1) for _ in range(4)]
        for tup in product(words,repeat=3):
            union = tup[0] | tup[1] | tup[2]
            actual[gf2_rank(tup)][union.bit_count()] += 1
        for h in range(1,4):
            count = 0
            for u in range(n+1):
                count += actual[h][u]
                upper = comb(n,u)*rank_total(caps[u],3,h) if caps[u] >= h else 0
                assert count <= upper
    print('Shortened dimensions and tuple counts passed exhaustive small-code checks', flush=True)


def main():
    self_test()
    prefix_tests()
    shortened_count_tests()
    spectrum = authenticated_caps()
    q = zero_feedback_probability()
    dimensions = dimension_caps()
    caps = support_caps(spectrum, dimensions=dimensions)
    groups = 8192 // 16
    print('Prefix diagnostic: 512 possible single active groups; untruncated log2 union bound', flush=True)
    for length in (180, 192, 204, 216, 224, 240, 255):
        probabilities = prefix_probabilities(256, length, q)
        bounds = [groups * weighted_cdf_bound(row, probabilities) for row in caps]
        total = sum(bounds, F(0))
        worst = max(range(16), key=lambda j: bounds[j])
        log = lambda value: log2(value.numerator) - log2(value.denominator)
        print(f'prefix={length}: all ranks={log(total):.3f}; rank1={log(bounds[0]):.3f}; '
              f'largest rank={worst+1} ({log(bounds[worst]):.3f})', flush=True)
        if length == 204:
            assert 0 < total < F(1, 1 << 41)
            print('EXACT CHECK: one-active-group zero-prefix union < 2^-41 at 204 regions', flush=True)
    print('No distance certificate: first activation does not control later cancellation/output weight.')


if __name__ == '__main__':
    main()
