"""Bounded exact checks for stratified and Toeplitz-batch proof interfaces.

These are finite combinatorial regressions, not a distance certificate or
benchmark. The default path uses only the standard library and writes no
files. --with-float additionally compares selected root diagnostic helpers.
Run serially with python -B; the root agent owns all executions.
Toeplitz fixtures use GF(2), at most three outputs and four state symbols.
Batch activity fixtures are tiny output patterns, not full-size code tests.
"""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def poisson_binomial(probabilities):
    law = [F(1)]
    for p in probabilities:
        nxt = [F(0)] * (len(law) + 1)
        for j, mass in enumerate(law):
            nxt[j] += mass * (1 - p)
            nxt[j + 1] += mass * p
        law = nxt
    return law


def binomial_moment(B, f, r, p, marker):
    return sum(
        F(choose(r, j)) * p**j * (1 - p)**(r - j)
        * marker**max(0, f + j - B // 2)
        for j in range(r + 1)
    )


@lru_cache(maxsize=4096)
def envelope(B, r, mu, marker):
    """Exact maximum in (8); bounded cache contains only tiny fixtures."""
    if r == 0:
        require(mu.denominator == 1, ("deterministic occupancy", mu))
        return marker**max(0, int(mu) - B // 2)
    values = [
        binomial_moment(B, f, r, (mu - f) / r, marker)
        for f in range(B - r + 1)
        if 0 <= mu - f <= r
    ]
    require(bool(values), ("empty envelope", B, r, mu))
    return max(values)


def remove_zero_bound(moment, marker):
    require(moment >= 1, ("moment below one", moment))
    return min(moment, (moment - 1) / (1 - 1 / marker))


def check_regional_moments():
    """Enumerate all tiny regional subsets, grouping by their fixed sizes."""
    placements = profiles = comparisons = positive_shift_zero_complements = 0
    markers = (F(3, 2), F(2), F(3))
    for M, B in ((3, 2), (4, 2), (3, 4), (4, 4)):
        k = B // 2
        histograms = defaultdict(Counter)
        for regions in product(range(1 << M), repeat=B):
            sizes = tuple(mask.bit_count() for mask in regions)
            occupancies = [
                sum((mask >> i) & 1 for mask in regions)
                for i in range(M)
            ]
            u = sum(sizes)
            d = sum(max(0, t - k) for t in occupancies)
            dc = sum(max(0, B - t - k) for t in occupancies)
            require(d == u - M * k + dc,
                    ("complement identity", M, B, regions, d, dc))
            if u > M * k and dc == 0:
                require(d > 0, ("positive shift lost", regions))
                positive_shift_zero_complements += 1
            histograms[sizes][d] += 1
            placements += 1

        for sizes, hist in histograms.items():
            total = sum(hist.values())
            require(total == product_of(choose(M, s) for s in sizes),
                    ("fixed-size normalization", M, B, sizes))
            law = poisson_binomial(tuple(F(s, M) for s in sizes))
            for marker in markers:
                actual = sum(F(n, total) * marker**d for d, n in hist.items())
                nonzero = sum(F(n, total) * marker**d
                              for d, n in hist.items() if d)
                independent = sum(mass * marker**max(0, t - k)
                                  for t, mass in enumerate(law))**M
                complementary = (
                    marker**(sum(sizes) - M * k)
                    * sum(mass * marker**max(0, k - t)
                          for t, mass in enumerate(law))**M
                )
                require(actual <= independent,
                        ("regional replacement", M, B, sizes, marker))
                require(actual <= complementary,
                        ("fixed-size complement bound", M, B, sizes, marker))
                require(nonzero <= remove_zero_bound(
                    min(independent, complementary), marker),
                    ("zero exclusion", M, B, sizes, marker))
                comparisons += 1
            profiles += 1
    require(positive_shift_zero_complements > 0,
            "Missing the complement-zero negative control")
    return dict(placements=placements, size_profiles=profiles,
                moment_comparisons=comparisons,
                positive_shift_zero_complements=positive_shift_zero_complements)


def product_of(values):
    answer = 1
    for value in values:
        answer *= value
    return answer


def check_convex_envelopes():
    equalizations = enlargements = 0
    for B, alphabet in (
        (2, (F(0), F(1, 3), F(2, 3), F(1))),
        (4, (F(0), F(1, 4), F(1, 2), F(3, 4), F(1))),
    ):
        for f in range(B + 1):
            for h in range(B - f + 1):
                for probabilities in product(alphabet, repeat=h):
                    law = poisson_binomial(probabilities)
                    mu = F(f) + sum(probabilities, F(0))
                    for marker in (F(3, 2), F(2)):
                        actual = sum(
                            mass * marker**max(0, f + j - B // 2)
                            for j, mass in enumerate(law)
                        )
                        if h:
                            equal = binomial_moment(
                                B, f, h, (mu - f) / h, marker)
                            require(actual <= equal,
                                    ("pair equalization", B, f, probabilities, marker))
                            equalizations += 1
                        for r in range(h, B + 1):
                            require(actual <= envelope(B, r, mu, marker),
                                    ("enlarged partial count", B, f,
                                     probabilities, r, marker))
                            enlargements += 1
    return dict(equalizations=equalizations, partial_enlargements=enlargements)


def interval_count(N, L, u, c, kind):
    free = u - L * (c if kind == "complete" else c - 1)
    if free < c or not 1 <= u <= N:
        return 0
    gaps = (choose(N - u + c, c) if kind == "complete"
            else choose(N - u + c - 1, c - 1))
    return choose(free - 1, c - 1) * gaps


def loose_interval_count(N, L, u, c):
    free = u - L * (c - 1)
    if free < c or not 1 <= u <= N:
        return 0
    return choose(free - 1, c - 1) * choose(N - u + c, c)


def endpoint_representations(N, L, c, kind):
    """Direct endpoint recursion; adjacent intervals are permitted."""
    def visit(index, cursor, intervals):
        if index == c:
            yield tuple(intervals)
            return
        last = index == c - 1
        minimum = 1 if last and kind != "complete" else L + 1
        for start in range(cursor, N + 1):
            ends = ((N,) if last and kind == "truncated"
                    else range(start + minimum - 1, N + 1))
            for end in ends:
                if end - start + 1 >= minimum:
                    yield from visit(index + 1, end + 1,
                                     intervals + [(start, end)])
    yield from visit(0, 1, [])


def volume(q, n, e):
    if n < 0 or e < 0:
        return 0
    return sum(choose(n, j) * (q - 1)**j for j in range(min(n, e) + 1))


def output_intervals(N, L, word):
    positions = [j + 1 for j in range(N) if (word >> j) & 1]
    clusters = [[positions[0]]]
    for position in positions[1:]:
        if position - clusters[-1][-1] > L:
            clusters.append([position])
        else:
            clusters[-1].append(position)
    intervals = tuple((part[0], min(N, part[-1] + L)) for part in clusters)
    kind = "complete" if positions[-1] <= N - L else "truncated"
    return intervals, kind


def check_intervals():
    syntactic_representations = binary_outputs = volume_checks = 0
    for N in range(1, 10):
        for L in range(1, min(3, N) + 1):
            for kind in ("complete", "truncated", "loose"):
                for c in range(1, N + 1):
                    counts = Counter()
                    for intervals in endpoint_representations(N, L, c, kind):
                        counts[sum(b - a + 1 for a, b in intervals)] += 1
                        syntactic_representations += 1
                    for u in range(1, N + 1):
                        expected = (loose_interval_count(N, L, u, c)
                                    if kind == "loose"
                                    else interval_count(N, L, u, c, kind))
                        require(counts[u] == expected,
                                ("syntactic interval count", N, L, kind, u, c))
                if kind == "loose":
                    for c in range(1, N + 1):
                        require(sum(loose_interval_count(N, L, u, c)
                                    for u in range(1, N + 1))
                                == choose(N + c - L * (c - 1), 2 * c),
                                ("total loose count", N, L, c))

            actual_representations = defaultdict(set)
            actual_weights = defaultdict(Counter)
            for word in range(1, 1 << N):
                intervals, kind = output_intervals(N, L, word)
                c = len(intervals)
                u = sum(b - a + 1 for a, b in intervals)
                U = {j for a, b in intervals for j in range(a, b + 1)}
                active = {j for j in range(1, N + 1)
                          if any((word >> (i - 1)) & 1
                                 for i in range(max(1, j - L), j))}
                support = {j + 1 for j in range(N) if (word >> j) & 1}
                require(U == active | support and len(active) == u - c,
                        ("activity sets", N, L, word))
                complete = intervals if kind == "complete" else intervals[:-1]
                forced = {j for _, end in complete
                          for j in range(end - L + 1, end + 1)}
                require(not forced & support and forced <= U,
                        ("forced tail", N, L, word))
                free = u - L * (c if kind == "complete" else c - 1)
                require(len(U - forced) == free,
                        ("free-coordinate count", N, L, word))
                actual_representations[kind, u, c].add(intervals)
                actual_weights[kind, intervals][word.bit_count()] += 1
                binary_outputs += 1
            for (kind, u, c), representations in actual_representations.items():
                require(len(representations) <= interval_count(N, L, u, c, kind),
                        ("actual representations only bounded", N, L, kind, u, c))
            for (kind, intervals), weights in actual_weights.items():
                c = len(intervals)
                u = sum(b - a + 1 for a, b in intervals)
                free = u - L * (c if kind == "complete" else c - 1)
                for e in range(N + 1):
                    actual = sum(count for weight, count in weights.items() if weight <= e)
                    require(actual <= volume(2, free - c, e - c),
                            ("activation output volume", N, L, kind, intervals, e))
                    volume_checks += 1
    return dict(syntactic_representations=syntactic_representations,
                binary_outputs=binary_outputs, output_volume_checks=volume_checks)


def check_qary_volumes():
    checks = 0
    for q in (2, 3):
        for free in range(6):
            for c in range(free + 1):
                weights = Counter(sum(x != 0 for x in word)
                                  for word in product(range(q), repeat=free)
                                  if all(word[j] != 0 for j in range(c)))
                for e in range(free + 1):
                    actual = sum(n for w, n in weights.items() if w <= e)
                    bound = (q - 1)**c * volume(q, free - c, e - c)
                    require(actual == bound, ("q-ary prescribed activations", q, free, c, e))
                    checks += 1
    return dict(qary_volume_checks=checks)


def check_entropy_volumes():
    """Exact rational form of the entropy Chernoff bound before logarithms."""
    checks = 0
    for q in (2, 4, 16, 256):
        for free in range(1, 13):
            for e in range(free):
                if e > q*(free-e):
                    continue
                cap = (F(1) if e == 0 else
                       F(q**e * free**free, e**e * (free-e)**(free-e)))
                require(volume(q,free,e) <= cap,
                        ('entropy volume',q,free,e))
                checks += 1
    return dict(entropy_volume_checks=checks)


def check_box_components():
    """Exact component checks for the coarse boxes, not a large-size cover."""
    comparisons = 0
    B = 4
    for N in (12, 16):
        M = N // B
        for L in (1, 2, 3):
            cmax = (N + L) // (L + 1)
            uranges = ((1, N), (1, N // 2), (N // 2, N),
                       (N // 2 + 1, N), (3, N - 2))
            cranges = ((1, cmax), (1, min(2, cmax)), (2, min(4, cmax)))
            for ul0, uh in uranges:
                for cl, ch0 in cranges:
                    ul = max(ul0, L * (cl - 1) + cl)
                    ch = min(ch0, (uh + L) // (L + 1))
                    if ul > uh or cl > ch:
                        continue
                    free_high = uh - L * (cl - 1)
                    n1, n2 = free_high - 1, N - ul + ch
                    count_cap = (choose(n1, min(ch - 1, n1 // 2))
                                 * choose(n2, min(ch, n2 // 2)))
                    closest = min(uh, N // 2) if ul <= N // 2 else N - ul
                    mu_high = F(closest, M)
                    for u in range(ul, uh + 1):
                        for c in range(cl, ch + 1):
                            count = loose_interval_count(N, L, u, c)
                            if not count:
                                continue
                            require(count <= count_cap,
                                    ("box count maximum", N, L, ul, uh, cl, ch, u, c))
                            free = u - L * (c - 1)
                            for e in (0, N // 3, N // 2):
                                w = min(F(free), F(e) + F(free, 128))
                                wcap = min(F(free_high), F(e) + F(free_high, 128))
                                require(c - min(u, N // 2) + w
                                        <= ch - min(ul, N // 2) + wcap,
                                        ("box coefficient maximum", N, L, u, c, e))
                            for marker in (F(5, 4), F(2)):
                                H = envelope(B, min(2 * c, B),
                                             F(min(u, N - u), M), marker)
                                Hcap = envelope(B, min(2 * ch, B), mu_high, marker)
                                require(H <= Hcap,
                                        ("box route maximum", N, L, ul, uh, cl, ch, u, c))
                                routed = (H**M if u > N // 2
                                          else remove_zero_bound(H**M, marker))
                                routecap = (Hcap**M if uh > N // 2
                                            else remove_zero_bound(Hcap**M, marker))
                                require(routed <= routecap,
                                        ("box crossing midpoint", N, L, u, c))
                                comparisons += 1
    return dict(box_component_comparisons=comparisons)


def toeplitz_binary_product(b, L, diagonals, state):
    """Direct matrix definition T[a,j]=d[a-j], independent of the proof."""
    result = 0
    for a in range(b):
        value = 0
        for j in range(L):
            value ^= ((diagonals >> (a - j + L - 1)) & 1) * ((state >> j) & 1)
        result |= value << a
    return result


def check_toeplitz_uniformity():
    nonzero_states = zero_states = products = 0
    for b in range(1, 4):
        for L in range(1, 5):
            for state in range(1 << L):
                outputs = Counter()
                for diagonals in range(1 << (b + L - 1)):
                    outputs[toeplitz_binary_product(b, L, diagonals, state)] += 1
                    products += 1
                if state:
                    require(len(outputs) == 1 << b
                            and all(outputs[value] == 1 << (L - 1)
                                    for value in range(1 << b)),
                            ("joint Toeplitz uniformity", b, L, state, outputs))
                    nonzero_states += 1
                else:
                    require(outputs == {0: 1 << (b + L - 1)},
                            ("quiet Toeplitz state", b, L))
                    zero_states += 1
    return dict(toeplitz_nonzero_states=nonzero_states,
                toeplitz_zero_states=zero_states, toeplitz_products=products)


def polynomial_product(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] += x * y
    return result


def batch_volume(q, b, free, c, e):
    if free < b * c or e < c:
        return 0
    activation = [0] + [choose(b, j) * (q - 1)**j for j in range(1, b + 1)]
    coefficients = [1]
    for _ in range(c):
        coefficients = polynomial_product(coefficients, activation)
    remaining = free - b * c
    coefficients = polynomial_product(
        coefficients, [choose(remaining, j) * (q - 1)**j
                       for j in range(remaining + 1)])
    return sum(coefficients[:e + 1])


def check_batch_activity():
    patterns = count_checks = volume_checks = 0
    for b in (1, 2, 3, 4):
        for n in range(1, 7):
            for m in range(1, min(3, n) + 1):
                L = m * b
                representations = defaultdict(set)
                for pattern in range(1, 1 << n):
                    intervals, _ = output_intervals(n, m, pattern)
                    c = len(intervals)
                    v = sum(end - start + 1 for start, end in intervals)
                    U = {j * b + a for start, end in intervals
                         for j in range(start - 1, end) for a in range(b)}
                    activations = {b * (start - 1) + a
                                   for start, _ in intervals for a in range(b)}
                    forced = {j * b + a for _, end in intervals[:-1]
                              for j in range(end - m, end) for a in range(b)}
                    free = b * v - L * (c - 1)
                    require(len(U - forced) == free and free >= b * c,
                            ("batch free coordinates", b, n, m, pattern))
                    for dense in ((False,) if b == 1 else (False, True)):
                        batches = [
                            ((1 << b) - 1 if dense else 1 << (j % b))
                            if (pattern >> j) & 1 else 0
                            for j in range(n)
                        ]
                        active = {j * b + a for j in range(n)
                                  if any(batches[max(0, j - m):j])
                                  for a in range(b)}
                        expanded_nonzero = {j * b + a for j in range(n)
                                            if batches[j] for a in range(b)}
                        scalar_support = {j * b + a for j in range(n)
                                          for a in range(b) if (batches[j] >> a) & 1}
                        require(U == active | expanded_nonzero,
                                ("expanded batch activity union", b, n, m, pattern))
                        require(active == U - activations
                                and len(active) == b * (v - c),
                                ("remove whole activation batches", b, n, m, pattern))
                        require(not forced & scalar_support and forced <= U,
                                ("batch forced zero tails", b, n, m, pattern))
                        components = sum(j - 1 not in active for j in active)
                        require(components <= c,
                                ("A remains an interval union", b, n, m, pattern))
                        patterns += 1
                    representations[v, c].add(intervals)
                for (v, c), actual in representations.items():
                    require(len(actual) <= loose_interval_count(n, m, v, c),
                            ("batch representation count", b, n, m, v, c))
                    count_checks += 1

    # Count the entire tiny free-coordinate assignment space independently
    # of the polynomial formula, including sparse nonzero activation batches.
    for q in (2, 3):
        for b in (1, 2, 3):
            for free in range(7):
                for c in range(min(2, free // b) + 1):
                    weights = Counter(
                        sum(symbol != 0 for symbol in word)
                        for word in product(range(q), repeat=free)
                        if all(any(word[j * b:(j + 1) * b]) for j in range(c))
                    )
                    for e in range(free + 1):
                        actual = sum(count for weight, count in weights.items() if weight <= e)
                        require(actual == batch_volume(q, b, free, c, e),
                                ("nonzero-batch volume", q, b, free, c, e))
                        volume_checks += 1
    return dict(batch_activity_patterns=patterns,
                batch_count_checks=count_checks, batch_volume_checks=volume_checks)


def batch_boxes():
    """Small clamped rectangles in the directed batch checker's coordinates."""
    B = 4
    for b in (2, 3, 4):
        for n in (4, 8):
            N = b * n
            for m in (1, 2, 3):
                L = b * m
                cmax = (n + m) // (m + 1)
                vranges = ((1, n), (1, n // 2), (n // 2, n),
                           (n // 2 + 1, n), (2, n - 1))
                cranges = ((1, cmax), (1, min(2, cmax)), (2, min(4, cmax)))
                for vl0, vh in vranges:
                    for cl, ch0 in cranges:
                        vl = max(vl0, m * (cl - 1) + cl)
                        ch = min(ch0, (vh + m) // (m + 1))
                        if vl <= vh and cl <= ch:
                            yield N, B, L, b, vl, vh, cl, ch


def check_batch_box_components():
    comparisons = zero_A_moments = 0
    for N, B, L, b, vl, vh, cl, ch in batch_boxes():
        n, m, M = N // b, L // b, N // B
        free_high = b * vh - L * (cl - 1)
        n1, n2 = vh - m * (cl - 1) - 1, n - vl + ch
        count_cap = (choose(n1, min(ch - 1, n1 // 2))
                     * choose(n2, min(ch, n2 // 2)))
        al, ah = max(0, b * (vl - ch)), b * (vh - cl)
        closest = min(ah, N // 2) if al <= N // 2 else N - al
        for v in range(vl, vh + 1):
            for c in range(cl, ch + 1):
                count = loose_interval_count(n, m, v, c)
                if not count:
                    continue
                active = b * (v - c)
                free = b * v - L * (c - 1)
                require(al <= active <= ah,
                        ("batch A size range", N, L, b, vl, vh, cl, ch, v, c))
                require(count <= count_cap,
                        ("batch box count", N, L, b, vl, vh, cl, ch, v, c))
                for e in (0, N // 3, N // 2):
                    w = min(F(free), F(e) + F(free, 128))
                    wcap = min(F(free_high), F(e) + F(free_high, 128))
                    require(w - min(active, N // 2) <= wcap - min(al, N // 2),
                            ("batch random-coordinate coefficient", N, L, b, v, c, e))
                for marker in (F(5, 4), F(2)):
                    H = envelope(B, min(2 * c, B),
                                 F(min(active, N - active), M), marker)
                    Hcap = envelope(B, min(2 * ch, B), F(closest, M), marker)
                    require(H <= Hcap,
                            ("batch box route", N, L, b, vl, vh, cl, ch, v, c))
                    # D(A)=0 need not exclude an affine inverse fiber. The
                    # batch certificate must retain its whole moment here.
                    require(1 <= H**M <= Hcap**M,
                            ("retain zero shortened A", N, L, b, v, c))
                    zero_A_moments += H == 1
                    comparisons += 1
    require(zero_A_moments > 0, "Missing the D(A)=0 moment regression")
    return dict(batch_box_component_comparisons=comparisons,
                batch_zero_A_moments=zero_A_moments)


def check_float_helpers():
    """Optional comparisons only; exact checks above have no NumPy dependency."""
    import math
    import stratified_fractional as diagnostic

    checks = 0

    def close_log(actual, expected, context):
        if expected == 0:
            require(actual == float('-inf'), context)
        else:
            reference = math.log(float(expected))
            require(abs(actual - reference) <= 2e-11 * max(1, abs(reference)), context)

    for B in (2, 4):
        for c in (1, 2, 3):
            for numerator in range(4 * B + 1):
                mu = F(numerator, 4)
                for marker in (F(5, 4), F(2), F(3)):
                    theta = math.log(float(marker))
                    H = envelope(B, min(2 * c, B), mu, marker)
                    excess = diagnostic.excess_log(B, B // 2, c, float(mu), theta)
                    close_log(excess, H - 1, ("float binomial excess", B, c, mu, marker))
                    for M in (3, 4):
                        close_log(diagnostic.route_log(M, excess, theta),
                                  remove_zero_bound(H**M, marker),
                                  ("float route moment", B, c, mu, marker, M))
                        checks += 1
    for N in range(1, 10):
        for L in (1, 2, 3):
            for u in range(1, N + 1):
                for c in range(1, N + 1):
                    close_log(diagnostic.representation_log(N, L, u, c),
                              F(loose_interval_count(N, L, u, c)),
                              ("float interval count", N, L, u, c))
                    checks += 1
    return dict(optional_float_comparisons=checks)


def check_float_batch_boxes():
    import math
    import stratified_batch_certificate as diagnostic

    comparisons = 0
    for N, B, L, b, vl, vh, cl, ch in batch_boxes():
        n, m, M = N // b, L // b, N // B
        e = N // 3
        bound, raw_marker = diagnostic.float_box(N, B, L, b, e, vl, vh, cl, ch)
        if b * vh * B + 2 * ch * N <= (B // 2) * N:
            require(bound == float('-inf'), ("batch zero guard", N, L, b, vl, vh, cl, ch))
            continue
        marker = F(*raw_marker)
        for v in range(vl, vh + 1):
            for c in range(cl, ch + 1):
                count = loose_interval_count(n, m, v, c)
                if not count:
                    continue
                active = b * (v - c)
                free = b * v - L * (c - 1)
                coefficient = min(F(free), F(e) + F(free, 128)) - min(active, N // 2)
                H = envelope(B, min(2 * c, B),
                             F(min(active, N - active), M), marker)
                point = (math.log2(count) + float(coefficient) * math.log2(float(marker))
                         + M * math.log2(float(H)))
                require(point <= bound + 2e-10 * max(1, abs(point), abs(bound)),
                        ("floating batch box above exact component point",
                         N, L, b, vl, vh, cl, ch, v, c, point, bound))
                comparisons += 1
    return dict(optional_batch_float_comparisons=comparisons)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--with-float', action='store_true',
                        help='also import the NumPy-based floating diagnostic helpers')
    args = parser.parse_args()
    summary = {}
    for check in (check_regional_moments, check_convex_envelopes,
                  check_intervals, check_qary_volumes, check_entropy_volumes, check_box_components,
                  check_toeplitz_uniformity, check_batch_activity,
                  check_batch_box_components):
        summary.update(check())
    if args.with_float:
        summary.update(check_float_helpers())
        summary.update(check_float_batch_boxes())
    print('PASS stratified fractional small checks', summary)


if __name__ == '__main__':
    main()
