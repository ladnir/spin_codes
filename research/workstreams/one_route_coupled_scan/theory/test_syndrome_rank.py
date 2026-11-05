#!/usr/bin/env python3
"""Bounded exact GF16 checks for SINGLE_DIAGONAL_RANK.md; no files written.

The root agent owns execution. Run with python -B. These are algebra and
counting checks, not a distance certificate or a benchmark.
"""

from collections import Counter
from itertools import combinations, product
from random import Random

from test_persistent_mds import INV, MUL, Q, dot, inverse, matmul, matvec, power, rank, require


def reduced(rows, width):
    data = [list(row) for row in rows]
    require(all(len(row) == width for row in data), "ragged subspace")
    pivots = []
    for column in range(width):
        pivot = next((i for i in range(len(pivots), len(data)) if data[i][column]), None)
        if pivot is None:
            continue
        used = len(pivots)
        data[used], data[pivot] = data[pivot], data[used]
        scale = INV[data[used][column]]
        data[used] = [MUL[scale][x] for x in data[used]]
        for i in range(len(data)):
            if i != used and data[i][column]:
                scale = data[i][column]
                data[i] = [x ^ MUL[scale][y] for x, y in zip(data[i], data[used])]
        pivots.append(column)
    return tuple(tuple(row) for row in data[:len(pivots)]), tuple(pivots)


def nullspace(rows, width):
    data, pivots = reduced(rows, width)
    basis = []
    for free in range(width):
        if free in pivots:
            continue
        vector = [0] * width
        vector[free] = 1
        for row, pivot in zip(data, pivots):
            vector[pivot] = row[free]  # Minus equals plus in characteristic two.
        basis.append(tuple(vector))
    return tuple(basis)


def span(rows, width):
    return reduced(rows, width)[0]


def intersection(left, right, width):
    columns = tuple(left) + tuple(right)
    equations = tuple(tuple(column[i] for column in columns) for i in range(width))
    relations = nullspace(equations, len(columns))
    return span((tuple(dot(relation[:len(left)], tuple(row[j] for row in left))
                       for j in range(width)) for relation in relations), width)


def indices(mask, width):
    return tuple(i for i in range(width) if mask >> i & 1)


def columns(matrix, mask):
    return tuple(tuple(row[j] for row in matrix) for j in indices(mask, len(matrix)))


def matrices(width, batches):
    require(width * batches <= Q, "toy global points repeat")
    return tuple(tuple(tuple(power(width * t + j, i) for j in range(width))
                       for i in range(width)) for t in range(batches))


def relaxed_step(states, h, w, overlap, width):
    result = set()
    for dimension, d, ell, core in states:
        upper = min(d, max(0, ell + h - width), d - min(core, width - h))
        for e in range(1 if h else 0, upper + 1):
            if not h and e:
                continue
            for shared in range(max(0, e + overlap - h), min(e, overlap) + 1):
                m_next = dimension + w - d + e
                d_next = e + w - shared
                ell_next = w if d_next == w else h + w - overlap
                require(m_next >= d_next >= 0 and d_next <= ell_next <= width,
                        "invalid relaxed dimension")
                result.add((m_next, d_next, ell_next, w))
    return result


def partition_rank(state_masks, output_masks, width, closed=False):
    """Best selected-time lower bound; tiny fixtures only, not a full-code DP."""
    count = len(state_masks)
    prefix = [0]
    for mask in output_masks:
        prefix.append(prefix[-1] + mask.bit_count())
    best = [0] * count
    for t in range(1, count):
        best[t] = max(best[s] + min(width - state_masks[t].bit_count(), prefix[t] - prefix[s])
                      for s in range(t))
    if closed:
        return max(best[s] + min(width, prefix[count] - prefix[s]) for s in range(count))
    return max(best)


def check_masks(width, state_masks, output_masks):
    batches = len(state_masks)
    require(batches == len(output_masks), "mask length mismatch")
    maps = matrices(width, batches)
    total = sum(mask.bit_count() for mask in output_masks)
    syndrome_rows = [[0] * total for _ in range(width)]
    constraints = []
    reachable = ()
    dimension = d = ell = core = cursor = shared_sum = 0
    retained = True
    relaxed = {(0, 0, 0, 0)}
    fractional_steps = 0
    for t, (hm, wm) in enumerate(zip(state_masks, output_masks)):
        h, w, overlap = hm.bit_count(), wm.bit_count(), (hm & wm).bit_count()
        observed = intersection(reachable, columns(maps[t], hm), width)
        e = len(observed)
        common = len(intersection(observed, columns(maps[t], wm), width))
        require(e <= min(d, max(0, ell + h - width)), "MDS intersection cap")
        require(e <= d - min(core, width - h), "known-column core cap")
        require(max(0, e + overlap - h) <= common <= min(e, overlap), "same-coset overlap cap")
        if h and not e:
            retained = False
        if not h:
            require(not e, "empty state mask must force zero intersection")
        fractional_steps += w - d + e < 0
        dimension += w - d + e
        shared_sum += common
        reachable = span(observed + columns(maps[t], wm), width)
        d = len(reachable)
        ell = w if d == w else h + w - overlap
        core = w
        require(dimension == shared_sum + d, "dimension telescoping identity")

        # Independent dense zero equations on every output variable.
        observations = matmul(inverse(maps[t]), syndrome_rows)
        constraints.extend(row for j, row in enumerate(observations) if not (hm >> j & 1))
        for j in indices(wm, width):
            for i in range(width):
                syndrome_rows[i][cursor] = maps[t][i][j]
            cursor += 1
        dense_rank = rank(constraints)
        require(dimension == cursor - dense_rank, "dense output dimension mismatch")
        image_dimension = rank(tuple(constraints) + tuple(syndrome_rows)) - dense_rank
        require(d == image_dimension, "dense reachable syndrome dimension mismatch")
        relaxed = relaxed_step(relaxed, h, w, overlap, width)
        if retained:
            require((dimension, d, ell, core) in relaxed, "dimension relaxation lost actual branch")

    closed_dimension = total - rank(tuple(constraints) + tuple(syndrome_rows))
    require(closed_dimension == dimension - d, "terminal rank mismatch")
    require(closed_dimension <= max(0, total - width), "global MDS reset rank")
    require(partition_rank(state_masks, output_masks, width) <= total - dimension,
            "selected-time free-terminal rank bound")
    require(partition_rank(state_masks, output_masks, width, True) <= total - closed_dimension,
            "selected-time closed-terminal rank bound")
    excursion = (retained and state_masks[0] == 0 and output_masks[0] != 0
                 and all(state_masks[t] for t in range(1, batches)))
    if excursion:
        h_total = sum(mask.bit_count() for mask in state_masks)
        rhs = total - width * (batches - 1)
        require(2 * (dimension - h_total) <= rhs + width, "free excursion amortized bound")
        require(2 * (closed_dimension - h_total) <= rhs, "closed excursion amortized bound")
    return dimension, closed_dimension, retained, fractional_steps, excursion


def check_dense_fixtures():
    fixtures = fractional = excursions = 0
    for masks in product(range(4), repeat=4):
        result = check_masks(2, masks[:2], masks[2:])
        fractional += result[3]
        excursions += result[4]
        fixtures += 1
    rng = Random(1907)
    for width, batches in ((2, 3), (4, 3), (4, 4)):
        for _ in range(256):
            hm = (0,) + tuple(rng.randrange(1, 1 << width) for _ in range(batches - 1))
            wm = tuple(rng.randrange(1, 1 << width) for _ in range(batches))
            result = check_masks(width, hm, wm)
            fractional += result[3]
            excursions += result[4]
            fixtures += 1
    # The fresh-previous-batch rank bound is zero, but delayed rank is b-1.
    for width in (2, 4):
        full = (1 << width) - 1
        dimension, _, retained, _, _ = check_masks(width, (0, full, full, 1), (full, full, 0, 0))
        require(retained and 2 * width - dimension == width - 1, "delayed-rank regression")
        fixtures += 1
    require(fractional > 0 and excursions > 0, "fixture classes missing")
    return fixtures, fractional, excursions


def check_cauchy_observation():
    entries = minors = 0
    for width in (2, 4):
        maps = matrices(width, 4)
        kappa = 1
        for point in range(1, width):
            kappa = MUL[kappa][point]

        def vanish(value):
            result = 1
            for point in range(width):
                result = MUL[result][value ^ point]
            return result

        for t in range(1, 4):
            blocks = tuple(matmul(inverse(maps[t]), maps[s]) for s in range(t))
            impulse = tuple(tuple(blocks[s][j][k] for s in range(t) for k in range(width))
                            for j in range(width))
            for j in range(width):
                for s in range(t):
                    numerator = vanish(width * t) ^ vanish(width * s)
                    require(numerator != 0, "source and observation cosets coincide")
                    for k in range(width):
                        denominator = MUL[kappa][(width * t + j) ^ (width * s + k)]
                        expected = MUL[numerator][INV[denominator]]
                        require(impulse[j][width * s + k] == expected, "scaled-Cauchy impulse formula")
                        entries += 1
            for size in range(1, width + 1):
                for rows in combinations(range(width), size):
                    for cols in combinations(range(t * width), size):
                        minor = tuple(tuple(impulse[i][j] for j in cols) for i in rows)
                        require(rank(minor) == size, "single-observation minor singular")
                        minors += 1
    return entries, minors


def support(vector):
    return sum(1 << j for j, value in enumerate(vector) if value)


def check_actual_outputs():
    width = 2
    maps = matrices(width, 2)
    transition = matmul(inverse(maps[1]), maps[0])
    counts = Counter()
    words = 0
    for first in product(range(Q), repeat=width):
        v = matvec(transition, first)
        syndrome = matvec(maps[0], first)
        for second in product(range(Q), repeat=width):
            contribution = matvec(maps[1], second)
            closed = not any(a ^ b for a, b in zip(syndrome, contribution))
            counts[((0, support(v)), (support(first), support(second)), closed)] += 1
            words += 1
    for (hm, wm, closed), count in counts.items():
        dimension, closed_dimension, retained, _, _ = check_masks(width, hm, wm)
        require(retained, "actual exact-support branch was discarded")
        bound = Q ** (closed_dimension if closed else dimension)
        require(count <= bound, "exact support count exceeds rank bound")
    require(words == Q ** 4, "output enumeration incomplete")
    return words, len(counts)


def main():
    fixtures, fractional, excursions = check_dense_fixtures()
    entries, minors = check_cauchy_observation()
    words, supports = check_actual_outputs()
    print("PASS syndrome-rank toy checks:", dict(
        dense_fixtures=fixtures, fractional_step_cases=fractional,
        eligible_excursion_fixtures=excursions,
        cauchy_entries=entries, single_observation_minors=minors,
        enumerated_outputs=words, exact_support_classes=supports))


if __name__ == "__main__":
    main()
