"""Scalar eight-state Cauchy repair; no kernel, benchmark, or distance claim."""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import combinations
from math import comb
from random import Random

from rate_design_coset import mul, power, xor_sum
from region64_reference import phase_geometry


H16, PHASES = phase_geometry()
W8 = H16[:8]
OTHER = next(x for x in H16 if x not in W8)
LABELS = tuple(tuple(p ^ c ^ w for c in (0, OTHER) for w in W8) for p in PHASES)
BETA = next(x for x in range(256) if x not in set(sum(LABELS, ())))
NODES = tuple(BETA ^ w for w in W8)
ETA = 2


def linearized(x):
    return power(x, 8) ^ power(x, 4) ^ power(x, 2) ^ x


def emit_coefficient(x):
    return power(x, 7) ^ power(x, 3) ^ x ^ 1


def coefficients(offset, emission):
    function = emit_coefficient if emission else lambda x: power(x, 7)
    return tuple(function(offset ^ w) for w in W8)


def zeta(values):
    result = list(values)
    for bit in (1, 2, 4):
        for mask in range(8):
            if not mask & bit:
                result[mask] ^= result[mask | bit]
    return tuple(result)


def circulant(values, fixed):
    """Nineteen products, gamma[0]=1; two input and two output blocks pack."""
    transformed, gamma = zeta(values), zeta(fixed)
    assert gamma[0] == 1
    result = list(transformed)
    for mask in range(1, 8):
        subset = mask
        while subset:
            result[mask] ^= mul(gamma[subset], transformed[mask ^ subset])
            subset = (subset - 1) & mask
    return zeta(result)


@lru_cache(maxsize=None)
def matrices(phase=0, symmetric=False):
    a = tuple(tuple(emit_coefficient(u ^ v) for v in NODES) for u in LABELS[phase])
    c = tuple(tuple(emit_coefficient(u ^ v) if symmetric else power(u ^ v, 7)
                    for u in LABELS[phase]) for v in NODES)
    return a, c


def matvec(matrix, vector):
    return tuple(xor_sum(mul(a, b) for a, b in zip(row, vector)) for row in matrix)


def emission(state, phase):
    return sum((circulant(state, coefficients(LABELS[phase][first] ^ BETA, True))
                for first in (0, 8)), ())


def injection(values, phase, symmetric=False):
    blocks = [circulant(values[first:first+8],
                        coefficients(LABELS[phase][first] ^ BETA, symmetric))
              for first in (0, 8)]
    return tuple(a ^ b for a, b in zip(*blocks))


def step(values, state, diagonal, phase, symmetric=False):
    refreshed = tuple(mul(a, b) for a, b in zip(state, diagonal))
    output = tuple(x ^ mul(ETA, a) for x, a in zip(values, emission(refreshed, phase)))
    following = tuple(a ^ b for a, b in zip(refreshed, injection(values, phase, symmetric)))
    return output, following


def adjoint(outputs, terminal, diagonal, phase, symmetric=False):
    a, c = matrices(phase, symmetric)
    inputs = tuple(x ^ y for x, y in zip(outputs, matvec(tuple(zip(*c)), terminal)))
    prior = tuple(mul(d, h ^ mul(ETA, y)) for d, h, y in
                  zip(diagonal, terminal, matvec(tuple(zip(*a)), outputs)))
    return inputs, prior


@lru_cache(maxsize=1)
def multiplication_table():
    return tuple(bytes(mul(a, b) for b in range(256)) for a in range(256))


def rank(matrix):
    if not matrix:
        return 0
    table = multiplication_table()
    work = [list(row) for row in matrix]
    found = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(found, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[found], work[pivot] = work[pivot], work[found]
        row = work[found]
        reciprocal = power(row[column], 254)
        products = table[reciprocal]
        row[column:] = [products[v] for v in row[column:]]
        for i in range(found+1, len(work)):
            coefficient = work[i][column]
            if coefficient:
                products = table[coefficient]
                work[i][column:] = [a ^ products[b] for a, b in zip(work[i][column:], row[column:])]
        found += 1
        if found == len(work):
            break
    return found


def subset_update_rank(indices, eta=ETA, symmetric=False):
    # det(I+C_J A_J)=det(I+A_J C_J). Pairing depends only on label differences.
    points = tuple(LABELS[0][i] for i in indices)
    if symmetric:
        pair = lambda x: emit_coefficient(x)
    else:
        pair = lambda x: power(x, 7)
    gram = tuple(tuple(int(i == j) ^ mul(eta, pair(u ^ v))
                       for j, v in enumerate(points)) for i, u in enumerate(points))
    return rank(gram)


def check_all_subsets(eta=ETA, symmetric=False):
    # CA=0 identifies complementary updates. Retain one of each pair.
    checked, singular = 0, []
    for size in range(9):
        for indices in combinations(range(16), size):
            if size == 8 and 0 not in indices:
                continue
            checked += 1
            if subset_update_rank(indices, eta, symmetric) != size:
                singular.append(indices)
    assert checked == 32768
    return {"eta": eta, "symmetric": symmetric, "complement_classes": checked,
            "singular_count": len(singular), "first_singular": singular[:3]}


def eta_obstruction_orbits():
    """Two exact bad masks and algebraic symmetries cover every eta outside H."""
    witnesses = ((2, (0, 1, 4, 10, 12)), (7, (0, 3, 8, 12)))
    orbits = []
    for eta, indices in witnesses:
        assert subset_update_rank(indices, eta) < len(indices)
        orbit = {mul(h, power(eta, 2**j)) for h in H16 if h for j in range(8)}
        orbits.append(orbit)
    assert len(orbits[0]) == len(orbits[1]) == 120
    assert orbits[0].isdisjoint(orbits[1])
    assert orbits[0] | orbits[1] == set(range(256)) - set(H16)
    # For eta in H*, choose label difference d with eta*d^7=1.
    for eta in H16:
        if eta:
            difference = next(d for d in H16 if mul(eta, power(d, 7)) == 1)
            assert subset_update_rank((0, LABELS[0].index(difference)), eta) == 1
    return tuple(tuple(sorted(orbit)) for orbit in orbits)


def binary_rank(values):
    pivots = {}
    for value in values:
        while value:
            bit = value.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
    return len(pivots)


def shared_product_counts():
    """Scalar fixed products in coefficient-basis circuits, before packing."""
    result = []
    for phase in range(4):
        counts = []
        for emit in (False, True):
            gamma = [zeta(coefficients(LABELS[phase][first] ^ BETA, emit)) for first in (0, 8)]
            counts.append(tuple(binary_rank([1] + [g[s] for g in gamma for s in range(1, 8)
                                                  if not s & t]) - 1 for t in range(8)))
        result.append(tuple(counts))
    return tuple(result)


def symmetric_return_envelope(weight, y, tau, incoming_size=None,
                              flat_cover=True, state_rank=8):
    """Pure-return alive-rank cover; no selected-zero factor is multiplied in.

    A known incoming actual support pays one branch. If its size is unknown,
    the positive cover sums all feasible nonempty actual support masks.
    """
    if not 0 <= weight <= 16:
        raise ValueError("packet weight is outside 0..16")
    if not 0 <= state_rank <= 8:
        raise ValueError("state-form rank is outside 0..8")
    if not weight or not state_rank:
        return y * 0
    lower = 9 if weight == 1 else max(1, 9-weight)
    if incoming_size is not None:
        if not 1 <= incoming_size <= 8:
            raise ValueError("incoming support size is outside 1..8")
        multiplicity = int(weight > 8-incoming_size)
    elif not flat_cover:
        multiplicity = sum(comb(8, h) for h in range(max(1, 9-weight), 9))
    else:
        multiplicity = sum(comb(8, k) / tau**k for k in range(min(weight-1, state_rank-1)+1))
    return multiplicity * y**lower / tau**8


def small_return_geometry():
    """All placements/phases; count zeros of a three-column torus combination.

    Every two-by-two minor is nonzero. Hence a combination with at least two
    zero residues lies on the unique torus line obtained from those two rows.
    Enumerating all pairs of residue rows exhausts every exceptional line.
    """
    table = multiplication_table()
    totals = []
    for phase in range(4):
        _, c = matrices(phase)
        maximum_zeros, min_exponent, line_count = 0, 8, 0
        output_min, witness = 16, None
        for inputs in combinations(range(16), 3):
            rows = [tuple(row[i] for i in inputs) for row in c]
            lines = set()
            for j, k in combinations(range(8), 2):
                a, b = rows[j], rows[k]
                vector = (table[a[1]][b[2]] ^ table[a[2]][b[1]],
                          table[a[2]][b[0]] ^ table[a[0]][b[2]],
                          table[a[0]][b[1]] ^ table[a[1]][b[0]])
                assert all(vector), "violates the seventh-power two-minor theorem"
                inverse = table[power(vector[0], 254)]
                lines.add(tuple(inverse[x] for x in vector))
            for vector in lines:
                zeros = sum(not xor_sum(table[a][b] for a, b in zip(row, vector)) for row in rows)
                assert zeros >= 2
                line_count += 1
                if zeros > maximum_zeros:
                    maximum_zeros, witness = zeros, (inputs, vector)
                min_exponent = min(min_exponent, 3 + 8 - zeros - 1)
                y = []
                for j, u in enumerate(LABELS[phase]):
                    emitted = xor_sum(table[power(u ^ LABELS[phase][i], 7)][x] for i, x in zip(inputs, vector))
                    direct = vector[inputs.index(j)] if j in inputs else 0
                    y.append(direct ^ table[ETA][emitted])
                output_min = min(output_min, sum(bool(x) for x in y))
        totals.append({"phase": phase, "placements": 560, "exceptional_lines": line_count,
                       "maximum_injection_zeros": maximum_zeros,
                       "minimum_diagonal_exponent": min_exponent,
                       "minimum_output_on_exceptional_lines": output_min,
                       "maximum_zero_witness": witness})
    return totals


def check_maps():
    assert len(set(W8)) == 8 and all(linearized(x) == 0 for x in W8)
    assert all(linearized(x) != 0 for x in H16 if x not in W8)
    rng = Random(20261022)
    checks = 0
    for phase in range(4):
        a, c = matrices(phase)
        for u, row in zip(LABELS[phase], a):
            assert linearized(u ^ BETA) != 0
            assert all(mul(value, u ^ v) == linearized(u ^ BETA)
                       for value, v in zip(row, NODES))
        for j in range(8):
            state = tuple(int(i == j) for i in range(8))
            assert emission(state, phase) == matvec(a, state)
            assert injection(emission(state, phase), phase) == (0,) * 8
            checks += 1
        for j in range(16):
            values = tuple(int(i == j) for i in range(16))
            assert injection(values, phase) == matvec(c, values)
            assert emission(injection(values, phase), phase) == tuple(power(u ^ LABELS[phase][j], 7) for u in LABELS[phase])
            checks += 1
        for _ in range(16):
            values = tuple(rng.randrange(256) for _ in range(16))
            state = tuple(rng.randrange(256) for _ in range(8))
            diagonal = tuple(rng.randrange(1, 256) for _ in range(8))
            outputs = tuple(rng.randrange(256) for _ in range(16))
            terminal = tuple(rng.randrange(256) for _ in range(8))
            y, following = step(values, state, diagonal, phase)
            dualx, duals = adjoint(outputs, terminal, diagonal, phase)
            dot = lambda x, y: xor_sum(mul(a, b) for a, b in zip(x, y))
            assert dot(y, outputs) ^ dot(following, terminal) == dot(values, dualx) ^ dot(state, duals)
            restored, repeat = step(y, state, diagonal, phase)
            assert restored == values and repeat == following
            checks += 1
    return checks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all-subsets", action="store_true")
    parser.add_argument("--symmetric", action="store_true")
    parser.add_argument("--small-returns", action="store_true")
    args = parser.parse_args()
    print({"map_checks": check_maps(), "W8": W8, "nodes": NODES,
           "packet_products_upper": 27, "full_volume_upper": "257/88"})
    if args.all_subsets:
        print(check_all_subsets(symmetric=args.symmetric))
    print({"eta_obstruction_orbit_sizes": tuple(map(len, eta_obstruction_orbits())),
           "shared_scalar_products_C_A": tuple(tuple(sum(x) for x in phase) for phase in shared_product_counts())})
    if args.small_returns:
        print(small_return_geometry())
