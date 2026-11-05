#!/usr/bin/env python3
"""Bounded exact checks for the persistent MDS-state proof interfaces.

The toy field is GF(16)=GF(2)[x]/(x^4+x+1), with primitive alpha=x.
For batch=state=2, use P=diag(1,alpha) and C[i,j]=beta_j^i, where
beta=(1,alpha^5) lies in GF(4). The recurrence is Y=X+T*S,
S_next=P*S+C*Y, with initial state zero and no terminal constraint.

These finite algebra and counting checks are not a distance certificate,
a large-field kernel test, or a performance result. This standard-library
script writes no files. Run with python -B; the root agent owns execution.
"""

from collections import Counter
from itertools import combinations, product
from math import comb
from random import Random


Q = 16
MODULUS = 0x13
ALPHA = 2


def require(condition, message):
    """Retain every check when Python is invoked with -O."""
    if not condition:
        raise AssertionError(message)


def multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & Q:
            a ^= MODULUS
    return result


MUL = tuple(tuple(multiply(a, b) for b in range(Q)) for a in range(Q))


def power(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = MUL[result][a]
        a = MUL[a][a]
        exponent >>= 1
    return result


INV = (0,) + tuple(power(a, Q - 2) for a in range(1, Q))
BETA = (1, power(ALPHA, 5))
POINTS = (0, 1, BETA[1], BETA[1] ^ 1)
P = ((1, 0), (0, ALPHA))
C = ((1, 1), BETA)


def dot(left, right):
    require(len(left) == len(right), "dot-product dimension mismatch")
    result = 0
    for a, b in zip(left, right):
        result ^= MUL[a][b]
    return result


def add(left, right):
    require(len(left) == len(right), "vector dimension mismatch")
    return tuple(a ^ b for a, b in zip(left, right))


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def matmul(left, right):
    columns = tuple(zip(*right))
    return tuple(tuple(dot(row, column) for column in columns) for row in left)


def rank(matrix):
    """Ordinary field elimination, independent of the Vandermonde argument."""
    rows = [list(row) for row in matrix]
    if not rows:
        return 0
    columns = len(rows[0])
    require(all(len(row) == columns for row in rows), "ragged matrix")
    used = 0
    for column in range(columns):
        pivot = next((i for i in range(used, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[used], rows[pivot] = rows[pivot], rows[used]
        scale = INV[rows[used][column]]
        rows[used] = [MUL[scale][x] for x in rows[used]]
        for i in range(len(rows)):
            if i != used and rows[i][column]:
                scale = rows[i][column]
                rows[i] = [x ^ MUL[scale][y] for x, y in zip(rows[i], rows[used])]
        used += 1
        if used == len(rows):
            break
    return used


def inverse(matrix):
    n = len(matrix)
    require(n and all(len(row) == n for row in matrix), "inverse requires a square matrix")
    rows = [list(row) + [int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for column in range(n):
        pivot = next((i for i in range(column, n) if rows[i][column]), None)
        require(pivot is not None, ("singular inverse fixture", matrix))
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = INV[rows[column][column]]
        rows[column] = [MUL[scale][x] for x in rows[column]]
        for i in range(n):
            if i != column and rows[i][column]:
                scale = rows[i][column]
                rows[i] = [x ^ MUL[scale][y] for x, y in zip(rows[i], rows[column])]
    return tuple(tuple(row[n:]) for row in rows)


def check_field():
    require({power(ALPHA, i) for i in range(Q - 1)} == set(range(1, Q)),
            "alpha must have multiplicative order 15")
    require(all(MUL[a][INV[a]] == 1 for a in range(1, Q)), "invalid field inverses")
    require(len(set(POINTS)) == 4 and all(power(x, 4) == x for x in POINTS),
            "evaluation points must form GF(4)")
    require(rank(C) == rank(P) == 2, "persistent update matrices must be invertible")


def reset_columns(batches):
    return tuple((1, MUL[power(ALPHA, batches - 1 - t)][beta])
                 for t in range(batches) for beta in BETA)


def check_reset_columns():
    cases = 0
    for batches in range(1, 5):
        columns = reset_columns(batches)
        require(len({column[1] for column in columns}) == 2 * batches,
                ("scaled points repeat inside the toy range", batches))
        for size in range(len(columns) + 1):
            for chosen in combinations(columns, size):
                require(rank(tuple(zip(*chosen))) == min(2, size),
                        ("reset-column rank", batches, chosen))
                cases += 1
    # This is a guard on the time-range hypothesis, not a code-distance claim.
    repeated = reset_columns(6)
    require(len(set(repeated)) < len(repeated), "expected quotient-order wraparound")
    return cases


def spectral_matrices(novel=False):
    vandermonde = tuple(tuple(power(x, degree) for degree in range(4)) for x in POINTS)
    if novel:
        # The columns describe a nontrivial monic degree-graded basis:
        # 1, X+1, X^2+X, X^3+alpha*X^2+X+1.
        basis = ((1, 1, 0, 1), (0, 1, 1, 1), (0, 0, 1, ALPHA), (0, 0, 0, 1))
        evaluation = matmul(vandermonde, basis)
    else:
        evaluation = vandermonde
    interpolation = inverse(evaluation)
    identity = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    require(matmul(interpolation, evaluation) == identity, "interpolation inverse")
    return tuple(row[:2] for row in evaluation), interpolation[2:]


def feedback_matrix(g, h, keys):
    return matmul(h, tuple(tuple(MUL[key][x] for x in row) for key, row in zip(keys, g)))


def check_spectral_uniformity():
    states = columns = 0
    for novel in (False, True):
        g, h = spectral_matrices(novel)
        for chosen in combinations(range(4), 2):
            require(rank(tuple(tuple(row[j] for j in chosen) for row in h)) == 2,
                    ("highest-coefficient columns", novel, chosen))
            columns += 1
        for state in product(range(Q), repeat=2):
            values = matvec(g, state)
            map_from_keys = tuple(tuple(MUL[row[j]][values[j]] for j in range(4)) for row in h)
            if not any(state):
                require(rank(map_from_keys) == 0, "zero state must have zero feedback")
                continue
            require(sum(x != 0 for x in values) >= 3, ("degree-one zero count", novel, state))
            require(rank(map_from_keys) == 2, ("joint feedback rank", novel, state))
            # A surjective linear map GF16^4 -> GF16^2 has 16^2 preimages
            # per output. Uniform independent keys therefore give joint uniformity.
            require(Q ** (4 - rank(map_from_keys)) == Q ** 2, "feedback fiber size")
            states += 1
    return states, columns


def encode(inputs, matrices):
    require(len(inputs) == len(matrices), "encoder setup length mismatch")
    state = (0, 0)
    outputs, states = [], [state]
    for batch, matrix in zip(inputs, matrices):
        output = add(batch, matvec(matrix, state))
        state = add(matvec(P, state), matvec(C, output))
        outputs.append(output)
        states.append(state)
    return tuple(outputs), tuple(states)


def decode(outputs, matrices):
    require(len(outputs) == len(matrices), "inverse setup length mismatch")
    state = (0, 0)
    inputs, states = [], [state]
    for output, matrix in zip(outputs, matrices):
        inputs.append(add(output, matvec(matrix, state)))
        state = add(matvec(P, state), matvec(C, output))
        states.append(state)
    return tuple(inputs), tuple(states)


def check_recurrence():
    rng = Random(0x504D4453)
    g, h = spectral_matrices(True)
    cases = prefixes = 0
    for batches in range(1, 5):
        for setup in range(4):
            keys = [tuple((0 if setup == 0 else 1 if setup == 1 else rng.randrange(Q))
                          for _ in range(4)) for _ in range(batches)]
            matrices = tuple(feedback_matrix(g, h, key) for key in keys)
            fixtures = [((0, 0),) * batches]
            for coordinate in range(2 * batches):
                flat = [0] * (2 * batches)
                flat[coordinate] = 1
                fixtures.append(tuple(tuple(flat[2 * t:2 * t + 2]) for t in range(batches)))
            fixtures.extend(tuple(tuple(rng.randrange(Q) for _ in range(2))
                                  for _ in range(batches)) for _ in range(16))
            for inputs in fixtures:
                outputs, states = encode(inputs, matrices)
                recovered, inverse_states = decode(outputs, matrices)
                require(recovered == inputs and inverse_states == states,
                        ("recurrence inverse", batches, setup, inputs))
                for end in range(1, batches + 1):
                    columns = reset_columns(end)
                    direct = tuple(dot(row, sum(outputs[:end], ())) for row in zip(*columns))
                    require(direct == states[end], ("unrolled reset map", end, outputs))
                    prefixes += 1
                cases += 1
            left, right = fixtures[-2:]
            summed = tuple(add(a, b) for a, b in zip(left, right))
            left_output, _ = encode(left, matrices)
            right_output, _ = encode(right, matrices)
            sum_output, _ = encode(summed, matrices)
            require(sum_output == tuple(add(a, b) for a, b in zip(left_output, right_output)),
                    "fixed-setup encoder linearity")
            # A final activation is emitted without a state flush.
            last_only = ((0, 0),) * (batches - 1) + ((1, 0),)
            outputs, states = encode(last_only, matrices)
            require(outputs == last_only and states[-1] != (0, 0), "free terminal convention")
    return cases, prefixes


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def interval_representations(n):
    """Enumerate intervals directly; adjacent intervals remain distinct."""
    counts = Counter()
    adjacent = 0

    def visit(cursor, total, prefix):
        nonlocal adjacent
        for start in range(cursor, n):
            for end in range(start + 1, n + 1):
                intervals = prefix + ((start, end),)
                counts[total + end - start, len(intervals)] += 1
                adjacent += int(any(a[1] == b[0] for a, b in zip(intervals, intervals[1:])))
                # This interval may precede another one only if its length is >=2.
                if end - start >= 2:
                    visit(end, total + end - start, intervals)

    visit(0, 0, ())
    return counts, adjacent


def check_interval_count():
    cases = representations = adjacent = 0
    for n in range(1, 8):
        counts, adjacent_here = interval_representations(n)
        adjacent += adjacent_here
        representations += sum(counts.values())
        for total in range(1, n + 1):
            for count in range(1, n + 1):
                expected = choose(total - count, count - 1) * choose(n - total + count, count)
                require(counts[total, count] == expected,
                        ("interval representation count", n, total, count,
                         counts[total, count], expected))
                cases += 1
    require(adjacent > 0, "adjacent intervals were not exercised")
    return cases, representations, adjacent


def check_quiet_kernel_trap():
    nonzero = tuple(state for state in product(range(Q), repeat=2) if any(state))
    for state in nonzero:
        require(any(matvec(P, state)), "zero output cannot kill a live state")
        require(any(matvec(C, state)), "square C cannot have a quiet activation")
    # Negative control: s=1 < b=2 gives C=(1,1). The nonzero input
    # (1,1) starts and ends at zero state, regardless of the feedback map.
    quiet = 0
    for coefficients in product(range(Q), repeat=2):
        entering = 0
        output = tuple(1 ^ MUL[coefficient][entering] for coefficient in coefficients)
        next_state = entering ^ output[0] ^ output[1]
        require(output == (1, 1) and next_state == 0, "missing s<b quiet-kernel trap")
        quiet += 1
    return len(nonzero), quiet


def main():
    check_field()
    reset = check_reset_columns()
    spectral_states, spectral_columns = check_spectral_uniformity()
    recurrence, prefixes = check_recurrence()
    interval_cases, representations, adjacent = check_interval_count()
    live_states, quiet_traps = check_quiet_kernel_trap()
    print("PASS persistent MDS toy interfaces:", dict(
        reset_subsets=reset, spectral_nonzero_states=spectral_states,
        spectral_column_pairs=spectral_columns, recurrence_cases=recurrence,
        unrolled_prefixes=prefixes, interval_count_cases=interval_cases,
        interval_representations=representations, adjacent_representations=adjacent,
        live_states=live_states, quiet_kernel_controls=quiet_traps))


if __name__ == "__main__":
    main()
