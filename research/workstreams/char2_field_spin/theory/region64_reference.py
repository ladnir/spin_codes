"""Exact four-state, period-64 regional inner and extended adjoint.

The width-16 base maps use the actual GF16 subfield of AES GF256.  Four
additive cosets provide 64 distinct injection labels.  A fixed invertible
orientation is supplied for each region; no matrix is sampled per step.
This is a scalar algebra reference, not a kernel, benchmark, or complete
distance certificate.  External initial/final states are retained so tests
can detect errors at the first and last basis conversion.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from random import Random

from rate_design_coset import geometry, mul, power, xor_sum


IDENTITY = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))


def transpose(matrix):
    return tuple(zip(*matrix))


def matvec(matrix, vector):
    return tuple(xor_sum(mul(a, b) for a, b in zip(row, vector)) for row in matrix)


def matmul(left, right):
    columns = transpose(right)
    return tuple(tuple(xor_sum(mul(a, b) for a, b in zip(row, column))
                       for column in columns) for row in left)


def inverse(matrix):
    work = [list(row)+list(unit) for row, unit in zip(matrix, IDENTITY)]
    for column in range(4):
        pivot = next((row for row in range(column, 4) if work[row][column]), None)
        if pivot is None:
            raise ValueError("singular regional orientation")
        work[column], work[pivot] = work[pivot], work[column]
        reciprocal = power(work[column][column], 254)
        work[column] = [mul(reciprocal, value) for value in work[column]]
        for row in range(4):
            if row != column:
                coefficient = work[row][column]
                work[row] = [a ^ mul(coefficient, b) for a, b in zip(work[row], work[column])]
    return tuple(tuple(row[4:]) for row in work)


@lru_cache(maxsize=1)
def phase_geometry():
    subfield, representatives = geometry()
    # Match Subfield16.h's abstract GF16 polynomial basis, rather than
    # interpreting abstract nibble labels as AES-field bytes.
    root = next(a for a in subfield if a > 1 and power(a, 4) ^ a ^ 1 == 0)
    basis = tuple(power(root, j) for j in range(4))
    points = tuple(xor_sum(basis[j] for j in range(4) if (label >> j) & 1)
                   for label in range(16))
    assert set(points) == set(subfield)
    a = representatives[1]
    first_two = set(points) | {a ^ u for u in points}
    b = next(beta for beta in representatives if beta not in first_two)
    phases = (0, a, a ^ b, b)
    assert len({beta ^ u for beta in phases for u in points}) == 64
    return points, phases


@lru_cache(maxsize=256)
def translation_coefficients(beta):
    beta2 = mul(beta, beta)
    return beta2, mul(beta2, beta)


def translation_matrix(beta):
    beta2, beta3 = translation_coefficients(beta)
    return ((1, 0, 0, 0), (beta, 1, 0, 0),
            (beta2, 0, 1, 0), (beta3, beta2, beta, 1))


def translation(state, beta, adjoint=False):
    """Four products/four XORs by two shears; coefficients are prepared once.

    The direct matrix has five products.  The factored update intentionally
    uses its new q1 term to generate beta^3*q0 without a fifth product.
    Transposing the two factors gives the adjoint update below.
    """
    q0, q1, q2, q3 = state
    beta2, _ = translation_coefficients(beta)
    if adjoint:
        a = q0 ^ mul(beta2, q2)
        b = q1 ^ mul(beta2, q3)
        c = q2 ^ mul(beta, q3)
        return (a ^ mul(beta, b), b, c, q3)
    b = q1 ^ mul(beta, q0)
    d = q3 ^ mul(beta, q2)
    c = q2 ^ mul(beta2, q0)
    return (q0, b, c, d ^ mul(beta2, b))


@lru_cache(maxsize=1)
def base_maps():
    points, _ = phase_geometry()
    emission = tuple((1, u, mul(u, u), power(u, 3)) for u in points)
    return emission, transpose(emission)


def local_step(values, state):
    emission, injection = base_maps()
    output = tuple(x ^ y for x, y in zip(values, matvec(emission, state)))
    following = tuple(q ^ c for q, c in zip(state, matvec(injection, values)))
    return output, following


def global_maps(orientation, beta):
    emission, injection = base_maps()
    change = translation_matrix(beta)
    # T_beta^{-1}=T_beta.  Emission and injection are generally not adjoints.
    return (matmul(matmul(emission, change), inverse(orientation)),
            matmul(matmul(orientation, change), injection))


def global_step(values, state, orientation, beta, adjoint=False):
    emission, injection = global_maps(orientation, beta)
    if adjoint:
        # State is the dual of the following global state; values dualize Y.
        output = tuple(v ^ c for v, c in zip(values, matvec(transpose(injection), state)))
        previous = tuple(h ^ a for h, a in zip(state, matvec(transpose(emission), values)))
        return output, previous
    output = tuple(x ^ a for x, a in zip(values, matvec(emission, state)))
    following = tuple(q ^ c for q, c in zip(state, matvec(injection, values)))
    return output, following


def validate(regions, orientations, diagonals):
    if not regions or len(regions) != len(orientations) or any(not region for region in regions):
        raise ValueError("one orientation per nonempty region required")
    if any(len(block) != 16 for region in regions for block in region):
        raise ValueError("every step has sixteen symbols")
    if any(len(matrix) != 4 or any(len(row) != 4 for row in matrix) for matrix in orientations):
        raise ValueError("regional orientations must be four-by-four")
    for matrix in orientations:
        inverse(matrix)
    if diagonals is not None:
        if (len(diagonals) != len(regions)
                or any(len(d) != len(r) for d, r in zip(diagonals, regions))
                or any(len(block) != 16 or any(not 0 < value < 256 for value in block)
                       for region in diagonals for block in region)):
            raise ValueError("diagonals must supply one nonzero GF256 scalar per input")


def scale(values, diagonals, region, step):
    return tuple(values) if diagonals is None else tuple(
        mul(value, coefficient) for value, coefficient in zip(values, diagonals[region][step]))


def encode_global(regions, orientations, initial=(0, 0, 0, 0), diagonals=None):
    validate(regions, orientations, diagonals)
    phases = phase_geometry()[1]
    state, result = tuple(initial), []
    for region, (blocks, orientation) in enumerate(zip(regions, orientations)):
        outputs = []
        for step, values in enumerate(blocks):
            output, state = global_step(scale(values, diagonals, region, step), state,
                                        orientation, phases[step % 4])
            outputs.append(output)
        result.append(tuple(outputs))
    return tuple(result), state


def entry_matrix(regions, orientations, region):
    """Map preceding local state to the first local state of this region."""
    current_inverse = inverse(orientations[region])
    if region == 0:
        return current_inverse  # First phase is zero; initial state is global.
    last_beta = phase_geometry()[1][(len(regions[region-1])-1) % 4]
    return matmul(matmul(current_inverse, orientations[region-1]), translation_matrix(last_beta))


def final_matrix(regions, orientations):
    """Convert the last local state back to the reported global state."""
    last_beta = phase_geometry()[1][(len(regions[-1])-1) % 4]
    return matmul(orientations[-1], translation_matrix(last_beta))


def encode_local(regions, orientations, initial=(0, 0, 0, 0), diagonals=None):
    validate(regions, orientations, diagonals)
    phases = phase_geometry()[1]
    state, result = tuple(initial), []
    for region, blocks in enumerate(regions):
        state = matvec(entry_matrix(regions, orientations, region), state)
        outputs = []
        for step, values in enumerate(blocks):
            if step:
                state = translation(state, phases[(step-1) % 4] ^ phases[step % 4])
            output, state = local_step(scale(values, diagonals, region, step), state)
            outputs.append(output)
        result.append(tuple(outputs))
    return tuple(result), matvec(final_matrix(regions, orientations), state)


def adjoint_global(regions, orientations, terminal=(0, 0, 0, 0), diagonals=None):
    validate(regions, orientations, diagonals)
    phases = phase_geometry()[1]
    state, result = tuple(terminal), [None]*len(regions)
    for region in reversed(range(len(regions))):
        outputs = [None]*len(regions[region])
        for step in reversed(range(len(regions[region]))):
            values, state = global_step(regions[region][step], state, orientations[region],
                                        phases[step % 4], adjoint=True)
            outputs[step] = scale(values, diagonals, region, step)
        result[region] = tuple(outputs)
    return tuple(result), state


def adjoint_local(regions, orientations, terminal=(0, 0, 0, 0), diagonals=None):
    validate(regions, orientations, diagonals)
    phases = phase_geometry()[1]
    state = matvec(transpose(final_matrix(regions, orientations)), terminal)
    result = [None]*len(regions)
    for region in reversed(range(len(regions))):
        outputs = [None]*len(regions[region])
        for step in reversed(range(len(regions[region]))):
            # The local map is symmetric because C0=A0^T.
            values, state = local_step(regions[region][step], state)
            outputs[step] = scale(values, diagonals, region, step)
            if step:
                state = translation(state, phases[(step-1) % 4] ^ phases[step % 4], adjoint=True)
        state = matvec(transpose(entry_matrix(regions, orientations, region)), state)
        result[region] = tuple(outputs)
    return tuple(result), state


def random_orientation(rng):
    while True:
        matrix = tuple(tuple(rng.randrange(256) for _ in range(4)) for _ in range(4))
        try:
            inverse(matrix)
            return matrix
        except ValueError:
            pass


def dot(left, right):
    return xor_sum(mul(a, b) for a, b in zip(left, right))


def flatten(regions):
    return tuple(value for region in regions for block in region for value in block)


def check():
    rng = Random(20261005)
    points, phases = phase_geometry()
    emission, injection = base_maps()
    assert matmul(injection, emission) == ((0, 0, 0, 0),)*4
    assert all(power(u, 16) == u for u in points)
    assert len({beta ^ u for beta in phases for u in points}) == 64
    assert {phases[i] ^ phases[(i+1) % 4] for i in range(4)} == {phases[1], phases[3]}
    for beta in range(256):
        change = translation_matrix(beta)
        assert matmul(change, change) == IDENTITY
        for u in points:
            vector = (1, u, mul(u, u), power(u, 3))
            shifted = u ^ beta
            assert matvec(change, vector) == (1, shifted, mul(shifted, shifted), power(shifted, 3))
        for coordinate in range(4):
            state = tuple(255 if j == coordinate else 0 for j in range(4))
            assert translation(state, beta) == matvec(change, state)
            assert translation(state, beta, True) == matvec(transpose(change), state)
        gamma = rng.randrange(256)
        assert matmul(change, translation_matrix(gamma)) == translation_matrix(beta ^ gamma)

    sequences = ((1,), (2,), (3,), (4,), (5,), (7,), (1, 2), (3, 4, 5), (7, 1, 6))
    cases = 0
    for lengths in sequences:
        for with_diagonal in (False, True):
            for trial in range(4):
                regions = tuple(tuple(tuple(rng.randrange(256) for _ in range(16))
                                      for _ in range(length)) for length in lengths)
                dual = tuple(tuple(tuple(rng.randrange(256) for _ in range(16))
                                   for _ in range(length)) for length in lengths)
                matrices = tuple(IDENTITY if trial == 0 else random_orientation(rng) for _ in lengths)
                diagonal = (tuple(tuple(tuple(rng.randrange(1, 256) for _ in range(16))
                                        for _ in range(length)) for length in lengths)
                            if with_diagonal else None)
                initial = tuple(rng.randrange(256) for _ in range(4))
                terminal = tuple(rng.randrange(256) for _ in range(4))
                direct = encode_global(regions, matrices, initial, diagonal)
                local = encode_local(regions, matrices, initial, diagonal)
                assert local == direct
                transposed = adjoint_global(dual, matrices, terminal, diagonal)
                assert adjoint_local(dual, matrices, terminal, diagonal) == transposed
                output, final = direct
                dual_input, dual_initial = transposed
                assert dot(flatten(output), flatten(dual)) ^ dot(final, terminal) == (
                    dot(flatten(regions), flatten(dual_input)) ^ dot(initial, dual_initial))
                # Also check the actual zero-initial, discarded-final encoder.
                output, _ = encode_local(regions, matrices, diagonals=diagonal)
                dual_input, _ = adjoint_local(dual, matrices, diagonals=diagonal)
                assert dot(flatten(output), flatten(dual)) == dot(flatten(regions), flatten(dual_input))
                cases += 1
    return cases


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        print(f"period64: translation identities and {check()} global/local/adjoint cases pass")
    else:
        points, phases = phase_geometry()
        print(f"embedded GF16 points: {points}")
        print(f"Gray-code phases: {phases}; 64 distinct injection labels")
