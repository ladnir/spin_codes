"""Exact four-state dual-MDS phase map and its extended adjoint.

This scalar reference preserves the archived period64 construction.  It
imports only its field, geometry, translation, and boundary utilities.
The new local map is Y=X+A0*J*S, S'=S+C0*X, where J reverses four state
coordinates.  Its adjoint is NOT the same local map.

For a label a, the global injection column is B*v(a), and the emission
row is v(a)^T*J*B^-1.  Any four distinct labels have nonzero Vandermonde
determinant, the product of their six pairwise differences.  Invertible
left/right factors preserve rank, proving the period64 MDS property.
The tests check the hypotheses and independent representative minors;
they do not substitute a sample of minors for that combinatorial proof.

No native implementation, benchmark, or whole-code certificate is supplied.
"""
from __future__ import annotations

import argparse
from functools import lru_cache

from rate_design_coset import mul, power
from region64_reference import (
    IDENTITY, base_maps as ordinary_base_maps, entry_matrix, final_matrix,
    inverse, matmul, matvec, phase_geometry, scale, translation, transpose,
    validate,
)


J = tuple(tuple(int(i+j == 3) for j in range(4)) for i in range(4))
ZERO_MATRIX = ((0, 0, 0, 0),)*4


def curve(label):
    square = mul(label, label)
    return 1, label, square, mul(square, label)


def period_labels():
    points, phases = phase_geometry()
    return tuple(beta ^ u for beta in phases for u in points)


@lru_cache(maxsize=1)
def base_maps():
    """Return new emission A0*J and unchanged injection C0."""
    ordinary_emission, injection = ordinary_base_maps()
    return matmul(ordinary_emission, J), injection


def local_step(values, state):
    """Forward step: reversal is on the emission input."""
    ordinary_emission, injection = ordinary_base_maps()
    output = tuple(x ^ y for x, y in zip(values, matvec(ordinary_emission, tuple(reversed(state)))))
    following = tuple(q ^ c for q, c in zip(state, matvec(injection, values)))
    return output, following


def adjoint_local_step(values, state):
    """True transpose: reversal is on the moment output, not emission."""
    ordinary_emission, injection = ordinary_base_maps()
    output = tuple(z ^ c for z, c in zip(values, matvec(ordinary_emission, state)))
    moment = tuple(reversed(matvec(injection, values)))
    previous = tuple(h ^ a for h, a in zip(state, moment))
    return output, previous


def global_maps(orientation, beta):
    """Build maps directly from translated labels, independently of T_beta."""
    points = phase_geometry()[0]
    rows = tuple(curve(u ^ beta) for u in points)
    emission = matmul(matmul(rows, J), inverse(orientation))
    injection = matmul(orientation, transpose(rows))
    return emission, injection


def period_maps(orientation=IDENTITY):
    rows = tuple(curve(label) for label in period_labels())
    return matmul(matmul(rows, J), inverse(orientation)), matmul(orientation, transpose(rows))


def global_step(values, state, orientation, beta, adjoint=False):
    emission, injection = global_maps(orientation, beta)
    if adjoint:
        output = tuple(z ^ c for z, c in zip(values, matvec(transpose(injection), state)))
        previous = tuple(h ^ a for h, a in zip(state, matvec(transpose(emission), values)))
        return output, previous
    output = tuple(x ^ a for x, a in zip(values, matvec(emission, state)))
    following = tuple(q ^ c for q, c in zip(state, matvec(injection, values)))
    return output, following


def _state(values):
    result = tuple(values)
    if len(result) != 4 or any(not isinstance(value, int) or not 0 <= value < 256 for value in result):
        raise ValueError("state must contain four GF256 bytes")
    return result


def encode_global(regions, orientations, initial=(0, 0, 0, 0), diagonals=None):
    validate(regions, orientations, diagonals)
    phases = phase_geometry()[1]
    state, result = _state(initial), []
    for region, (blocks, orientation) in enumerate(zip(regions, orientations)):
        outputs = []
        for step, values in enumerate(blocks):
            output, state = global_step(scale(values, diagonals, region, step), state,
                                        orientation, phases[step % 4])
            outputs.append(output)
        result.append(tuple(outputs))
    return tuple(result), state


def encode_local(regions, orientations, initial=(0, 0, 0, 0), diagonals=None):
    validate(regions, orientations, diagonals)
    phases = phase_geometry()[1]
    state, result = _state(initial), []
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
    state, result = _state(terminal), [None]*len(regions)
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
    state = matvec(transpose(final_matrix(regions, orientations)), _state(terminal))
    result = [None]*len(regions)
    for region in reversed(range(len(regions))):
        outputs = [None]*len(regions[region])
        for step in reversed(range(len(regions[region]))):
            values, state = adjoint_local_step(regions[region][step], state)
            outputs[step] = scale(values, diagonals, region, step)
            if step:
                state = translation(state, phases[(step-1) % 4] ^ phases[step % 4], adjoint=True)
        state = matvec(transpose(entry_matrix(regions, orientations, region)), state)
        result[region] = tuple(outputs)
    return tuple(result), state


def vandermonde_determinant(labels):
    """Determinant of the square leading-power Vandermonde minor."""
    result = 1
    for index, label in enumerate(labels):
        for previous in labels[:index]:
            result = mul(result, label ^ previous)
    return result


def determinant(matrix):
    """Independent elimination oracle for the small test minors."""
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("square matrix required")
    work = [list(row) for row in matrix]
    result = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return 0
        work[column], work[pivot] = work[pivot], work[column]
        diagonal = work[column][column]
        result = mul(result, diagonal)
        reciprocal = power(diagonal, 254)
        work[column] = [mul(reciprocal, value) for value in work[column]]
        for row in range(column+1, size):
            factor = work[row][column]
            work[row] = [x ^ mul(factor, y) for x, y in zip(work[row], work[column])]
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        import unittest
        suite = unittest.defaultTestLoader.loadTestsFromName("test_dual_mds_phase")
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        raise SystemExit(0 if result.wasSuccessful() else 1)
    print("Dual MDS period64: forward emission reverses state; adjoint reverses moment output.")
