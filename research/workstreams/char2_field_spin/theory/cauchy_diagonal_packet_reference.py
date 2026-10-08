"""Five-product GF4-circulant maps for a separate diagonal-refresh packet.

No native encoder, benchmark, or whole-code certificate.  The field is AES
GF256, with actual embedded GF4/GF16.  Both maps use the same fixed tables.
See CAUCHY_DIAGONAL_PACKET_COST.md for the construction and operation count.
"""
from __future__ import annotations

import argparse
from random import Random

from rate_design_coset import mul, power, xor_sum
from region64_reference import phase_geometry


def geometry():
    h16, phases = phase_geometry()
    h4 = tuple(x for x in h16 if power(x, 4) == x)
    omega = next(x for x in h4 if x not in (0, 1))
    h4 = (0, 1, omega, omega ^ 1)
    alpha = next(x for x in h16 if x not in h4)
    representatives = (0, alpha, mul(alpha, omega), mul(alpha, omega ^ 1))
    labels = tuple(tuple(phase ^ representative ^ x
                         for representative in representatives for x in h4)
                   for phase in phases)
    union = set(sum(labels, ()))
    beta = next(x for x in range(256) if x not in union)
    nodes = tuple(beta ^ x for x in h4)
    assert len(union) == 64 and len(set(nodes)) == 4 and union.isdisjoint(nodes)
    return h4, representatives, labels, beta, nodes


H4, REPRESENTATIVES, LABELS, BETA, NODES = geometry()


def add(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def zeta4(values):
    # The superset zeta transform is its own inverse in characteristic two.
    x0, x1, x2, x3 = values
    a, b = x0 ^ x1, x2 ^ x3
    return a ^ b, x1 ^ x3, b, x3


def circulant_coefficients(offset):
    return tuple(power(offset ^ x, 3) for x in H4)


def circulant4(coefficients, values):
    """XOR convolution by four fixed coefficients; exactly five products.

    The special cubic coefficients have zeta[0]=1.  The five multiplication
    sites below are scalar GF256 products.  Four independent cosets can fill
    the four 128-bit-record lanes of each corresponding GFNI512 instruction.
    """
    one, g1, g2, g3 = zeta4(coefficients)
    assert one == 1
    z0, z1, z2, z3 = zeta4(values)
    w0 = z0
    w1 = z1 ^ mul(g1, z0)
    w2 = z2 ^ mul(g2, z0)
    w3 = z3 ^ mul(g1, z2) ^ mul(g2, z1) ^ mul(g3, z0)
    return zeta4((w0, w1, w2, w3))


def dense_circulant(coefficients, values):
    return tuple(xor_sum(mul(coefficients[i ^ j], values[j]) for j in range(4))
                 for i in range(4))


def emission(state, phase):
    """A[u,j]=(u+v_j)^3+1, before the separate eta emission scale."""
    total = xor_sum(state)
    result = []
    for first in range(0, 16, 4):
        offset = LABELS[phase][first] ^ BETA
        block = circulant4(circulant_coefficients(offset), state)
        result.extend(x ^ total for x in block)
    return tuple(result)


def injection(values, phase):
    """C[j,u]=(u+v_j)^3; accumulate four symmetric circulant blocks."""
    result = (0, 0, 0, 0)
    for first in range(0, 16, 4):
        offset = LABELS[phase][first] ^ BETA
        result = add(result, circulant4(circulant_coefficients(offset), values[first:first+4]))
    return result


def direct_emission(state, phase):
    return tuple(xor_sum(mul(power(u ^ v, 3) ^ 1, s) for v, s in zip(NODES, state))
                 for u in LABELS[phase])


def direct_injection(values, phase):
    return tuple(xor_sum(mul(power(u ^ v, 3), x) for u, x in zip(LABELS[phase], values))
                 for v in NODES)


def step(values, state, diagonal, phase, eta=2):
    refreshed = tuple(mul(d, s) for d, s in zip(diagonal, state))
    output = add(values, tuple(mul(eta, value) for value in emission(refreshed, phase)))
    following = add(refreshed, injection(values, phase))
    return output, following


def adjoint_step(outputs, terminal, diagonal, phase, eta=2):
    # C^T h differs from A h by sum(h) in each output coordinate.
    state_sum = xor_sum(terminal)
    input_dual = add(outputs, tuple(value ^ state_sum for value in emission(terminal, phase)))
    # A^T z differs from C z by sum(z) in each state coordinate.
    output_sum = xor_sum(outputs)
    moment = tuple(value ^ output_sum for value in injection(outputs, phase))
    prior_dual = tuple(mul(d, h ^ mul(eta, value))
                       for d, h, value in zip(diagonal, terminal, moment))
    return input_dual, prior_dual


def inner_product(a, b):
    return xor_sum(mul(x, y) for x, y in zip(a, b))


def encode_inner(packets, state_diagonals):
    """Zero-initial-state inner, without outer/input coordinate diagonals."""
    assert len(packets) == len(state_diagonals)
    state, outputs = (0, 0, 0, 0), []
    for index, (values, diagonal) in enumerate(zip(packets, state_diagonals)):
        output, state = step(values, state, diagonal, index % 4)
        outputs.append(output)
    return tuple(outputs), state


def check():
    cases = 0
    for phase in range(4):
        for j in range(4):
            for scalar in range(256):
                state = tuple(scalar if i == j else 0 for i in range(4))
                value = emission(state, phase)
                assert value == direct_emission(state, phase)
                assert injection(value, phase) == (0, 0, 0, 0)
                cases += 1
        for j in range(16):
            for scalar in range(256):
                values = tuple(scalar if i == j else 0 for i in range(16))
                assert injection(values, phase) == direct_injection(values, phase)
                cases += 1
        for u in LABELS[phase]:
            denominator = power(u, 4) ^ u ^ power(BETA, 4) ^ BETA
            assert denominator != 0
            for v in NODES:
                # Row-scaled Cauchy identity, including every actual phase.
                assert mul(power(u ^ v, 3) ^ 1, u ^ v) == denominator
            for w in LABELS[phase]:
                # Exact pairing needed for the all-J similarity argument.
                pairing = xor_sum(mul(power(u ^ v, 3) ^ 1, power(w ^ v, 3)) for v in NODES)
                assert pairing == power(u ^ w, 3)
    rng = Random(20261014)
    for case in range(128):
        phase = case % 4
        values = tuple(rng.randrange(256) for _ in range(16))
        state = tuple(rng.randrange(256) for _ in range(4))
        outputs = tuple(rng.randrange(256) for _ in range(16))
        terminal = tuple(rng.randrange(256) for _ in range(4))
        diagonal = tuple(rng.randrange(1, 256) for _ in range(4))
        output, following = step(values, state, diagonal, phase)
        input_dual, prior_dual = adjoint_step(outputs, terminal, diagonal, phase)
        assert inner_product(output, outputs) ^ inner_product(following, terminal) == (
            inner_product(values, input_dual) ^ inner_product(state, prior_dual))
        # Apply the same local inner a second time, with the same incoming
        # state: CA=0 makes both state trajectories equal.
        restored, repeated_following = step(output, state, diagonal, phase)
        assert restored == values and repeated_following == following
    for length in (1, 3, 4, 5, 8, 17):
        packets = tuple(tuple(rng.randrange(256) for _ in range(16)) for _ in range(length))
        diagonals = tuple(tuple(rng.randrange(1, 256) for _ in range(4)) for _ in range(length))
        output, final = encode_inner(packets, diagonals)
        recovered, repeated_final = encode_inner(output, diagonals)
        assert recovered == packets and repeated_final == final
    return cases


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        print(f"passed {check()} exact map basis/scalar checks, 128 extended adjoints, pairing identities, and six stream involutions")
    print(f"AES embedded GF4={H4}; residue nodes={NODES}; beta={BETA}")
    print("Each map: four 5-product circulants =20 scalar products =5 GFNI512 record lanes.")
    print("Conservative packet:10 maps +1 eta +1 state diagonal +4 input diagonal =16 GFNI512 equivalents.")
    print("Complete [32,22] arithmetic volume:41/88 +16/11 =169/88 per source record; not a runtime claim.")
