"""Separate width32/state8 symmetric Cauchy scalar reference.

The two label sets combine the four old width16 phases in pairs. No old
construction is changed; no native kernel, timing, or distance certificate.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from random import Random

from state_capacity_reference import (BETA, ETA, NODES, LABELS as LABELS16,
                                     coefficients, circulant, emit_coefficient,
                                     linearized, matvec, mul, power, xor_sum)


LABELS = (LABELS16[0]+LABELS16[1], LABELS16[2]+LABELS16[3])


@lru_cache(maxsize=2)
def matrices(phase=0):
    a = tuple(tuple(emit_coefficient(u ^ v) for v in NODES) for u in LABELS[phase])
    return a, tuple(zip(*a))


def emission(state, phase=0):
    return sum((circulant(state, coefficients(LABELS[phase][first] ^ BETA, True))
                for first in (0, 8, 16, 24)), ())


def injection(values, phase=0):
    blocks = [circulant(values[first:first+8],
                        coefficients(LABELS[phase][first] ^ BETA, True))
              for first in (0, 8, 16, 24)]
    return tuple(xor_sum(block[j] for block in blocks) for j in range(8))


def step(values, state, diagonal, phase=0):
    refreshed = tuple(mul(d, x) for d, x in zip(diagonal, state))
    output = tuple(x ^ mul(ETA, y) for x, y in zip(values, emission(refreshed, phase)))
    following = tuple(x ^ y for x, y in zip(refreshed, injection(values, phase)))
    return output, following


def adjoint(outputs, terminal, diagonal, phase=0):
    values = tuple(x ^ y for x, y in zip(outputs, emission(terminal, phase)))
    prior = tuple(mul(d, h ^ mul(ETA, y)) for d, h, y in
                  zip(diagonal, terminal, injection(outputs, phase)))
    return values, prior


def return_minimum_weight(weight):
    if not 1 <= weight <= 32:
        raise ValueError("nonempty packet weight must lie in1..32")
    return 25 if weight == 1 else max(1, 25-weight)


def check():
    counts = dict(map_basis=0, singleton_returns=0, pairings=0,
                  extended_adjoints=0, power_sums=0, stream_involutions=0)
    for phase, labels in enumerate(LABELS):
        assert len(labels) == len(set(labels)) == 32
        assert set(labels).isdisjoint(NODES)
        origin = labels[0]
        differences = {u ^ origin for u in labels}
        assert all(a ^ b in differences for a in differences for b in differences)
        for exponent in range(31):
            assert xor_sum(power(u, exponent) for u in labels) == 0
            counts['power_sums'] += 1
        a, c = matrices(phase)
        for j in range(8):
            state = tuple(int(i == j) for i in range(8))
            out = emission(state, phase)
            assert out == matvec(a, state)
            assert injection(out, phase) == (0,)*8
            counts['map_basis'] += 1
        for j in range(32):
            values = tuple(int(i == j) for i in range(32))
            state = injection(values, phase)
            assert state == matvec(c, values)
            emitted = emission(state, phase)
            returned = tuple(x ^ mul(ETA, y) for x, y in zip(values, emitted))
            assert sum(bool(x) for x in returned) == 25
            assert sum(x == 0 for i, x in enumerate(emitted) if i != j) == 7
            counts['map_basis'] += 1
            counts['singleton_returns'] += 1
            for i, u in enumerate(labels):
                expected = 1 if i == j else mul(linearized(u ^ labels[j]), power(u ^ labels[j], 254))
                assert emitted[i] == expected
                counts['pairings'] += 1
    rng = Random(832256176)
    dot = lambda x, y: xor_sum(mul(a, b) for a, b in zip(x, y))
    for index in range(32):
        phase = index % 2
        values = tuple(rng.randrange(256) for _ in range(32))
        state = tuple(rng.randrange(256) for _ in range(8))
        diagonal = tuple(rng.randrange(1, 256) for _ in range(8))
        outputs = tuple(rng.randrange(256) for _ in range(32))
        terminal = tuple(rng.randrange(256) for _ in range(8))
        y, following = step(values, state, diagonal, phase)
        dx, ds = adjoint(outputs, terminal, diagonal, phase)
        assert dot(y, outputs) ^ dot(following, terminal) == dot(values, dx) ^ dot(state, ds)
        restored, repeated = step(y, state, diagonal, phase)
        assert restored == values and repeated == following
        counts['extended_adjoints'] += 1
    for length in (1, 2, 3, 7):
        packets = [tuple(rng.randrange(256) for _ in range(32)) for _ in range(length)]
        diagonals = [tuple(rng.randrange(1, 256) for _ in range(8)) for _ in range(length)]
        def encode(source):
            state, output = (0,)*8, []
            for j, values in enumerate(source):
                y, state = step(values, state, diagonals[j], j % 2)
                output.append(y)
            return output, state
        output, final = encode(packets)
        restored, repeated = encode(output)
        assert restored == packets and repeated == final
        counts['stream_involutions'] += 1
    return counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(check())
