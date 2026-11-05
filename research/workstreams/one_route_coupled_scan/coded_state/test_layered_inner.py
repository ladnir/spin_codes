"""Bounded exact checks for the actual pair-shear family, not a full certificate.

No files are generated. The tiny-field cycle census is exhaustive; the GF16
affine-fiber checks exercise all incoming pair values and all line directions.
"""
from collections import Counter
from itertools import product
from fractions import Fraction
from math import comb
import time


def field(bits, polynomial):
    q = 1 << bits
    table = [[0]*q for _ in range(q)]
    for a, b in product(range(q), repeat=2):
        x, y, z = a, b, 0
        while y:
            if y & 1:
                z ^= x
            y >>= 1
            x <<= 1
            if x & q:
                x ^= polynomial
        table[a][b] = z
    return table


def setup(table, delta=2):
    q = len(table)
    inverse = next(x for x in range(1, q) if table[delta][x] == 1)
    excluded = {0, 1, 1 ^ delta, delta ^ inverse}
    keys = [(a, b) for a in range(1, q) for b in range(1, q)
            if table[a][b] not in excluded]
    return keys, q-len(excluded)


def pair(table, key, x, y):
    a, b = key
    p = x ^ table[a][y]
    return p, y ^ table[b][p]


def cycle_census():
    table = field(2, 7)
    keys, allowed = setup(table)
    assert len(keys) == 3 and allowed == 1
    choices = list(product(keys, repeat=4))
    maximum = 0
    worst = None
    for state in product(range(4), repeat=4):
        if not any(state):
            continue
        counts = Counter()
        for chosen in choices:
            out = list(state)
            for k, (i, j) in enumerate(((0, 1), (2, 3), (0, 2), (1, 3))):
                out[i], out[j] = pair(table, chosen[k], out[i], out[j])
            assert any(out)
            counts[tuple(out)] += 1
        if max(counts.values()) > maximum:
            maximum = max(counts.values())
            worst = state
        if state == (1, 0, 0, 0):
            assert len(counts) == 27 and set(counts.values()) == {3}
            assert all(out[0] == 1 for out in counts)
    return dict(field=4, states=255, setups_per_state=81,
                max_point_numerator=maximum, denominator=81, witness=worst)


def affine_fibers():
    table = field(4, 19)
    keys, allowed = setup(table)
    q = len(table)
    # All one-dimensional subspaces, then {0} and the whole plane.
    codes = [{(a, table[slope][a]) for a in range(q)} for slope in range(q)]
    codes += [{(0, a) for a in range(q)}, {(0, 0)}, set(product(range(q), repeat=2))]
    checks = 0
    for x, y in product(range(q), repeat=2):
        distribution = Counter()
        for key in keys:
            p, r = pair(table, key, x, y)
            distribution[(p ^ table[2][x], r ^ table[2][y])] += 1
        for code in codes:
            for offset in ((0, 0), (1, 3)):
                affine = {(a ^ offset[0], b ^ offset[1]) for a, b in code}
                actual = Fraction(sum(distribution[z] for z in affine), len(keys))
                if not (x or y):
                    bound = Fraction((0, 0) in affine)
                elif y == 0:
                    compatible = [z for z in affine if z[0] == table[3][x]]
                    shortened = sum(a == 0 for a, b in code)
                    assert not compatible or len(compatible) == shortened
                    bound = Fraction(shortened if compatible else 0, q-1)
                else:
                    # Regular branch has two free coordinates. The exceptional
                    # branch has none, and is possible only when x != 0.
                    bound = Fraction(len(code), (q-1)*allowed)
                    exceptional = (table[2][x], table[3][y])
                    if x and exceptional in affine:
                        bound += Fraction(1, q-1)
                assert actual <= bound
                checks += 1
    return checks


def sparse_paths():
    table = field(4, 19)
    keys, _ = setup(table)
    # A three-symbol input can close a one-coordinate state in two batches.
    u, v = 1, 5
    reset_keys = sum(pair(table, key, u, 0) == (u, v) for key in keys)
    assert Fraction(reset_keys, len(keys)) == Fraction(1, 15)
    # Sustained one-coordinate states exist for every setup if arbitrary input
    # is allowed. This is not an outer-codeword witness.
    state = [0]*512
    state[0] = 1
    total_x = total_y = 0
    for t in range(27):
        mapped = list(state)
        stage = t % 9
        for i in range(512):
            if i & (1 << stage):
                continue
            j = i ^ (1 << stage)
            mapped[i], mapped[j] = pair(table, keys[(i+t) % len(keys)], state[i], state[j])
        # Specify Q_next = Q; compute the input that maintains this state.
        x = [a ^ b for a, b in zip(state, mapped)]
        y = [a ^ table[2][a] for a in state]
        assert sum(a != 0 for a in x) == 1
        assert sum(a != 0 for a in y) == 1
        total_x += sum(a != 0 for a in x)
        total_y += sum(a != 0 for a in y)
    return dict(reset_probability=str(Fraction(reset_keys, len(keys))),
                sustained_batches=27, sustained_input_weight=total_x,
                sustained_output_weight=total_y)


def dense_receipt():
    N, W, zero_budget, bits = 2**21, (45*2**21)//100, 8192, 32
    # Union over <=2^N supports, <=q^W amplitudes, <=2^(N/2)
    # exceptional-branch masks; point-density inflation is strictly <2.
    log2_upper = N + bits*W + N//2 - (bits//2)*(N-zero_budget) + 1
    assert log2_upper == -78655
    assert 4*N < (1 << bits)//2  # Bernoulli makes the inflation <2.
    return dict(distance_permille=450, state_zero_budget=zero_budget,
                log2_failure_strict_upper=log2_upper, whole_code_certificate=False)


def fixed_coordinate_receipt():
    # Padding includes every output whose spatial support has size <=23.
    g = 23
    numerator = comb(512, g)*8192*comb(256, 129)*g**129
    denominator = 256**129
    assert numerator*(1 << 51) < denominator
    # Negative control: the same coarse bound does not certify g=24 at 40 bits.
    assert comb(512, 24)*8192*comb(256, 129)*24**129*(1 << 40) > denominator
    return dict(coordinates_at_most=g, failure_below_power_of_two=-51,
                arbitrary_inner_keys=True, payload_field_bits=128)


def main():
    start = time.monotonic()
    print(dict(cycle_census=cycle_census(), affine_fiber_checks=affine_fibers(),
               sparse_paths=sparse_paths(), dense_receipt=dense_receipt(),
               fixed_coordinate_receipt=fixed_coordinate_receipt()))
    assert time.monotonic()-start < 30


if __name__ == '__main__':
    main()
