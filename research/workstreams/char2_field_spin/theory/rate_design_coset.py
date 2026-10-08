"""Exact scalar reference for a GF256 coset-state inner and its adjoint.

This is a mathematical reference, not a native kernel or a benchmark.
The within-region maps use the GF16 subfield of the AES GF256 field.
Only boundary changes use GF256 coefficients outside that subfield.
"""

from functools import lru_cache


def mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 256:
            a ^= 0x11b
    return result


def power(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = mul(result, a)
        a = mul(a, a)
        exponent >>= 1
    return result


def xor_sum(values):
    result = 0
    for value in values:
        result ^= value
    return result


@lru_cache(maxsize=None)
def geometry():
    subfield = tuple(a for a in range(256) if power(a, 16) == a)
    covered, representatives = set(), []
    for beta in range(256):
        if beta not in covered:
            representatives.append(beta)
            covered.update(beta ^ u for u in subfield)
    return subfield, tuple(representatives)


def basis(state, delta, transpose=False):
    a, b, c = state
    d2 = mul(delta, delta)
    if transpose:
        return [a ^ mul(delta, b) ^ mul(d2, c), b, c]
    return [a, b ^ mul(delta, a), c ^ mul(d2, a)]


def emit(state):
    return [state[0] ^ mul(u, state[1]) ^ mul(mul(u, u), state[2]) for u in geometry()[0]]


def moments(values):
    points = geometry()[0]
    return [xor_sum(values), xor_sum(mul(u, x) for u, x in zip(points, values)),
            xor_sum(mul(mul(u, u), x) for u, x in zip(points, values))]


def local_step(values, state):
    output = [x ^ y for x, y in zip(values, emit(state))]
    updated = [a ^ b for a, b in zip(state, moments(values))]
    return output, updated


def global_step(values, state, beta):
    points = geometry()[0]
    output = [x ^ mul(1 ^ mul(u, beta) ^ mul(mul(u, u), mul(beta, beta)), state[0])
              ^ mul(u, state[1]) ^ mul(mul(u, u), state[2]) for u, x in zip(points, values)]
    injected = [xor_sum(values), xor_sum(mul(u ^ beta, x) for u, x in zip(points, values)),
                xor_sum(mul(mul(u ^ beta, u ^ beta), x) for u, x in zip(points, values))]
    return output, [a ^ b for a, b in zip(state, injected)]


def encode(regions, betas):
    state, previous, output = [0, 0, 0], 0, []
    for blocks, beta in zip(regions, betas):
        state = basis(state, previous ^ beta)
        encoded = []
        for block in blocks:
            y, state = local_step(block, state)
            encoded.append(y)
        output.append(encoded)
        previous = beta
    return output


def adjoint(regions, betas):
    state, result = [0, 0, 0], []
    for index in reversed(range(len(regions))):
        values = []
        for block in reversed(regions[index]):
            x = [a ^ b for a, b in zip(block, emit(state))]
            state = [a ^ b for a, b in zip(state, moments(block))]
            values.append(x)
        result.append(list(reversed(values)))
        previous = betas[index-1] if index else 0
        state = basis(state, previous ^ betas[index], transpose=True)
    return list(reversed(result))


def one_row_lower(region_steps, outer_distance=5):
    if region_steps < 1 or not 3 <= outer_distance <= 16:
        raise ValueError('require positive regional steps and outer distance at least three')
    # A first return needs four distinct-region injections, giving two whole
    # interior regions. Without a return, at least d_outer-1 whole regions
    # follow the first injection. One active packet can cost one symbol.
    return min(2, outer_distance-1)*(14*region_steps-1)
