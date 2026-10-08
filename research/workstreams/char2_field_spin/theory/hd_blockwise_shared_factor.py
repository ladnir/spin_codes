"""Exact low-rank and cube circuits for the shared eight-point map.

This is a bounded scalar algebra prototype, not native code or a performance
measurement. Matrices are over AES GF256; they extend linearly to GF65536.
"""

from functools import lru_cache

from hd_blockwise_maps import matmul, matvec, mul, power, xor_sum
from hd_blockwise_shared_maps import (
    basis, evaluation, identity, inverse, rank, subspace_polynomial,
)


MOMENT_MASKS = tuple(i for i in range(64) if i.bit_count() <= 2)
MOMENT_TO_HIGH = (
    (1, 1, 1, 1),
    (1, 0, 1, 0),
    (1, 1, 0, 0),
    (1, 0, 0, 0),
)


def matrix_xor(a, b):
    return tuple(tuple(x ^ y for x, y in zip(ar, br)) for ar, br in zip(a, b))


def evaluation_at(offset, size=8):
    """Unlike the earlier helper, also allow nonaligned translated nodes."""
    return tuple(tuple(basis(j, offset ^ i) for j in range(size))
                 for i in range(size))


@lru_cache(None)
def tile_matrix(beta=64):
    if not 8 <= beta < 256:
        raise ValueError('translation must be outside the input eight-space')
    return matmul(evaluation_at(beta), inverse(evaluation(0)))


@lru_cache(None)
def pair_factor(beta=64):
    """F=I+UV, U=[D;BD], V=[B,I], using four-coordinate blocks."""
    n = matrix_xor(tile_matrix(beta), identity(8))
    e = tuple(tuple(row[:4]) for row in n[:4])
    d = tuple(tuple(row[4:]) for row in n[:4])
    b = matmul(inverse(d), e)
    u = d + matmul(b, d)
    v = tuple(tuple(b[i]) + identity(4)[i] for i in range(4))
    return u, v, b, d


@lru_cache(None)
def polynomial_factor(beta=64):
    """Use p(x)=a(z)+x*b(z), z=x*x+beta*x; U evaluates beta*b(z)."""
    change = []
    for x in range(8):
        z = mul(x, x) ^ mul(beta, x)
        powers = tuple(power(z, j) for j in range(4))
        change.append(powers + tuple(mul(x, value) for value in powers))
    inv = inverse(tuple(change))
    u = tuple(tuple(mul(beta, change[i][j]) for j in range(4)) for i in range(8))
    return u, inv[4:]


@lru_cache(None)
def cube_coefficients(beta=64):
    """F=sum_S a[S] product_{i in S}(I+P_(2**i)); a[0]=1."""
    a = list(tile_matrix(beta)[0])
    for bit in range(3):
        for mask in range(8):
            if not mask & (1 << bit):
                a[mask] ^= a[mask | (1 << bit)]
    return tuple(a)


def cube_differences(values):
    """Keep one copy of each repeated face/edge/whole-cube sum (19 records)."""
    if len(values) != 8:
        raise ValueError('eight input coordinates required')
    return {(mask, base): xor_sum(values[base ^ sub] for sub in range(8)
                                 if sub & ~mask == 0)
            for mask in range(1, 8) for base in range(8) if base & mask == 0}


def cube_feedthrough(values, beta=64, differences=None):
    differences = cube_differences(values) if differences is None else differences
    a = cube_coefficients(beta)
    products = {(mask, base): mul(a[mask], value)
                for (mask, base), value in differences.items()}
    return tuple(values[j] ^ xor_sum(products[mask, j & ~mask]
                                    for mask in range(1, 8)) for j in range(8))


def anf(values):
    """Multilinear Boolean coefficients, with field addition equal to XOR."""
    result = list(values)
    if not result or len(result) & (len(result) - 1):
        raise ValueError('power-of-two table required')
    for bit in range(len(result).bit_length() - 1):
        for mask in range(len(result)):
            if mask & (1 << bit):
                result[mask] ^= result[mask ^ (1 << bit)]
    return tuple(result)


@lru_cache(None)
def moment_coefficients():
    return tuple(anf(tuple(basis(j, x) for x in range(64))) for j in range(4))


@lru_cache(None)
def emission_coefficients(beta=64):
    return tuple(anf(tuple(basis(j, beta ^ x) for x in range(64))) for j in range(4))


def accumulate_moments(tile_differences):
    if len(tile_differences) != 8:
        raise ValueError('eight tiles required')
    # Every requested low-bit mask has size <=2, so its complementary cube
    # difference is nonempty and was already computed for feedthrough.
    return {mask: xor_sum(differences[7 ^ (mask & 7), mask & 7]
                         for tile, differences in enumerate(tile_differences)
                         if tile & (mask >> 3) == mask >> 3)
            for mask in MOMENT_MASKS}


def feedback_from_moments(moments):
    coefficients = moment_coefficients()
    low = tuple(xor_sum(mul(coefficients[j][mask], moments[mask])
                        for mask in MOMENT_MASKS) for j in range(4))
    return matvec(MOMENT_TO_HIGH, low)


def emission_from_anf(state, beta=64):
    if len(state) != 4:
        raise ValueError('four state coordinates required')
    coefficients = emission_coefficients(beta)
    # Reuse each identical (state coordinate, constant) multiplication.
    bank = {(j, constant): mul(constant, state[j]) for j in range(4)
            for constant in set(coefficients[j]) if constant not in (0, 1)}
    polynomial = {}
    for mask in MOMENT_MASKS:
        value = 0
        for j in range(4):
            c = coefficients[j][mask]
            value ^= state[j] if c == 1 else bank[j, c] if c else 0
        polynomial[mask] = value
    return tuple(xor_sum(value for mask, value in polynomial.items()
                         if x & mask == mask) for x in range(64))


def packet_circuit(values, state, beta=64):
    """Return exact emitted packet and feedback; M*state is deliberately external."""
    if len(values) != 64 or len(state) != 4:
        raise ValueError('64 inputs and four state coordinates required')
    differences = tuple(cube_differences(values[8*i:8*i+8]) for i in range(8))
    feedthrough = tuple(value for i in range(8)
                        for value in cube_feedthrough(values[8*i:8*i+8], beta,
                                                      differences[i]))
    emission = emission_from_anf(state, beta)
    return tuple(x ^ y for x, y in zip(feedthrough, emission)), \
        feedback_from_moments(accumulate_moments(differences))


def reference_packet(values, state, beta=64):
    out = []
    feedback = [0] * 4
    for alpha in range(0, 64, 8):
        coefficients = list(matvec(inverse(evaluation(alpha)), values[alpha:alpha+8]))
        for j in range(4):
            feedback[j] ^= coefficients[j+4]
            coefficients[j] ^= state[j]
        out.extend(matvec(evaluation_at(alpha ^ beta), tuple(coefficients)))
    return tuple(out), tuple(feedback)


def product_receipt(beta=64):
    a = cube_coefficients(beta)
    feedthrough = sum(8 >> mask.bit_count() for mask in range(1, 8)
                      if a[mask] not in (0, 1))
    feedback = sum(c not in (0, 1) for row in moment_coefficients() for c in row)
    emission = sum(len(set(row) - {0, 1}) for row in emission_coefficients(beta))
    refresh = 4  # Two genuine-GF65536 packed products, two GFNI each.
    return dict(beta=beta, tile_cube_coefficients=a,
                feedthrough_record_products_per_tile=feedthrough,
                feedthrough_GFNI_per_packet=8*((feedthrough+3)//4),
                feedback_record_products_per_packet=feedback,
                feedback_GFNI_per_packet=(feedback+3)//4,
                emission_record_products_per_packet=emission,
                emission_GFNI_per_packet=(emission+3)//4,
                rankone_refresh_GFNI_per_packet=refresh,
                packed_arithmetic_GFNI_per_packet=(8*((feedthrough+3)//4)
                    +(feedback+3)//4+(emission+3)//4+refresh),
                moment_records=22, native_instruction_count=False,
                native_timing=False, whole_distance_certificate=False)


def receipt():
    data = product_receipt(64)
    u, v, _, _ = pair_factor(64)
    c = inverse(evaluation(0))[4:]
    data.update(rank_local_residual=rank(v),
                rank_residual_and_feedback=rank(v+c),
                VU_zero=not any(x for row in matmul(v,u) for x in row),
                packet_beta8_count=product_receipt(8)['packed_arithmetic_GFNI_per_packet'],
                packet_beta128_count=product_receipt(128)['packed_arithmetic_GFNI_per_packet'])
    return data


if __name__ == '__main__':
    import json
    print(json.dumps(receipt(), indent=2))
