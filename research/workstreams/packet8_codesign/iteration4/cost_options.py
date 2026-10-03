"""Exact small algebra and source-level operation counts; no timing model.

All field bytes use the AES polynomial basis. Scalar adjoints are binary
transposes under the coordinate dot product, not ordinary multiplication.
This file implements no encoder kernel and changes no construction default.
"""
from __future__ import annotations


def mul8(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a = ((a << 1) ^ (0x11b if a & 128 else 0)) & 255
    return result


def adj8(coefficient, value):
    return sum(((mul8(coefficient, 1 << bit) & value).bit_count() & 1) << bit
        for bit in range(8))


def bytes_of(value, count):
    return tuple((value >> (8*i)) & 255 for i in range(count))


def word_of(values):
    return sum(value << (8*i) for i, value in enumerate(values))


def expansion24(state, scales=(1,)*8):
    if len(scales) != 8 or any(not 1 <= x <= 255 for x in scales):
        raise ValueError('eight nonzero byte scales required')
    a, b, c = bytes_of(state, 3)
    return word_of(mul8(scales[h], a ^ mul8(h, b) ^ mul8(mul8(h, h), c))
        for h in range(8))


def feedback24(word):
    x = bytes_of(word, 8)
    result = [0, 0, 0]
    for h, value in enumerate(x):
        result[0] ^= value
        result[1] ^= mul8(h, value)
        result[2] ^= mul8(mul8(h, h), value)
    return word_of(result)


def literal_expansion_adjoint(word, scales=(1,)*8):
    return sum(((expansion24(1 << bit, scales) & word).bit_count() & 1) << bit
        for bit in range(24))


def two_band_expansion_adjoint(word, d):
    """A'^T via three scaled high-band Boolean moments, then four maps."""
    if not 1 <= d <= 255:
        raise ValueError('nonzero high-band scale required')
    x = bytes_of(word, 8)
    low0, high0 = x[0]^x[1]^x[2]^x[3], x[4]^x[5]^x[6]^x[7]
    low_odd, high_odd = x[1]^x[3], x[5]^x[7]
    low_bit1, high_bit1 = x[2]^x[3], x[6]^x[7]
    # Three new byte affine applications; each is applied to both payload halves.
    scaled0, scaled_odd, scaled_bit1 = adj8(d, high0), adj8(d, high_odd), adj8(d, high_bit1)
    total, odd, bit1, bit2 = low0^scaled0, low_odd^scaled_odd, low_bit1^scaled_bit1, scaled0
    return word_of((total, odd^adj8(2, bit1)^adj8(4, bit2),
        odd^adj8(4, bit1)^adj8(16, bit2)))


def mul24_reference(left, right):
    """Ordinary polynomial multiplication modulo z^3+z+1 over GF256."""
    a, b = bytes_of(left, 3), bytes_of(right, 3)
    coefficients = [0]*5
    for i in range(3):
        for j in range(3):
            coefficients[i+j] ^= mul8(a[i], b[j])
    for degree in (4, 3):
        coefficients[degree-3] ^= coefficients[degree]
        coefficients[degree-2] ^= coefficients[degree]
    return word_of(coefficients[:3])


def mul24_six(left, right):
    a, b, c = bytes_of(left, 3)
    r0, r1, r2 = bytes_of(right, 3)
    p0, p1, p2 = mul8(a, r0), mul8(b, r1), mul8(c, r2)
    p01, p02, p12 = mul8(a^b, r0^r1), mul8(a^c, r0^r2), mul8(b^c, r1^r2)
    return word_of((p0^p12^p1^p2, p01^p0^p12, p02^p0^p1))


def adj24_six(value, scalar):
    """Binary adjoint of the six-product circuit; also six byte maps."""
    t0, t1, t2 = bytes_of(value, 3)
    r0, r1, r2 = bytes_of(scalar, 3)
    q0 = adj8(r0, t0^t1^t2)
    q1 = adj8(r1, t0^t2)
    q2 = adj8(r2, t0)
    q01, q02, q12 = adj8(r0^r1, t1), adj8(r0^r2, t2), adj8(r1^r2, t0^t1)
    return word_of((q0^q01^q02, q1^q01^q12, q2^q02^q12))


def pow24(value, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = mul24_six(result, value)
        exponent >>= 1
        value = mul24_six(value, value)
    return result


def mul16_reference(left, right):
    """GF256[u]/(u^2+u+0x20), matching the existing outer tower."""
    a, b = bytes_of(left, 2)
    c, d = bytes_of(right, 2)
    return word_of((mul8(a, c)^mul8(0x20, mul8(b, d)),
        mul8(a, d)^mul8(b, c)^mul8(b, d)))


def mul16_three(left, right):
    a, b = bytes_of(left, 2)
    c, d = bytes_of(right, 2)
    p0, p1, p2 = mul8(a, c), mul8(b, mul8(0x20, d)), mul8(a^b, c^d)
    return word_of((p0^p1, p0^p2))


def adj16_three(value, scalar):
    y0, y1 = bytes_of(value, 2)
    c, d = bytes_of(scalar, 2)
    q0, q1, q2 = adj8(c, y0^y1), adj8(mul8(0x20, d), y0), adj8(c^d, y1)
    return word_of((q0^q2, q1^q2))


def gfni_counts(k=65536):
    """Static instruction applications, not cycles; 128-bit payloads only."""
    if type(k) is not int or k < 256 or k % 256:
        raise ValueError('K must be a positive multiple of256')
    steps, small_groups, wide_groups = 2*k//64, k//128, k//256
    small_outer = small_groups*(16*2*3 + 4*15 + 4*8)
    wide_outer = wide_groups*(16*2*9 + 8*15 + 8*8)
    kernels = {
        'baseline16_GL2': 20*steps + 12*(steps-1),
        'state24_GL3': 24*steps + 26*(steps-1),
        'state24_scalar': 24*steps + 20*(steps-1),
        'state24_scalar_two_band': 24*steps + 26*(steps-1),
        'state24_scalar_seven_scales': 24*steps + 34*(steps-1),
        'width32_state16_GL2': 20*(2*steps)-10,
        'width32_state16_scalar': 18*(2*steps)-8,
    }
    baseline = kernels['baseline16_GL2']+small_outer
    return dict(physical_steps=steps, update_steps=steps-1, inner_gfni=kernels,
        small_outer_gfni=small_outer, wide_outer_shared_gfni=wide_outer,
        wide_outer_old_parity_penalty=wide_groups*8*3,
        outer_randomizer_extra_xors=wide_groups*16*2*15-small_groups*16*2*3,
        baseline_total_gfni=baseline,
        wide_plus_scalar24_total_gfni=wide_outer+kernels['state24_scalar'],
        wide_plus_scalar24_extra_gfni=wide_outer+kernels['state24_scalar']-baseline,
        changed_outer_needs_native_byte_input=True,
        instruction_counts_are_not_timing_predictions=True)
