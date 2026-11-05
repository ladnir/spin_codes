"""One-fixed-product affine-syndrome mixer with exact a1/a2 raw counts.

H_lambda=I+lambda*11^T+sum_k b_k*b_k^T, lambda outside{0,1}.
The first implementation fixes the AES field for q256; q4 uses x^2+x+1
for exhaustive toy validation. Every source ordering is caller supplied.
No native code, sampler normalization, regional provider, or certificate.
"""
from collections import Counter
from functools import lru_cache

from hd_outer_pair_layers import field_multiply

CONSTRUCTION = 'G16-RS16x8-D2-cross-pair-Hlambda-D1-research'
FULL = (1 << 16)-1
MODULI = {4: 0x7, 256: 0x11b}


def _parameters(lambda_value, q):
    if type(q) is not int or q not in MODULI:
        raise ValueError('Only the explicit GF4 and AES-GF256 fields are implemented')
    if type(lambda_value) is not int or not 1 < lambda_value < q:
        raise ValueError('lambda must be a field element outside{0,1}')
    return MODULI[q]


def _multiply(a, b, q):
    return field_multiply(a, b, modulus=MODULI[q])


def _inverse(a, q):
    if not a:
        raise ValueError('Cannot invert zero')
    result, exponent = 1, q-2
    while exponent:
        if exponent & 1:
            result = _multiply(result, a, q)
        a = _multiply(a, a, q)
        exponent >>= 1
    return result


@lru_cache(maxsize=8)
def matrix(lambda_value=2, q=256):
    _parameters(lambda_value, q)
    return tuple(tuple(lambda_value ^ int(i == j) ^ ((i & j).bit_count() & 1)
                       for j in range(16)) for i in range(16))


def direct(values, lambda_value=2, q=256):
    _parameters(lambda_value, q)
    values = tuple(values)
    if len(values) != 16 or any(type(x) is not int or not 0 <= x < q for x in values):
        raise ValueError('Sixteen values in the declared field are required')
    output = []
    for row in matrix(lambda_value, q):
        value = 0
        for coefficient, x in zip(row, values):
            if x:
                value ^= _multiply(coefficient, x, q)
        output.append(value)
    return tuple(output)


def _circuit(x, add, multiply, lambda_value):
    p = tuple(add(x[2*i], x[2*i+1]) for i in range(8))
    a0 = add(add(add(x[1], x[3]), add(x[5], x[7])), add(add(x[9], x[11]), add(x[13], x[15])))
    a1 = add(add(p[1], p[3]), add(p[5], p[7]))
    q = tuple(add(p[2*i], p[2*i+1]) for i in range(4))
    a2, a3 = add(q[1], q[3]), add(q[2], q[3])
    total = multiply(lambda_value, add(add(q[0], q[1]), a3))
    syndromes = (a0, a1, a2, a3)
    correction = [total]*16
    for i in range(1, 16):
        low = i & -i
        correction[i] = add(correction[i ^ low], syndromes[low.bit_length()-1])
    return tuple(add(x[i], correction[i]) for i in range(16))


def circuit(values, lambda_value=2, q=256):
    _parameters(lambda_value, q)
    values = tuple(values)
    if len(values) != 16 or any(type(x) is not int or not 0 <= x < q for x in values):
        raise ValueError('Sixteen values in the declared field are required')
    return _circuit(values, lambda a, b: a ^ b,
                    lambda a, b: _multiply(a, b, q), lambda_value)


def primitive_counts():
    counts = Counter()
    def add(a, b):
        counts['binary_xors'] += 1
        return a ^ b
    def multiply(a, b):
        counts['fixed_products'] += 1
        return _multiply(a, b, 256)
    _circuit((0,)*16, add, multiply, 2)
    return dict(binary_xors=counts['binary_xors'], fixed_products=counts['fixed_products'],
                variable_products=0, includes_diagonals=False,
                includes_pair_layer=False, native_performance_claim=False)


def exceptional_ratios(lambda_value=2, q=256):
    """Fixed input ratio x_first/x_second, in controller branch order."""
    _parameters(lambda_value, q)
    return (1, _multiply(lambda_value, _inverse(lambda_value ^ 1, q), q),
            _multiply(lambda_value ^ 1, _inverse(lambda_value, q), q))


def column_branches(source, lambda_value=2, q=256):
    """Unmerged ordered (original_input_weight, output_mask, raw_count).

    a1 has branches zero,single. a2 has zero,left-single,right-single,
    ratio1,ratio(lambda/(lambda+1)),reciprocal,remaining. Coincident masks
    and zero-count branches are deliberately retained. Input source order
    is never sorted, so a caller can apply a column-dependent permutation
    without exchanging the original left/right row's activity markers.
    """
    _parameters(lambda_value, q)
    source = tuple(source)
    if (len(source) not in (1, 2) or len(set(source)) != len(source)
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('One or two ordered distinct source labels are required')
    d = q-1
    if len(source) == 1:
        return ((0, 0, 1), (1, FULL, d))
    first, second = source
    columns = matrix(lambda_value, q)
    masks = tuple(sum(bool(_multiply(row[first], ratio, q) ^ row[second]) << i
                      for i, row in enumerate(columns))
                  for ratio in exceptional_ratios(lambda_value, q))
    result = ((0, 0, 1), (1, FULL, d), (1, FULL, d),
              *((2, mask, d) for mask in masks), (2, FULL, d*(q-4)))
    if sum(count for _, _, count in result) != q*q or any(weight and not mask for weight, mask, count in result if count):
        raise AssertionError('Allowed-input mass or invertibility failed')
    return result


def joint_counts(source, lambda_value=2, q=256):
    """Merged raw counts; no probability, RS, row-choice, or projective factor."""
    counts = Counter()
    for weight, mask, count in column_branches(source, lambda_value, q):
        if count:
            counts[weight, mask] += count
    return dict(sorted(counts.items()))
