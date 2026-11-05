"""Exact binary rank-five mixer and two-source joint counting interfaces.

Hbin=I+11^T+sum_k b_k b_k^T on the sixteen binary4 labels.
It is not H*: the field-weighted u*u^T correction has been removed.
No regional provider, native code, or distance claim is included.
"""
from collections import Counter
from functools import lru_cache

CONSTRUCTION = 'G16-RS16x8-D2-cross-pair-binary-rank5-D1-research'


@lru_cache(maxsize=2)
def matrix(mixer='binary_syndrome'):
    if mixer not in ('binary_syndrome', 'parity'):
        raise ValueError('Choose binary_syndrome or parity')
    return tuple(tuple(int(i == j) ^ 1 ^ (((i & j).bit_count() & 1) if mixer == 'binary_syndrome' else 0)
                       for j in range(16)) for i in range(16))


def column_mask(index, *, mixer='binary_syndrome'):
    if type(index) is not int or not 0 <= index < 16:
        raise ValueError('A row index0..15 is required')
    return sum(row[index] << i for i, row in enumerate(matrix(mixer)))


def direct(values, *, mixer='binary_syndrome'):
    values = tuple(values)
    if len(values) != 16 or any(type(x) is not int or x < 0 for x in values):
        raise ValueError('Sixteen binary-field representations required')
    output = []
    for row in matrix(mixer):
        value = 0
        for coefficient, x in zip(row, values):
            if coefficient:
                value ^= x
        output.append(value)
    return tuple(output)


def _circuit(x, add):
    p = tuple(add(x[2*i], x[2*i+1]) for i in range(8))
    a0 = add(add(add(x[1], x[3]), add(x[5], x[7])), add(add(x[9], x[11]), add(x[13], x[15])))
    a1 = add(add(p[1], p[3]), add(p[5], p[7]))
    q = tuple(add(p[2*i], p[2*i+1]) for i in range(4))
    a2, a3 = add(q[1], q[3]), add(q[2], q[3])
    total = add(add(q[0], q[1]), a3)
    syndromes = (a0, a1, a2, a3)
    correction = [total]*16
    for i in range(1, 16):
        low = i & -i
        correction[i] = add(correction[i ^ low], syndromes[low.bit_length()-1])
    return tuple(add(x[i], correction[i]) for i in range(16))


def circuit(values):
    values = tuple(values)
    if len(values) != 16 or any(type(x) is not int or x < 0 for x in values):
        raise ValueError('Sixteen binary-field representations required')
    return _circuit(values, lambda a, b: a ^ b)


def primitive_counts():
    counts = Counter()
    def add(a, b):
        counts['binary_xors'] += 1
        return a ^ b
    _circuit((0,)*16, add)
    return dict(binary_xors=counts['binary_xors'], variable_products=0, fixed_products=0,
                includes_diagonals=False, includes_pair_layer=False, native_performance_claim=False)


def two_row_joint_counts(first, second, *, mixer='binary_syndrome', q=256):
    """Raw allowed-input (input_weight,output_mask)->count for two sources.

    Counts include the unique zero, two singleton branches, equal-NZ and
    distinct-NZ branches.  The actual row pair stays fixed across columns.
    Input weight is pre-mixer; sum it across columns to retain the two RS
    rows' aggregate floor18.  No probability or RS normalizer is inserted.
    """
    if (type(first) is not int or type(second) is not int or not 0 <= first < second < 16
            or type(q) is not int or q < 2 or q & (q-1)):
        raise ValueError('Ordered distinct row indices and a binary-field order required')
    a, b = column_mask(first, mixer=mixer), column_mask(second, mixer=mixer)
    d = q-1
    result = Counter({(0, 0): 1})
    for key, count in (((1, a), d), ((1, b), d), ((2, a ^ b), d), ((2, a | b), d*(q-2))):
        if count:
            result[key] += count
    if sum(result.values()) != q*q or any(not mask for (weight, mask) in result if weight):
        raise AssertionError('Two-source mass or injectivity changed')
    return dict(sorted(result.items()))


def two_row_support_counts(first, second, *, mixer='binary_syndrome', q=256):
    result = Counter()
    for (_, mask), count in two_row_joint_counts(first, second, mixer=mixer, q=q).items():
        result[mask] += count
    return dict(sorted(result.items()))
