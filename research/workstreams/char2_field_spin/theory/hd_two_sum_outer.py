"""Exact local algebra for a two-sum alternative to the XOR2 outer mixer.

H = I + 11^T + uu^T, u=(0,...,7), over AES GF256. The prospective outer
keeps eight RS[32,16] rows and independent nonzero D1/D2 in every column.
These checks concern one column, not a whole-code distance certificate.
The four fixed products in the explicit circuit are not timing results.
"""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from time import perf_counter

from gf8_packet_maps import MUL, matvec, matmul, rank, xor_sum
from hd_xor_column_probe import mds_prescribed_support


LABELS = tuple(range(8))


def matrix():
    return tuple(tuple(int(i == j) ^ 1 ^ MUL[i][j] for j in LABELS) for i in LABELS)


def direct(values):
    if len(values) != 8 or any(type(v) is not int or not 0 <= v < 256 for v in values):
        raise ValueError('Eight GF256 bytes required')
    s = xor_sum(values)
    t = xor_sum(MUL[i][v] for i, v in enumerate(values))
    return tuple(v ^ s ^ MUL[i][t] for i, v in enumerate(values))


def circuit(values):
    """Four fixed products and28 binary XORs, no extra sampled coefficients."""
    if len(values) != 8 or any(type(v) is not int or not 0 <= v < 256 for v in values):
        raise ValueError('Eight GF256 bytes required')
    x0, x1, x2, x3, x4, x5, x6, x7 = values
    p01, p23, p45, p67 = x0 ^ x1, x2 ^ x3, x4 ^ x5, x6 ^ x7
    b = p23 ^ p67
    s = (p01 ^ p45) ^ b
    c = p45 ^ p67
    a = (x1 ^ x3) ^ (x5 ^ x7)
    t = a ^ MUL[2][b] ^ MUL[4][c]
    t2, t4 = MUL[2][t], MUL[4][t]
    c0, c1, c2 = s, s ^ t, s ^ t2
    c3, c4, c5 = c1 ^ t2, s ^ t4, c1 ^ t4
    c6, c7 = c2 ^ t4, c3 ^ t4
    return tuple(x ^ c for x, c in zip(values, (c0, c1, c2, c3, c4, c5, c6, c7)))


@lru_cache(None)
def uniform_mask_counts(active):
    """Exact output-mask counts for a uniform vector on coordinate space F^A.

    For T, q^(|A|-rank(H_[outside T,A])) counts outputs supported inside T.
    Subset Mobius inversion then gives exact support, with integer arithmetic.
    Zero input entries are allowed; this is not the D1 nonzero-input kernel.
    """
    if type(active) is not int or not 0 <= active < 256:
        raise ValueError('Eight-bit active coordinate mask required')
    columns = [i for i in LABELS if active >> i & 1]
    h = matrix()
    counts = [256**(len(columns)-rank([[h[i][j] for j in columns]
              for i in LABELS if not target >> i & 1])) for target in range(256)]
    for bit in LABELS:
        for target in range(256):
            if target >> bit & 1:
                counts[target] -= counts[target ^ (1 << bit)]
    if min(counts) < 0 or sum(counts) != 256**len(columns):
        raise ArithmeticError('Exact uniform mask inversion failed')
    return tuple(counts)


@lru_cache(None)
def nonzero_mask_counts(active):
    """Exact output-mask counts when every selected input coordinate is nonzero."""
    if type(active) is not int or not 0 <= active < 256:
        raise ValueError('Eight-bit active coordinate mask required')
    result = [0]*256
    subset = active
    while True:
        sign = -1 if (active.bit_count()-subset.bit_count()) % 2 else 1
        for target, count in enumerate(uniform_mask_counts(subset)):
            result[target] += sign*count
        if not subset:
            break
        subset = (subset-1) & active
    if min(result) < 0 or sum(result) != 255**active.bit_count():
        raise ArithmeticError('Exact nonzero input inversion failed')
    return tuple(result)


def silent_moment_count(arity):
    """NZ assignments with both sum x and sum u_i*x_i zero, for any fixed A."""
    if type(arity) is not int or not 0 <= arity <= 8:
        raise ValueError('Arity0..8 required')
    if arity == 0:
        return 1
    if arity < 3:
        return 0
    return mds_prescribed_support(256, arity, arity-2, arity)


def summarize():
    started = perf_counter()
    h = matrix()
    identity = tuple(tuple(int(i == j) for j in LABELS) for i in LABELS)
    if matmul(h, h) != identity or rank(h) != 8:
        raise ArithmeticError('Expected involutory rank-eight mixer')
    rows = []
    for arity in range(1, 9):
        weights, by_mask, lowest = [0]*9, [], 9
        for positions in combinations(LABELS, arity):
            active = sum(1 << i for i in positions)
            counts = nonzero_mask_counts(active)
            local = [0]*9
            for target, count in enumerate(counts):
                local[target.bit_count()] += count
                weights[target.bit_count()] += count
            minimum = next(w for w, count in enumerate(local) if count)
            lowest = min(lowest, minimum)
            by_mask.append(dict(active_mask=active, minimum_output_weight=minimum, weight_counts=local))
        rows.append(dict(active_inputs=arity, minimum_output_weight=lowest,
            minimum_branch_weight=arity+lowest, summed_weight_counts=weights,
            silent_two_moment_count_per_support=silent_moment_count(arity), supports=by_mask))
    return dict(construction='prospective-two-sum-eight-row-outer', field_polynomial='0x11b',
        labels=list(LABELS), involution=True, matrix_rank=8, matrix=h,
        branch_number=min(row['minimum_branch_weight'] for row in rows), entries=rows,
        local_mask_counts_exact=True, fixed_products_per_eight=4, binary_xors_per_eight=28,
        old_xor2_fixed_products_per_eight=0, old_xor2_binary_xors_per_eight=15,
        same_normalized_variable_products_per_eight=15, extra_sampled_coefficients=0,
        native_implemented=False, benchmarked=False, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a local algebra artifact')
    result = summarize()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf8') as target:
        json.dump(result, target, indent=2)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('matrix', 'entries')}))
    for row in result['entries']:
        print(json.dumps({key: value for key, value in row.items() if key != 'supports'}))
