"""Whole-code K20 certificate for GRS[128,83] and repaired cubic packets.

Replays all three disjoint/overlapping covers and checks their geometry
and joins. The guarantee is over once-sampled setup: a GF256 code of
minimum symbol distance at least162202, except with the reported bound.
No timing result or native implementation is claimed.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
import json

from gf8_packet_maps import (CUBIC_RETURN_MINIMA, cubic_matrices,
    cubic_profile, cubic_returned_matrix, determinant, exact_return_minima, matmul)
from gf8_packet_middle_certificate import metadata, certify as middle_certificate, certify_dense
from gf8_packet_sparse_certificate import certify as sparse_certificate


def validate_maps():
    """Finite exact check of every needed minor and returned support bound."""
    actual, expected = cubic_profile(), metadata()
    for key in ('construction', 'width', 'state', 'eta', 'injection_constant'):
        if actual[key] != expected[key]:
            raise ArithmeticError(f'local map/receipt mismatch: {key}')
    if actual['exact_return_minima'] != expected['return_minima']:
        raise ArithmeticError('local return minima do not match the receipt')
    returned = cubic_returned_matrix(0)
    minors = 0
    for phase in range(8):
        a, c = cubic_matrices(phase)
        if any(v for row in matmul(c, a) for v in row):
            raise ArithmeticError('packet maps do not satisfy CA=0')
        if cubic_returned_matrix(phase) != returned:
            raise ArithmeticError('phase return matrix changed')
        for size in range(1, 5):
            for rows in combinations(range(8), size):
                for columns in combinations(range(4), size):
                    if not determinant([[a[i][j] for j in columns] for i in rows]):
                        raise ArithmeticError('emission is not superregular')
                    if not determinant([[c[j][i] for j in columns] for i in rows]):
                        raise ArithmeticError('injection is not superregular')
                    minors += 2
    if exact_return_minima(returned) != CUBIC_RETURN_MINIMA:
        raise ArithmeticError('return minima do not match the finite map')
    return dict(phases=8, checked_nonzero_minors=minors,
                exact_return_minima=list(CUBIC_RETURN_MINIMA))


def combine(sparse, middle, dense):
    expected = metadata()
    keys = ('construction', 'outer', 'q', 'width', 'state', 'k_log',
            'requested_k', 'padded_k', 'rows', 'n', 'output_threshold', 'return_minima',
            'eta', 'injection_constant')
    for receipt in (sparse, middle, dense):
        if receipt.get('complete') is not True or receipt.get('rigorous') is not True:
            raise ValueError('all components must be complete directed certificates')
        for key in keys:
            if receipt.get(key) != expected[key]:
                raise ValueError(f'component geometry mismatch: {key}')
        if receipt['margin_bits_lower'] < 40:
            raise ArithmeticError('component margin below40')
    first_active, last_active = sparse['active_range']
    first_middle, last_middle = middle['weights']
    first_dense, last_dense = dense['weights']
    outer_distance = 128-83+1
    if first_active != 1 or not 1 <= last_active <= expected['rows']:
        raise ValueError('sparse cover must start at one')
    if not 1 <= first_middle <= last_middle < first_dense <= last_dense == expected['n']:
        raise ValueError('invalid middle/dense cover ordering')
    if first_dense != last_middle+1:
        raise ValueError('gap in middle/dense join')
    if outer_distance*(last_active+1) < first_middle:
        raise ValueError('gap in sparse/middle join')
    p, q = (int(value, 16) for value in sparse['failure_upper_dyadic_hex'])
    sparse_upper = F(p, q)
    if sparse_upper <= 0 or sparse_upper > F(1, 1 << sparse['margin_bits_lower']):
        raise ValueError('invalid sparse failure upper bound')
    total = sparse_upper + sum((F(1, 1 << r['margin_bits_lower']) for r in (middle, dense)), F(0))
    bits = min(r['margin_bits_lower'] for r in (sparse, middle, dense))
    while total > F(1, 1 << bits):
        bits -= 1
    if bits < 40:
        raise ArithmeticError('whole-code setup-failure margin below40')
    return dict(**expected, complete=True, rigorous=True, whole_code=True,
                minimum_distance=expected['output_threshold']+1,
                margin_bits_lower=bits, failure_upper_dyadic_exponent=bits,
                component_margins_bits=[r['margin_bits_lower'] for r in (sparse, middle, dense)],
                sparse_active_range=sparse['active_range'], middle_weights=middle['weights'],
                dense_weights=dense['weights'], minimum_outer_distance=outer_distance,
                first_remaining_outer_weight=outer_distance*(last_active+1),
                setup='independent uniform regional permutations and independent nonzero GF256 diagonals',
                scalar_extension='same scalar code on all16 byte lanes; no16-fold failure union',
                native_timing_claim=False)


def certify():
    maps = validate_maps()
    sparse = sparse_certificate()
    middle = middle_certificate()
    dense = certify_dense()
    result = combine(sparse, middle, dense)
    result['local_map_checks'] = maps
    result['component_partition_sha256'] = [sparse['partition_sha256'],
        middle['partition_sha256'], dense['witness_sha256']]
    return result


if __name__ == '__main__':
    print(json.dumps(certify()), flush=True)
