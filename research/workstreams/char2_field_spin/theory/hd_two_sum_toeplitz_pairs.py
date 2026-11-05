"""Exact Toeplitz2 support law and positive fixed-A column contractions.

The exact-support mode does NOT assert independent nonzero amplitudes:
its lazy line branch preserves an input pair's ratio.  The iid-majorant
mode is a pointwise measure bound compatible with an iid-NZ comparator.
This module does not build regional operators or certify distance.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, log2

from hd_two_sum_input_weight import marked_count


MODES = ('exact-support', 'iid-majorant', 'diagonal', 'uniform-cap')


def checked_q(q):
    if type(q) is not int or not 2 <= q <= 256 or q & (q-1):
        raise ValueError('A binary field order2..256 is required')
    return q


def checked_matching(matching, group):
    result = tuple(tuple(pair) for pair in matching)
    if (type(group) is not int or group < 2 or group % 2
            or any(len(pair) != 2 for pair in result)
            or any(type(i) is not int for pair in result for i in pair)
            or sorted(i for pair in result for i in pair) != list(range(group))):
        raise ValueError('An exact perfect matching is required')
    return result


def xor_matching(direction, group=16):
    if (type(group) is not int or group < 2 or group & (group-1)
            or type(direction) is not int or not 1 <= direction < group):
        raise ValueError('A nonzero XOR direction and power-of-two group are required')
    return tuple((i, i ^ direction) for i in range(group) if i < (i ^ direction))


def xor_schedule(group=16):
    """All nonzero XOR matchings, then direction1 once more."""
    return tuple(xor_matching(d, group) for d in range(1, group)) + (xor_matching(1, group),)


@lru_cache(None)
def pair_kernel(q=256, mode='exact-support'):
    """Rows/columns are labelled support masks0,1,2,3; entries are rational.

    Exact support: q^-2 lazy support plus q^-2 times each nonzero vector.
    iid majorant: the lazy double-support line costs an additional q-1.
    uniform cap: the earlier rho=1/[q(q-1)] bound on every nonzero vector.
    """
    checked_q(q)
    if mode not in MODES:
        raise ValueError('Unknown pair-law mode')
    rows = []
    for incoming in range(4):
        if incoming == 0 or mode == 'diagonal':
            rows.append(tuple(F(int(outgoing == incoming)) for outgoing in range(4)))
            continue
        if mode == 'uniform-cap':
            rows.append((F(0), F(1, q), F(1, q), F(q-1, q)))
            continue
        row = [F(0), F(q-1, q*q), F(q-1, q*q), F((q-1)**2, q*q)]
        lazy = q-1 if mode == 'iid-majorant' and incoming == 3 else 1
        row[incoming] += F(lazy, q*q)
        rows.append(tuple(row))
    return tuple(rows)


def pair_mark(incoming, left, right, *, q=256, mode='exact-support'):
    if type(incoming) is not int or not 0 <= incoming <= 3:
        raise ValueError('A labelled two-coordinate support is required')
    left, right = F(left), F(right)
    if min(left, right) < 0:
        raise ValueError('Nonnegative output markers are required')
    row = pair_kernel(q, mode)[incoming]
    return row[0] + row[1]*left + row[2]*right + row[3]*left*right


def support_factor(support, matching, markers, *, q=256, mode='exact-support'):
    markers = tuple(F(x) for x in markers)
    matching = checked_matching(matching, len(markers))
    if type(support) is not int or not 0 <= support < 1 << len(markers):
        raise ValueError('A coordinate support mask is required')
    result = F(1)
    for left, right in matching:
        incoming = ((support >> left) & 1) | (((support >> right) & 1) << 1)
        result *= pair_mark(incoming, markers[left], markers[right], q=q, mode=mode)
    return result


def _coefficient_table(active, matching, markers, mode, pair_mask=None):
    group = len(markers)
    if (group not in (8, 16) or type(active) is not int
            or not 0 <= active < 1 << group):
        raise ValueError('A fixed g8/g16 input coordinate space is required')
    matching = checked_matching(matching, group)
    if pair_mask is not None and (type(pair_mask) is not int or not 0 <= pair_mask < 1 << len(matching)):
        raise ValueError('A matching-pair support mask is required')
    table = {(0, 0): F(1)}
    for index, (left, right) in enumerate(matching):
        branches = []
        for local in range(4):
            if pair_mask is not None and bool(local) != bool(pair_mask >> index & 1):
                continue
            bits = ((local & 1) << left) | (((local >> 1) & 1) << right)
            mass = (F(1) if pair_mask is not None else
                    pair_mark(local, markers[left], markers[right], mode=mode))
            if mass:
                branches.append((local.bit_count(), (active & bits).bit_count(), mass))
        next_table = {}
        for (weight, overlap), coefficient in table.items():
            for dw, db, mass in branches:
                key = weight+dw, overlap+db
                next_table[key] = next_table.get(key, F(0)) + coefficient*mass
        table = next_table
    return table


def column_contraction(active, matching, markers, *, input_marker=F(1), mode='exact-support'):
    """Sum_{x in F256^A} v^wt(x) E[prod_i marker_i^1_{output_i!=0}].

    Exact-support gives equality. iid-majorant/uniform-cap give upper bounds.
    The fixed-A word normalization255^(-8|A|), row floors and alpha are NOT
    applied here.  Keep A and all markers fixed through their claimed word.
    """
    markers = tuple(F(x) for x in markers)
    if any(x < 0 for x in markers) or mode not in MODES:
        raise ValueError('Nonnegative markers and a declared law are required')
    table = _coefficient_table(active, matching, markers, mode)
    arity, group = active.bit_count(), len(markers)
    return sum((coefficient*marked_count(arity, overlap, weight, F(input_marker), group)
                for (weight, overlap), coefficient in table.items()), F(0))


def pair_mask_count(active, matching, selected, *, group=16, input_marker=F(1)):
    """Exact F_A(B;v), with B a fixed set of nonzero matching pairs after H."""
    table = _coefficient_table(active, matching, (F(1),)*group, 'exact-support', selected)
    return sum((coefficient*marked_count(active.bit_count(), overlap, weight, F(input_marker), group)
                for (weight, overlap), coefficient in table.items()), F(0))


def quiet_pair_floor(arity):
    """Necessary active-pair floor for a word quiet in all16 columns.

    Every active RS[16,8] row has>=9 nonzeros.  The matching schedule is all
    15 XOR directions plus one repeat.  This is NOT a nonquiet-word floor.
    """
    if type(arity) is not int or not 0 <= arity <= 16:
        raise ValueError('Arity0..16 is required')
    return max((9*arity+1)//2, 9*arity-comb(arity, 2)-arity//2)


def diagnostics():
    """Bounded local gate only: no regional provider, optimization or cover."""
    q, active = 256, 255
    rows = []
    for z in (F(1, 2), F(49, 50)):
        factors = {}
        for mode in MODES:
            value = F(1)
            for matching in xor_schedule():
                value *= support_factor(active, matching, (z,)*16, mode=mode)
            factors[mode] = log2(value)
        rows.append({'output_marker': str(z), 'fixed_pre_toeplitz_support': active,
                     'columns': 16, 'log2_word_support_moment': factors,
                     'exact_improvement_over_diagonal_bits': factors['diagonal']-factors['exact-support']})
    contraction_rows = []
    for active in (1, 3, 15, 255, 65535):
        matching, marks = xor_matching(8), (F(1, 2),)*16
        masses = {mode: column_contraction(active, matching, marks, mode=mode) for mode in MODES}
        contraction_rows.append({'active': active, 'arity': active.bit_count(),
                                 'matching_direction': 8, 'output_marker': '1/2',
                                 'log2_uniform_input_moment': {mode: log2(value)-8*active.bit_count()
                                                             for mode, value in masses.items()}})
    return {'method': 'exact-local-toeplitz2-mixture-and-positive-fixed-A-contraction',
            'exact_support_kernel_normalized': True,
            'exact_support_kernel_implies_iid_amplitudes': False,
            'iid_majorant_available': True, 'regional_banks_built': 0,
            'witness_searches': 0, 'directed_distance_receipt': False,
            'whole_code_certificate': False, 'fixed_mask_word_controls': rows,
            'uniform_input_column_controls': contraction_rows,
            'quiet_word_pair_floors': {str(a): quiet_pair_floor(a) for a in (3, 4, 5, 8, 16)}}


if __name__ == '__main__':
    import argparse
    from hashlib import sha256
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new artifact')
    result = diagnostics()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    data = json.dumps(result, indent=2, allow_nan=False)+'\n'
    with args.output.open('x', encoding='utf-8', newline='\n') as destination:
        destination.write(data)
    print(json.dumps({'output': str(args.output), 'sha256': sha256(args.output.read_bytes()).hexdigest(),
                      'results': result}, allow_nan=False))
