"""Exact small-field checks and certificate transport for byte MDS sandwiches."""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, log2
from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'iteration5'))
import scaled_positive as sp
import outward_positive as op


def mul(a, b, modulus):
    result = 0
    degree = modulus.bit_length()-1
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & (1 << degree):
            a ^= modulus
    return result


def inverse(a, modulus):
    if not a:
        raise ZeroDivisionError('nonzero field element required')
    q = 1 << (modulus.bit_length()-1)
    result = 1
    for _ in range(q-2):
        result = mul(result, a, modulus)
    return result


def rank(matrix, modulus):
    a = [list(row) for row in matrix]
    rows, columns = len(a), len(a[0]) if a else 0
    r = 0
    for c in range(columns):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = inverse(a[r][c], modulus)
        a[r] = [mul(x, scale, modulus) for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                scale = a[i][c]
                a[i] = [x ^ mul(scale, y, modulus) for x, y in zip(a[i], a[r])]
        r += 1
        if r == rows:
            break
    return r


def mds_minors(matrix, modulus):
    n = len(matrix)
    checked = 0
    for k in range(1, n+1):
        for rows in combinations(range(n), k):
            for columns in combinations(range(n), k):
                checked += 1
                minor = [[matrix[r][c] for c in columns] for r in rows]
                if rank(minor, modulus) != k:
                    return False, dict(checked=checked, rows=rows, columns=columns)
    return True, dict(checked=checked)


def matvec(matrix, x, modulus):
    out = []
    for row in matrix:
        value = 0
        for a, b in zip(row, x):
            value ^= mul(a, b, modulus)
        out.append(value)
    return tuple(out)


def tower256_mul(a, b):
    """GF16[t]/(t^2+t+8), with GF16 modulus0x13 and low nibble first."""
    a0,a1,b0,b1 = a&15,a>>4,b&15,b>>4
    cross = mul(a0,b1,0x13) ^ mul(a1,b0,0x13)
    high = mul(a1,b1,0x13)
    low = mul(a0,b0,0x13) ^ mul(8,high,0x13)
    return low | ((cross ^ high) << 4)


def support(x):
    return sum(bool(value) << i for i, value in enumerate(x))


def sandwich_point_masses(matrix, modulus):
    """Enumerate D_in's exact support law; D_out labels are counted analytically."""
    n = len(matrix)
    q = 1 << (modulus.bit_length()-1)
    counts = Counter()
    for x in product(range(q), repeat=n):
        if any(x):
            counts[support(x), support(matvec(matrix, x, modulus))] += 1
    maximum = F(0)
    for (before, after), count in counts.items():
        # Each fixed input with this support becomes uniform on nonzero
        # labels there. Each fixed output label then has this exact mass.
        denominator = (q-1)**(before.bit_count()+after.bit_count())
        maximum = max(maximum, F(count, denominator))
    return maximum, counts


def normalized_diagonal_ensemble(matrix, modulus):
    """Tiny exact check: fixing D_in[0]=1 leaves the full matrix distribution."""
    n = len(matrix)
    q = 1 << (modulus.bit_length()-1)
    original, normalized = Counter(), Counter()
    for left in product(range(1, q), repeat=n):
        for right in product(range(1, q), repeat=n):
            key = tuple(mul(mul(left[i], matrix[i][j], modulus), right[j], modulus)
                        for i in range(n) for j in range(n))
            original[key] += 1
        for tail in product(range(1, q), repeat=n-1):
            right = (1,)+tail
            key = tuple(mul(mul(left[i], matrix[i][j], modulus), right[j], modulus)
                        for i in range(n) for j in range(n))
            normalized[key] += 1
    return original, normalized


def rank_two_search():
    """Bounded search H=D+2J+4uv^T with u=v=(0,0,1,1), D binary."""
    table = np.array([[mul(a, b, 0x11b) for b in range(256)] for a in range(256)], dtype=np.uint8)
    numbers = np.arange(1 << 16, dtype=np.uint32)
    bits = ((numbers[:, None] >> np.arange(16)) & 1).astype(np.uint8).reshape(-1, 4, 4)
    base = np.full((4, 4), 2, dtype=np.uint8)
    base[2:, 2:] ^= 4
    matrices = bits ^ base
    indices = np.arange(len(numbers))
    for rows in combinations(range(4), 2):
        for cols in combinations(range(4), 2):
            a, b, c, d = (matrices[indices, rows[0], cols[0]], matrices[indices, rows[0], cols[1]],
                           matrices[indices, rows[1], cols[0]], matrices[indices, rows[1], cols[1]])
            indices = indices[(table[a, d] ^ table[b, c]) != 0]
            if not len(indices):
                return dict(two_by_two_survivors=0, mds_survivors=0)
    survivors = []
    for i in indices:
        matrix = matrices[i].tolist()
        if mds_minors(matrix, 0x11b)[0]:
            binary = bits[i].tolist()
            # Three XOR/XOR3 ops for the two common sums, then ceil(m/2)
            # per output for its binary row plus the shared product(s).
            alu = 3+sum((sum(row)+(r >= 2)+1)//2 for r, row in enumerate(binary))
            survivors.append((alu, int(i), matrix, binary))
    survivors.sort()
    return dict(two_by_two_survivors=len(indices), mds_survivors=len(survivors),
                best=survivors[:4])


def decode(record):
    return F(int(record['numerator_hex'], 16))*F(2)**record['exponent']


def encode(value):
    value = F(value)
    if value.denominator & (value.denominator-1):
        raise ValueError('dyadic endpoint required')
    return dict(numerator_hex=hex(value.numerator), exponent=1-value.denominator.bit_length())


def margin(value):
    value = F(value)
    return log2(value.denominator)-log2(value.numerator) if value else float('inf')


def correct_receipt(path, gamma, *, symbols=16):
    """Transport authenticated selected endpoints; no inner or placement replay."""
    gamma = F(gamma)
    if gamma < 1:
        raise ValueError('domination factor must be at least1')
    op.check_runtime()
    saved = json.loads(Path(path).read_text(encoding='utf-8'))
    if (saved.get('schema') != 'packet8-wider24-whole-code-outward-1'
            or not saved.get('whole_code_certificate')
            or not saved.get('source_pins_verified_at_finish')):
        raise ValueError('completed whole-code certificate required')
    pins = dict(saved['source_sha256'])
    pins[str(Path(path).resolve())] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    for name, digest in pins.items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ArithmeticError(f'pin mismatch: {name}')
    old = {int(q):decode(item) for q, item in saved['selected_endpoints'].items()}
    if set(old) != set(range(1, 257)) or sum(old.values(), F(0)) != decode(saved['exact_union']):
        raise ArithmeticError('complete exact original union required')
    trials = {t['macro_receipt']:F(t['alpha']) for t in saved['trials']}
    changed = {}
    factors = {}
    for q in range(1, 257):
        alpha = F(1) if q == 1 else trials[saved['selected_witnesses'][str(q)]]
        factor = sp.integer_power(gamma, symbols*q)
        factor = sp.fractional_power(factor, alpha)
        factor = sp.scalar_fraction(factor)
        endpoint = min(F(1), old[q]*factor)
        changed[q], factors[q] = endpoint, factor
    total = sum(changed.values(), F(0))
    for p in (Path(__file__).resolve(), HERE/'test_structured_randomizers.py',
              Path(sp.__file__).resolve(), Path(op.__file__).resolve()):
        pins[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return dict(schema='packet8-wider24-mds-sandwich-transport-1',
        original_receipt=str(Path(path).resolve()), gamma=str(gamma),
        gamma_numerator=gamma.numerator,gamma_denominator=gamma.denominator,
        independent_randomized_symbols_per_group=symbols,
        proof_premise='each fixed nonzero symbol image dominated by gamma times uniform nonzero; independent symbols/groups',
        q1_correction='gamma^symbols', tail_correction='gamma^(symbols*q*alpha)',
        inner_and_routing_unchanged=True, source_sha256=pins,
        source_pins_verified_at_finish=True, numerical_correction='positive outward scaled arithmetic',
        endpoints={str(q):encode(v) for q,v in changed.items()},
        factors={str(q):encode(v) for q,v in factors.items()},
        exact_union=encode(total), exact_union_at_most_2_minus40=total <= F(1,1 << 40),
        combined_margin_bits_display=margin(total),
        weakest_q=max(changed, key=changed.get), q1_margin_bits_display=margin(changed[1]),
        whole_code_certificate_under_stated_randomizer_premise=total <= F(1,1 << 40))


def main(args):
    if args.output.exists():
        raise ValueError('fresh output path required')
    aes = [[2,3,1,1],[1,2,3,1],[1,1,2,3],[3,1,1,2]]
    passed, details = mds_minors(aes, 0x13)
    if not passed:
        raise ArithmeticError('proposed GF16 AES matrix is not MDS')
    result = correct_receipt(args.receipt, F((1 << 32)-1,255**4))
    result['H_GF16_polynomial'] = 'x^4+x+1'
    result['H_rows'] = aes
    result['H_all_minors_checked'] = details['checked']
    result['diagonal_GF256_polynomial'] = 't^2+t+8 over GF16(0x13), low nibble constant coefficient'
    result['rank_two_search'] = rank_two_search()
    for name,digest in result['source_sha256'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ArithmeticError('source changed during transport')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('gamma','combined_margin_bits_display','q1_margin_bits_display',
        'weakest_q','exact_union_at_most_2_minus40','H_all_minors_checked','rank_two_search')},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,default=HERE.parent/'iteration5'/'whole_wider24_v1.json')
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())
