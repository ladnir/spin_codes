"""Independent exact-rational check of the saved certificate corrections.

This verifies authentication, every correction factor and the final union.
It does not rerun the frozen local census or establish the randomizer lemma.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json


def dyadic(record):
    return F(int(record['numerator_hex'], 16)) * F(2)**record['exponent']


def authenticate(record):
    for name, expected in record['source_sha256'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
            raise ArithmeticError(f'source pin mismatch: {name}')


def verify(path):
    saved = json.loads(path.read_text(encoding='utf-8'))
    expected_schema = {
        'packet8-wider24-mds-sandwich-transport-1': 16,
        'packet8-wider24-one-identity-mds-sandwich-transport-1': 15,
    }
    symbols = expected_schema[saved['schema']]
    if saved['independent_randomized_symbols_per_group'] != symbols:
        raise ArithmeticError('wrong symbol count')
    gamma = F((1 << 32)-1, 255**4)
    if F(saved['gamma']) != gamma:
        raise ArithmeticError('wrong pointwise envelope')
    authenticate(saved)
    original = json.loads(Path(saved['original_receipt']).read_text(encoding='utf-8'))
    authenticate(original)
    if not original['whole_code_certificate']:
        raise ArithmeticError('original certificate incomplete')
    alphas = {t['macro_receipt']: F(t['alpha']) for t in original['trials']}
    keys = {str(q) for q in range(1, 257)}
    if (set(saved['endpoints']) != keys or set(saved['factors']) != keys
            or set(original['selected_endpoints']) != keys):
        raise ArithmeticError('incomplete occupancy range')
    total, old_total = F(0), F(0)
    for q in range(1, 257):
        key = str(q)
        alpha = F(1) if q == 1 else alphas[original['selected_witnesses'][key]]
        exponent = symbols*q*alpha
        factor = dyadic(saved['factors'][key])
        # Raise both sides to the positive integer denominator. All arithmetic
        # is rational: neither floating logs nor the producer's tangent routine.
        if factor**exponent.denominator < gamma**exponent.numerator:
            raise ArithmeticError(f'inward correction at q={q}')
        old = dyadic(original['selected_endpoints'][key])
        changed = dyadic(saved['endpoints'][key])
        if changed != min(F(1), old*factor):
            raise ArithmeticError(f'endpoint mismatch at q={q}')
        old_total += old
        total += changed
    if old_total != dyadic(original['exact_union']):
        raise ArithmeticError('original union mismatch')
    if total != dyadic(saved['exact_union']) or total > F(1, 1 << 40):
        raise ArithmeticError('corrected union mismatch or exceeds target')
    return {'receipt': str(path), 'endpoints': 256,
            'pins': len(saved['source_sha256']), 'exact_comparison_passes': True}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipts', type=Path, nargs='+')
    for receipt in parser.parse_args().receipts:
        print(json.dumps(verify(receipt)))
