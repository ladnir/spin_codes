"""Directed middle and dense certificates for repaired cubic width8/state4.

Frozen K20 GRS[128,83], one regional route, independent nonzero input and
state diagonals. Only explicitly stated integer weight classes are covered.
The middle replay uses exact rational packet matrices before outward
decimal rounding. No optimizer or floating objective decides acceptance.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb

from cauchy_firstmoment_screen import geometry
from gf8_outer_rate_middle import conditioning
from gf8_packet_width_screen import local_matrix
from return_cover_certificate import Witness, interval, row_upper, terminal_power
from state8_density_certificate import integer_span, refined_cells, retarget_witness
from subfield_certificate import HIGH, LOW, ZERO, integer_margin


LOW.Emax = HIGH.Emax = 1_000_000_000
LOW.Emin = HIGH.Emin = -1_000_000_000
CONSTRUCTION = 'gf256-repaired-cubic-w8-s4-eta2-c3-regional-nonzero-diagonals'
RETURN_MINIMA = (8, 6, 5, 4, 3, 2, 2, 1)
WITNESSES = tuple((F(density), Witness(*(F(v) for v in values)))
    for density, values in (
        ('3/100', ('2342/7955', '6508/8863', '82937/9552', '3289/809')),
        ('1/25', ('3659/8876', '6478/9805', '3485/384', '20186/6759')),
        ('3/50', ('250/379', '3792/7001', '65614/6563', '18218/9607')),
        ('2/25', ('4091/4474', '2215/4898', '37077/3358', '8237/6112')),
        ('1/10', ('11317/9540', '3447/9064', '95784/7853', '6871/6748')),
        ('3/25', ('14104/9429', '755/2376', '48293/3579', '7425/9296')),
        ('7/50', ('11467/6112', '2570/9901', '125249/8359', '5916/9217')),
        ('4/25', ('55495/4292', '112/3861', '114211/2997', '3433/9730')),
        ('9/50', ('323063/8704', '37/4067', '507207/6650', '505/2384')),
        ('1/5', ('911509/9790', '14/3943', '1483345/9011', '1080/8959')),
        ('6/25', ('283572/1711', '9/5078', '255', '617/8331'))))


def metadata():
    g = geometry(20, 128, 83)
    return dict(construction=CONSTRUCTION, outer=[128, 83], q=256, width=8,
                state=4, k_log=20, requested_k=2**20, padded_k=g['k'],
                rows=g['rows'], n=g['n'], output_threshold=g['threshold'],
                eta=2, injection_constant=3, return_minima=list(RETURN_MINIMA),
                distance='1/10', nominal_rate='83/128',
                requested_rate=str(F(2**20, g['n'])))


def exact_transfer(witness):
    """Aggregate the proved local covers without floating arithmetic."""
    if any(not isinstance(getattr(witness, name), F) for name in ('x', 'y', 'tau', 'tilt')):
        raise ValueError('exact rational packet witnesses required')
    if not witness.x > 0 or not 0 < witness.y <= 1 or not 1 <= witness.tau <= 255 or not witness.tilt > 0:
        raise ValueError('invalid exact witness')
    result = [[F(0), F(0)], [F(0), F(0)]]
    for w in range(9):
        matrix = local_matrix(w, witness.y, witness.tau, 8, 4, RETURN_MINIMA)
        factor = comb(8, w)*witness.x**w
        for i in range(2):
            for j in range(2):
                result[i][j] += factor*matrix[i][j]
    return result


def prefactor(witness, g):
    matrix = [[interval(x) for x in row] for row in exact_transfer(witness)]
    return (interval(conditioning(g['rows']))**128*terminal_power(matrix, g['n']//8)
            /(interval(witness.tau)*interval(witness.y)**g['output_threshold']))


def cell_upper(lo, hi, g, witness, factor):
    if not 0 < lo <= hi < g['n']:
        raise ValueError('invalid integer middle cell')
    source = row_upper(F(lo, g['n']), F(hi, g['n']), witness, length=128, parity=45)
    base = witness.x*witness.tilt
    exponent = lo if base >= 1 else hi
    return (hi-lo+1)*factor*source**g['rows']/interval(base)**exponent


def certify(first=F(3, 100), last=F(6, 25), initial_grid=1000,
            cell_bits=128, witnesses=WITNESSES):
    if initial_grid < 1 or cell_bits < 40 or not witnesses:
        raise ValueError('invalid middle replay parameters')
    if any(not isinstance(anchor, F) or not 0 < anchor < 1 for anchor, _ in witnesses):
        raise ValueError('exact rational witness anchors required')
    g = metadata()
    first_weight, last_weight = integer_span(g['n'], first, last)
    factors = tuple(prefactor(w, g) for _, w in witnesses)

    def evaluate(lo, hi):
        middle = F(lo+hi, 2*g['n'])
        candidates = sorted(range(len(witnesses)),
                            key=lambda i: abs(witnesses[i][0]-middle))[:3]
        return min(((i, cell_upper(lo, hi, g,
                                  retarget_witness(*witnesses[i], middle), factors[i]))
                    for i in candidates), key=lambda item: item[1].hi)

    maximum = interval(F(1, 2**cell_bits)).lo
    span = (g['n']+initial_grid-1)//initial_grid
    next_weight, cells, total, largest = first_weight, 0, ZERO, ZERO.hi
    digest, uses = sha256(), [0]*len(witnesses)
    for lo, hi, index, value in refined_cells(first_weight, last_weight, span, evaluate, maximum):
        if lo != next_weight or hi < lo:
            raise ArithmeticError('middle coverage is not contiguous')
        next_weight = hi+1
        cells += 1
        uses[index] += 1
        total += value
        largest = max(largest, value.hi)
        digest.update(f'{lo}:{hi}:{index}\n'.encode('ascii'))
    if next_weight != last_weight+1:
        raise ArithmeticError('middle coverage ends early')
    bits = integer_margin(total.hi)
    if bits < 40:
        raise ArithmeticError(f'middle component below40 bits: {bits}')
    return dict(**g, weights=[first_weight, last_weight],
                covered_weights=last_weight-first_weight+1, checked_cells=cells,
                density_min=str(first), density_max_exclusive=str(last),
                failure_upper=str(total.hi), margin_bits_lower=bits,
                class_margin_bits_lower=bits,
                minimum_cell_margin_bits_lower=integer_margin(largest),
                partition_sha256=digest.hexdigest(), witness_cell_counts=uses,
                conditioning='central-mode-reciprocal^128', precision_digits=80,
                complete=True, rigorous=True, whole_code=False,
                scope='complete middle inverse-input-weight interval only')


def certify_dense(first=F(6, 25), cell_bits=128):
    """Apply the inner-independent dense theorem with fresh packet metadata.

    The theorem uses the invertible causal feedthrough map and independent
    input diagonals. It depends on outer geometry, not packet width/state.
    Reuse its directed arithmetic, never its old construction metadata.
    """
    from gf8_outer_rate_dense import certify as dense_arithmetic
    g = metadata()
    old = dense_arithmetic(83, 20, first, cell_bits=cell_bits)
    for key in ('outer', 'q', 'k_log', 'requested_k', 'padded_k', 'rows', 'n', 'output_threshold'):
        if old[key] != g[key]:
            raise ArithmeticError(f'dense arithmetic geometry mismatch: {key}')
    keys = ('weights', 'first_density', 'covered_weights', 'checked_cells',
            'precision_bits', 'initial_grid', 'target_cell_bits', 'ball_exponent',
            'conditioning', 'margin_bits_lower', 'failure_upper_dyadic_exponent',
            'witness_sha256')
    return dict(**g, **{key: old[key] for key in keys}, complete=True,
                rigorous=True, whole_code=False,
                scope='complete dense inverse-input-weight interval only',
                dense_inner_hypothesis='causal identity-feedthrough inner and independent nonzero input diagonals')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dense', action='store_true')
    parser.add_argument('--cell-bits', type=int, default=128)
    args = parser.parse_args()
    result = certify_dense(cell_bits=args.cell_bits) if args.dense else certify(cell_bits=args.cell_bits)
    print(json.dumps(result), flush=True)
