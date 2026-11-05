"""Fresh directed inverse-image terminal receipt for a full-RS state code.

For fixed setup, identity feedthrough gives a unique inverse: recover each
X=Y-AQ, then advance the known state recurrence. The inner map therefore
preserves the target-counting argument used by the frozen terminal helper.
The outer counting law and dimensions are unchanged. Each call recomputes
the directed ball and tail bound; no earlier numerical receipt is reused.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import hd_field128_full_rs_state_certificate as certificate
import hd_field128_lift_tail_certificate as terminal


base, local = certificate.base, certificate.local


def source_pins():
    pins = dict(certificate.source_pins(), **terminal.source_pins())
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def terminal_receipt(*, q=base.FIELD_ORDER, state=24, distance=F(2, 5), first=None,
                     failure_bits=40, precision=256):
    name = local.construction(state)
    distance = base.exact(distance)
    pins = source_pins()
    automatic = first is None
    # Fresh directed evaluation of the construction-independent inverse-image formula.
    result = terminal.terminal_receipt(q=q, distance=distance, first=first,
                                       failure_bits=failure_bits, precision=precision)
    # Exact-threshold rounding can conservatively cost one boundary weight.
    while automatic and not result['meets_requested_range_target']:
        next_first = result['weight_range'][0]+1
        if next_first > base.GROUPS:
            raise ArithmeticError('No terminal suffix meets the requested directed target')
        result = terminal.terminal_receipt(q=q, distance=distance, first=next_first,
                                           failure_bits=failure_bits, precision=precision)
    if source_pins() != pins:
        raise RuntimeError('A pinned helper changed during the terminal replay')
    return dict(result, construction=name, state_dimension=state, inner_width=64,
                inner_emission=f'single-RS[64,{state}]', emission_minimum=65-state,
                feedback='A transpose', CA_zero=True, fixed_map_field_order=256,
                fixed_nodes='GF256 additive subspace0..63', state_nodes=f'128..{127+state}',
                refresh=f'D2 H_{state} D1 with independent full-field nonzero diagonals',
                pair_swaps=True, state_persists_across_regions=True, full_field_random_coefficients=True,
                outer_groups=base.GROUPS,
                outer_randomizer='independent full-field invertible4x4 Toeplitz per outer group and output column',
                outer_randomizers_per_group=64, automatic_terminal_cutoff=automatic,
                covers_all_weights_in_range_simultaneously=True, no_union_over_tail_weights=True,
                causal_inverse='For fixed setup, recover X=Y-AQ and advance Q_next=D2 H_state D1 (Q+C X)',
                inner_formula_dependency='identity feedthrough and setup independence; no state-specific local cap',
                outer_point_counting_cap='q^(-96)*(q-1)^(-32) per nonzero outer group',
                old_numerical_bounds_reused=False,
                source_sha256=pins, pin_scope='complete recursive static local-import closure, including lazy imports',
                arithmetic='directed Arb with exact rational Chernoff mark',
                scope='one inverse-image event for every source weight in this terminal suffix; not a whole-code certificate')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', type=int, choices=local.STATES, default=24)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--first', type=int)
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--failure-bits', type=int, default=40)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new file')
    if not 8 <= args.field_bits <= 128 or args.field_bits%8:
        parser.error('Field bits must be a multiple of8in8..128')
    result = terminal_receipt(q=1 << args.field_bits, state=args.state, distance=F(args.distance),
                              first=args.first, failure_bits=args.failure_bits, precision=args.precision)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('construction', 'weight_range', 'target_distance',
                                                  'log2_upper_rational', 'meets_requested_range_target')}, allow_nan=False))


if __name__ == '__main__':
    main()
