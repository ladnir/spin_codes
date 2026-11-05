"""Fresh directed terminal suffix for the fixed-gamma full-RS construction.

Fix all route and inner setup independently of the outer. Identity
feedthrough recovers X=Y-AQ, then advances Q with the known gamma. Thus
every target has one inverse image, regardless of the new local caps.
The frozen inverse-image formula is freshly evaluated, not copied from
a receipt for another state or construction.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import hd_field128_full_rs_gamma_certificate as gamma
import hd_field128_lift_tail_certificate as terminal


base, local = gamma.base, gamma.local


def source_pins():
    pins = dict(gamma.source_pins(), **terminal.source_pins())
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def terminal_receipt(*, state=30, distance=F(2, 5), first=None,
                     failure_bits=41, precision=256):
    name, distance = local.construction(state), base.exact(distance)
    pins, automatic = source_pins(), first is None
    result = terminal.terminal_receipt(q=base.FIELD_ORDER, distance=distance,
                                      first=first, failure_bits=failure_bits,
                                      precision=precision)
    while automatic and not result['meets_requested_range_target']:
        first = result['weight_range'][0]+1
        if first > base.GROUPS:
            raise ArithmeticError('No directed terminal suffix meets the target')
        result = terminal.terminal_receipt(q=base.FIELD_ORDER, distance=distance,
                                          first=first, failure_bits=failure_bits,
                                          precision=precision)
    if source_pins() != pins:
        raise RuntimeError('A pinned helper changed during terminal replay')
    return dict(result, construction=name, state_dimension=state, inner_width=64,
                inner_emission=f'single-RS[64,{state}]', emission_minimum=65-state,
                feedback='gamma A transpose', gamma_extension_degree=16,
                gamma_base_field_order=256, gamma_choice='fixed degree16 element over GF256',
                CA_zero=True, fixed_map_field_order=256,
                fixed_nodes='GF256 additive subspace0..63', state_nodes=f'128..{127+state}',
                refresh=f'D2 H_{state} D1 with independent full-field nonzero diagonals',
                initial_state='zero', final_state='unrestricted',
                pair_swaps=True, state_persists_across_regions=True,
                full_field_random_coefficients=True, outer_independent_of_route_and_inner=True,
                outer_groups=base.GROUPS, outer_randomizers_per_group=64,
                outer_randomizer='independent uniform invertible4x4 Toeplitz per group and column',
                outer_point_counting_cap='q^(-96)*(q-1)^(-32) per nonzero outer group',
                causal_inverse='Recover X=Y-AQ, then Q_next=D2 H_state D1 (Q+gamma C X)',
                inner_formula_dependency='causal identity feedthrough, not the state-specific local caps',
                automatic_terminal_cutoff=automatic,
                covers_all_weights_in_range_simultaneously=True, no_union_over_tail_weights=True,
                old_numerical_bounds_reused=False, source_sha256=pins,
                arithmetic='fresh directed Arb evaluation with exact rational Chernoff mark',
                scope='one event covering this terminal suffix only; not a whole-code certificate')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', type=int, choices=local.STATES, default=30)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--first', type=int)
    parser.add_argument('--failure-bits', type=int, default=41)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new file')
    result = terminal_receipt(state=args.state, distance=F(args.distance), first=args.first,
                              failure_bits=args.failure_bits, precision=args.precision)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('construction', 'weight_range',
                     'log2_upper_rational', 'meets_requested_range_target')}))


if __name__ == '__main__':
    main()
