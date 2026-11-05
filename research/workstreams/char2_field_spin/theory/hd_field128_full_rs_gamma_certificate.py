"""Directed selected-weight verifier for the new gamma-feedback construction.

This verifier requires GF(2^128), width64, state16/24/28/30, and a fixed
gamma of degree16 over GF256. Exact local caps use deterministic principal
nonsingularity and the joint active/inactive return count. Every requested
coefficient and even/odd matching profile is retained at one weight r<=64.
The receipt is not a whole-code certificate or an implementation timing.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
import sys
from time import monotonic

from flint import arb, arb_poly, ctx

import hd_field128_full_rs_gamma_screen as local
import hd_field128_lift_certificate as base


STATES, FIELD_ORDER, MAXIMUM = local.STATES, base.FIELD_ORDER, 64


def packet_table(z, maximum, *, state=24, bits=80, budget=None):
    z = base.exact(z)
    local.construction(state)
    if (type(maximum) is not int or not 0 <= maximum <= MAXIMUM or not 0 < z < 1
            or type(bits) is not int or not 32 <= bits <= 512):
        raise ValueError('Invalid bounded packet scope, exact mark, or mantissa budget')
    budget = budget or base.Budget()
    table = {}
    for j in range(maximum+1):
        budget.check()
        table[j] = tuple(base.rational(base.round_up(value, bits))
                         for value in local.exact_packet(j, z, FIELD_ORDER, 64, state))
    budget.check()
    return table


def directed_parts(z, alpha, maximum, state, bits, budget, paired=True):
    z, alpha = base.exact(z), base.exact(alpha)
    if not 0 < alpha <= 1 or type(paired) is not bool:
        raise ValueError('Expected exact0<alpha<=1 and boolean paired flag')
    table = packet_table(z, maximum, state=state, bits=bits, budget=budget)
    exponent = base.rational(alpha)
    powered = {j: tuple(base.upper(value**exponent) if not value.is_zero() else arb(0) for value in row)
               for j, row in table.items()}
    geometry = local.frozen.fixed16.geometry
    single = {(s, f): tuple(count*value for value in powered[s+2*f])
              for s, f, count in geometry(maximum)}
    if not paired or alpha == 1:
        return single, None
    good, composed = geometry(maximum, good=True), {}
    for s, f, count in good:
        budget.check()
        for t, g, other in good:
            if s+t+2*(f+g) > maximum:
                continue
            target = composed.setdefault((s+t, f+g), [arb(0) for unused in range(4)])
            # Raw chronological multiplication precedes alpha; counts remain outside it.
            for entry, value in enumerate(base.matrix_product(table[s+2*f], table[t+2*g])):
                if not value.is_zero():
                    target[entry] += count*other*value**exponent
    budget.check()
    return single, composed


def regional_directed(z, alpha, maximum, *, state=24, bits=80, paired=True, budget=None):
    z, alpha = base.exact(z), base.exact(alpha)
    local.construction(state)
    if (type(maximum) is not int or not 1 <= maximum <= MAXIMUM or not 0 < z < 1
            or not 0 < alpha <= 1 or type(paired) is not bool):
        raise ValueError('Invalid bounded regional scope or exact witness')
    budget = budget or base.Budget()
    marker = F(maximum, base.PAIRS)
    single, composed = directed_parts(z, alpha, maximum, state, bits, budget, paired)
    if composed is None:
        run = base.polynomial_power(base.encode(single, maximum, marker), 128, maximum, budget)
    else:
        good = base.encode({key: row for key, row in single.items() if key[0] <= 1}, maximum, marker)
        dirty = base.encode({key: row for key, row in single.items() if key[0] >= 2}, maximum, marker)
        parts = (base.encode(composed, maximum, marker),
                 base.polynomial_product(good, dirty, maximum, budget),
                 base.polynomial_product(dirty, good, maximum, budget),
                 base.polynomial_product(dirty, dirty, maximum, budget))
        operator = tuple(sum((part[i] for part in parts), arb_poly()) for i in range(4))
        run = base.polynomial_power(operator, 64, maximum, budget)
    stride, regions = 2*maximum+1, {}
    for s, f in base.valid_degrees(maximum):
        normalizer = base.comb(base.PAIRS, s)*base.comb(base.PAIRS-s, f)*base.rational(marker)**(s+f)
        regions[s+2*f, f] = tuple(base.upper(run[i][s+stride*f]/normalizer) for i in range(4))
    budget.check()
    return regions


def source_pins():
    directory, pending, seen = Path(__file__).resolve().parent, [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory/name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def weight_receipt(r, z, alpha, *, state=24, distance=F(2, 5), precision=256, bits=80,
                   paired=True, failure_bits=40, max_seconds=60., _pins=None):
    z, alpha, distance = base.exact(z), base.exact(alpha), base.exact(distance)
    name = local.construction(state)
    if (type(r) is not int or not 1 <= r <= MAXIMUM or not 0 < distance < F(1, 2)
            or type(precision) is not int or not 128 <= precision <= 4096
            or type(failure_bits) is not int or not 0 <= failure_bits <= 4096):
        raise ValueError('Invalid selected weight, target, or precision')
    pins = source_pins() if _pins is None else _pins
    budget = base.Budget(max_seconds)
    forbidden = int(distance*base.OUTPUT_LENGTH)
    with ctx.workprec(precision):
        regions = regional_directed(z, alpha, r, state=state, bits=bits, paired=paired, budget=budget)
        union, profiles = arb(0), 0
        for even in range(r//2+1):
            budget.check()
            for odd in range(r//2+1):
                count = base.profile_count_upper(base.PAIRS, r, even, odd)
                if count:
                    product = base.matrix_power(base.matrix_product(regions[r, even], regions[r, odd]), base.REGIONS//2)
                    union += count*(product[0]+product[1])
                    profiles += 1
        log2arb = lambda value: value.log()/arb(2).log()
        prefix = base.rational(alpha)*(r*log2arb(base.rational(base.outer_density(FIELD_ORDER)))
                                        -log2arb(arb(FIELD_ORDER-1))-forbidden*log2arb(base.rational(z)))
        value = prefix+log2arb(union)
        ceiling, enclosure = base.rational_ceiling(value), str(value)
    budget.check()
    if source_pins() != pins:
        raise RuntimeError('A pinned transitive helper changed during replay')
    return dict(construction=name, inner_width=64, state_dimension=state, inner_emission=f'single-RS[64,{state}]',
                emission_minimum=65-state, feedback='gamma A transpose', CA_zero=True,
                gamma_extension_degree=16, gamma_base_field_order=256, gamma_choice='any fixed element of degree16 over GF256',
                principal_nonsingularity='deterministic symmetric unit-diagonal determinant degree-halving lemma',
                fixed_map_field_order=256, fixed_nodes='GF256 additive subspace0..63', state_nodes=f'128..{127+state}',
                refresh=f'D2 H_{state} D1; fixed MDS Lagrange from0..{state-1} to{state}..{2*state-1}',
                field_order=FIELD_ORDER, field_bits=128, active_groups=r, profiles=profiles,
                outer_groups=base.GROUPS, source_dimension=base.DIMENSION, output_length=base.OUTPUT_LENGTH,
                forbidden_weight=forbidden, target_distance=str(distance), outer_measure='invertible-toeplitz',
                projective_divisor=str(FIELD_ORDER-1), pair_swaps=True, state_persists_across_regions=True,
                full_field_random_coefficients=True, covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
                log2_upper_rational=str(ceiling), log2_upper_enclosure=enclosure, requested_failure_bits=failure_bits,
                meets_requested_selected_weight_target=ceiling <= -failure_bits,
                witness=dict(z=str(z), alpha=str(alpha), precision=precision, packet_mantissa_bits=bits,
                             paired=paired, coefficient_maximum=r, coefficient_marker=str(F(r, base.PAIRS))),
                coefficient_method='gamma-full-RS-exact-packets-positive-carry-free-Arb-truncated-bivariate',
                local_caps='exact rational joint active/inactive return count; upward dyadic packet rounding',
                composition_steps=2 if paired and alpha != 1 else 1,
                source_sha256=pins, pin_scope='complete recursive static local-import closure, including lazy imports',
                runtime=dict(python=sys.version, python_flint=version('python-flint'), numpy=version('numpy')),
                elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                arithmetic='directed Arb with exact rational local caps', directed=True, whole_code_certificate=False,
                scope='new gamma construction; selected weight only; no encoder or timing claim')


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', required=True, type=int)
    parser.add_argument('--z', required=True)
    parser.add_argument('--alpha', required=True)
    parser.add_argument('--state', type=int, choices=STATES, default=24)
    parser.add_argument('--distance', required=True)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--max-seconds', type=float, default=60.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output already exists')
    result = weight_receipt(args.r, F(args.z), F(args.alpha), state=args.state, distance=F(args.distance),
                            precision=args.precision, bits=args.bits, max_seconds=args.max_seconds)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: result[key] for key in ('construction', 'active_groups', 'target_distance',
                                                  'log2_upper_rational', 'directed', 'whole_code_certificate')}))


if __name__ == '__main__':
    main()
