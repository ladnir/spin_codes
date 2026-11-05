"""Directed selected-weight verifier for full-RS[width,state] candidates.

Public widths are 64 and 128; 1 <= state <= width/2. The fixed maps and
exact rational packet caps come from the frozen wide provider. The outer
Toeplitz setup, 4096 pair slots per region, 256 regions, and persistent
state remain unchanged. Width 128 uses 64 pair slots per packet, 64 packets
per region, and 32 composed packet pairs per region.

Every local cap is rounded upward to a dyadic mantissa before Arb work.
Good-good raw matrices are multiplied before entrywise alpha; arrangement
counts stay outside alpha. Frozen positive carry-free polynomial helpers
retain every requested histogram coefficient. Each receipt covers every
support and feasible even/odd profile at one selected weight, at most 64.
It is not a whole-code certificate or an encoder performance claim.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
import sys
from time import monotonic

from flint import arb, arb_poly, ctx

import hd_field128_wide_rs_state_screen as local
import hd_field128_lift_certificate as base


FIELD_ORDER = base.FIELD_ORDER
WIDTHS = local.WIDTHS
DEFAULT_WIDTH, DEFAULT_STATE = 128, 48
MAXIMUM = 64


def packet_table(z, maximum, *, q=FIELD_ORDER, width=DEFAULT_WIDTH, state=DEFAULT_STATE,
                 bits=80, joint=True, budget=None):
    """At most 65 relevant occupancy caps; exact rational arithmetic only."""
    z = base.exact(z)
    local.construction(width, state)
    base.checked_q(q)
    if (type(maximum) is not int or not 0 <= maximum <= MAXIMUM or not 0 < z < 1
            or type(bits) is not int or not 32 <= bits <= 512 or type(joint) is not bool):
        raise ValueError('Invalid bounded packet scope, exact mark, or mantissa budget')
    budget = budget or base.Budget()
    table = {}
    for occupancy in range(maximum + 1):
        budget.check()
        exact = local.exact_packet(occupancy, z, q, n=width, state=state, joint=joint)
        table[occupancy] = tuple(base.rational(base.round_up(value, bits)) for value in exact)
    budget.check()
    return table


def directed_parts(z, alpha, maximum, q, width, state, bits, budget, paired=True, joint=True):
    """Swap-averaged single coefficients and optional raw-GG fractional caps."""
    z, alpha = base.exact(z), base.exact(alpha)
    if not 0 < alpha <= 1 or type(paired) is not bool:
        raise ValueError('Expected exact 0<alpha<=1 and a boolean paired flag')
    table = packet_table(z, maximum, q=q, width=width, state=state, bits=bits, joint=joint, budget=budget)
    exponent = base.rational(alpha)
    powered = {j: tuple(base.upper(value**exponent) if not value.is_zero() else arb(0) for value in values)
               for j, values in table.items()}
    # Fair singleton swaps cancel 2**s. The unpowered multiplicity is
    # multinomial(width/2;s,f), not the old width-64 packet multiplicity.
    single = {(s, f): tuple(count * value for value in powered[s + 2*f])
              for s, f, count in local.geometry(maximum, width=width)}
    if not paired or alpha == 1:
        return single, None
    good = local.geometry(maximum, width=width, good=True)
    composed = {}
    for s, f, count in good:
        budget.check()
        for t, g, other in good:
            if s + t + 2*(f + g) > maximum:
                continue
            values = composed.setdefault((s + t, f + g), [arb(0) for _ in range(4)])
            # Chronological raw product sums both intermediate-state paths
            # before alpha. Arrangement multiplicities remain outside alpha.
            raw = base.matrix_product(table[s + 2*f], table[t + 2*g])
            for entry, value in enumerate(raw):
                if not value.is_zero():
                    values[entry] += count * other * value**exponent
    budget.check()
    return single, composed


def regional_directed(z, alpha, maximum, *, q=FIELD_ORDER, width=DEFAULT_WIDTH,
                      state=DEFAULT_STATE, bits=80, paired=True, joint=True, budget=None):
    """Normalized regional matrices for every degree r<=maximum<=64."""
    z, alpha = base.exact(z), base.exact(alpha)
    local.construction(width, state)
    base.checked_q(q)
    if (type(maximum) is not int or not 1 <= maximum <= MAXIMUM
            or not 0 < z < 1 or not 0 < alpha <= 1 or type(paired) is not bool
            or type(joint) is not bool or type(bits) is not int or not 32 <= bits <= 512):
        raise ValueError('Invalid selected-weight cutoff, exact witness, or packet parameters')
    budget = budget or base.Budget()
    route = local.routing(width)
    marker = F(maximum, base.PAIRS)
    single, composed = directed_parts(z, alpha, maximum, q, width, state, bits, budget, paired, joint)
    if composed is None:
        run = base.polynomial_power(base.encode(single, maximum, marker),
                                    route['packets_per_region'], maximum, budget)
    else:
        good = base.encode({key: values for key, values in single.items() if key[0] <= 1}, maximum, marker)
        dirty = base.encode({key: values for key, values in single.items() if key[0] >= 2}, maximum, marker)
        parts = (base.encode(composed, maximum, marker),
                 base.polynomial_product(good, dirty, maximum, budget),
                 base.polynomial_product(dirty, good, maximum, budget),
                 base.polynomial_product(dirty, dirty, maximum, budget))
        operator = tuple(sum((part[i] for part in parts), arb_poly()) for i in range(4))
        run = base.polynomial_power(operator, route['composed_pairs_per_region'], maximum, budget)
    stride, result = 2*maximum + 1, {}
    for s, f in base.valid_degrees(maximum):
        normalizer = base.comb(base.PAIRS, s) * base.comb(base.PAIRS - s, f) * base.rational(marker)**(s + f)
        result[s + 2*f, f] = tuple(base.upper(run[i][s + stride*f] / normalizer) for i in range(4))
    budget.check()
    return result


def source_pins():
    """Complete recursive static local-import closure, including lazy imports."""
    directory = Path(__file__).resolve().parent
    pending, seen = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory / name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0] + '.py'
                if (directory / target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def weight_receipt(r, z, alpha, *, q=FIELD_ORDER, width=DEFAULT_WIDTH, state=DEFAULT_STATE,
                   distance=F(12, 25), precision=256, bits=80, paired=True, joint=True,
                   outer_measure='invertible-toeplitz', failure_bits=40, max_seconds=60., _pins=None):
    """Bound every support and feasible matching profile at selected weight r."""
    z, alpha, distance = base.exact(z), base.exact(alpha), base.exact(distance)
    details, field_bits = local.metadata(width, state), base.checked_q(q)
    if (type(r) is not int or not 1 <= r <= MAXIMUM or not 0 < distance < F(1, 2)
            or type(precision) is not int or not 128 <= precision <= 4096
            or type(failure_bits) is not int or not 0 <= failure_bits <= 4096):
        raise ValueError('Invalid selected-weight scope, target, or precision')
    pins = source_pins() if _pins is None else _pins
    budget = base.Budget(max_seconds)
    forbidden = int(distance * base.OUTPUT_LENGTH)
    with ctx.workprec(precision):
        regions = regional_directed(z, alpha, r, q=q, width=width, state=state, bits=bits,
                                    paired=paired, joint=joint, budget=budget)
        union, profiles = arb(0), 0
        for even in range(r//2 + 1):
            budget.check()
            for odd in range(r//2 + 1):
                count = base.profile_count_upper(base.PAIRS, r, even, odd)
                if count:
                    run = base.matrix_power(base.matrix_product(regions[r, even], regions[r, odd]), base.REGIONS // 2)
                    union += count * (run[0] + run[1])
                    profiles += 1
        log2arb = lambda value: value.log() / arb(2).log()
        prefix = base.rational(alpha) * (r * log2arb(base.rational(base.outer_density(q, outer_measure)))
                                        - log2arb(arb(q - 1)) - forbidden * log2arb(base.rational(z)))
        value = prefix + log2arb(union)
        ceiling, enclosure = base.rational_ceiling(value), str(value)
    budget.check()
    if source_pins() != pins:
        raise RuntimeError('A pinned transitive helper changed during the replay')
    return dict(**details, fixed_map_field_order=256,
                field_order=q, field_bits=field_bits, active_groups=r, profiles=profiles,
                covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
                forbidden_weight=forbidden, target_distance=str(distance),
                outer_measure=outer_measure, projective_divisor=str(q - 1), joint_rank_return_cap=joint,
                log2_upper_rational=str(ceiling), log2_upper_enclosure=enclosure,
                requested_failure_bits=failure_bits, meets_requested_selected_weight_target=ceiling <= -failure_bits,
                meets_40_bit_selected_weight_target=ceiling <= -40,
                witness=dict(z=str(z), alpha=str(alpha), precision=precision, packet_mantissa_bits=bits,
                             paired=paired, joint_rank_return_cap=joint,
                             coefficient_maximum=r, coefficient_marker=str(F(r, base.PAIRS))),
                coefficient_method='wide-full-RS-state-exact-packets-positive-carry-free-Arb-truncated-bivariate',
                local_caps='exact rational generic MDS rank/fiber return bound; upward dyadic packet rounding',
                composition_steps=2 if paired and alpha != 1 else 1,
                source_sha256=pins, pin_scope='complete recursive static local-import closure, including lazy imports',
                runtime=dict(python=sys.version, python_flint=version('python-flint'), numpy=version('numpy')),
                elapsed_seconds=monotonic() - budget.started, max_seconds=max_seconds,
                arithmetic='directed Arb with exact rational local caps', directed=True, whole_code_certificate=False,
                scope='new width/state construction; every support at this selected weight only; no encoder or timing claim')


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--z', required=True, help='Exact rational or decimal witness')
    parser.add_argument('--alpha', required=True, help='Exact rational or decimal witness')
    parser.add_argument('--width', type=int, choices=WIDTHS, default=DEFAULT_WIDTH)
    parser.add_argument('--state', type=int, required=True)
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--distance', required=True, help='Exact target distance for this construction')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--failure-bits', type=int, default=40)
    parser.add_argument('--max-seconds', type=float, default=60.)
    parser.add_argument('--unpaired', action='store_true')
    parser.add_argument('--simple-return', action='store_true')
    parser.add_argument('--outer-measure', choices=base.OUTER_MEASURES, default='invertible-toeplitz')
    parser.add_argument('--output', type=Path, required=True, help='New output path; exclusive create')
    args = parser.parse_args()
    if not 8 <= args.field_bits <= 128 or args.field_bits % 8:
        parser.error('Field bits must be a multiple of 8 in 8..128')
    if args.output.exists():
        parser.error('Output already exists')
    row = weight_receipt(args.r, F(args.z), F(args.alpha), q=1 << args.field_bits,
                         width=args.width, state=args.state, distance=F(args.distance),
                         precision=args.precision, bits=args.bits, paired=not args.unpaired,
                         joint=not args.simple_return, outer_measure=args.outer_measure,
                         failure_bits=args.failure_bits, max_seconds=args.max_seconds)
    with args.output.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(row, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: row[key] for key in ('construction', 'active_groups', 'target_distance',
                                               'log2_upper_rational', 'profiles', 'directed',
                                               'meets_requested_selected_weight_target', 'whole_code_certificate')}))


if __name__ == '__main__':
    main()
