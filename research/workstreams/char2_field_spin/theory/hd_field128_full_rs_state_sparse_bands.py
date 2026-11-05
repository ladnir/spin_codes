"""Directed full-RS state sparse bands from one exact regional down-set.

Every requested integer weight has its own outer prefix and complete
even/odd profile sum. A common coefficient cutoff and marker retain all
lower coefficients exactly; no interpolation between weights is used.
Completed rows are checkpointed before the next row, including when a
later table or profile calculation reaches its explicit resource limit.
"""

from fractions import Fraction as F
from hashlib import sha256
from importlib.metadata import version
import json
from pathlib import Path
import sys
from time import monotonic

from flint import arb, ctx

import hd_field128_full_rs_state_certificate as certificate


base, local = certificate.base, certificate.local


def source_pins():
    pins = certificate.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def write_new(path, value):
    raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
    with Path(path).open('xb') as stream:
        stream.write(raw)
    return dict(path=str(Path(path).resolve()), sha256=sha256(raw).hexdigest())


def _union_ceiling(rows, precision):
    if not rows:
        return None
    with ctx.workprec(precision):
        mass = sum((arb(2)**base.rational(F(row['log2_upper_rational'])) for row in rows), arb(0))
        return str(base.rational_ceiling(mass.log()/arb(2).log()))


def _weight_from_regions(r, regions, z, alpha, *, q, distance, outer_measure, budget):
    """Called within the same Arb precision context as the shared table."""
    budget.check()
    union, profiles = arb(0), 0
    for even in range(r//2+1):
        budget.check()
        for odd in range(r//2+1):
            count = base.profile_count_upper(base.PAIRS, r, even, odd)
            if count:
                run = base.matrix_power(base.matrix_product(regions[r, even], regions[r, odd]), base.REGIONS//2)
                union += count*(run[0]+run[1])
                profiles += 1
    log = lambda value: value.log()/arb(2).log()
    forbidden = int(distance*base.OUTPUT_LENGTH)
    prefix = base.rational(alpha)*(r*log(base.rational(base.outer_density(q, outer_measure)))
                                   -log(arb(q-1))-forbidden*log(base.rational(z)))
    value = prefix+log(union)
    result = dict(active_groups=r, profiles=profiles,
                  log2_upper_rational=str(base.rational_ceiling(value)), log2_upper_enclosure=str(value),
                  covers_all_supports_at_weight=True, covers_all_even_odd_splits=True)
    budget.check()
    return result


def batch_receipt(first, last, z, alpha, *, maximum=None, q=base.FIELD_ORDER, state=24,
                  distance=F(2, 5), precision=256, bits=80, paired=True, joint=True,
                  outer_measure='invertible-toeplitz', failure_bits=40, max_seconds=600.,
                  witness_source=None, checkpoint=None):
    z, alpha, distance = map(base.exact, (z, alpha, distance))
    name, field_bits = local.construction(state), base.checked_q(q)
    maximum = last if maximum is None else maximum
    if (any(type(value) is not int for value in (first, last, maximum, precision, bits, failure_bits))
            or not 1 <= first <= last <= maximum <= certificate.MAXIMUM or not 128 <= precision <= 4096
            or not 32 <= bits <= 512 or not 0 <= failure_bits <= 4096
            or not 0 < z < 1 or not 0 < alpha <= 1 or not 0 < distance < F(1, 2)
            or type(paired) is not bool or type(joint) is not bool or outer_measure not in base.OUTER_MEASURES):
        raise ValueError('Invalid exact sparse band, witness, construction, or precision')
    pins = source_pins()
    budget = base.Budget(max_seconds)
    witness = dict(z=str(z), alpha=str(alpha), precision=precision, packet_mantissa_bits=bits,
                   paired=paired, joint_rank_return_cap=joint, coefficient_maximum=maximum,
                   coefficient_marker=str(F(maximum, base.PAIRS)))
    metadata = dict(construction=name, state_dimension=state, inner_width=64, emission_minimum=65-state,
                    inner_emission=f'single-RS[64,{state}]', feedback='A transpose', CA_zero=True,
                    fixed_map_field_order=256, fixed_nodes='GF256 additive subspace0..63',
                    state_nodes=f'128..{127+state}',
                    refresh=f'D2 H_{state} D1 with independent full-field nonzero diagonals',
                    field_order=q, field_bits=field_bits, target_distance=str(distance),
                    forbidden_weight=int(distance*base.OUTPUT_LENGTH), source_dimension=base.DIMENSION,
                    output_length=base.OUTPUT_LENGTH, outer_groups=base.GROUPS,
                    outer_measure=outer_measure, projective_divisor=str(q-1), pair_swaps=True,
                    state_persists_across_regions=True, full_field_random_coefficients=True,
                    joint_rank_return_cap=joint, composition_steps=2 if paired and alpha != 1 else 1,
                    witness=witness, witness_source=witness_source, old_numerical_bounds_reused=False,
                    source_sha256=pins, pin_scope='complete recursive static local-import closure, including lazy imports',
                    coefficient_method='shared-full-RS-state-positive-carry-free-Arb-truncated-bivariate',
                    local_caps='exact rational generic MDS rank/fiber bound; upward dyadic packet rounding',
                    requested_failure_bits=failure_bits,
                    runtime=dict(python=sys.version, python_flint=version('python-flint')),
                    arithmetic='directed Arb with exact rational local caps', directed=True, whole_code_certificate=False)
    rows, failure, tables = [], None, 0
    phase, active = 'regional-table', None
    with ctx.workprec(precision):
        try:
            budget.check()
            regions = certificate.regional_directed(z, alpha, maximum, q=q, state=state, bits=bits,
                                                    paired=paired, joint=joint, budget=budget)
            tables = 1
            phase = 'weight'
            for active in range(first, last+1):
                result = _weight_from_regions(active, regions, z, alpha, q=q, distance=distance,
                                              outer_measure=outer_measure, budget=budget)
                if source_pins() != pins:
                    raise RuntimeError('A pinned helper changed during the sparse band replay')
                row = dict(metadata, **result,
                           meets_requested_selected_weight_target=F(result['log2_upper_rational']) <= -failure_bits,
                           meets_40_bit_selected_weight_target=F(result['log2_upper_rational']) <= -40,
                           elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                           scope='all supports at this exact selected weight; shared down-set, no interpolation')
                rows.append(row)
                if checkpoint:
                    checkpoint(row)
        except (TimeoutError, ArithmeticError) as error:
            failure = dict(phase=phase, active_groups=active, error_type=type(error).__name__, message=str(error))
    if source_pins() != pins:
        raise RuntimeError('A pinned helper changed during the sparse band replay')
    complete = len(rows) == last-first+1
    union = _union_ceiling(rows, precision)
    return dict(metadata, requested_weight_range=[first, last], weight_range=[first, last],
                coefficient_maximum=maximum, coefficient_marker=str(F(maximum, base.PAIRS)),
                regional_tables_built=tables, weight_receipts=rows, weights_evaluated=len(rows),
                every_weight_evaluated=complete, complete_requested_range=complete,
                completed_weight_range=[first, first+len(rows)-1] if rows else None,
                pending_weight_ranges=[] if complete else [[first+len(rows), last]],
                completed_union_log2_upper_rational=union, log2_upper_rational=union if complete else None,
                meets_requested_band_target=complete and F(union) <= -failure_bits,
                covers_all_supports_in_range=complete, covers_all_even_odd_splits=True,
                no_weight_interpolation=True, failure=failure, source_pins_unchanged_at_exit=True,
                elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                scope='only completed integer weights in the requested sparse band; not a whole-code certificate')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=int, required=True)
    parser.add_argument('--last', type=int, required=True)
    parser.add_argument('--maximum', type=int)
    parser.add_argument('--z', required=True)
    parser.add_argument('--alpha', required=True)
    parser.add_argument('--state', type=int, choices=local.STATES, default=24)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--failure-bits', type=int, default=40)
    parser.add_argument('--max-seconds', type=float, default=600.)
    parser.add_argument('--unpaired', action='store_true')
    parser.add_argument('--simple-return', action='store_true')
    parser.add_argument('--witness-source', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('Output directory must be unused')
    if not 8 <= args.field_bits <= 128 or args.field_bits%8:
        parser.error('Field bits must be a multiple of8in8..128')
    provenance = None
    if args.witness_source:
        raw = args.witness_source.read_bytes()
        provenance = dict(path=str(args.witness_source.resolve()), sha256=sha256(raw).hexdigest(),
                          role='exact mark provenance only; no numerical bound reused')
    args.output_dir.mkdir(exist_ok=False)

    def checkpoint(row):
        write_new(args.output_dir/f"weight_{row['active_groups']:04d}.json", row)
        print(json.dumps({key: row[key] for key in ('active_groups', 'log2_upper_rational', 'profiles',
                                                   'meets_requested_selected_weight_target')}, allow_nan=False), flush=True)

    result = batch_receipt(args.first, args.last, F(args.z), F(args.alpha), maximum=args.maximum,
                           q=1 << args.field_bits, state=args.state, distance=F(args.distance),
                           precision=args.precision, bits=args.bits, paired=not args.unpaired,
                           joint=not args.simple_return, failure_bits=args.failure_bits,
                           max_seconds=args.max_seconds, witness_source=provenance, checkpoint=checkpoint)
    write_new(args.output_dir/'summary.json', result)
    print(json.dumps({key: result[key] for key in ('weights_evaluated', 'pending_weight_ranges',
                                                  'completed_union_log2_upper_rational', 'failure')}, allow_nan=False))


if __name__ == '__main__':
    main()
