"""Serial directed sparse bands with one regional table per exact witness.

The positive polynomial down-set s+2*f<=maximum retains every coefficient
needed at a smaller source weight. Its common marker maximum/4096 cancels
in each target normalizer. Every row uses its own source-weight prefix and
all feasible alternating matching profiles. The driver covers only the
explicitly requested band, never the whole code, and preserves unfinished
work when a resource bound ends a run. Frozen certificate helpers are read
without modification.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from math import isfinite, log2
from pathlib import Path
from time import monotonic

from flint import arb, ctx

import hd_field128_lift_certificate as cert


DEFAULT_SCHEDULE = ((1, 1, 1), (2, 2, 2), (4, 3, 5), (8, 6, 10),
                    (16, 11, 20), (24, 21, 28), (32, 29, 40),
                    (48, 41, 56), (64, 57, 64))


def source_pins():
    """Complete recursive static local-import closure, including optional FLOAT imports."""
    directory = Path(__file__).resolve().parent
    pending, seen = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        path = directory/name
        raw = path.read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw, filename=str(path))):
            modules = ([entry.name for entry in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def default_witness(center, q=cert.FIELD_ORDER):
    bits = cert.checked_q(q)
    if type(center) is not int or not 1 <= center <= cert.MAXIMUM:
        raise ValueError('Sparse witness center must lie in 1..64')
    alpha = F(1) if center <= 8 else F(2, 5) if center <= 17 else F(1, 5) if center <= 32 else F(1, 10)
    z = (1-F(center, 1200 if center <= 17 else 1400))**(bits//8)
    return z, alpha*F(8, bits)


def scheduled_bands(first=1, last=64, schedule=DEFAULT_SCHEDULE):
    if type(first) is not int or type(last) is not int or not 1 <= first <= last <= cert.MAXIMUM:
        raise ValueError('Requested sparse interval must lie in 1..64')
    result = []
    for center, lo, hi in schedule:
        if (any(type(value) is not int for value in (center, lo, hi))
                or not 1 <= center <= cert.MAXIMUM or not 1 <= lo <= hi <= cert.MAXIMUM):
            raise ValueError('Invalid bounded schedule row')
        if max(first, lo) <= min(last, hi):
            result.append((center, max(first, lo), min(last, hi)))
    covered = [weight for _, lo, hi in result for weight in range(lo, hi+1)]
    if covered != list(range(first, last+1)):
        raise ValueError('Schedule must cover each requested integer exactly once, in order')
    return tuple(result)


def _weight_from_regions(r, regions, z, alpha, q, distance, outer_measure, precision, bits,
                         paired, maximum, budget):
    """Caller keeps the regional precision context alive through aggregation."""
    budget.check()
    union, profiles = arb(0), 0
    for a in range(r//2+1):
        budget.check()
        for b in range(r//2+1):
            count = cert.profile_count_upper(cert.PAIRS, r, a, b)
            if count:
                run = cert.matrix_power(cert.matrix_product(regions[r, a], regions[r, b]), cert.REGIONS//2)
                union += count*(run[0]+run[1])
                profiles += 1
    log = lambda value: value.log()/arb(2).log()
    forbidden = int(distance*cert.OUTPUT_LENGTH)
    prefix = cert.rational(alpha)*(r*log(cert.rational(cert.outer_density(q, outer_measure)))
                                   -log(arb(q-1))-forbidden*log(cert.rational(z)))
    value = prefix+log(union)
    ceiling = cert.rational_ceiling(value)
    return dict(active_groups=r, profiles=profiles, log2_upper_rational=str(ceiling),
                log2_upper_enclosure=str(value), covers_all_supports_at_weight=True,
                covers_all_even_odd_splits=True,
                witness=dict(z=str(z), alpha=str(alpha), precision=precision,
                             packet_mantissa_bits=bits, paired=paired,
                             coefficient_maximum=maximum, coefficient_marker=str(F(maximum, cert.PAIRS))),
                directed=True, whole_code_certificate=False)


def _union_ceiling(rows, precision):
    if not rows:
        return None
    with ctx.workprec(precision):
        mass = sum((arb(2)**cert.rational(F(row['log2_upper_rational'])) for row in rows), arb(0))
        return str(cert.rational_ceiling(mass.log()/arb(2).log()))


def batch_receipt(first, last, z, alpha, *, maximum=None, q=cert.FIELD_ORDER,
                  distance=F(7, 20), precision=256, bits=80, paired=True,
                  outer_measure='invertible-toeplitz', max_seconds=60., _pins=None):
    """One directed regional table, then every weight first..last separately."""
    z, alpha, distance = map(cert.exact, (z, alpha, distance))
    field_bits = cert.checked_q(q)
    maximum = last if maximum is None else maximum
    if (any(type(value) is not int for value in (first, last, maximum, precision))
            or not 1 <= first <= last <= maximum <= cert.MAXIMUM
            or not 128 <= precision <= 4096 or not 0 < distance < F(1, 2)):
        raise ValueError('Invalid directed sparse band, maximum, or precision')
    pins = source_pins() if _pins is None else _pins
    if source_pins() != pins:
        raise RuntimeError('A pinned source changed before the band replay')
    budget = cert.Budget(max_seconds)
    rows = []
    with ctx.workprec(precision):
        regions = cert.regional_directed(z, alpha, maximum, q=q, bits=bits,
                                         paired=paired, budget=budget)
        for weight in range(first, last+1):
            rows.append(_weight_from_regions(weight, regions, z, alpha, q, distance, outer_measure,
                                              precision, bits, paired, maximum, budget))
        ceiling = _union_ceiling(rows, precision)
    budget.check()
    if source_pins() != pins:
        raise RuntimeError('A pinned transitive helper changed during the band replay')
    return dict(construction=cert.CONSTRUCTION, field_order=q, field_bits=field_bits,
                weight_range=[first, last], coefficient_maximum=maximum,
                coefficient_marker=str(F(maximum, cert.PAIRS)), regional_tables_built=1,
                weight_receipts=rows, weights_evaluated=len(rows), every_weight_evaluated=True,
                log2_upper_rational=ceiling, meets_40_bit_band_target=F(ceiling) <= -40,
                covers_all_supports_in_range=True, covers_all_even_odd_splits=True,
                target_distance=str(distance), forbidden_weight=int(distance*cert.OUTPUT_LENGTH),
                source_dimension=cert.DIMENSION, output_length=cert.OUTPUT_LENGTH,
                outer_groups=cert.GROUPS, outer_measure=outer_measure, projective_divisor=str(q-1),
                pair_swaps=True, state_persists_across_regions=True, full_field_random_coefficients=True,
                witness=dict(z=str(z), alpha=str(alpha), precision=precision,
                             packet_mantissa_bits=bits, paired=paired),
                source_sha256=pins, pin_scope='complete recursive static local-import closure, including lazy imports',
                coefficient_method='shared-positive-carry-free-Arb-truncated-bivariate',
                elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                directed=True, whole_code_certificate=False,
                scope='all indexed nonempty supports at every integer weight in the stated band')


def float_band_screen(first, last, z, alpha, *, maximum=None, q=cert.FIELD_ORDER,
                      distance=F(7, 20), paired=True, outer_measure='invertible-toeplitz', max_seconds=60.):
    """Optional bounded proposal check sharing one FLOAT regional table."""
    import hd_field128_lift_screen as screen
    maximum = last if maximum is None else maximum
    if any(type(n) is not int for n in (first, last, maximum)) or not 1 <= first <= last <= maximum <= 64:
        raise ValueError('Invalid FLOAT sparse band')
    distance, z, alpha = cert.exact(distance), float(cert.exact(z)), float(cert.exact(alpha))
    if not 0 < distance < F(1, 2):
        raise ValueError('Invalid target distance')
    budget = screen.Budget(max_seconds)
    regions = screen.regional_screen(z, alpha, maximum, q=q, paired=paired, _budget=budget)
    forbidden = int(distance*screen.OUTPUT_LENGTH)
    density = screen.outer_log_density(q, outer_measure)
    bits = screen.checked_q(q, construction=True)
    rows = []
    for weight in range(first, last+1):
        budget.check()
        prefix = alpha*(weight*density-(bits+screen.log1mexp2(-bits))-forbidden*log2(z))
        values = []
        for a in range(weight//2+1):
            budget.check()
            for b in range(weight//2+1):
                count = cert.profile_count_upper(cert.PAIRS, weight, a, b)
                if count:
                    run = screen.power(screen.multiply(regions[weight, a], regions[weight, b]), 128)
                    values.append(log2(count)+prefix+screen.total(run[0]))
        rows.append(dict(active_groups=weight, log2_upper_screen=float(screen.total(values)),
                         profiles=len(values), directed=False, whole_code_certificate=False))
    return dict(weight_range=[first, last], coefficient_maximum=maximum, regional_tables_built=1,
                rows=rows, log2_upper_screen=float(screen.total(row['log2_upper_screen'] for row in rows)),
                every_weight_evaluated=True, covers_all_supports_in_range=True,
                covers_all_even_odd_splits=True, directed=False, whole_code_certificate=False)


def run_schedule(*, first=1, last=64, schedule=DEFAULT_SCHEDULE, witnesses=None,
                 q=cert.FIELD_ORDER, distance=F(7, 20), precision=256, bits=80,
                 paired=True, outer_measure='invertible-toeplitz', max_seconds=3600.,
                 max_band_seconds=600., failure_bits=40, checkpoint=None):
    """Serial resource-bounded run; checkpoint receives immutable progress snapshots."""
    bands = scheduled_bands(first, last, schedule)
    cert.checked_q(q)
    distance = cert.exact(distance)
    if (not isfinite(max_band_seconds) or not 0 < max_band_seconds <= 3600
            or not isfinite(failure_bits) or failure_bits < 0
            or not 0 < distance < F(1, 2) or type(precision) is not int or not 128 <= precision <= 4096
            or type(bits) is not int or not 32 <= bits <= 512 or outer_measure not in cert.OUTER_MEASURES):
        raise ValueError('Invalid per-band resource limit or failure target')
    pins = source_pins()
    budget = cert.Budget(max_seconds)
    completed, failure, attempted = [], None, 0

    def snapshot():
        next_weight = first if not completed else completed[-1]['weight_range'][1]+1
        complete = next_weight > last
        union = _union_ceiling(completed, precision)
        return dict(construction=cert.CONSTRUCTION, field_order=q, field_bits=q.bit_length()-1,
                    requested_weight_range=[first, last],
                    planned_bands=[dict(center=c, weight_range=[lo, hi], coefficient_maximum=hi)
                                   for c, lo, hi in bands],
                    completed_bands=completed.copy(), completed_band_count=len(completed),
                    completed_weight_range=[first, next_weight-1] if completed else None,
                    pending_weight_ranges=[] if complete else [[next_weight, last]],
                    complete_requested_range=complete,
                    completed_union_log2_upper_rational=union,
                    requested_range_meets_target=complete and F(union) <= -failure_bits,
                    failure_bits=failure_bits, target_distance=str(distance),
                    forbidden_weight=int(distance*cert.OUTPUT_LENGTH),
                    precision=precision, packet_mantissa_bits=bits, outer_measure=outer_measure,
                    max_seconds=max_seconds, max_band_seconds=max_band_seconds,
                    elapsed_seconds=monotonic()-budget.started, failure=failure,
                    source_sha256=pins, regional_tables_in_completed_bands=len(completed),
                    regional_table_builds_attempted=attempted,
                    directed=True, whole_code_certificate=False,
                    scope='completed listed bands only; every pending weight remains unclaimed')

    if checkpoint:
        checkpoint(snapshot())
    for center, lo, hi in bands:
        try:
            budget.check()
            z, alpha = default_witness(center, q) if witnesses is None or center not in witnesses else witnesses[center]
            remaining = min(max_band_seconds, max_seconds-(monotonic()-budget.started))
            if remaining <= 0:
                raise TimeoutError('Sparse-band schedule exhausted max_seconds before the next table')
            attempted += 1
            row = batch_receipt(lo, hi, z, alpha, maximum=hi, q=q, distance=distance,
                                precision=precision, bits=bits, paired=paired, outer_measure=outer_measure,
                                max_seconds=remaining, _pins=pins)
            row['witness_center'] = center
            completed.append(row)
            if checkpoint:
                checkpoint(snapshot())
        except (TimeoutError, ArithmeticError, RuntimeError) as error:
            failure = dict(center=center, weight_range=[lo, hi], error_type=type(error).__name__, message=str(error))
            break
    result = snapshot()
    result['source_pins_unchanged_at_exit'] = source_pins() == pins
    if checkpoint:
        checkpoint(result)
    return result


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first', type=int, default=1)
    parser.add_argument('--last', type=int, default=64)
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--distance', default='7/20')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--max-seconds', type=float, default=3600.)
    parser.add_argument('--max-band-seconds', type=float, default=600.)
    parser.add_argument('--failure-bits', type=float, default=40.)
    parser.add_argument('--checkpoint-dir', type=Path, help='New directory for immutable progress snapshots')
    parser.add_argument('--output', type=Path, help='Unused final JSON output path')
    args = parser.parse_args()
    if not 8 <= args.field_bits <= 128 or args.field_bits % 8:
        parser.error('--field-bits must be a multiple of eight in 8..128')
    if args.output and args.output.exists():
        parser.error('Output already exists')
    if args.checkpoint_dir and args.checkpoint_dir.exists():
        parser.error('Checkpoint directory already exists')
    sequence = 0
    if args.checkpoint_dir:
        args.checkpoint_dir.mkdir(parents=False, exist_ok=False)

    def save(snapshot):
        nonlocal sequence
        if args.checkpoint_dir:
            path = args.checkpoint_dir/f'snapshot_{sequence:03d}.json'
            with path.open('x', encoding='utf-8') as stream:
                json.dump(snapshot, stream, indent=2, allow_nan=False)
                stream.write('\n')
            sequence += 1
        print(json.dumps(dict(completed_bands=snapshot['completed_band_count'],
                              pending_weight_ranges=snapshot['pending_weight_ranges'],
                              completed_union_log2_upper_rational=snapshot['completed_union_log2_upper_rational'],
                              failure=snapshot['failure'])), flush=True)

    result = run_schedule(first=args.first, last=args.last, q=1 << args.field_bits, distance=F(args.distance),
                          precision=args.precision, bits=args.bits, max_seconds=args.max_seconds,
                          max_band_seconds=args.max_band_seconds, failure_bits=args.failure_bits, checkpoint=save)
    encoded = json.dumps(result, indent=2, allow_nan=False)
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(encoded+'\n')
    print(encoded)


if __name__ == '__main__':
    main()
