"""Directed sparse coefficients through weight256, without changing frozen helpers.

The local packet has at most64 occupied coordinates, and a composed packet
pair at most128. Its complete local coefficients therefore suffice for any
regional truncation through256. Positive down-set truncation and the marker
normalizer are unchanged. Every retained result covers its explicit weights
only; an interrupted run leaves all uncomputed weights pending.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, isfinite
from pathlib import Path
from time import monotonic

from flint import arb_poly, ctx

import hd_field128_lift_certificate as base
import hd_field128_lift_sparse_bands as bands


MAXIMUM, LOCAL_MAXIMUM = 256, 128
DEFAULT_RANGES = ((65, 70), (194, 196), (210, 256))
OLD_NAME = 'hd_shifted_sparse_checkpoint_directed256.json'
OLD_SHA256 = '765fc1ae5086d91031692338e4930266333a15c60f838ae92731a75b4088e5d8'
OLD_WITNESSES = ((65, 72, '333/350', 128), (73, 90, '659/700', 128),
                 (91, 110, '13/14', 128), (111, 128, '32/35', 128),
                 (129, 160, '157/175', 256), (161, 192, '153/175', 256),
                 (193, 225, '149/175', 256), (226, 256, '29/35', 256))


def source_pins():
    # bands recursively pins both local imports above and their dependencies.
    pins = bands.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def checked_ranges(ranges, *, lower=1):
    result, previous = [], lower-1
    for row in ranges:
        if (len(row) != 2 or any(type(x) is not int for x in row)
                or not lower <= row[0] <= row[1] <= MAXIMUM or row[0] <= previous):
            raise ValueError('Weight ranges must be ordered, nonempty, disjoint, and bounded')
        result.append(tuple(row))
        previous = row[1]
    if not result:
        raise ValueError('At least one requested range is required')
    return tuple(result)


def pending_ranges(requested, completed):
    missing, start, previous = [], None, None
    for lo, hi in requested:
        for r in range(lo, hi+1):
            if r in completed:
                if start is not None:
                    missing.append([start, previous])
                    start = None
            elif start is None:
                start = r
            elif previous != r-1:
                missing.append([start, previous])
                start = r
            previous = r
        if start is not None:
            missing.append([start, hi])
            start = None
    return missing


def regional_from_parts(single, composed, maximum, *, budget, return_maximum=None):
    """Unchanged positive regional algebra; optional extraction of a prefix.

    Inputs are local caps, not receipt data. A test may restrict their local
    degrees when comparing only a prefix that cannot depend on omitted terms.
    Production calls supply every local degree up to min(maximum,128).
    """
    if type(maximum) is not int or not 1 <= maximum <= MAXIMUM:
        raise ValueError('Regional coefficient cutoff must lie in1..256')
    return_maximum = maximum if return_maximum is None else return_maximum
    if type(return_maximum) is not int or not 0 <= return_maximum <= maximum:
        raise ValueError('Invalid regional extraction prefix')
    marker = F(maximum, base.PAIRS)
    if composed is None:
        run = base.polynomial_power(base.encode(single, maximum, marker), 128, maximum, budget)
    else:
        good = base.encode({k: v for k, v in single.items() if k[0] <= 1}, maximum, marker)
        dirty = base.encode({k: v for k, v in single.items() if k[0] >= 2}, maximum, marker)
        parts = (base.encode(composed, maximum, marker),
                 base.polynomial_product(good, dirty, maximum, budget),
                 base.polynomial_product(dirty, good, maximum, budget),
                 base.polynomial_product(dirty, dirty, maximum, budget))
        polynomial = tuple(sum((part[i] for part in parts), arb_poly()) for i in range(4))
        run = base.polynomial_power(polynomial, 64, maximum, budget)
    stride, regions = 2*maximum+1, {}
    for s, f in base.valid_degrees(return_maximum):
        if not s:
            budget.check()
        normalizer = comb(base.PAIRS, s)*comb(base.PAIRS-s, f)*base.rational(marker)**(s+f)
        regions[s+2*f, f] = tuple(base.upper(run[i][s+stride*f]/normalizer) for i in range(4))
    budget.check()
    return regions


def regional_directed(z, alpha, maximum, *, q=base.FIELD_ORDER, bits=80, paired=True,
                      budget=None, return_maximum=None):
    z, alpha = map(base.exact, (z, alpha))
    base.checked_q(q)
    if (type(maximum) is not int or not 1 <= maximum <= MAXIMUM or not 0 < z < 1
            or not 0 < alpha <= 1 or type(bits) is not int or not 32 <= bits <= 512):
        raise ValueError('Invalid wide sparse exact witness or resource scope')
    budget = budget or base.Budget()
    single, composed = base.directed_parts(z, alpha, min(maximum, LOCAL_MAXIMUM), q, bits, budget, paired)
    return regional_from_parts(single, composed, maximum, budget=budget, return_maximum=return_maximum)


def band_receipt(ranges, z, alpha, *, maximum, q=base.FIELD_ORDER, distance=F(7, 20),
                 precision=256, bits=80, paired=True, outer_measure='invertible-toeplitz',
                 max_seconds=600., failure_bits=40, progress=None):
    ranges = checked_ranges(ranges)
    z, alpha, distance = map(base.exact, (z, alpha, distance))
    field_bits = base.checked_q(q)
    if (type(maximum) is not int or not ranges[-1][1] <= maximum <= MAXIMUM
            or type(precision) is not int or not 128 <= precision <= 4096
            or type(failure_bits) is not int or failure_bits < 0 or not 0 < distance < F(1, 2)
            or outer_measure not in base.OUTER_MEASURES):
        raise ValueError('Invalid wide band witness or target')
    pins, budget = source_pins(), base.Budget(max_seconds)
    rows, failure, table_built = [], None, False
    with ctx.workprec(precision):
        try:
            regions = regional_directed(z, alpha, maximum, q=q, bits=bits, paired=paired,
                                         budget=budget, return_maximum=ranges[-1][1])
            table_built = True
            for lo, hi in ranges:
                for r in range(lo, hi+1):
                    row = bands._weight_from_regions(r, regions, z, alpha, q, distance, outer_measure,
                                                     precision, bits, paired, maximum, budget)
                    budget.check()
                    rows.append(row)
                    if progress:
                        progress(dict(stage='complete-weight', active_groups=r,
                                      log2_upper_rational=row['log2_upper_rational']))
        except (TimeoutError, ArithmeticError) as error:
            failure = dict(error_type=type(error).__name__, message=str(error))
        ceiling = bands._union_ceiling(rows, precision)
    if source_pins() != pins:
        raise RuntimeError('A pinned wide sparse helper changed during replay')
    pending = pending_ranges(ranges, {row['active_groups'] for row in rows})
    complete = not pending
    return dict(construction=base.CONSTRUCTION, field_order=q, field_bits=field_bits,
                requested_weight_ranges=[list(row) for row in ranges], pending_weight_ranges=pending,
                weight_receipts=rows, weights_evaluated=len(rows), every_requested_weight_evaluated=complete,
                complete_requested_ranges=complete, coefficient_maximum=maximum,
                coefficient_marker=str(F(maximum, base.PAIRS)), local_coefficient_maximum=min(maximum, LOCAL_MAXIMUM),
                regional_tables_built=int(table_built), completed_union_log2_upper_rational=ceiling,
                log2_upper_rational=ceiling if complete else None,
                requested_ranges_meet_target=complete and F(ceiling) <= -failure_bits,
                covers_all_supports_at_completed_weights=True, covers_all_even_odd_splits=True,
                target_distance=str(distance), forbidden_weight=int(distance*base.OUTPUT_LENGTH),
                source_dimension=base.DIMENSION, output_length=base.OUTPUT_LENGTH, outer_groups=base.GROUPS,
                outer_measure=outer_measure, projective_divisor=str(q-1), pair_swaps=True,
                state_persists_across_regions=True, full_field_random_coefficients=True,
                witness=dict(z=str(z), alpha=str(alpha), precision=precision, packet_mantissa_bits=bits, paired=paired),
                coefficient_method='positive-carry-free-Arb-truncated-bivariate-cutoff256',
                source_sha256=pins, elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                failure=failure, failure_bits=failure_bits, directed=True, whole_code_certificate=False,
                scope='completed explicit weights only; no claim for gaps or pending weights')


def old_witness_plan(path=None, ranges=DEFAULT_RANGES):
    ranges = checked_ranges(ranges, lower=65)
    path = Path(path or Path(__file__).with_name(OLD_NAME)).resolve()
    raw = path.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != OLD_SHA256:
        raise ValueError('Pinned old sparse checkpoint changed')
    receipt = json.loads(raw)
    accepted = {tuple(row['active_group_range']): row for row in receipt['accepted']}
    jobs = []
    for lo, hi, old_z, maximum in OLD_WITNESSES:
        row = accepted[lo, hi]
        witness = row['witness']
        if (row.get('construction') != 'toeplitz-gf256-h4-alternating-swapped-pairs-k20-rate-half'
                or row.get('field_order') != 256 or row.get('pair_swaps') is not True
                or witness['z'] != old_z or witness['alpha'] != '1/10'
                or witness['coefficient_maximum'] != maximum):
            raise ValueError('Old sparse witness semantics differ from the recorded proposal')
        clipped = [list((max(lo, a), min(hi, b))) for a, b in ranges if max(lo, a) <= min(hi, b)]
        if clipped:
            z = F(str(float(F(old_z)**16)))
            jobs.append(dict(old_weight_range=[lo, hi], requested_weight_ranges=clipped,
                             coefficient_maximum=maximum, z=str(z), alpha='1/160',
                             old_z=old_z, old_alpha='1/10', old_numerical_ceilings_reused=False))
    covered = [r for job in jobs for a, b in job['requested_weight_ranges'] for r in range(a, b+1)]
    expected = [r for a, b in ranges for r in range(a, b+1)]
    if covered != expected:
        raise ValueError('Old witness plan leaves a gap or duplicate in requested weights')
    return dict(input_provenance=dict(path=str(path), sha256=digest), requested_weight_ranges=list(map(list, ranges)),
                jobs=jobs, transform='z=Fraction(str(float(old_z**16))); alpha=old_alpha/16; new exact proposals',
                old_numerical_ceilings_reused=False, directed=False, whole_code_certificate=False)


def run_plan(*, ranges=DEFAULT_RANGES, old_checkpoint=None, precision=256, bits=80,
             distance=F(7, 20), max_seconds=3600., max_job_seconds=600., failure_bits=40, checkpoint=None):
    distance = base.exact(distance)
    if (not isfinite(max_job_seconds) or not 0 < max_job_seconds <= 3600
            or type(precision) is not int or not 128 <= precision <= 4096
            or type(bits) is not int or not 32 <= bits <= 512
            or type(failure_bits) is not int or failure_bits < 0 or not 0 < distance < F(1, 2)):
        raise ValueError('Invalid wide plan parameters or per-job deadline')
    plan = old_witness_plan(old_checkpoint, ranges)
    budget, pins = base.Budget(max_seconds), source_pins()
    receipts, rows, failure = [], [], None

    def snapshot():
        ceiling = bands._union_ceiling(rows, precision)
        covered = {row['active_groups'] for row in rows}
        pending = pending_ranges(plan['requested_weight_ranges'], covered)
        return dict(construction=base.CONSTRUCTION, field_order=base.FIELD_ORDER, field_bits=128,
                    target_distance=str(distance), forbidden_weight=int(distance*base.OUTPUT_LENGTH),
                    source_dimension=base.DIMENSION, output_length=base.OUTPUT_LENGTH, outer_groups=base.GROUPS,
                    outer_measure='invertible-toeplitz', projective_divisor=str(base.FIELD_ORDER-1), pair_swaps=True,
                    state_persists_across_regions=True, full_field_random_coefficients=True,
                    plan=plan, band_receipts=receipts.copy(), weight_receipts=rows.copy(),
                    requested_weight_ranges=plan['requested_weight_ranges'], pending_weight_ranges=pending,
                    completed_weight_count=len(rows), duplicate_weight_count=len(rows)-len(covered),
                    complete_requested_ranges=not pending, completed_union_log2_upper_rational=ceiling,
                    log2_upper_rational=ceiling if not pending else None,
                    requested_ranges_meet_target=not pending and F(ceiling) <= -failure_bits,
                    failure_bits=failure_bits, source_sha256=pins, precision=precision, packet_mantissa_bits=bits,
                    elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds,
                    max_job_seconds=max_job_seconds, failure=failure,
                    directed=True, whole_code_certificate=False,
                    scope='completed explicit weights only; no claim for omitted or pending weights')

    if checkpoint:
        checkpoint(snapshot())
    for index, job in enumerate(plan['jobs']):
        try:
            budget.check()
            remaining = min(max_job_seconds, max_seconds-(monotonic()-budget.started))
            if remaining <= 0:
                raise TimeoutError('Wide sparse plan exhausted its time budget')
            receipt = band_receipt(job['requested_weight_ranges'], job['z'], job['alpha'],
                                   maximum=job['coefficient_maximum'], precision=precision, bits=bits,
                                   distance=distance, max_seconds=remaining, failure_bits=failure_bits)
            receipt['old_witness_job'] = index
            receipts.append(receipt)
            rows.extend(dict(row, job_index=index) for row in receipt['weight_receipts'])
            if not receipt['complete_requested_ranges']:
                failure = dict(job_index=index, detail=receipt['failure'])
            if checkpoint:
                checkpoint(snapshot())
            if failure:
                break
        except (TimeoutError, ArithmeticError) as error:
            failure = dict(job_index=index, error_type=type(error).__name__, message=str(error))
            break
    if source_pins() != pins:
        raise RuntimeError('Pinned wide sparse sources changed during the plan')
    result = snapshot()
    if checkpoint:
        checkpoint(result)
    return result


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--range', type=int, nargs=2, action='append', dest='ranges')
    parser.add_argument('--old-checkpoint', type=Path)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--distance', default='7/20')
    parser.add_argument('--max-seconds', type=float, default=3600.)
    parser.add_argument('--max-job-seconds', type=float, default=600.)
    parser.add_argument('--failure-bits', type=int, default=40)
    parser.add_argument('--checkpoint-dir', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    if args.output.exists() or (args.checkpoint_dir and args.checkpoint_dir.exists()):
        parser.error('Output/checkpoint paths must be unused')
    if args.checkpoint_dir:
        args.checkpoint_dir.mkdir(parents=False, exist_ok=False)
    sequence = 0

    def save(result):
        nonlocal sequence
        if args.checkpoint_dir:
            with (args.checkpoint_dir/f'snapshot_{sequence:03d}.json').open('x', encoding='utf-8') as stream:
                json.dump(result, stream, indent=2, allow_nan=False)
                stream.write('\n')
            sequence += 1
        print(json.dumps(dict(completed_weight_count=result['completed_weight_count'],
                              pending_weight_ranges=result['pending_weight_ranges'],
                              completed_union_log2_upper_rational=result['completed_union_log2_upper_rational'],
                              failure=result['failure'])), flush=True)

    ranges = args.ranges or DEFAULT_RANGES
    result = (old_witness_plan(args.old_checkpoint, ranges) if args.plan_only else
              run_plan(ranges=ranges, old_checkpoint=args.old_checkpoint, precision=args.precision, bits=args.bits,
                       distance=F(args.distance), max_seconds=args.max_seconds, max_job_seconds=args.max_job_seconds,
                       failure_bits=args.failure_bits, checkpoint=save))
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=str(args.output))))


if __name__ == '__main__':
    main()
