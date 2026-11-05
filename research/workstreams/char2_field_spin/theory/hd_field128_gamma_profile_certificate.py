"""Directed all-profile chronological transfer for a frozen gamma alias bank.

Every even/odd full-pair profile receives the existing exact joint-profile
count bound. Positive Arb arithmetic evaluates (R_even R_odd)^128 with initial
state zero and unrestricted final state. Completed even-profile rows have
immutable checkpoints; a partial run never supplies a selected-weight bound.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

from flint import arb, ctx

import hd_field128_full_rs_gamma_alias_certificate as provider
import hd_field128_lift_alias_cover as cover
from hd_shifted_pairs_count import profile_count_upper


base, tail = provider.base, provider.tail


def source_pins():
    pins = dict(provider.source_pins(), **cover.source_pins())
    directory = Path(__file__).resolve().parent
    for name in (Path(__file__).name, 'hd_shifted_pairs_count.py'):
        pins[name] = sha256((directory/name).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_source_pins(pins):
    # The complete closure is resolved once, then only its bytes are rehashed.
    directory = Path(__file__).resolve().parent
    if any(sha256((directory/name).read_bytes()).hexdigest() != digest for name, digest in pins.items()):
        raise RuntimeError('A pinned profile-transfer helper changed during replay')


def write_new(path, value):
    raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
    with Path(path).open('xb') as stream:
        stream.write(raw)
    return dict(path=str(Path(path).resolve()), sha256=sha256(raw).hexdigest())


class ProfileCounts:
    """Cached exact equivalent of profile_count_upper on feasible marginals."""

    def __init__(self, pairs, weight, budget):
        self.full = tail.full_counts(pairs, weight)
        self.weight = weight
        budget.check()
        self.runs = tail.circular_run_counts(2*pairs, weight)
        self.contracted = tail.circular_run_counts(pairs, weight//2) if weight % 2 == 0 else None
        first = self.full.start
        marginal = tail.one_pairing_count(pairs, weight, first)
        self.marginals = {}
        for full in self.full:
            if not full % 128:
                budget.check()
            self.marginals[full] = marginal
            if full+1 < self.full.stop:
                single, empty = weight-2*full, pairs-weight+full
                marginal, remainder = divmod(marginal*single*(single-1), 4*(full+1)*(empty+1))
                if remainder:
                    raise ArithmeticError('Nonintegral exact matching-marginal recurrence')
        budget.check()

    def count(self, first, second):
        runs = self.weight-first-second
        if 2*first == self.weight or 2*second == self.weight:
            return self.contracted[runs] if runs < len(self.contracted) else 0
        return min(self.runs[runs], self.marginals[first], self.marginals[second])


def pair_row_sum128(left, right):
    """Row-zero sum of the chronological product to power128, directed upward."""
    upper = base.upper
    a, b, c, d = left
    e, f, g, h = right
    a, b, c, d = (upper(a*e+b*g), upper(a*f+b*h),
                  upper(c*e+d*g), upper(c*f+d*h))
    # Seven fixed 2x2 squarings. Reuse bc and a+d in each positive square.
    for unused in range(7):
        bc, trace = b*c, a+d
        a, b, c, d = upper(a*a+bc), upper(b*trace), upper(c*trace), upper(d*d+bc)
    return upper(a+b)


def support_log(value):
    if value.is_zero():
        return dict(zero=True, log2_upper_rational=None, log2_upper_enclosure=None)
    logarithm = tail.log2arb(value)
    return dict(zero=False, log2_upper_rational=str(base.rational_ceiling(logarithm)),
                log2_upper_enclosure=str(logarithm))


class IncompleteTransfer(RuntimeError):
    def __init__(self, cause, partial):
        super().__init__(str(cause))
        self.cause, self.partial = cause, partial


def transfer_sum(pairs, weight, regions, *, repetitions=256, budget=None, checkpoint=None):
    """Generic small-instance arithmetic interface; only complete runs return a sum."""
    full = tail.full_counts(pairs, weight)
    if (type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or set(regions) != set(full)
            or any(len(row) != 4 for row in regions.values())):
        raise ValueError('Complete regional matrices and a positive even repetition count are required')
    budget = budget or base.Budget()
    # Freeze positive endpoints before multiplication; signed intervals are not matrices of caps.
    regions = {f: tuple(base.upper(value) for value in regions[f]) for f in full}
    total, rows, positive_profiles = arb(0), [], 0
    try:
        counts = ProfileCounts(pairs, weight, budget)
        for first in full:
            budget.check()
            row_sum, nonzero = arb(0), 0
            left = regions[first]
            for index, second in enumerate(full):
                if not index % 32:
                    budget.check()
                count = counts.count(first, second)
                if not count:
                    continue
                if repetitions == 256:
                    value = pair_row_sum128(left, regions[second])
                else:
                    powered = base.matrix_power(base.matrix_product(left, regions[second]), repetitions//2)
                    value = base.upper(powered[0]+powered[1])
                row_sum += count*value
                nonzero += 1
            row_sum = base.upper(row_sum)
            total = base.upper(total+row_sum)
            budget.check()
            row = dict(even_full=first, odd_full_range=[full.start, full.stop-1],
                       profile_pairs_evaluated=len(full), profiles_with_positive_count=nonzero,
                       row_support=support_log(row_sum), completed_support=support_log(total),
                       complete_even_profile_row=True, covers_all_odd_profiles_for_this_even_profile=True,
                       complete_selected_weight=False, covers_all_supports_at_weight=False,
                       common_potential_used=False, elapsed_seconds=monotonic()-budget.started)
            rows.append(row)
            positive_profiles += nonzero
            if checkpoint:
                checkpoint(row)
        budget.check()
    except (TimeoutError, ArithmeticError) as error:
        raise IncompleteTransfer(error, dict(rows=rows, completed_even_rows=len(rows),
                                completed_even_full_range=[full.start, full.start+len(rows)-1] if rows else None,
                                pending_even_full_range=([full.start+len(rows), full.stop-1]
                                                         if len(rows) < len(full) else None),
                                profile_pairs_evaluated=len(rows)*len(full),
                                profiles_with_positive_count=positive_profiles)) from error
    return total, dict(rows=rows, completed_even_rows=len(rows),
                       completed_even_full_range=[full.start, full.stop-1], pending_even_full_range=None,
                       profile_pairs_evaluated=len(full)**2, profiles_with_positive_count=positive_profiles)


def regional_matrices(bank, weight, budget):
    return cover.RegionalCache(bank).regional(weight, budget)


def weight_receipt(bank, weight, *, max_seconds=600., failure_bits=40,
                   witness_source=None, checkpoint=None, progress=None):
    if (not isinstance(bank, provider.AliasBank) or type(weight) is not int
            or not 1 <= weight <= base.GROUPS or type(failure_bits) is not int
            or not 0 <= failure_bits <= 4096):
        raise ValueError('Expected the frozen gamma AliasBank, weight in1..8192, and a valid failure-bit target')
    bank.check_pins()
    pins, budget = source_pins(), base.Budget(max_seconds)
    full = tail.full_counts(base.PAIRS, weight)
    metadata = dict(bank.metadata(), source_sha256=pins, coverage_driver_sha256=pins,
                    witness_source=witness_source, active_groups=weight, common_potential_used=False,
                    support_count_method='exact profile_count_upper via cached run and marginal counts',
                    transfer_method='positive Arb row-zero sum of (R_even R_odd)^128',
                    repetitions=256, paired_region_power=128, requested_failure_bits=failure_bits,
                    expected_even_rows=len(full), expected_profile_pairs=len(full)**2,
                    max_seconds=max_seconds, old_numerical_bounds_reused=False,
                    directed=True, whole_code_certificate=False)
    partial = dict(rows=[], completed_even_rows=0, completed_even_full_range=None,
                   pending_even_full_range=[full.start, full.stop-1], profile_pairs_evaluated=0,
                   profiles_with_positive_count=0)
    failure, ceiling, enclosure = None, None, None

    def save(row):
        check_source_pins(pins)
        record = dict(metadata, **row,
                      scope='completed even-profile row only; not a selected-weight bound')
        if checkpoint:
            checkpoint(record)
        if progress:
            progress(dict(stage='completed-even-profile-row', even_full=row['even_full'],
                          completed_even_rows=row['even_full']-full.start+1, expected_even_rows=len(full),
                          elapsed_seconds=row['elapsed_seconds']))

    with ctx.workprec(bank.precision):
        try:
            budget.check()
            if progress:
                progress(dict(stage='regional-normalization', active_groups=weight))
            regions = regional_matrices(bank, weight, budget)
            support, partial = transfer_sum(base.PAIRS, weight, regions, budget=budget, checkpoint=save)
            value = bank.constant+bank.slope*weight+tail.log2arb(support)
            ceiling, enclosure = str(base.rational_ceiling(value)), str(value)
            budget.check()
        except IncompleteTransfer as error:
            partial = error.partial
            failure = dict(error_type=type(error.cause).__name__, message=str(error.cause))
        except (TimeoutError, ArithmeticError) as error:
            failure = dict(error_type=type(error).__name__, message=str(error))
    bank.check_pins()
    check_source_pins(pins)
    complete = failure is None and partial['completed_even_rows'] == len(full)
    return dict(metadata, **partial, status='complete' if complete else 'incomplete', failure=failure,
                complete_selected_weight=complete, covers_all_supports_at_weight=complete,
                covers_all_even_odd_splits=complete,
                log2_upper_rational=ceiling if complete else None,
                log2_upper_enclosure=enclosure if complete else None,
                meets_requested_selected_weight_target=complete and F(ceiling) <= -failure_bits,
                elapsed_seconds=monotonic()-budget.started, source_pins_unchanged_at_exit=True,
                scope='every matching profile at this selected weight only, if complete; no partial full-bound claim')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--state', type=int, choices=provider.STATES, default=24)
    parser.add_argument('--z', required=True)
    parser.add_argument('--alpha', required=True)
    parser.add_argument('--marker', nargs=2, action='append', required=True)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--grid', type=int, default=128)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--failure-bits', type=int, default=40)
    parser.add_argument('--build-seconds', type=float, default=600.)
    parser.add_argument('--max-seconds', type=float, default=600.)
    parser.add_argument('--max-table-bytes', type=int, default=384 << 20)
    parser.add_argument('--witness-source', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() or not 1 <= args.r <= base.GROUPS:
        parser.error('Expected an unused output directory and weight in1..8192')
    base.Budget(args.build_seconds)
    base.Budget(args.max_seconds)
    raw = args.witness_source.read_bytes() if args.witness_source else None
    provenance = (dict(path=str(args.witness_source.resolve()), sha256=sha256(raw).hexdigest(),
                       role='witness provenance only; no numerical bounds reused') if raw is not None else None)
    args.output_dir.mkdir(exist_ok=False)
    progress = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
    try:
        bank = provider.AliasBank(args.z, args.alpha, args.marker, state=args.state, distance=args.distance,
                                 grid=args.grid, precision=args.precision, bits=args.bits,
                                 max_seconds=args.build_seconds, max_table_bytes=args.max_table_bytes,
                                 progress=progress)
    except (TimeoutError, ArithmeticError) as error:
        write_new(args.output_dir/'summary.json', dict(status='bank-build-incomplete',
                  construction=provider.local.construction(args.state), active_groups=args.r,
                  complete_selected_weight=False, covers_all_supports_at_weight=False,
                  log2_upper_rational=None, directed=False, whole_code_certificate=False,
                  failure=dict(error_type=type(error).__name__, message=str(error)),
                  source_sha256=source_pins(), witness_source=provenance))
        raise SystemExit(2)

    def checkpoint(row):
        if raw is not None and args.witness_source.read_bytes() != raw:
            raise RuntimeError('Pinned witness source changed')
        write_new(args.output_dir/f"even_full_{row['even_full']:04d}.json", row)

    result = weight_receipt(bank, args.r, max_seconds=args.max_seconds, failure_bits=args.failure_bits,
                            witness_source=provenance, checkpoint=checkpoint, progress=progress)
    if raw is not None and args.witness_source.read_bytes() != raw:
        raise RuntimeError('Pinned witness source changed')
    result.update(build_seconds=args.build_seconds, max_table_bytes=args.max_table_bytes)
    write_new(args.output_dir/'summary.json', result)
    progress({key: result[key] for key in ('status', 'completed_even_rows', 'pending_even_full_range',
                                          'log2_upper_rational', 'failure')})
    if not result['complete_selected_weight']:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
