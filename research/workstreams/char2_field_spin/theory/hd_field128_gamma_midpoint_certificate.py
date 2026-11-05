"""Directed all-profile transfer using the proved midpoint count upper bound.

Only support counts change. The frozen all-pair bank, normalized regional
matrices, chronological transfer, and outer prefix are reused. An incomplete
run retains completed rows but supplies no selected-weight bound.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

from flint import arb, ctx

import hd_field128_gamma_midpoint_counts as midpoint
import hd_field128_gamma_near_face_certificate as frozen


provider, base, tail = frozen.provider, frozen.base, frozen.tail
write_new, IncompleteTransfer = frozen.write_new, frozen.IncompleteTransfer
check_source_pins = frozen.check_source_pins


def source_pins():
    pins = {**frozen.source_pins(), **midpoint.source_pins()}
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def get_counts(pairs, weight, counts, max_singletons, budget):
    if counts is None:
        return midpoint.MidpointCounts(pairs, weight,
            max_singletons=64 if max_singletons is None else max_singletons, budget=budget)
    if (type(counts) is not midpoint.MidpointCounts or (counts.pairs, counts.weight) != (pairs, weight)
            or (max_singletons is not None and
                (type(max_singletons) is not int or max_singletons != counts.base.max_singletons))):
        raise ValueError('Expected matching midpoint counts without conflicting cutoff')
    budget.check()
    return counts


def transfer_sum(pairs, weight, regions, *, repetitions=256, counts=None,
                 max_singletons=None, budget=None, checkpoint=None):
    """Directed complete sum; generic dimensions serve tiny exact tests."""
    full = tail.full_counts(pairs, weight)
    if (type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or set(regions) != set(full) or any(len(row) != 4 for row in regions.values())):
        raise ValueError('Complete regional matrices and positive even repetitions required')
    budget = budget or base.Budget()
    regions = {f: tuple(base.upper(value) for value in regions[f]) for f in full}
    total, rows, positive_profiles = arb(0), [], 0
    provenance = None
    try:
        counts = get_counts(pairs, weight, counts, max_singletons, budget)
        provenance = counts.diagnostics()
        for first in full:
            budget.check()
            row_sum, nonzero = arb(0), 0
            left = regions[first]
            for position, second in enumerate(full):
                if not position % 32:
                    budget.check()
                count = counts.count(first, second)
                if not count:
                    continue
                if repetitions == 256:
                    value = frozen.frozen.pair_row_sum128(left, regions[second])
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
                       row_support=frozen.frozen.support_log(row_sum),
                       completed_support=frozen.frozen.support_log(total),
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
                  pending_even_full_range=([full.start+len(rows), full.stop-1] if len(rows) < len(full) else None),
                  profile_pairs_evaluated=len(rows)*len(full), profiles_with_positive_count=positive_profiles,
                  count_provenance=provenance)) from error
    return total, dict(rows=rows, completed_even_rows=len(rows),
                       completed_even_full_range=[full.start, full.stop-1], pending_even_full_range=None,
                       profile_pairs_evaluated=len(full)**2, profiles_with_positive_count=positive_profiles,
                       count_provenance=provenance)


def weight_receipt(bank, weight, *, counts=None, max_singletons=None,
                   max_seconds=600., failure_bits=40, witness_source=None,
                   checkpoint=None, context_checkpoint=None, progress=None):
    """Selected-weight certificate from an actual frozen gamma all-pair bank.

    Count objects can be reused across marks at the same weight. Every
    profile receives its midpoint cap. Context-linked row checkpoints make
    no full-weight claim. Use write_new for immutable files.
    """
    if (type(bank) is not provider.AliasBank or type(weight) is not int
            or not 1 <= weight <= base.GROUPS or type(failure_bits) is not int
            or not 0 <= failure_bits <= 4096):
        raise ValueError('Expected the actual gamma all-pair bank, weight 1..8192, and a valid failure target')
    bank.check_pins()
    pins, budget = source_pins(), base.Budget(max_seconds)
    full = tail.full_counts(base.PAIRS, weight)
    metadata = dict(bank.metadata(), source_sha256=pins, coverage_driver_sha256=pins,
                    witness_source=witness_source, active_groups=weight, common_potential_used=False,
                    transfer_method='frozen positive Arb row-zero sum of (R_even R_odd)^128',
                    regional_normalization='frozen RegionalCache, entrywise minimum over all bank markers',
                    repetitions=256, paired_region_power=128, requested_failure_bits=failure_bits,
                    expected_even_rows=len(full), expected_profile_pairs=len(full)**2,
                    max_seconds=max_seconds, old_numerical_bounds_reused=False,
                    directed=True, whole_code_certificate=False)
    partial = dict(rows=[], completed_even_rows=0, completed_even_full_range=None,
                   pending_even_full_range=[full.start, full.stop-1], profile_pairs_evaluated=0,
                   profiles_with_positive_count=0)
    failure = ceiling = enclosure = support_summary = context_hash = None

    def save(row):
        check_source_pins(pins)
        if checkpoint:
            checkpoint(dict(row, context_sha256=context_hash, active_groups=weight,
                            directed=True, whole_code_certificate=False,
                            scope='completed even-profile row only; requires midpoint-count metadata context'))
        if progress:
            progress(dict(stage='completed-even-profile-row', even_full=row['even_full'],
                          expected_even_rows=len(full), elapsed_seconds=row['elapsed_seconds']))

    with ctx.workprec(bank.precision):
        try:
            budget.check()
            counts = get_counts(base.PAIRS, weight, counts, max_singletons, budget)
            metadata['count_provenance'] = counts.diagnostics()
            metadata['count_provenance']['float_log_count_path_used'] = False
            context_hash = sha256(json.dumps(metadata, sort_keys=True, separators=(',', ':'),
                                            allow_nan=False).encode('utf-8')).hexdigest()
            check_source_pins(pins)
            if context_checkpoint:
                context_checkpoint(dict(context=metadata, context_sha256=context_hash,
                                        complete_selected_weight=False, covers_all_supports_at_weight=False,
                                        whole_code_certificate=False, scope='metadata only; no numerical bound'))
            if progress:
                progress(dict(stage='regional-normalization', active_groups=weight))
            regions = frozen.regional_matrices(bank, weight, budget)
            support, partial = transfer_sum(base.PAIRS, weight, regions, counts=counts,
                                            budget=budget, checkpoint=save)
            # The context records the integer midpoint policy used by every row.
            partial.pop('count_provenance')
            support_summary = frozen.frozen.support_log(support)
            value = bank.constant+bank.slope*weight+tail.log2arb(support)
            ceiling, enclosure = str(base.rational_ceiling(value)), str(value)
            budget.check()
        except IncompleteTransfer as error:
            partial = dict(error.partial)
            partial.pop('count_provenance', None)
            failure = dict(error_type=type(error.cause).__name__, message=str(error.cause))
        except (TimeoutError, ArithmeticError) as error:
            failure = dict(error_type=type(error).__name__, message=str(error))
    bank.check_pins()
    check_source_pins(pins)
    complete = failure is None and partial['completed_even_rows'] == len(full)
    return dict(metadata, **partial, context_sha256=context_hash,
                status='complete' if complete else 'incomplete', failure=failure,
                complete_selected_weight=complete, covers_all_supports_at_weight=complete,
                covers_all_even_odd_splits=complete, support_sum=support_summary if complete else None,
                log2_upper_rational=ceiling if complete else None,
                log2_upper_enclosure=enclosure if complete else None,
                meets_requested_selected_weight_target=complete and F(ceiling) <= -failure_bits,
                elapsed_seconds=monotonic()-budget.started, source_pins_unchanged_at_exit=True,
                scope='all matching profiles at one selected weight if complete; no partial or whole-code claim')
