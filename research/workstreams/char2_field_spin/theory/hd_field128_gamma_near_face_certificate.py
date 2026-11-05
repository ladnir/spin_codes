"""Directed all-profile transfer with exact singleton-deletion count bounds.

Only the support count changes from the frozen profile certificate. The
actual all-pair gamma bank, normalized regional matrices, chronological
power128, and outer prefix are reused. A partial run has no weight bound.
Compact row checkpoints refer to one exact metadata context by SHA256;
call context_checkpoint to persist that context before arithmetic starts.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

from flint import arb, ctx

import hd_field128_gamma_all_pair_certificate as provider
import hd_field128_gamma_near_face_counts as near
import hd_field128_gamma_profile_certificate as frozen


base, tail = provider.base, provider.tail
write_new = frozen.write_new
check_source_pins = frozen.check_source_pins
IncompleteTransfer = frozen.IncompleteTransfer


def source_pins():
    pins = {**near.source_pins(), **provider.source_pins(), **frozen.source_pins()}
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def count_provenance(max_singletons):
    if type(max_singletons) is not int or not 0 <= max_singletons <= 8192:
        raise ValueError('Expected an integer optional singleton cutoff in0..8192')
    return dict(method='minimum of frozen profile count and both matching deletion bounds',
                exact_integer_counts=True, maximum_refined_singletons=max_singletons,
                cutoff_role='omits only optional improvements; every profile is retained',
                matching_bound='C(P,s)*2^s*sum(A(P-s,f,t),t=max(0,k-s)..k)',
                variables='s=r-2f, k=r-a-b; indexed circular words, no rotation factor',
                residual_constants='A(n,0,0)=A(n,n,0)=A(0,0,0)=1; impossible cases zero',
                count_function='hd_field128_gamma_near_face_counts.NearFaceCounts.count',
                float_log_count_path_used=False)


def transfer_sum(pairs, weight, regions, *, repetitions=256, max_singletons=32,
                 budget=None, checkpoint=None):
    """Complete directed transfer; the generic dimensions are for tiny tests."""
    count_provenance(max_singletons)
    full = tail.full_counts(pairs, weight)
    if (type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or set(regions) != set(full) or any(len(row) != 4 for row in regions.values())):
        raise ValueError('Complete regional matrices and a positive even repetition count are required')
    budget = budget or base.Budget()
    regions = {f: tuple(base.upper(value) for value in regions[f]) for f in full}
    total, rows, positive_profiles = arb(0), [], 0
    try:
        counts = near.NearFaceCounts(pairs, weight, max_singletons=max_singletons, budget=budget)
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
                    value = frozen.pair_row_sum128(left, regions[second])
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
                       row_support=frozen.support_log(row_sum), completed_support=frozen.support_log(total),
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
                  profile_pairs_evaluated=len(rows)*len(full), profiles_with_positive_count=positive_profiles)) from error
    return total, dict(rows=rows, completed_even_rows=len(rows),
                       completed_even_full_range=[full.start, full.stop-1], pending_even_full_range=None,
                       profile_pairs_evaluated=len(full)**2, profiles_with_positive_count=positive_profiles)


def regional_matrices(bank, weight, budget):
    return frozen.cover.RegionalCache(bank).regional(weight, budget)


def weight_receipt(bank, weight, *, max_singletons=32, max_seconds=600., failure_bits=40,
                   witness_source=None, checkpoint=None, context_checkpoint=None, progress=None):
    """One selected weight; only a complete source-checked result is a bound.

    Row checkpoints carry no full-weight claim. Their context digest refers
    to the single metadata object returned here and optionally checkpointed
    first. Persist callbacks with write_new for immutable on-disk artifacts.
    """
    provenance = count_provenance(max_singletons)
    if (type(bank) is not provider.AliasBank or type(weight) is not int
            or not 1 <= weight <= base.GROUPS or type(failure_bits) is not int
            or not 0 <= failure_bits <= 4096):
        raise ValueError('Expected the actual all-pair gamma AliasBank, weight1..8192, and a valid failure target')
    bank.check_pins()
    pins, budget = source_pins(), base.Budget(max_seconds)
    full = tail.full_counts(base.PAIRS, weight)
    metadata = dict(bank.metadata(), source_sha256=pins, coverage_driver_sha256=pins,
                    witness_source=witness_source, active_groups=weight, common_potential_used=False,
                    support_count_method=provenance['method'], count_provenance=provenance,
                    transfer_method='frozen positive Arb row-zero sum of (R_even R_odd)^128',
                    regional_normalization='frozen RegionalCache, entrywise minimum over all bank markers',
                    repetitions=256, paired_region_power=128, requested_failure_bits=failure_bits,
                    expected_even_rows=len(full), expected_profile_pairs=len(full)**2,
                    max_seconds=max_seconds, old_numerical_bounds_reused=False,
                    directed=True, whole_code_certificate=False)
    # One context, rather than repeating the complete bank/marker metadata in every row.
    context_hash = sha256(json.dumps(metadata, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode('utf-8')).hexdigest()
    check_source_pins(pins)
    if context_checkpoint:
        context_checkpoint(dict(context=metadata, context_sha256=context_hash,
                                complete_selected_weight=False, covers_all_supports_at_weight=False,
                                whole_code_certificate=False, scope='metadata only; no numerical bound'))
    partial = dict(rows=[], completed_even_rows=0, completed_even_full_range=None,
                   pending_even_full_range=[full.start, full.stop-1], profile_pairs_evaluated=0,
                   profiles_with_positive_count=0)
    failure, ceiling, enclosure, support_summary = None, None, None, None

    def save(row):
        check_source_pins(pins)
        record = dict(row, context_sha256=context_hash, active_groups=weight,
                      maximum_refined_singletons=max_singletons, directed=True,
                      whole_code_certificate=False,
                      scope='completed even-profile row only; requires pinned metadata context; not a weight bound')
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
            support, partial = transfer_sum(base.PAIRS, weight, regions, max_singletons=max_singletons,
                                            budget=budget, checkpoint=save)
            support_summary = frozen.support_log(support)
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
    return dict(metadata, **partial, context_sha256=context_hash,
                status='complete' if complete else 'incomplete', failure=failure,
                complete_selected_weight=complete, covers_all_supports_at_weight=complete,
                covers_all_even_odd_splits=complete,
                support_sum=support_summary if complete else None,
                log2_upper_rational=ceiling if complete else None,
                log2_upper_enclosure=enclosure if complete else None,
                meets_requested_selected_weight_target=complete and F(ceiling) <= -failure_bits,
                elapsed_seconds=monotonic()-budget.started, source_pins_unchanged_at_exit=True,
                scope='every matching profile at this selected weight, if complete; no partial or whole-code claim')
