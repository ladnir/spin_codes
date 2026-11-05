"""Compact directed bands: one actual bank, every integer weight, no interpolation.

One frozen RegionalCache serves the complete all-profile transfer at each
requested weight. The context is saved once. Each weight artifact contains
all completed even-row sums and links the context by SHA256; the summary
links these artifacts rather than duplicating their rows or bank metadata.
A timeout saves the partial weight but never includes it in the union.
"""

from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from time import monotonic

from flint import ctx

import hd_field128_gamma_near_face_certificate as selected


base, tail, provider = selected.base, selected.tail, selected.provider
cover = selected.frozen.cover
RegionalCache = cover.RegionalCache
write_new = selected.write_new


def source_pins():
    pins = selected.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def _empty_progress(weight):
    full = tail.full_counts(base.PAIRS, weight)
    return dict(rows=[], completed_even_rows=0, completed_even_full_range=None,
                pending_even_full_range=[full.start, full.stop-1], profile_pairs_evaluated=0,
                profiles_with_positive_count=0)


def _check_links(links):
    if any(sha256(Path(link['path']).read_bytes()).hexdigest() != link['sha256'] for link in links):
        raise RuntimeError('An immutable band artifact changed during the replay')


def _weight(bank, cache, weight, *, budget, pins, context, max_singletons, failure_bits, progress):
    """Only receipt plumbing is new; count and transfer arithmetic are frozen."""
    full = tail.full_counts(base.PAIRS, weight)
    partial = _empty_progress(weight)
    failure, ceiling, enclosure, support_log, prefix = None, None, None, None, None

    def row_checkpoint(row):
        selected.check_source_pins(pins)
        if progress:
            progress(dict(stage='completed-even-profile-row', active_groups=weight,
                          even_full=row['even_full'], completed_even_rows=row['even_full']-full.start+1,
                          expected_even_rows=len(full), elapsed_seconds=monotonic()-budget.started))

    try:
        budget.check()
        if progress:
            progress(dict(stage='regional-normalization', active_groups=weight))
        regions = cache.regional(weight, budget)
        support, partial = selected.transfer_sum(base.PAIRS, weight, regions, budget=budget,
                            max_singletons=max_singletons, checkpoint=row_checkpoint)
        support_log = selected.frozen.support_log(support)
        prefix = bank.constant+bank.slope*weight
        value = prefix+tail.log2arb(support)
        ceiling, enclosure = str(base.rational_ceiling(value)), str(value)
        budget.check()
    except selected.IncompleteTransfer as error:
        partial = error.partial
        failure = dict(error_type=type(error.cause).__name__, message=str(error.cause))
    except (TimeoutError, ArithmeticError) as error:
        failure = dict(error_type=type(error).__name__, message=str(error))
    bank.check_pins()
    selected.check_source_pins(pins)
    complete = failure is None and partial['completed_even_rows'] == len(full)
    return dict(partial, context=context, active_groups=weight,
                count_provenance=selected.count_provenance(max_singletons),
                expected_even_rows=len(full), expected_profile_pairs=len(full)**2,
                status='complete' if complete else 'incomplete', failure=failure,
                complete_selected_weight=complete, covers_all_supports_at_weight=complete,
                covers_all_even_odd_splits=complete, common_potential_used=False,
                support_sum=support_log if complete else None,
                log2_outer_prefix_enclosure=str(prefix) if complete else None,
                log2_upper_rational=ceiling if complete else None,
                log2_upper_enclosure=enclosure if complete else None,
                meets_requested_selected_weight_target=complete and F(ceiling) <= -failure_bits,
                requested_failure_bits=failure_bits, elapsed_seconds=monotonic()-budget.started,
                directed=True, whole_code_certificate=False, source_pins_checked_before_publication=True,
                requires_linked_context=True,
                scope='this exact weight only if complete; context supplies bank witness and source closure')


def band_receipt(bank, first, last, *, output_dir, max_singletons=32, max_weights=32,
                 max_seconds=600., failure_bits=64, witness_source=None, progress=None):
    """Save immutable context/weight/summary artifacts and return the summary.

    Every completed integer weight is evaluated from this same actual bank.
    The global wall limit applies to normalization and all profile sums.
    Final pin checks and small union/receipt finalization follow the limit.
    A partial selected weight stays pending even if all its rows completed
    before the final wall check. No prior numerical receipt is imported.
    """
    cover.checked_range(first, last, max_weights)
    provenance = selected.count_provenance(max_singletons)
    if (type(bank) is not provider.AliasBank or type(failure_bits) is not int
            or not 0 <= failure_bits <= 4096):
        raise ValueError('Expected the actual all-pair gamma bank and an integer failure target in0..4096')
    budget = base.Budget(max_seconds)
    bank.check_pins()
    pins = source_pins()
    directory = Path(output_dir).resolve()
    if directory.exists():
        raise FileExistsError('Band output directory must be unused')
    metadata = dict(bank.metadata(), schema='gamma-near-face-band-context-v1',
                    source_sha256=pins, coverage_driver_sha256=pins, witness_source=witness_source,
                    requested_weight_range=[first, last], count_provenance=provenance,
                    requested_failure_bits=failure_bits, max_seconds=max_seconds, max_weights=max_weights,
                    repetitions=256, paired_region_power=128, common_potential_used=False,
                    regional_normalization='one frozen RegionalCache; fresh normalization at every integer weight',
                    support_transfer='frozen near-face transfer_sum; every feasible even/odd profile',
                    old_numerical_bounds_reused=False, no_weight_interpolation=True,
                    complete_requested_range=False, covers_all_supports_in_range=False,
                    directed=True, whole_code_certificate=False, scope='metadata only; no probability bound')
    selected.check_source_pins(pins)
    directory.mkdir(parents=False, exist_ok=False)
    context = write_new(directory/'context.json', metadata)
    cache, completed, links, partial_link, failure = RegionalCache(bank), [], [], None, None
    with ctx.workprec(bank.precision):
        for weight in range(first, last+1):
            _check_links((context,))
            row = _weight(bank, cache, weight, budget=budget, pins=pins, context=context,
                          max_singletons=max_singletons, failure_bits=failure_bits, progress=progress)
            artifact = write_new(directory/f'weight_{weight:04d}.json', row)
            link = dict(artifact, active_groups=weight, status=row['status'],
                        complete_selected_weight=row['complete_selected_weight'],
                        log2_upper_rational=row['log2_upper_rational'],
                        completed_even_rows=row['completed_even_rows'],
                        pending_even_full_range=row['pending_even_full_range'])
            links.append(link)
            if not row['complete_selected_weight']:
                partial_link = link
                failure = dict(active_groups=weight, **row['failure'])
                break
            completed.append(link)
            if progress:
                progress(dict(stage='complete-weight', active_groups=weight,
                              log2_upper_rational=row['log2_upper_rational'],
                              completed_weights=len(completed), requested_weights=last-first+1))
        # Reuse directed union arithmetic. Partial profile rows never enter it.
        union = cover.union_ceiling(completed, bank.precision)
    bank.check_pins()
    selected.check_source_pins(pins)
    _check_links((context, *links))
    next_weight = first+len(completed)
    complete = next_weight > last
    result = dict(schema='gamma-near-face-band-summary-v1', context=context,
                  construction=bank.construction if hasattr(bank, 'construction') else metadata['construction'],
                  requested_weight_range=[first, last], weight_receipts=links,
                  completed_weight_range=[first, next_weight-1] if completed else None,
                  weights_evaluated=len(completed), every_integer_weight_evaluated=complete,
                  pending_weight_ranges=[] if complete else [[next_weight, last]],
                  partial_weight_receipt=partial_link, complete_requested_range=complete,
                  covers_all_supports_in_range=complete, covers_all_even_odd_splits=complete,
                  completed_union_log2_upper_rational=union,
                  log2_upper_rational=union if complete else None,
                  requested_range_meets_target=complete and F(union) <= -failure_bits,
                  requested_failure_bits=failure_bits,
                  weights_above_target=[row['active_groups'] for row in completed
                                        if F(row['log2_upper_rational']) > -failure_bits],
                  regional_caches_built=1, common_actual_bank=True, no_weight_interpolation=True,
                  source_pins_unchanged_at_exit=True, failure=failure,
                  elapsed_seconds=monotonic()-budget.started, max_seconds=max_seconds, max_weights=max_weights,
                  directed=True, whole_code_certificate=False,
                  scope='union over completed integer weights only; no pending-weight or whole-code claim')
    write_new(directory/'summary.json', result)
    return result
