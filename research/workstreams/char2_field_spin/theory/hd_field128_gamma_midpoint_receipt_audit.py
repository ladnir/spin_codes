"""Saved-row audit of midpoint receipts; no coefficient or FFT replay.

The midpoint integer is independently obtained from four binomials in the
central fixed-flow formula, then intersected with independent near-face
counts. Every profile is compared with the producer. The mathematical
midpoint-maximum lemma is inherited from its reviewed, pinned proof note.
"""

from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path

from flint import arb, ctx

import hd_field128_gamma_midpoint_certificate as producer
import hd_field128_gamma_near_face_band_audit as previous


base, tail, index = previous.base, previous.tail, previous.index
require, same, fraction = previous.require, previous.same, previous.fraction
DIRECTORY = Path(__file__).resolve().parent


def source_pins():
    return dict(sorted({**index.source_closure(Path(__file__).name), **producer.source_pins()}.items()))


class DirectCounts:
    def __init__(self, pairs, weight, cutoff):
        self.pairs, self.weight = pairs, weight
        self.base = previous.IndependentCounts(pairs, weight, cutoff)
        self.full = self.base.full
        # Central degrees are independent of the profile. Cache only scalar
        # binomials, not a quadratic table of profile counts.
        self.choose = lru_cache(maxsize=16384)(comb)

    def midpoint(self, first, second):
        p, r, a, b = self.pairs, self.weight, first, second
        if r > p:
            a, b, r = p-r+a, p-r+b, 2*p-r
        if a < b:
            a, b = b, a
        if r == p and a == b == 0:
            return 2
        singles = r-2*a
        u, v = singles//2, singles-singles//2
        e, f = a+u-b, a+v-b
        e0, e1, o0, o1 = p-a-u, a+u, p-a-v, a+v
        if not e0 or not o0:
            central = self.choose(p, a)
        elif not e1:
            central = self.choose(p, v)
        elif not o1:
            central = self.choose(p, u)
        else:
            tree = u*(p-r+b)*o1+a*f*o0
            numerator = p*tree*self.choose(e0, v)*self.choose(e1, a)*self.choose(o0, e)*self.choose(o1, b)
            central, remainder = divmod(numerator, e0*e1*o0*o1)
            require(not remainder, 'Nonintegral independent central fixed-flow count')
        return (singles+1)*central

    def count(self, first, second):
        return min(self.base.count(first, second), self.midpoint(first, second))


def audit_record(record, *, pairs=4096, precision=384):
    """Generic tiny-pair arithmetic core; public receipts always use P4096."""
    weight = record['active_groups']
    cutoff = record['count_provenance']['maximum_refined_singletons']
    independent = DirectCounts(pairs, weight, cutoff)
    reference = producer.midpoint.MidpointCounts(pairs, weight, max_singletons=cutoff)
    provenance = dict(reference.diagnostics(), float_log_count_path_used=False)
    require(record['count_provenance'] == provenance, 'Wrong midpoint count provenance')
    full, rows = independent.full, record['rows']
    same(record, dict(completed_even_rows=len(full), completed_even_full_range=[full.start, full.stop-1],
                      pending_even_full_range=None, expected_even_rows=len(full), expected_profile_pairs=len(full)**2,
                      profile_pairs_evaluated=len(full)**2), 'selected-weight coverage')
    require(len(rows) == len(full), 'Missing or duplicated rows')
    total, rounded, positives, digest = arb(0), arb(0), 0, sha256()
    with ctx.workprec(precision):
        for first, saved in zip(full, rows):
            nonzero = 0
            for second in full:
                central = independent.midpoint(first, second)
                require(central == reference.midpoint_upper(first, second), 'Midpoint differs from independent fixed flow')
                count = min(independent.base.count(first, second), central)
                require(count == reference.count(first, second), 'Final midpoint count differs from independent minimum')
                nonzero += bool(count)
                encoded = count.to_bytes((count.bit_length()+7)//8, 'big')
                digest.update(len(encoded).to_bytes(4, 'big'))
                digest.update(encoded)
            same(saved, dict(even_full=first, odd_full_range=[full.start, full.stop-1], profile_pairs_evaluated=len(full),
                             profiles_with_positive_count=nonzero, complete_even_profile_row=True,
                             covers_all_odd_profiles_for_this_even_profile=True, complete_selected_weight=False,
                             covers_all_supports_at_weight=False, common_potential_used=False), 'row count/coverage')
            logarithm = previous._log_record(saved['row_support'], precision)
            if logarithm is None:
                require(nonzero == 0, 'Zero row with positive-count profiles')
            else:
                total += arb(2)**logarithm.upper()
                rounded += arb(2)**base.rational(fraction(saved['row_support']['log2_upper_rational']))
            cumulative = previous._log_record(saved['completed_support'], precision)
            if total.is_zero():
                require(cumulative is None, 'Nonzero cumulative support for zero rows')
            else:
                require(cumulative is not None and base.rational(fraction(saved['completed_support']['log2_upper_rational']))
                        >= tail.log2arb(total).upper(), 'Cumulative ceiling misses saved rows')
            positives += nonzero
        same(record, dict(profiles_with_positive_count=positives), 'positive profile total')
        support = previous._log_record(record['support_sum'], precision)
        require(support is not None and not total.is_zero(), 'Expected positive complete support')
        require(base.rational(fraction(record['support_sum']['log2_upper_rational'])) >= tail.log2arb(total).upper(),
                'Support ceiling misses saved rows')
        prefix = previous.outer_prefix(record, weight)
        value = prefix+tail.log2arb(total)
        bound = index.scalar_ceiling(record, precision)
        require(base.rational(bound) >= value.upper(), 'Selected ceiling misses fresh prefix plus rows')
        require(arb(record['log2_upper_enclosure']).overlaps(prefix+support), 'Summary and support/prefix disagree')
        return dict(active_groups=weight, even_rows=len(rows), profile_pairs=len(full)**2,
                    profiles_with_positive_count=positives, exact_count_table_sha256=digest.hexdigest(),
                    count_hash_encoding='row-major unsigned big-endian counts with4-byte big-endian byte lengths',
                    midpoint_formula='(s+1)*independent four-binomial fixed_flow_count at floor(s/2), with constant/alternating exceptions',
                    all_producer_midpoint_and_final_counts_match=True, fresh_outer_prefix_enclosure=str(prefix),
                    saved_log2_upper_rational=str(bound),
                    saved_row_reaggregation_log2_upper_rational=str(base.rational_ceiling(value)),
                    rounded_row_reaggregation_log2_upper_rational=str(base.rational_ceiling(prefix+tail.log2arb(rounded))),
                    saved_ceiling_dominates_row_enclosures=True, binomial_cache_entries=independent.choose.cache_info().currsize)


def audit_receipts(paths):
    require(isinstance(paths, (list, tuple)) and 1 <= len(paths) <= 8, 'Expected1..8 explicit selected receipts')
    pins, loader, reports, seen = source_pins(), index.Validator(F(2, 5)), [], set()
    for path in paths:
        record, pin = loader.load(path)
        same(record, dict(previous.EXPECTED, status='complete', failure=None, complete_selected_weight=True,
                          covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
                          repetitions=256, paired_region_power=128, common_potential_used=False,
                          source_pins_unchanged_at_exit=True, source_sha256=producer.source_pins(),
                          coverage_driver_sha256=producer.source_pins(), bank_source_sha256=producer.provider.source_pins()),
             'complete gamma30/d40 midpoint receipt')
        weight = record['active_groups']
        require(type(weight) is int and 1 <= weight <= 8192 and weight not in seen, 'Invalid or duplicate selected weight')
        seen.add(weight)
        witness = record['witness']
        require(0 < fraction(witness['z']) < 1 and 0 < fraction(witness['alpha']) <= 1, 'Invalid exact marks')
        require(1 <= len(witness['markers']) <= 128 and all(len(row) == 2 and all(fraction(x) > 0 for x in row)
                for row in witness['markers']), 'Invalid positive marker bank')
        same(witness, dict(state_dimension=30, inner_width=64, gamma_extension_degree=16,
                           gamma_base_field_order=256, joint_rank_return_cap=True, all_pair_before_alpha=True), 'witness scope')
        require(type(witness['precision']) is int and 128 <= witness['precision'] <= 4096
                and witness['grid'] in (64, 128, 256) and 32 <= witness['packet_mantissa_bits'] <= 512, 'Invalid bank arithmetic parameters')
        proposal, proposal_pin = loader.dependency(record['witness_source'], pin, 'literal-proposal')
        same(witness, proposal['witness'], 'exact proposal bank literals')
        source = proposal.get('source')
        if source:
            target = Path(source['path'])
            loader.load(target if target.is_absolute() else DIRECTORY/target, expected=source['sha256'], role='proposal-source-literals')
        report = audit_record(record)
        target = record['requested_failure_bits']
        require(type(target) is int and 0 <= target <= 4096, 'Invalid failure target')
        same(record, dict(meets_requested_selected_weight_target=fraction(record['log2_upper_rational']) <= -target), 'selected target flag')
        reports.append(dict(report, receipt=pin, proposal=proposal_pin, requested_failure_bits=target,
                            meets_requested_target=record['meets_requested_selected_weight_target']))
    band_union = index.union_ceiling([fraction(report['saved_log2_upper_rational']) for report in reports], 384)
    require(source_pins() == pins, 'Audit closure changed')
    producer.check_source_pins(pins)
    previous.band._check_links(tuple(loader.files.values()))
    return dict(schema='gamma-midpoint-selected-receipts-audit-v1', audit_result='PASS',
                construction=previous.EXPECTED['construction'], state_dimension=30, target_distance='2/5',
                selected_weights=sorted(seen), weight_audits=reports, precision_bits=384,
                selected_events_union_log2_upper_rational=str(band_union), source_sha256=pins,
                checked_file_pins=list(loader.files.values()), coefficient_replay=False,
                whole_code_certificate=False, coverage_index_modified=False,
                audit_limit='Counts and saved row/prefix arithmetic only; no transfer coefficient replay. '
                            'Inline complete rows are pinned by receipt bytes; no independent checkpoint context reconstruction. '
                            'The midpoint-maximum theorem is inherited from its reviewed pinned note, not reproved here.')
