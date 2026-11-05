"""Two fixed FLOAT cap banks shared across all sparse r17 cohort splits.

Default scope is the single ordered source class (8,0).  The explicit
all-sources option covers all81 ordered source classes without new fits.
Every group keeps its constituent row set through all32 columns.  The
same universal c[s,d] bank bounds each split n1+n2=17; only the positive
support polynomial and once-only class multiplicities change.  No fit,
directed replay, or whole-code claim is performed by this module.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_geometry_strata import MultiplicityCapsFloat, geometry_stratum_counts
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_gate import ordered_region_tables, sparse_log_count
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION


GROUPS = 17
SOURCE = (8, 0)
CANONICAL_BINS = ((1, 8), (9, 14), (15, 15), (16, 16), (17, 17))
WITNESSES = (
    ('inherited', F('996248217/1000000000'), F('-127/2')),
    ('sparse-mixed-refit', F('4487948476156805/4503599627370496'), F('-8814329/131072')),
)


def cohort_bins(n2):
    """Disjoint exhaustive multiplicity bins; zero pairs is one empty class."""
    if type(n2) is not int or not 0 <= n2 <= GROUPS:
        raise ValueError('An integer two-row cohort size from0 through17 is required')
    if n2 == 0:
        return ((0, 0),)
    return tuple((low, min(high, n2)) for low, high in CANONICAL_BINS if low <= n2)


def geometry_bins(n2):
    bins = cohort_bins(n2)
    counts = geometry_stratum_counts(n2) if n2 else {0: 1}
    rows = tuple((low, high, sum(counts[mu] for mu in range(low, high+1)))
                 for low, high in bins)
    if sum(count for _, _, count in rows) != 28**n2:
        raise ArithmeticError('Multiplicity bins must count each fixed assignment once')
    return rows


def class_row(n1, bounds, geometries, source_count, caps, label, z, logp):
    """Evaluate one fixed cohort/bin using a previously computed universal cap."""
    if type(n1) is not int or not 0 <= n1 <= GROUPS:
        raise ValueError('An integer one-row cohort size from0 through17 is required')
    n2 = GROUPS-n1
    if tuple(bounds) not in cohort_bins(n2):
        raise ValueError('Multiplicity bin is outside this cohort partition')
    rows = n1+2*n2
    coloring = comb(GROUPS, n1)
    value = (log2(source_count)+log2(coloring)+3*n1+log2(geometries)
             -log2(Q-1)-LIMIT*log2(float(z))-float(logp)
             +sparse_log_count(caps, n1, n2))
    return dict(cohort_counts=[n1, n2, 0], mu_range=list(bounds),
        cap_mu=max(1, bounds[1]), source_supports=source_count,
        cohort_colorings=coloring, one_row_choices=8**n1,
        two_row_geometries=geometries, sparse_active_rows=rows,
        sparse_support_floor=17*rows, sparse_support_normalizer_exponent=16*rows,
        witness=label, mark=str(z), log2_potential=str(logp),
        failure_log2_screen=float(value))


def _source_screen(source, source_count, partitions, raw_banks, *, progress=None):
    """Share each source-conditioned cap across every applicable cohort."""
    started = perf_counter()
    required = sorted({max(1, high) for bins in partitions.values() for _, high, _ in bins})
    banks, selected = [], {}
    for label, z, logp in WITNESSES:
        bank_start = perf_counter()
        tables = ordered_region_tables(raw_banks[z], GROUPS, *source)
        evaluator = MultiplicityCapsFloat(tables, float(logp), GROUPS)
        caps = {}
        for upper in required:
            caps[upper] = evaluator(upper)
            if progress:
                progress(dict(stage='shared-cap', witness=label, cap_mu=upper,
                              elapsed_seconds=perf_counter()-bank_start))
        entries = []
        for n1 in range(GROUPS+1):
            n2 = GROUPS-n1
            for low, high, geometries in partitions[n2]:
                row = class_row(n1, (low, high), geometries, source_count,
                    caps[max(1, high)], label, z, logp)
                row['source_matching_counts'] = list(source)
                entries.append(row)
                key = n1, low, high
                if key not in selected or row['failure_log2_screen'] < selected[key]['failure_log2_screen']:
                    selected[key] = row
            if progress:
                value = float(np.logaddexp2.reduce([
                    row['failure_log2_screen'] for row in entries if row['cohort_counts'][0] == n1]))
                progress(dict(stage='cohort', witness=label, cohort_counts=[n1, n2, 0],
                              failure_log2_screen=value))
        banks.append(dict(witness=label, mark=str(z), log2_potential=str(logp), entries=entries,
            cap_mu_values=required, resource_tables=evaluator._envelope.cache_info().currsize,
            failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])),
            elapsed_seconds=perf_counter()-bank_start))
    cohorts = []
    for n1 in range(GROUPS+1):
        entries = [selected[n1, low, high] for low, high, _ in partitions[GROUPS-n1]]
        value = float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries]))
        per_bank = {bank['witness']: float(np.logaddexp2.reduce([
            row['failure_log2_screen'] for row in bank['entries'] if row['cohort_counts'][0] == n1]))
            for bank in banks}
        cohorts.append(dict(cohort_counts=[n1, GROUPS-n1, 0], entries=entries,
            failure_log2_screen=value, bank_log2_screens=per_bank, selected_cohort_reaches_40_bits=value < -40))
    return dict(construction=CONSTRUCTION, field_order=Q, forbidden_weight=LIMIT,
        active_groups=GROUPS, source_matching_counts=list(source), source_supports=source_count,
        source_classes=1, all_source_classes=81, complete_source_coverage=False,
        cohort_order=['a1', 'a2', 'a3-through-a8'], complete_sparse_cohorts_for_selected_source=True,
        cohort_count=len(cohorts), pair_swaps=True, fixed_geometry_across_columns=True,
        projective_normalizations=1, endpoint_potential_factors=1,
        witness_banks=len(banks), maximum_witness_banks=2, fitted_new_witnesses=0,
        directed=False, proposal_requires_fresh_directed_replay=True, whole_code_certificate=False,
        banks=banks, cohorts=cohorts,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in cohorts])),
        failed_cohorts_at_40_bits=[row['cohort_counts'] for row in cohorts if not row['selected_cohort_reaches_40_bits']],
        positive_cohorts=[row['cohort_counts'] for row in cohorts if row['failure_log2_screen'] > 0],
        elapsed_seconds=perf_counter()-started)


def run(*, all_sources=False, progress=None):
    """Two fixed witnesses,18 cohorts and either one or all81 source classes."""
    if type(all_sources) is not bool:
        raise ValueError('all_sources must be an explicit boolean')
    started = perf_counter()
    counts = matching_joint_counts(REGION_LENGTH, GROUPS)
    if all_sources:
        expected = {(a, b) for a in range(9) for b in range(9)}
        if (set(counts) != expected or any(type(value) is not int or value <= 0 for value in counts.values())
                or sum(counts.values()) != comb(REGION_LENGTH, GROUPS)):
            raise ArithmeticError('The complete81-class joint source partition is required')
    partitions = {n2: geometry_bins(n2) for n2 in range(GROUPS+1)}
    raw_banks = {z: pair_regions_screen(float(z), GROUPS) for _, z, _ in WITNESSES}
    sources = []
    for source in sorted(counts) if all_sources else [SOURCE]:
        result = _source_screen(source, counts[source], partitions, raw_banks,
                                progress=None if all_sources else progress)
        sources.append(result)
        if all_sources and progress:
            progress(dict(stage='source', completed=len(sources), total=len(counts),
                source_matching_counts=list(source), failure_log2_screen=result['failure_log2_screen'],
                failed_cohorts_at_40_bits=result['failed_cohorts_at_40_bits'],
                elapsed_seconds=perf_counter()-started))
    if not all_sources:
        result['elapsed_seconds'] = perf_counter()-started
        return result
    cohorts = []
    for n1 in range(GROUPS+1):
        values = [source['cohorts'][n1]['failure_log2_screen'] for source in sources]
        value = float(np.logaddexp2.reduce(values))
        cohorts.append(dict(cohort_counts=[n1, GROUPS-n1, 0], failure_log2_screen=value,
            selected_cohort_reaches_40_bits=value < -40,
            worst_source_matching_counts=sources[int(np.argmax(values))]['source_matching_counts']))
    entries = [row for source in sources for cohort in source['cohorts'] for row in cohort['entries']]
    return dict(construction=CONSTRUCTION, field_order=Q, forbidden_weight=LIMIT,
        active_groups=GROUPS, source_classes=len(sources), all_source_classes=len(counts),
        source_supports=sum(counts.values()), complete_source_coverage=True,
        cohort_order=['a1', 'a2', 'a3-through-a8'], complete_sparse_cohort_coverage=True,
        cohort_count=len(cohorts), pair_swaps=True, fixed_geometry_across_columns=True,
        projective_normalizations=1, endpoint_potential_factors=1,
        witness_banks=2, maximum_witness_banks=2, fitted_new_witnesses=0,
        source_conditioned_cap_banks=2*len(sources), raw_regional_banks=2,
        directed=False, proposal_requires_fresh_directed_replay=True, whole_code_certificate=False,
        sources=sources, cohorts=cohorts, selected_cells=len(entries),
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in cohorts])),
        failed_cohorts_at_40_bits=[row['cohort_counts'] for row in cohorts if not row['selected_cohort_reaches_40_bits']],
        positive_cohorts=[row['cohort_counts'] for row in cohorts if row['failure_log2_screen'] > 0],
        failed_cells_at_40_bits=[dict(source_matching_counts=row['source_matching_counts'],
            cohort_counts=row['cohort_counts'], mu_range=row['mu_range'],
            failure_log2_screen=row['failure_log2_screen']) for row in entries if row['failure_log2_screen'] >= -40],
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all-sources', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite frozen FLOAT artifacts')
    result = run(all_sources=args.all_sources,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('banks', 'cohorts', 'sources', 'failed_cells_at_40_bits')},
                     allow_nan=False), flush=True)
