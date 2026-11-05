"""One inherited-witness FLOAT screen of all sparse r17 affine source classes.

The new u=8,...,15 two-sum outer has no pair-multiplicity strata.  This
wrapper reuses one raw regional bank and one category-cap bank per ordered
source class across all18 sparse cohorts.  The regional matrices still
envelope the feasible matching counts of each regional subset; only the
COLUMN occupancy resource maximum has been removed by the exact joint
hole mixture.  Fixed constituent row sets persist through all32 columns.

This is a bounded proposal, not a directed or whole-code certificate.
The previously frozen affine helper is imported without modification.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_affine_outer as affine
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import ordered_region_tables, sparse_log_count
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


GROUPS, LENGTH = 17, 8192
MARK, LOG_POTENTIAL = F('996248217/1000000000'), F('-127/2')
MAX_SECONDS = 120
AFFINE_SHA256 = '52d9c3ad412e79d4ac424270734490d3831fc4bbab02e238fedf7ce24fa3d2af'


def checked_source_counts():
    counts = matching_joint_counts(LENGTH, GROUPS)
    expected = {(a, b) for a in range(9) for b in range(9)}
    if (set(counts) != expected
            or any(type(count) is not int or count <= 0 for count in counts.values())
            or sum(counts.values()) != comb(LENGTH, GROUPS)):
        raise ArithmeticError('Exactly81 positive source classes must partition C(8192,17)')
    return counts


def screen(*, budget_seconds=MAX_SECONDS, progress=None):
    if (isinstance(budget_seconds, bool) or not isinstance(budget_seconds, (int, float))
            or not 0 < budget_seconds <= MAX_SECONDS):
        raise ValueError('The total screen budget must be positive and at most120 seconds')
    started = perf_counter()
    deadline = started+budget_seconds

    def check_budget():
        if perf_counter() >= deadline:
            raise TimeoutError('The fixed total FLOAT screen budget was exhausted')

    counts = checked_source_counts()
    sources, attempted, raw_banks = [], 0, 0
    timed_out = False
    logp = float(LOG_POTENTIAL)
    common = -log2(255)-LIMIT*log2(float(MARK))-logp
    try:
        check_budget()
        raw = pair_regions_screen(float(MARK), GROUPS)
        raw_banks = 1
        check_budget()
        for source in sorted(counts):
            check_budget()
            attempted += 1
            source_started = perf_counter()
            regions = ordered_region_tables(raw, GROUPS, *source)
            caps = {}
            for single, double in categories(GROUPS):
                check_budget()
                vector = affine.column_action_float(regions, logp, single, double)
                caps[single, double] = float(np.max(vector-[0., logp]))
            entries = []
            for n1 in range(GROUPS+1):
                check_budget()
                n2 = GROUPS-n1
                value = (log2(counts[source])+log2(comb(GROUPS, n1))
                         +3*n1+n2*log2(28)+common+sparse_log_count(caps, n1, n2))
                entries.append(dict(cohort_counts=[n1, n2, 0],
                    source_supports=counts[source], source_colorings=comb(GROUPS, n1),
                    singleton_row_choices=8**n1, pair_row_choices=28**n2,
                    sparse_support_floor=17*(n1+2*n2),
                    sparse_normalizer_exponent=16*(n1+2*n2),
                    failure_log2_screen=float(value)))
            check_budget()
            result = dict(source_matching_counts=list(source), source_supports=counts[source],
                entries=entries,
                failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])),
                elapsed_seconds=perf_counter()-source_started)
            sources.append(result)
            if progress:
                progress(dict(stage='source', completed=len(sources), total=len(counts),
                    source_matching_counts=list(source), failure_log2_screen=result['failure_log2_screen'],
                    worst_cohort=entries[int(np.argmax([row['failure_log2_screen'] for row in entries]))]['cohort_counts'],
                    elapsed_seconds=perf_counter()-started))
    except TimeoutError:
        timed_out = True

    complete = len(sources) == len(counts) and not timed_out
    cohorts = []
    if complete:
        for n1 in range(GROUPS+1):
            values = [source['entries'][n1]['failure_log2_screen'] for source in sources]
            cohorts.append(dict(cohort_counts=[n1, GROUPS-n1, 0],
                failure_log2_screen=float(np.logaddexp2.reduce(values)),
                worst_source_matching_counts=sources[int(np.argmax(values))]['source_matching_counts']))
    return dict(construction=affine.CONSTRUCTION, labels=list(affine.LABELS),
        field_order=256, field_polynomial='0x11b', forbidden_weight=LIMIT,
        active_groups=GROUPS, total_columns=32,
        witness_label='inherited', mark=str(MARK), log2_potential=str(LOG_POTENTIAL),
        witness_banks=1, fitted_new_witnesses=0,
        source_classes=len(sources), expected_source_classes=len(counts),
        source_count_partition=comb(LENGTH, GROUPS),
        completed_source_supports=sum(row['source_supports'] for row in sources),
        source_conditioned_cap_banks_attempted=attempted, raw_regional_banks=raw_banks,
        cohort_order=['a1', 'a2', 'a3-through-a8'], cohort_count=18,
        complete_source_coverage=complete, complete_sparse_cohort_coverage=complete,
        selected_cells=18*len(sources), pair_swaps=True,
        fixed_source_geometry=True, fixed_rows_across_columns=True,
        local_mask_mixture_exact=True, regional_matrices_are_upper_bounds=True,
        regional_matching_count_envelopes_used=True,
        column_occupancy_resource_maximum_used=False, pair_multiplicity_strata_used=False,
        potential_quotient_after_hole_mixture=True, ordered_regional_matrices=True,
        projective_normalizations=1, endpoint_potential_factors=1,
        directed=False, old_certificates_transferred=False, whole_code_certificate=False,
        proposal_requires_fresh_directed_replay=True,
        budget_seconds=budget_seconds, budget_exhausted=timed_out,
        elapsed_seconds=perf_counter()-started, sources=sources, cohorts=cohorts,
        failure_log2_screen=(float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in cohorts]))
                             if complete else None),
        completed_sources_log2_screen=(float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in sources]))
                                       if sources else None),
        failed_cohorts_at_40_bits=([row['cohort_counts'] for row in cohorts if row['failure_log2_screen'] >= -40]
                                  if complete else None),
        failed_cells_at_40_bits=[dict(source_matching_counts=source['source_matching_counts'],
                                     cohort_counts=row['cohort_counts'], failure_log2_screen=row['failure_log2_screen'])
            for source in sources for row in source['entries'] if row['failure_log2_screen'] >= -40])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite an existing research artifact')
    if sha256(Path(affine.__file__).read_bytes()).hexdigest() != AFFINE_SHA256:
        parser.error('The frozen affine helper pin changed')
    result = screen(progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    dependencies = (Path(__file__), Path(affine.__file__),
        Path(__file__).with_name('hd_xor_pair_route.py'),
        Path(__file__).with_name('hd_xor_mixed_gate.py'),
        Path(__file__).with_name('hd_xor_group_geometry.py'))
    result['source_sha256'] = {path.name: sha256(path.read_bytes()).hexdigest() for path in dependencies}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('sources', 'cohorts', 'failed_cells_at_40_bits')}, allow_nan=False), flush=True)
