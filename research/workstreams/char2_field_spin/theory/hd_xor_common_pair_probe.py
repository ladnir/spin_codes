"""Fixed two-mark FLOAT diagnostic for the paired-route common-pair stratum.

The only constituent geometries covered have all 17 groups on one common
unordered pair. The count28 occurs once, outside the32-column support
polynomial. Source support geometry is conditioned on its exact ordered
matching counts (a,b); each supplied regional table dominates its subsets.
These are proposals requiring fresh directed replay, not certificates.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_geometry_strata import common_pair_caps_float
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import joint_log_count
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_joint_gate import CONSTRUCTION, ordered_tables
from hd_xor_multigroup import LIMIT, Q, REGION_LENGTH


MARKS = ('996548217/1000000000', '996598217/1000000000')
LOG_POTENTIAL = '-127/2'
GROUPS = 17


def run(all_sources=False, progress=None):
    started = perf_counter()
    source_counts = matching_joint_counts(REGION_LENGTH, GROUPS)
    selected = sorted(source_counts) if all_sources else [(0, 0), (8, 8)]
    potential = float(Fraction(LOG_POTENTIAL))
    rows, summaries = [], []
    for mark in MARKS:
        mark_started = perf_counter()
        z = float(Fraction(mark))
        raw = pair_regions_screen(z, GROUPS)
        mark_rows = []
        for a, b in selected:
            row_started = perf_counter()
            caps = common_pair_caps_float(ordered_tables(raw, a, b), potential, GROUPS)
            value = (log2(source_counts[a, b])+log2(28)-log2(Q-1)
                     -LIMIT*log2(z)-potential+joint_log_count(caps, GROUPS))
            row = dict(source_matching_counts=[a, b], source_supports=source_counts[a, b],
                       constituent_pair_choices=28, maximum_pair_multiplicity=17,
                       mark=mark, log2_potential=LOG_POTENTIAL,
                       failure_log2_screen=value, fully_active_column_log2_cap=caps[0, 17],
                       elapsed_seconds=perf_counter()-row_started)
            rows.append(row)
            mark_rows.append(row)
            if progress:
                progress(row)
        summaries.append(dict(mark=mark, source_classes=len(selected),
            failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in mark_rows])),
            largest_class=max(mark_rows, key=lambda row: row['failure_log2_screen'])['source_matching_counts'],
            largest_class_log2_screen=max(row['failure_log2_screen'] for row in mark_rows),
            elapsed_seconds=perf_counter()-mark_started))
    return dict(construction=CONSTRUCTION, directed=False, proposal_requires_fresh_directed_replay=True,
        active_groups=17, active_rows_per_group=2, maximum_pair_multiplicity=17,
        constituent_geometry='one shared unordered constituent pair',
        complete_source_geometry=bool(all_sources), complete_constituent_geometry=False,
        geometry_fixed_across_columns=True, exact_ordered_matching_counts=True,
        labelled_region_order=True, cancellation_mixture_before_potential=True,
        exact_common_pair_column_law=True, individual_support_floors_relaxed=True,
        projective_normalizations=1, whole_code_certificate=False,
        field_order=Q, forbidden_weight=LIMIT, source_classes=len(selected),
        summaries=summaries, entries=rows, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all-sources', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a previous probe artifact')
    result = run(args.all_sources, progress=lambda row: print(json.dumps(row), flush=True))
    result['probe_source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(output=str(args.output.resolve()), summaries=result['summaries'],
                          elapsed_seconds=result['elapsed_seconds'])), flush=True)
