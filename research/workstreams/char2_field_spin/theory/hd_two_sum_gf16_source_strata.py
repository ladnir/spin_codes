"""Bounded FLOAT control retaining the source's maximum matching count.

Fix one indexed source set G before all columns and regional placements.
Let F be the larger full-pair count in the two adjacent matchings. Every
regional subset U of G has f(U)<=F. Thus max_{f<=F} R[m,f] is a uniform
regional cap for this source stratum, without resampling G or its row sets.

The exact indexed-source count N(r,F) is charged once, outside alpha.
The three frozen marks, numeric slopes, and eight fixed type templates
are unchanged. Each stratum selects one deterministic template before
the positive union over F. No fitting or directed claim occurs here.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_gf16_emptycut_control as control
import hd_two_sum_gf16_paired_nz as paired
import hd_xor_group_geometry as geometry


DIRECTORY = Path(__file__).resolve().parent
low, cut, state = control.low, control.cut, control.low.state
MAXIMUM = 32
PINS = {
    'hd_two_sum_gf16_emptycut_control.py': 'cd1990bc0c5ae23b76de4077fd6d6a46c42fda915c40c0a8cd0234d1f6d1b8b7',
    'hd_two_sum_gf16_state_emptycut.py': '5cca0a50e5f0517a0a3066e09e50a98833c8deb6bb8e31925c4c966051f9c84f',
    'hd_two_sum_gf16_paired_nz.py': 'e6938a76e750924d03329e091a21582d5686b606e41517d43da6118c3995e048',
    'hd_xor_group_geometry.py': '880b51f316341ed4d77d5020694627805cfd708cdf6278a95e142b8ddf5aff17',
}


def load_reference():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Frozen dependency changed: {name}')
    raw = (DIRECTORY/control.PAIRED_NAME).read_bytes()
    if sha256(raw).hexdigest() != control.PAIRED_SHA:
        raise ValueError('Frozen paired regional reference changed')
    source = json.loads(raw)
    low.checked_dependencies()
    if (source.get('directed') is not False or source.get('complete_requested_anchors') is not True
            or len(source.get('banks', [])) != 3 or len(source.get('witnesses', [])) != 45):
        raise ValueError('The complete three-bank, unrefined FLOAT reference is required')
    for index, (r, _, mark) in enumerate(low.ANCHORS):
        bank, rows = source['banks'][index], source['witnesses'][15*index:15*(index+1)]
        if (bank['active_groups'], bank['mark'], bank['alpha']) != (r, str(mark), str(low.ALPHA)):
            raise ValueError('Frozen anchor mismatch')
        if any(row['log2_slopes'] != rows[0]['log2_slopes'] for row in rows):
            raise ValueError('All original numeric slopes must agree within each anchor')
    return source


def restricted_curves(regions, maximum, full_max):
    """Four log curves, with only the sound upper constraint f(U)<=F."""
    if (type(maximum) is not int or not 1 <= maximum <= MAXIMUM
            or type(full_max) is not int or not 0 <= full_max <= maximum//2):
        raise ValueError('Occupancy1..32 and a feasible matching maximum required')
    rows = []
    for m in range(maximum+1):
        values = np.asarray([regions[m, f] for f in range(min(full_max, m//2)+1)])
        if values.shape != (min(full_max, m//2)+1, 2, 2) or np.any(np.isnan(values)) or np.any(np.isposinf(values)):
            raise ValueError('Complete finite-or-zero regional histograms required')
        rows.append(np.max(values, axis=0).ravel())
    return np.asarray(rows).T


def checked_source_counts(size, groups):
    counts = geometry.matching_max_counts(size, groups, max_reduced_weight=MAXIMUM)
    if (set(counts) != set(range(groups//2+1)) or any(type(c) is not int or c <= 0 for c in counts.values())
            or sum(counts.values()) != comb(size, groups)):
        raise ArithmeticError('Matching strata must partition every indexed source exactly once')
    return counts


def evaluate_stratum(prepared, groups, full_max, source_count, mark):
    """Use N(r,F), not N(r,F)**alpha, in this full-cohort stratum."""
    if (type(groups) is not int or not 1 <= groups <= MAXIMUM
            or type(full_max) is not int or not 0 <= full_max <= groups//2
            or type(source_count) is not int or not 1 <= source_count <= comb(8192, groups)
            or not 0 < F(mark) < 1 or prepared.get('regions') != 16
            or np.shape(prepared['typed'])[1:] != (16,)):
        raise ValueError('Explicit G16 source stratum and complete fixed type laws required')
    kernel = state.kernel_from_types(prepared['classes'], prepared['path_a'], prepared['typed'], groups)
    power = state.matrix_power_log(kernel, 16)
    prefix = float(low.ALPHA)*(-low.fractional.LIMIT*log2(float(mark))-log2(255))
    terms = (prepared['path_a']+groups/16*state.logsum(16*prepared['typed'], axis=1)
             +np.log2([row[3] for row in prepared['classes']]))
    candidates = [i for i, row in enumerate(prepared['classes']) if row[:2] == (0, 0)]
    dominant = max(candidates, key=lambda i: terms[i])
    return dict(active_groups=groups, source_full_pair_max=full_max, source_supports=source_count,
        source_count_outside_alpha=True, mark=str(mark), alpha=str(low.ALPHA),
        failure_log2_screen=float(prefix+log2(source_count)+state.logsum(power[0])),
        log2_nonempty_00_intercept=prepared['log_c'], log2_slopes=prepared['slopes'].tolist(),
        log2_intercepts=prepared['intercepts'].tolist(), envelope_minimum_slack=float(np.min(prepared['slack'])),
        dominant_00_path_transition_counts=list(prepared['classes'][dominant][2]),
        dominant_00_path_allowed_counts=prepared['allowed_counts'][dominant].tolist(),
        dominant_type=int(np.argmax(prepared['typed'][dominant]))+1,
        dominant_00_column_log_term=float(terms[dominant]),
        all_cohorts=True, all_constituent_row_sets=True, fixed_source_across_all_columns=True,
        fixed_types_across_all_columns=True, exact_total_arity_restriction=False,
        directed=False, whole_code_certificate=False, requires_fresh_directed_replay=True)


def saved_matrix(matrix):
    return [None if np.isneginf(x) else float(x) for x in np.asarray(matrix).ravel()]


def run(max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not np.isfinite(max_seconds) or not 0 < max_seconds <= 30:
        raise ValueError('An explicit at-most30s budget is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Fixed-source stratum control exhausted its budget')
    result = dict(construction=low.fractional.CONSTRUCTION, field_order=256,
        forbidden_weight=low.fractional.LIMIT, source_size=8192,
        anchors=[dict(active_groups=r, mark=str(mark), alpha=str(low.ALPHA)) for r, _, mark in low.ANCHORS],
        source_partition='F=max(full even-matching pairs,full odd-matching pairs)',
        source_count_outside_alpha=True, fixed_source_across_all_columns=True,
        all_source_supports_per_anchor=True, all_cohorts_per_anchor=True,
        all_constituent_row_sets_per_anchor=True, fixed_type_weights_across_paths_and_columns=True,
        regional_subset_constraint='0<=f<=min(F,floor(m/2))', regional_maximum=MAXIMUM,
        frozen_reference=control.PAIRED_NAME, frozen_reference_sha256=control.PAIRED_SHA,
        dependency_pins=dict(PINS), provider=paired.provider_metadata(),
        unchanged_numeric_slopes=True, unchanged_eight_type_templates=True, fitted_parameters=0,
        templates_per_stratum=8, requested_strata=sum(r//2+1 for r, _, _ in low.ANCHORS),
        exact_total_arity_restriction=False, continuous_weight_band_covered=False,
        directed=False, whole_code_certificate=False, requires_fresh_directed_replay=True,
        compute_budget_seconds=max_seconds, regional_banks_rebuilt=0,
        complete_requested_anchors=False, banks=[], witnesses=[], strata=[], selected=[])
    try:
        check_time()
        reference = load_reference()
        for index, (r, _, mark) in enumerate(low.ANCHORS):
            check_time()
            regions = paired.paired_nz_regions(mark, low.ALPHA, MAXIMUM)
            result['regional_banks_rebuilt'] += 1
            required = {(m, f) for m in range(MAXIMUM+1) for f in range(m//2+1)}
            if set(regions) != required:
                raise ArithmeticError('Fresh bank must retain every R[m,f] through32')
            curves_all = state.regional_curves(regions, MAXIMUM)
            old = np.asarray([[(-np.inf if x is None else x) for x in row]
                              for row in reference['banks'][index]['curves']])
            if not np.allclose(curves_all, old, atol=1e-8, rtol=0):
                raise ArithmeticError('Fresh unrestricted FLOAT curves changed from the pinned provider')
            slopes = np.asarray(reference['witnesses'][15*index]['log2_slopes'])
            templates = list(control.type_templates(r, mark))
            if len(templates) != 8:
                raise ArithmeticError('Exactly the eight frozen type templates are required')
            counts = checked_source_counts(8192, r)
            result['banks'].append(dict(active_groups=r, mark=str(mark), alpha=str(low.ALPHA),
                log2_slopes=slopes.tolist(), maximum=MAXIMUM,
                regions=[dict(occupancy=m, full_pairs=f, log2_matrix=saved_matrix(regions[m, f]))
                         for m, f in sorted(regions)], zero_entries_are_none=True,
                source_strata=[dict(full_max=k, count=c) for k, c in counts.items()],
                type_templates=[metadata for metadata, _ in templates]))
            anchor_strata = []
            for full_max, count in counts.items():
                check_time()
                curves = restricted_curves(regions, r, full_max)
                if not np.all(curves <= curves_all[:, :r+1]):
                    raise ArithmeticError('A source restriction cannot enlarge a regional cap')
                indices = []
                for template_index, (metadata, weights) in enumerate(templates):
                    check_time()
                    prepared = cut.prepared_envelope(curves, slopes, weights)
                    row = evaluate_stratum(prepared, r, full_max, count, mark)
                    row.update(template=metadata, template_index=template_index)
                    indices.append(len(result['witnesses']))
                    result['witnesses'].append(row)
                winner = min(indices, key=lambda i: result['witnesses'][i]['failure_log2_screen'])
                selected = dict(result['witnesses'][winner], witness_index=winner)
                result['strata'].append(selected)
                anchor_strata.append(selected)
                if progress:
                    progress(dict(active_groups=r, source_full_pair_max=full_max,
                        failure_log2_screen=selected['failure_log2_screen'],
                        template=selected['template']['name'], elapsed_seconds=perf_counter()-started))
            result['selected'].append(dict(active_groups=r,
                failure_log2_screen=float(state.logsum([row['failure_log2_screen'] for row in anchor_strata])),
                exact_source_support_count=sum(counts.values()), complete_source_partition=True,
                selected_strata=len(anchor_strata), deterministic_template_per_stratum=True,
                extra_template_factor=False))
        check_time()
        result.update(status='completed', complete_requested_anchors=True,
            failure_log2_screen=float(state.logsum([row['failure_log2_screen'] for row in result['selected']])))
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen artifact')
    result = run(args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('banks', 'witnesses', 'strata')}))
