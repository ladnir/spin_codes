"""Six saved-bank FLOAT controls, ONLY the F=1 indexed-source stratum.

Weights18/24/32, old-zero and joint-empty-floor1, frozen numeric slopes.
No regional provider or optimizer is called. The source count replaces
C(8192,r) explicitly, outside alpha. A30-second cap preserves completed
rows but never presents a partial run as a complete requested control.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_gf16_state_cover as cover
import hd_two_sum_gf16_low_arity as low
import hd_two_sum_input_weight as joint


DIRECTORY = Path(__file__).resolve().parent
BANK_NAME = 'hd_two_sum_gf16_source_strata_anchors_float.json'
BANK_SHA = '5050033694bb6a1ef7852aa309ee3abca252c9740fc713a83aa8945a380fb1cf'
TEMPLATES = ('old-all-zero', 'joint-empty-floor1')
FULL_MAX = 1
PINS = {
    'hd_two_sum_gf16_state_cover.py': '789ab7ab0ad1f4ea4b3e2874a2c7393ff00aac4f016a77bdeada3cc5c8dfbda5',
    'hd_two_sum_gf16_state_emptycut.py': '5cca0a50e5f0517a0a3066e09e50a98833c8deb6bb8e31925c4c966051f9c84f',
    'hd_two_sum_input_weight.py': '5b13c75a6bcc1a47308dd634ed6b3382f67d1a009f6c7ed419881f21aeb7d0cc',
    'hd_two_sum_gf16_low_arity.py': '943dd058318f0f4c5d5a371fa3a84210646e71ca8c01327399d9a60ce7b9064a',
}


def load_source():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Checked coverage dependency changed: {name}')
    raw = (DIRECTORY/BANK_NAME).read_bytes()
    if sha256(raw).hexdigest() != BANK_SHA:
        raise ValueError('Frozen source-stratum regional bank changed')
    low.checked_dependencies()
    source = json.loads(raw)
    if (source.get('construction') != low.fractional.CONSTRUCTION
            or source.get('field_order') != 256 or source.get('source_size') != 8192
            or source.get('forbidden_weight') != low.fractional.LIMIT
            or source.get('directed') is not False
            or source.get('complete_requested_anchors') is not True
            or len(source.get('banks', [])) != 3):
        raise ValueError('The complete three-anchor FLOAT source-stratum artifact is required')
    for bank, (r, _, mark) in zip(source['banks'], low.ANCHORS):
        if (bank.get('active_groups'), bank.get('mark'), bank.get('alpha'), bank.get('maximum')) != (
                r, str(mark), str(low.ALPHA), 32):
            raise ValueError('Frozen anchor metadata changed')
        counts = {row['full_max']: row['count'] for row in bank['source_strata']}
        if (len(counts) != len(bank['source_strata']) or set(counts) != set(range(r//2+1))
                or any(type(n) is not int or n <= 0 for n in counts.values())
                or sum(counts.values()) != comb(8192, r)):
            raise ValueError('Pinned source counts must form the full indexed partition')
        if np.shape(bank['log2_slopes']) != (4,) or np.any(~np.isfinite(bank['log2_slopes'])):
            raise ValueError('Four frozen numeric slopes required')
        keys = [(row['occupancy'], row['full_pairs']) for row in bank['regions']]
        required = {(m, f) for m in range(33) for f in range(m//2+1)}
        if len(keys) != len(required) or set(keys) != required:
            raise ValueError('Complete saved R[m,f] tables through32 required')
    return source


def restricted_curves(bank):
    """Use only the retained f<=1 tables, never reconstruct a regional cap."""
    regions = {}
    for row in bank['regions']:
        value = np.array([-np.inf if x is None else x for x in row['log2_matrix']])
        if value.shape != (4,) or np.any(np.isnan(value)) or np.any(np.isposinf(value)):
            raise ValueError('Four finite-or-zero saved entries required')
        regions[row['occupancy'], row['full_pairs']] = value
    r = bank['active_groups']
    return np.stack([np.max([regions[m, f] for f in range(min(FULL_MAX, m//2)+1)], axis=0)
                     for m in range(r+1)], axis=1)


def fixed_templates(bank):
    """Reconstruct exactly two stored marker templates; do not recompute them."""
    templates = {row['name']: row for row in bank['type_templates']}
    if len(templates) != len(bank['type_templates']) or not set(TEMPLATES) <= set(templates):
        raise ValueError('The two specified frozen templates must occur exactly once')
    for name in TEMPLATES:
        metadata = templates[name]
        markers = [F(value) for value in metadata['markers']]
        if len(markers) != 16 or any(x < 0 for x in markers):
            raise ValueError('Exactly16 fixed nonnegative marker logs required')
        if name == 'old-all-zero':
            if any(markers) or metadata['family'] != 'nonempty-column-floor':
                raise ValueError('Old-zero marker template changed')
            dense = low.fractional.fractional_dense_weights(low.ALPHA)
            weights = cover.emptycut.fractional_cardinality_weights(low.ALPHA, dense, markers)
        else:
            if metadata['family'] != 'pre-H-total-support-floor':
                raise ValueError('Joint input-weight marker template changed')
            cache = {value: joint.type_weights(low.ALPHA, float(value)) for value in set(markers)}
            weights = {a: cache[markers[a-1]][a] for a in range(1, 17)}
        yield metadata.copy(), weights


def restrict_count(row, r, source_count):
    """Replace, do not multiply or alpha-power, the indexed-source count."""
    if type(source_count) is not int or not 0 < source_count <= comb(8192, r):
        raise ValueError('Positive exact source-stratum count required')
    result = dict(row)
    result['failure_log2_screen'] += log2(source_count)-log2(comb(8192, r))
    result.update(active_groups=r, source_counts='N(r,F=1)', source_supports=source_count,
        source_full_pair_max=FULL_MAX, source_count_outside_alpha=True,
        source_stratum_only=True, source_partition_complete=False,
        all_source_supports=False, all_source_supports_in_declared_stratum=True,
        fixed_source_across_all_columns=True)
    return result


def run(max_seconds=30, progress=None):
    if (type(max_seconds) not in (int, float) or not np.isfinite(max_seconds)
            or not 0 < max_seconds <= 30):
        raise ValueError('An explicit at-most30s control budget is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Six fixed-source coverage controls exhausted their budget')
    result = dict(construction=low.fractional.CONSTRUCTION, field_order=256,
        forbidden_weight=low.fractional.LIMIT, source_size=8192,
        source_bank=BANK_NAME, source_bank_sha256=BANK_SHA, dependency_pins=dict(PINS),
        source_full_pair_max=FULL_MAX, source_stratum_only=True,
        all_source_supports=False, source_partition_complete=False,
        all_source_supports_in_declared_stratum=True, source_count_outside_alpha=True,
        all_cohorts_per_anchor=True, all_constituent_row_sets_per_anchor=True,
        regional_banks_rebuilt=0, fitted_parameters=0, unchanged_numeric_slopes=True,
        unchanged_literal_markers=True, requested_controls=6, templates=list(TEMPLATES),
        exact_arity_restriction_used=False, continuous_weight_band_covered=False,
        directed=False, whole_code_certificate=False, requires_fresh_directed_replay=True,
        compute_budget_seconds=max_seconds, witnesses=[], complete_requested_controls=False)
    try:
        check_time()
        source = load_source()
        for bank, (r, _, mark) in zip(source['banks'], low.ANCHORS):
            check_time()
            curves = restricted_curves(bank)
            count = next(row['count'] for row in bank['source_strata'] if row['full_max'] == FULL_MAX)
            slopes = np.asarray(bank['log2_slopes'])
            for metadata, weights in fixed_templates(bank):
                check_time()
                result['in_progress'] = dict(active_groups=r, template=metadata['name'])
                prepared = cover.prepared_envelope(curves, slopes, weights)
                old = restrict_count(cover.emptycut.evaluate_prepared(prepared, r, mark, low.ALPHA), r, count)
                reference = [row for row in source['witnesses'] if row['active_groups'] == r
                    and row['source_full_pair_max'] == FULL_MAX and row['template']['name'] == metadata['name']]
                if len(reference) != 1 or abs(old['failure_log2_screen']-reference[0]['failure_log2_screen']) > 1e-7:
                    raise ArithmeticError('Frozen F=1 empty-cut control failed its replay identity')
                row = restrict_count(cover.evaluate_prepared(prepared, r, mark, low.ALPHA,
                    check_time=check_time), r, count)
                if row['failure_log2_screen'] > old['failure_log2_screen']+1e-7:
                    raise ArithmeticError('A pathwise minimum cannot enlarge the previous bound')
                row.update(template=metadata, log2_slopes=slopes.tolist(),
                    previous_emptycut_failure_log2_screen=old['failure_log2_screen'],
                    improvement_bits=old['failure_log2_screen']-row['failure_log2_screen'])
                result['witnesses'].append(row)
                result.pop('in_progress', None)
                if progress:
                    progress(dict(active_groups=r, source_full_pair_max=FULL_MAX,
                        template=metadata['name'], failure_log2_screen=row['failure_log2_screen'],
                        improvement_bits=row['improvement_bits'], elapsed_seconds=perf_counter()-started))
        check_time()
        selected = [min((row for row in result['witnesses'] if row['active_groups'] == r),
                        key=lambda row: row['failure_log2_screen']) for r, _, _ in low.ANCHORS]
        result.update(status='completed', complete_requested_controls=True, selected=selected,
            selected_three_strata_log2_screen=float(low.state.logsum([row['failure_log2_screen'] for row in selected])))
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
        parser.error('Do not overwrite frozen research artifacts')
    result = run(args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('witnesses', 'selected')}))
