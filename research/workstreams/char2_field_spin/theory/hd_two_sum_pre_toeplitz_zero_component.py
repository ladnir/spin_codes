"""Frozen all00 component diagnostic; NOT a probability lower bound.

Restore the positive empty00 split inside the already-selected comparator.
One subset-zeta per old/new law suffices because all00 is matching-invariant.
Other state paths are omitted from this component, not bounded by it.
Also report the largest single state-word component of the ordinary cap.
No new counts, regional laws, parameters, or optimizer are constructed.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_pre_toeplitz_schedule as schedule


DIRECTORY = Path(__file__).resolve().parent
SOURCE = 'hd_two_sum_pre_toeplitz_schedule_gate_float.json'
SOURCE_SHA = '64b5aafb91ca7d2b11584ad9896836d9bd3e3239b5c4a36c82bb8921688d7e35'
SCHEDULE_SHA = '2c40606f0efc615a9fc8003d839a0bcb61c4050bfb3e1a6188d2d7ddf6199305'


def subset_zeta_log(log_weights):
    """At output U, return log2 sum_{T subset U} weight(T), positively."""
    values = np.asarray(log_weights, float).copy()
    size = len(values)
    group = size.bit_length()-1
    if (values.shape != (size,) or not 1 <= group <= 16 or size != 1 << group
            or np.any(np.isnan(values)) or np.any(np.isposinf(values))):
        raise ValueError('A complete finite-or-zero binary mask table required')
    for bit in range(group):
        view = values.reshape(-1, 2, 1 << bit)
        view[:, 1, :] = np.logaddexp2(view[:, 1, :], view[:, 0, :])
    return values


def cardinalities(group):
    index = np.arange(1 << group, dtype=np.uint32)
    counts = np.zeros(len(index), dtype=np.uint8)
    for bit in range(group):
        counts += ((index >> bit) & 1).astype(np.uint8)
    return counts


def log_complement(log_c):
    if np.isnan(log_c) or log_c > 0:
        raise ValueError('A probability c in[0,1] required')
    if np.isneginf(log_c):
        return 0.
    if log_c == 0:
        return -np.inf
    return float(np.log(-np.expm1(log_c*np.log(2.)))/np.log(2.))


def count_times(counts, value):
    # Define 0*log0=0, preserving exact zero-probability cut endpoints.
    return np.where(np.asarray(counts) == 0, 0., value) if np.isneginf(value) else counts*value


def prepare_zero_component(mask_weights, log_b):
    weights = np.asarray(mask_weights, float)
    group = len(weights).bit_length()-1
    if not np.isfinite(log_b) or weights.shape != (1 << group,) or not 1 <= group <= 16:
        raise ValueError('Complete weights and one fixed finite slope required')
    counts = cardinalities(group)
    empty = float(weights[0])
    nonempty = weights+counts*log_b
    nonempty[0] = -np.inf
    return dict(group=group, empty=empty, counts=counts, nonempty_zeta=subset_zeta_log(nonempty))


def zero_component(prepared, log_c, groups, log_activity=0., *, columns=16, minimum=9, ordinary_intercept=0.):
    """One all00 column term under refined/ordinary affine majorants."""
    if (type(groups) is not int or groups < 1 or not np.isfinite(log_activity) or log_activity < 0
            or not np.isfinite(ordinary_intercept) or not 1 <= minimum <= columns):
        raise ValueError('Fixed positive group count, activity mark, and floor required')
    group, counts = prepared['group'], prepared['counts']
    complement = log_complement(log_c)
    log_g = np.logaddexp2(prepared['empty'], prepared['nonempty_zeta']+log_activity)-minimum*log_activity/columns
    cut_weights = count_times(counts, log_c)+count_times(group-counts, complement)
    if not np.isclose(np.logaddexp2.reduce(cut_weights), 0., atol=1e-12, rtol=0):
        raise ArithmeticError('The positive cut weights must sum to1')
    terms = cut_weights+groups*log_g
    largest = int(np.argmax(terms))
    return dict(refined_column_log2=float(np.logaddexp2.reduce(terms)),
                ordinary_column_log2=float(group*ordinary_intercept+groups*log_g[-1]),
                dominant_allowed_mask=largest, dominant_allowed_size=int(counts[largest]),
                largest_cut_term_log2=float(terms[largest]))


def best_paths(nonempty, empty, intercepts, groups, log_activity=0., *, columns=16, minimum=9):
    """Maximum positive ordinary path component for each ordered endpoint."""
    group = np.ndim(nonempty)-1
    initial, final, counts = schedule.path_geometry(group)
    terms = (counts @ intercepts+groups*(np.logaddexp2(empty, np.asarray(nonempty).ravel()+log_activity)
                                        -minimum*log_activity/columns))
    matrices, indices = np.empty((2, 2)), np.empty((2, 2), dtype=int)
    for left in range(2):
        for right in range(2):
            selected = np.flatnonzero((initial == left) & (final == right))
            winner = int(selected[np.argmax(terms[selected])])
            matrices[left, right], indices[left, right] = terms[winner], winner
    return matrices, indices


def largest_state_word(columns):
    """Max-plus product of per-column path terms, starting in state0."""
    values, traces = np.array([0., -np.inf]), [[], []]
    for matrix, indices in columns:
        candidates = values[:, None]+matrix
        parents = np.argmax(candidates, axis=0)
        values = candidates[parents, np.arange(2)]
        traces = [traces[int(parents[final])]+[int(indices[parents[final], final])] for final in range(2)]
    final = int(np.argmax(values))
    return float(values[final]), traces[final], final


def load_saved():
    if sha256(Path(schedule.__file__).read_bytes()).hexdigest() != SCHEDULE_SHA:
        raise ValueError('Reviewed schedule source changed')
    raw = (DIRECTORY/SOURCE).read_bytes()
    if sha256(raw).hexdigest() != SOURCE_SHA:
        raise ValueError('Saved fixed-schedule artifact changed')
    source = json.loads(raw)
    if (source.get('construction') != schedule.CONSTRUCTION or source.get('active_groups') != 18
            or source.get('source_full_pair_max') != 1 or source.get('homogeneous_constituent_arity') != 8
            or source.get('mark') != str(schedule.MARK) or source.get('alpha') != str(schedule.ALPHA)
            or source.get('status') != 'completed' or source.get('complete_requested_controls') is not True
            or source.get('completed_matching_columns') != 16 or source.get('directed') is not False):
        raise ValueError('Complete frozen single-class FLOAT schedule required')
    weights = {name: np.array([[-np.inf if value is None else value for value in row] for row in table], float)
               for name, table in source['log2_matching_weights'].items()}
    if set(weights) != {'old', 'new'} or any(table.shape != (17, 9) for table in weights.values()):
        raise ValueError('Both saved incidence tables are required')
    return source, weights, schedule.load_frozen_anchor()


def run(max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not np.isfinite(max_seconds) or not 0 < max_seconds <= 30:
        raise ValueError('An explicit at-most30s budget is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Frozen component diagnostic exhausted its budget')
    result = dict(diagnostic='positive-components-of-the-selected-majorant',
        probability_lower_bound=False, whole_probability_upper_bound=False,
        whole_code_certificate=False, directed=False, active_groups=18,
        homogeneous_constituent_arity=8, source_full_pair_max=1,
        all_cohorts=False, regional_banks_built=0, new_mask_contractions=0, witness_fits=0,
        old_new_weight_tables_from_frozen_artifact=True, refined_all00_matching_invariant=True,
        source_artifact=SOURCE, source_artifact_sha256=SOURCE_SHA,
        compute_budget_seconds=max_seconds, complete_requested_controls=False,
        completed_matching_columns=0, components=[], ordinary_largest_words=[])
    try:
        source, weights, anchor = load_saved()
        prefix = log2(anchor['count'])-float(schedule.ALPHA)*(schedule.LIMIT*log2(float(schedule.MARK))+log2(255))
        log_c = float(np.max(anchor['curves'][0, 1:]-anchor['slopes'][0]*np.arange(1, 19)))
        if log_c > 0:
            raise ValueError('The fixed00 slope does not admit the positive empty split')
        result.update(log2_nonempty00_intercept=log_c, log2_complement=log_complement(log_c),
                      source_supports=anchor['count'], source_count_outside_alpha=True,
                      mark=str(schedule.MARK), alpha=str(schedule.ALPHA))
        matching = schedule.xor_schedule()[0]
        for mode in ('old', 'new'):
            check_time()
            table = schedule.mask_weight_table(weights[mode], matching)
            prepared = prepare_zero_component(table, anchor['slopes'][0])
            for name, marker in anchor['markers']:
                value = zero_component(prepared, log_c, 18, float(marker), ordinary_intercept=anchor['intercepts'][0])
                value.update(mode=mode, template=name, log2_activity=str(marker),
                    refined_all00_word_component_log2=prefix+16*value['refined_column_log2'],
                    ordinary_all00_word_component_log2=prefix+16*value['ordinary_column_log2'])
                result['components'].append(value)
            if progress:
                progress(dict(stage='all00-zeta', mode=mode, elapsed_seconds=perf_counter()-started))
        cache, collected = {}, {(mode, name): [] for mode in ('old', 'new') for name, _ in anchor['markers']}
        for index, matching in enumerate(schedule.xor_schedule()):
            check_time()
            for mode in ('old', 'new'):
                key = (mode, matching if mode == 'new' else ())
                if key not in cache:
                    table = schedule.mask_weight_table(weights[mode], matching)
                    empty, table[0] = table[0], -np.inf
                    nonempty = schedule.path_transform(table, anchor['slopes'].reshape(2, 2))
                    cache[key] = {name: best_paths(nonempty, empty, anchor['intercepts'], 18, float(marker))
                                  for name, marker in anchor['markers']}
                for name, _ in anchor['markers']:
                    collected[mode, name].append(cache[key][name])
            result['completed_matching_columns'] = index+1
        for mode in ('old', 'new'):
            for name, marker in anchor['markers']:
                value, traces, final = largest_state_word(collected[mode, name])
                _, _, counts = schedule.path_geometry(16)
                result['ordinary_largest_words'].append(dict(mode=mode, template=name,
                    log2_activity=str(marker), component_log2=prefix+value, terminal_state=final,
                    paths=[format(path, '017b') for path in traces],
                    transition_counts=[counts[path].tolist() for path in traces]))
        check_time()
        result.update(status='completed', complete_requested_controls=True)
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
        parser.error('Output must be a new artifact')
    result = run(args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8', newline='\n') as destination:
        destination.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(args.output), sha256=sha256(args.output.read_bytes()).hexdigest(),
                         status=result['status'], components=result['components'], elapsed_seconds=result['elapsed_seconds'])))
