"""Recover one large term of the pinned a6 proof majorant, not of failure.

The row-floor control failed although quiet input degrees were exact. This
diagnostic restricts each column to all18 moment replicas choosing the SAME
quiet-mask branch or the SAME selected nonquiet-norm branch. Their coefficient
is the branch norm raised to18; no multinomial coefficient is needed.

For A=0..5, the quiet input masks recovered through P_j/Q_j are0, A, or
{2,3,4,5}. The nonquiet surrogate is A. Six separate original-row floors
therefore reduce exactly to two counters, for{0,1} and{2,3,4,5}. Counters
are capped at9. A max-plus DP keeps these counters and the ordered two-state
path, then backtracks one term meeting EVERY original-row floor in each
replica. It does not sum terms or rebuild a regional provider.

A positive normalized term is a LOWER bound on this proof expression only.
It can rule out fixing this expression solely by individual row-weight floors
while preserving its nonquiet full-A surrogate and coefficients. Nonquiet
input masks need not actually equal A. Exact nonquiet counts, stricter RS
relations, or tighter common operators may still improve the bound. No
construction failure or probability lower bound follows.
"""
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path

import numpy as np

import hd_lambda_double_rowfloor_control as previous
from hd_pair_identity_kernel import entry_log_matrices


DIRECTORY = Path(__file__).resolve().parent
RECEIPT = 'hd_lambda_double_rowfloor_9927_float.json'
RECEIPT_SHA = '53d11e69d6a9375816db9eccc81babbe7ea62de1589d21e1e7e1358400a736a5'
PREVIOUS_SHA = 'de46c03f2633a35de3fd12ed2f3c608d5b5001d77c3c2c2fb2401bed00c96da4'
SOURCE_MASK = 63
ROW_CLASSES = (3, 60)


def dependencies():
    if sha256((DIRECTORY/'hd_lambda_double_rowfloor_control.py').read_bytes()).hexdigest() != PREVIOUS_SHA:
        raise ValueError('Frozen row-floor controller changed')
    return {**previous.dependencies(), Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def _array(value):
    result = np.asarray([[(-np.inf if v is None else v) for v in row] for row in value], dtype=float)
    if result.shape != (2, 2) or np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Two-state log matrix with exact zero entries required')
    return result


def original_mask(physical_mask, column):
    """Undo Q then P on a pointwise-quiet support; never infer a label orbit."""
    before, after = previous.base.column_permutations(previous.base.fiber.SCHEDULE, column)
    inside = sum(1 << after[before[i]] for i in previous.SOURCE)
    if type(physical_mask) is not int or physical_mask < 0 or physical_mask & ~inside:
        raise ValueError('Quiet physical support must lie in this scheduled A')
    return sum(bool(physical_mask & (1 << after[before[i]])) << i for i in previous.SOURCE)


def mask_moment(mask, zero, nonzero, *, rows=16):
    """One ordered common-measure moment, without constructing all2^16 masks."""
    if type(mask) is not int or not 0 <= mask < 1 << rows:
        raise ValueError('Mask exceeds the explicit chronological row count')
    word = np.array(((0., -np.inf), (-np.inf, 0.)))
    for i in range(rows):
        word = previous.base.log_matmul(word, nonzero if mask & (1 << i) else zero)
    return word


def load_branches():
    raw = (DIRECTORY/RECEIPT).read_bytes()
    if sha256(raw).hexdigest() != RECEIPT_SHA:
        raise ValueError('Pinned row-floor receipt changed')
    saved = json.loads(raw)
    if (saved.get('construction') != previous.base.fiber.outer.CONSTRUCTION
            or saved.get('schedule') != previous.base.fiber.SCHEDULE
            or saved.get('original_source') != list(previous.SOURCE)
            or saved.get('actual_source_arity') != 6 or saved.get('groups') != 18
            or saved.get('source_class') != 'Fmax<=1' or saved.get('mark') != '9927/10000'
            or saved.get('complete') is not True or saved.get('directed') is not False
            or saved.get('source_pins_unchanged') is not True
            or saved.get('reviewed_source_sha256') != previous.dependencies()
            or saved.get('input_sha256') != previous.base.INPUT_SHA):
        raise ValueError('Receipt scope, completeness, or source pins differ')
    previous.load_baseline()
    zero, nonzero = entry_log_matrices(saved['common_measure'], 18)
    details = saved['column_diagnostics']
    if [d['column'] for d in details] != list(range(16)):
        raise ValueError('All sixteen ordered columns are required')
    columns = []
    for detail in details:
        j = detail['column']
        if not detail['all_zero_fiber_coordinates_protected']:
            raise ValueError('Quiet support must equal input support before P/Q')
        branches = []
        for item in detail['quiet_scores']:
            physical, count = item['physical_output_mask'], item['exact_count']
            mask = original_mask(physical, j)
            if mask not in (0, SOURCE_MASK, 60) or count != (255 if mask else 1):
                raise ValueError('Unexpected exact quiet fibre branch')
            if int(bool(mask & 1)) != item['original_row0_activity']:
                raise ValueError('Saved row activity disagrees with P/Q inversion')
            branches.append(dict(kind='zero' if not mask else 'quiet', original_mask=mask,
                physical_output_mask=physical, exact_count=count,
                norm_log2=log2(count)+mask_moment(physical, zero, nonzero)/18))
        if sum(branch['kind'] == 'zero' for branch in branches) != 1:
            raise ValueError('Exactly one zero branch is required')
        branches.append(dict(kind='nonquiet-surrogate', original_mask=SOURCE_MASK,
            physical_output_mask=None, exact_count=None,
            norm_log2=np.minimum(_array(detail['factor_norm_log2']), _array(detail['direct_norm_log2']))))
        columns.append(branches)
    return saved, columns


def best_repeated_branch_term(columns, *, groups=18, floor=9, row_classes=ROW_CLASSES):
    """Exact max-plus DP on FLOAT branch logs; choose one path, do not sum."""
    if (type(groups) is not int or not 1 <= groups <= 18
            or type(floor) is not int or not 0 <= floor <= len(columns) <= 16
            or len(row_classes) != 2 or any(type(mask) is not int or mask <= 0 for mask in row_classes)
            or row_classes[0] & row_classes[1]):
        raise ValueError('At most16 columns,18 replicas and two disjoint row classes required')
    source = row_classes[0] | row_classes[1]
    current = {(0, 0, 0): 0.}
    histories = []
    for column in columns:
        following, history = {}, {}
        for index, branch in enumerate(column):
            mask = branch['original_mask']
            if (type(mask) is not int or mask < 0 or mask & ~source
                    or any(mask & cls not in (0, cls) for cls in row_classes)):
                raise ValueError('Each branch must be constant on each row class')
            activity = tuple(int(bool(mask & cls)) for cls in row_classes)
            matrix = _array(branch['norm_log2'])
            for old_key, score in current.items():
                old_state, x, y = old_key
                counts = (min(floor, x+activity[0]), min(floor, y+activity[1]))
                for new_state in (0, 1):
                    value = score+groups*matrix[old_state, new_state]
                    key = (new_state, *counts)
                    if isfinite(value) and (key not in following or value > following[key]):
                        following[key] = float(value)
                        history[key] = (old_key, index)
        current, histories = following, histories+[history]
    terminals = [key for key in current if key[1:] == (floor, floor)]
    if not terminals:
        return dict(found=False, term_log2=None, path=[], original_row_activity_counts={})
    final = max(terminals, key=lambda key: current[key])
    term, key, path = current[final], final, []
    for j in reversed(range(len(columns))):
        old, index = histories[j][key]
        branch = columns[j][index]
        path.append(dict(column=j, branch_index=index, kind=branch['kind'],
            original_surrogate_mask=branch['original_mask'],
            physical_output_mask=branch.get('physical_output_mask'), exact_count=branch.get('exact_count'),
            state_in=old[0], state_out=key[0],
            term_log2=groups*float(_array(branch['norm_log2'])[old[0], key[0]])))
        key = old
    path.reverse()
    counts = {i: sum(bool(row['original_surrogate_mask'] & (1 << i)) for row in path)
              for i in range(source.bit_length()) if source & (1 << i)}
    if min(counts.values()) < floor or abs(sum(row['term_log2'] for row in path)-term) > 1e-9:
        raise ArithmeticError('Recovered term does not match the DP or original-row floors')
    return dict(found=True, term_log2=term, path=path, original_row_activity_counts=counts,
        initial_state=0, terminal_state=final[0], each_replica_obeys_surrogate_row_floors=True,
        identical_branch_in_all_replicas=True, multinomial_lower_factor=1,
        counter_classes=row_classes, counter_states_per_boundary=2*(floor+1)**2)


def run():
    pins = dependencies()
    saved, columns = load_branches()
    witness = best_repeated_branch_term(columns)
    normalized = None
    if witness['found']:
        normalized = saved['common_prefix_log2']+witness['term_log2']-864*log2(255)
        if normalized > saved['result']['fixed_A_repeated18_bound_log2']+1e-8:
            raise ArithmeticError('Term exceeds its pinned positive row-floor majorant')
    if dependencies() != pins or sha256((DIRECTORY/RECEIPT).read_bytes()).hexdigest() != RECEIPT_SHA:
        raise ValueError('Sources or receipt changed during the saved-only diagnostic')
    return previous.base._jsonable(dict(construction=saved['construction'], original_source=previous.SOURCE,
        schedule=saved['schedule'], groups=18, source_class=saved['source_class'],
        mark=saved['mark'], common_prefix_log2=saved['common_prefix_log2'],
        input_receipt=RECEIPT, input_receipt_sha256=RECEIPT_SHA,
        method_lower_term_log2=normalized,
        pinned_rowfloor_majorant_log2=saved['result']['fixed_A_repeated18_bound_log2'],
        gap_to_pinned_majorant_bits=(saved['result']['fixed_A_repeated18_bound_log2']-normalized
                                    if normalized is not None else None),
        witness=witness, RS_factor_whole_word='255^(-864)',
        failure_probability_lower_bound=False, construction_failure_claim=False,
        nonquiet_surrogate_is_not_actual_input_mask=True,
        lower_bound_on_proof_majorant_only=True,
        only_surrogate_row_weight_floor_refinements_tested=True,
        stronger_RS_relations_not_ruled_out=True,
        directed=False, regional_provider_runs=0, full_word_tail_runs=0,
        reviewed_source_sha256=pins, source_pins_unchanged=True))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    result = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
