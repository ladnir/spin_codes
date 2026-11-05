"""One coupled five-cone FLOAT gate with a log-domain potential LP.

The scalar moment law is the pinned old F1 R00 law, not a fitted new
shape. A finite edge-budget LP proposes a common positive potential.
Every complete weighted row sum is recomputed independently afterward;
its largest surplus supplies lambda**256. No LP residual is a proof.
Full fixed-A Holder weights, whole-word floor1296 and once-only endpoint
normalization are inherited algebraically, not as old-construction receipts.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import linprog

from hd_two_sum_gf16_mean_column import column_polynomial, raw_holder_weights
from hd_two_sum_late_alpha_quiet import log_power
from hd_two_sum_matrix_moment import ARTIFACT as MOMENT_ARTIFACT, ARTIFACT_SHA as MOMENT_SHA, load_saved
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
STATES = ('Z', 'L', 'R', 'B', 'P')
GROUPS, ARITY, COLUMNS, FLOOR = 18, 8, 16, 1296
MARK = Fraction(691, 700)
LOG_PROBABILITY = -24.529436473994224
LOG_ATOM = -8.907630419774677
LOG_BOUND = 512.
SOURCE_PINS = {
    'hd_coupled_state_packet.py': 'faab3a80e582f0c1864b79f63b256788d99a5034f15bfbd9e781940b43bc2a43',
    'hd_coupled_state_regions.py': 'd843965b9e74d1aa86b80d09c7ab4e3a0afa594dc2c91057a7fa7497c97fa01f',
    'hd_two_sum_gf16_mean_column.py': 'fa6377002fa52fbfcb980ddd4b72d30df16b6c12c7c6db24eed54554e61a1593',
    'hd_two_sum_matrix_moment.py': 'bc192ebdafc9abec59f7977e660c5d4741737eb62cf479c920fdbb83a060ce67',
}


def logs(values):
    values = np.asarray(values, dtype=float)
    if np.any(np.isnan(values)) or np.any(np.isposinf(values)):
        raise ValueError('Finite-or-zero logarithms required')
    return values


def checked_bank(bank, maximum, full_max=1, dimension=5):
    if (type(maximum) is not int or maximum < 1 or type(full_max) is not int
            or full_max < 0 or type(dimension) is not int or dimension < 1):
        raise ValueError('Positive finite coefficient scope required')
    expected = {(m, f) for m in range(maximum+1)
                for f in range(min(full_max, m//2)+1)}
    if set(bank) != expected:
        raise ValueError('Complete and exact requested histogram scope required')
    result = {}
    for key in sorted(expected):
        value = logs(bank[key])
        if value.shape != (dimension, dimension):
            raise ValueError('Every regional matrix must use the declared state dimension')
        result[key] = value
    return result


def shape(maximum):
    result = LOG_PROBABILITY+LOG_ATOM*np.arange(maximum+1)
    result[0] = 0.
    return result


def edge_budget_program(bank, log_moments, *, log_bound=LOG_BOUND):
    """Build a small LP in LOG potentials, so tiny probabilities are never LP inputs.

    Variables are u0..u(d-1),t,s, with u0=0. Each finite edge in a row
    with k nonzero entries gets at most 1/k of that row's allowed budget.
    This is sufficient, but deliberately more conservative than logsumexp.
    """
    log_moments = logs(log_moments)
    if (log_moments.ndim != 1 or len(log_moments) < 2
            or np.any(~np.isfinite(log_moments)) or log_moments[0] != 0.
            or not np.isfinite(log_bound) or log_bound <= 0):
        raise ValueError('Positive finite moment law and explicit log-potential box required')
    first = next(iter(bank.values()), None)
    if first is None or np.ndim(first) != 2:
        raise ValueError('Nonempty regional matrix bank required')
    dimension = np.shape(first)[0]
    bank = checked_bank(bank, len(log_moments)-1, dimension=dimension)
    rows, rhs = [], []
    for (weight, _), matrix in bank.items():
        for i in range(dimension):
            targets = np.flatnonzero(np.isfinite(matrix[i]))
            if not len(targets):
                continue
            allocation = log2(len(targets))
            for j in targets:
                row = np.zeros(dimension+2)
                row[j] += 1.
                row[i] -= 1.
                row[dimension] = -1.
                rows.append(row)
                rhs.append(log_moments[weight]-matrix[i, j]-allocation)
    for i in range(dimension):
        row = np.zeros(dimension+2)
        row[i] = -1.
        row[dimension+1] = -1.
        rows.append(row)
        rhs.append(0.)
    objective = np.zeros(dimension+2)
    objective[dimension], objective[dimension+1] = 1., 1./256
    bounds = [(0., 0.)]+[(-log_bound, log_bound)]*(dimension-1)+[(0., None), (0., log_bound)]
    return objective, np.asarray(rows), np.asarray(rhs), bounds


def recompute_ratios(bank, log_moments, log_potential):
    """Compute actual full row ratios; ignore any LP objective/residual claim."""
    log_moments, potential = logs(log_moments), logs(log_potential)
    if (potential.ndim != 1 or np.any(~np.isfinite(potential)) or potential[0] != 0.
            or log_moments.ndim != 1 or np.any(~np.isfinite(log_moments))):
        raise ValueError('Finite moment and normalized positive potential required')
    bank = checked_bank(bank, len(log_moments)-1, dimension=len(potential))
    ratios = {}
    for (weight, full), matrix in bank.items():
        values = np.logaddexp2.reduce(matrix+potential[None, :], axis=1)-potential
        ratios[weight, full] = values-log_moments[weight]
    maximum = max(float(np.max(value)) for value in ratios.values())
    surplus = max(0., maximum)
    contacts = [dict(active_weight=m, full_pairs=f, state=int(i))
                for (m, f), values in ratios.items() for i, value in enumerate(values)
                if abs(value-maximum) <= 1e-9]
    return dict(log2_lambda=surplus, unclipped_maximum_log2_ratio=maximum,
        endpoint_log2=max(0., -float(np.min(potential))), contacts=contacts,
        row_log2_ratios=[dict(active_weight=m, full_pairs=f,
            rows=[None if value == -np.inf else float(value) for value in values])
            for (m, f), values in sorted(ratios.items())])


def propose_potential(bank, log_moments, *, log_bound=LOG_BOUND):
    objective, matrix, rhs, bounds = edge_budget_program(bank, log_moments, log_bound=log_bound)
    solution = linprog(objective, A_ub=matrix, b_ub=rhs, bounds=bounds, method='highs')
    if not solution.success or np.any(~np.isfinite(solution.x)):
        raise ArithmeticError('Log-domain potential LP failed: '+str(solution.message))
    dimension = len(solution.x)-2
    potential = solution.x[:dimension].copy()
    potential[0] = 0.
    actual = recompute_ratios(bank, log_moments, potential)
    actual.update(log2_potential=potential.tolist(),
        exact_log2_potential=[str(Fraction(float(value))) for value in potential],
        lp_status=int(solution.status), lp_objective=float(solution.fun),
        lp_t=float(solution.x[-2]), lp_s=float(solution.x[-1]),
        lp_maximum_edge_violation=float(np.max(matrix@solution.x-rhs)),
        lp_edge_constraints=len(rhs), log_potential_bounds=[-log_bound, log_bound],
        actual_ratios_recomputed=True, raw_tiny_potentials_used_in_lp=False,
        lp_is_only_a_candidate=True)
    return actual


def serialize_bank(bank):
    return [dict(active_weight=m, full_pairs=f,
        log2_matrix=[[None if value == -np.inf else float(value) for value in row]
                     for row in matrix]) for (m, f), matrix in sorted(bank.items())]


def run(*, max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit positive budget at most30 seconds required')
    from hd_coupled_state_regions import raw_region_bank
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: '+name)
    data, _ = load_saved()
    moment = data['two_atom_envelope']
    if (moment['log2_probability'] != LOG_PROBABILITY
            or moment['log2_nonzero_atom'] != LOG_ATOM):
        raise ValueError('Frozen old M00 moment shape changed')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Coupled mean control exhausted its budget')
    result = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        source_count=data['source_count'], mark=str(MARK), alpha=1,
        construction='proposed-g16-coupled-second-order-inner',
        scope='complete homogeneous r18/arity8/Fmax1 class only',
        state_order=list(STATES), regional_alpha=1, raw_packet_banks_built=0,
        matrix_rows_are_probability_laws=False, point_cone_is_arbitrary_mass=True,
        old_construction_receipts_transferred=False, scalar_shape='pinned old F1 R00 two-atom law',
        scalar_shape_log2_probability=LOG_PROBABILITY, scalar_shape_log2_atom=LOG_ATOM,
        scalar_shape_artifact=MOMENT_ARTIFACT, scalar_shape_artifact_sha256=MOMENT_SHA,
        extra_rowset_factor=False, source_count_once=True, projective_normalizations=1,
        rs_normalizer_exponent=1152, aggregate_surrogate_floor=FLOOR,
        all_outer_masks=True, all_boundary_state_paths=True, raw_holder_exponent=16,
        fixed_rowsets_across_columns=True, terminal_factor_once=True,
        positive_lp_potential_entries_may_exceed_one=True,
        reviewed_source_sha256=SOURCE_PINS,
        new_weight_marks=0, moment_shape_fits=0, potential_initializations=1,
        zero_matrix_entries_encoded_as_none=True, directed=False,
        whole_code_certificate=False, compute_budget_seconds=max_seconds,
        complete_requested_control=False)
    try:
        check_time()
        provided = raw_region_bank(float(MARK), GROUPS, full_max=1)
        bank = checked_bank(provided['regions'], GROUPS)
        result.update(raw_packet_banks_built=1, provider_metadata=provided['metadata'],
                      regional_log2_matrices=serialize_bank(bank), regional_histogram_count=len(bank))
        if progress:
            progress(dict(stage='raw_bank_ready', histograms=len(bank), metadata=provided['metadata']))
        check_time()
        potential = propose_potential(bank, shape(GROUPS))
        check_time()
        column = column_polynomial(raw_holder_weights(), GROUPS, LOG_PROBABILITY,
                                    LOG_ATOM, arity=ARITY, check_time=check_time)
        word = log_power(column, COLUMNS, check_time)
        tail = float(np.logaddexp2.reduce(word[FLOOR:]))
        prefix = log2(int(data['source_count']))-LIMIT*log2(float(MARK))-log2(255)-1152*log2(255)
        comparison = prefix+tail
        surplus = 256*potential['log2_lambda']+potential['endpoint_log2']
        result.update(potential=potential, frozen_scalar_comparator_log2=comparison,
            full_word_tail_log2=tail, prefix_log2=prefix, potential_penalty_log2=surplus,
            log2_screen=comparison+surplus, complete_requested_control=True)
        check_time()
        if progress:
            progress(dict(stage='gate_complete', log2_screen=result['log2_screen'],
                penalty=surplus, log2_potential=potential['log2_potential']))
    except TimeoutError as error:
        result['complete_requested_control'] = False
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=30)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    answer = run(max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
