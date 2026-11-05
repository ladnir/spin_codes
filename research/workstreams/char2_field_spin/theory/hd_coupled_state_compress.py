"""Saved-bank five-to-two cone compression with ONE fixed live mixture.

The live reference is sum_i w_i mu_i for i=L,R,B,P, with sumw=1.
The P summand is an arbitrary probability measure, not flat density.
All weights are fixed across occupancies/regions/columns. A positive
coefficient embedding proves chronological domination without assuming
actual state independence. This helper builds no regional bank or gate.
"""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from hd_coupled_state_gate import checked_bank


DIRECTORY = Path(__file__).resolve().parent
ARTIFACT = 'hd_coupled_state_gate_float.json'
ARTIFACT_SHA = '0c385a1a2bb3a04d3991cdefb6a312fe59e13add4b2b2afecfcf60a8173b5da0'
SOURCE_SHA = '974db96f33705f1d58f83fece7da5e1bd408e082946115bd9e865565eaa4efe0'
FIXED_LOG_WEIGHTS = (-64., -128., 0., -172.)
LIVE_ORDER = ('L', 'R', 'B', 'P')


def normalize_live_weights(log_weights=FIXED_LOG_WEIGHTS):
    values = np.asarray(log_weights, dtype=float)
    if values.shape != (4,) or np.any(~np.isfinite(values)):
        raise ValueError('Four strictly positive finite live weights required')
    total = float(np.logaddexp2.reduce(values))
    normalized = values-total
    return normalized, total


def compress_matrix(log_matrix, normalized_log_weights):
    """Return log2x2 S such that E T <= S E for E=(eZ; normalizedw).

    Inputs/outputs are coefficient caps, not stochastic matrices. The
    maximum over output live cones precedes no hidden summation or
    normalization. Empty rows/columns may retain exact zero (-inf).
    """
    matrix, weights = np.asarray(log_matrix, dtype=float), np.asarray(normalized_log_weights, dtype=float)
    if (matrix.shape != (5, 5) or np.any(np.isnan(matrix)) or np.any(np.isposinf(matrix))
            or weights.shape != (4,) or np.any(~np.isfinite(weights))
            or abs(float(np.logaddexp2.reduce(weights))) > 1e-12):
        raise ValueError('One finite-or-zero5x5 cap and normalized positive live weights required')
    incoming_live = np.logaddexp2.reduce(weights[:, None]+matrix[1:, :], axis=0)
    return np.array([[matrix[0, 0], np.max(matrix[0, 1:]-weights)],
                     [incoming_live[0], np.max(incoming_live[1:]-weights)]])


def compress_bank(bank, log_weights=FIXED_LOG_WEIGHTS, *, maximum=18):
    bank = checked_bank(bank, maximum)
    weights, normalization = normalize_live_weights(log_weights)
    regions = {key: compress_matrix(matrix, weights) for key, matrix in bank.items()}
    return dict(regions=regions, normalized_log_weights=weights,
        metadata=dict(state_order=['zero', 'weighted-live'], original_live_order=list(LIVE_ORDER),
            input_log_weights=list(map(float, log_weights)), normalized_log_weights=weights.tolist(),
            normalization_log2_total_mass=normalization, terminal_vector=[1, 1],
            extra_terminal_factor_log2=0., common_weights_across_regions=True,
            point_cone_is_arbitrary_probability_measure=True,
            arbitrary_point_measure_can_change_between_steps=True,
            probability_rows_assumed=False, raw_regional_banks_built=0,
            parameter_fits=0, directed=False, whole_code_certificate=False))


def entrywise_curves(regions, *, maximum=18):
    """Universal f<=1 entry envelope AFTER using one common live reference."""
    regions = checked_bank(regions, maximum, dimension=2)
    return np.stack([np.maximum.reduce([regions[m, f]
        for f in range(min(1, m//2)+1)]).reshape(4) for m in range(maximum+1)], axis=1)


def load_saved_bank():
    raw = (DIRECTORY/ARTIFACT).read_bytes()
    if sha256(raw).hexdigest() != ARTIFACT_SHA:
        raise ValueError('Frozen coupled regional artifact changed')
    source = json.loads(raw)
    required = dict(active_groups=18, arity=8, source_full_pair_max=1, mark='691/700', alpha=1,
        state_order=['Z', 'L', 'R', 'B', 'P'], regional_alpha=1,
        construction='proposed-g16-coupled-second-order-inner',
        point_cone_is_arbitrary_mass=True, complete_requested_control=True,
        source_sha256=SOURCE_SHA, zero_matrix_entries_encoded_as_none=True,
        directed=False, whole_code_certificate=False)
    if any(source.get(key) != value for key, value in required.items()):
        raise ValueError('Exact completed coupled r18/a8/Fmax1 artifact required')
    records = source['regional_log2_matrices']
    bank = {}
    for row in records:
        key = row['active_weight'], row['full_pairs']
        if key in bank:
            raise ValueError('Duplicate saved histogram')
        bank[key] = np.array([[-np.inf if value is None else value for value in line]
                              for line in row['log2_matrix']], dtype=float)
    bank = checked_bank(bank, 18)
    if source.get('regional_histogram_count') != len(bank):
        raise ValueError('Saved histogram-count metadata disagrees')
    return source, bank
