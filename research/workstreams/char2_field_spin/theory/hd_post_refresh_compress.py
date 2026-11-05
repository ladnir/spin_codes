"""One fixed normalized live reference for the post-refresh eight cones.

This is coefficient-cone domination, not stochastic lumping. P and G are
arbitrary probability families; F is an arbitrary first-coordinate family
with conditional flat nonzero second coordinate. Positive mixtures remain
in those families. The reference weights are fixed before any raw bank is
seen and are identical for every occupancy, region, and outer column.
No regional bank, optimizer, or quantitative gate runs in this module.
"""
import numpy as np


STATE_ORDER = ('Z', 'L', 'R', 'B', 'D', 'P', 'F', 'G')
LIVE_ORDER = STATE_ORDER[1:]
FIXED_LOG_WEIGHTS = (0., -128., 0., 0., -128., -128., -256.)


def checked_bank(bank, maximum=18, *, dimension=8):
    if type(maximum) is not int or maximum < 1 or type(dimension) is not int or dimension < 1:
        raise ValueError('Positive occupancy maximum and state dimension required')
    expected = {(m, f) for m in range(maximum+1) for f in range(min(1, m//2)+1)}
    if set(bank) != expected:
        raise ValueError('Complete exact f<=1 histogram scope required')
    result = {}
    for key in sorted(expected):
        matrix = np.asarray(bank[key], dtype=float)
        if (matrix.shape != (dimension, dimension) or np.any(np.isnan(matrix))
                or np.any(np.isposinf(matrix))):
            raise ValueError('Finite-or-zero matrix with the declared dimension required')
        result[key] = matrix
    return result


def normalize_live_weights(log_weights=FIXED_LOG_WEIGHTS):
    values = np.asarray(log_weights, dtype=float)
    if values.shape != (7,) or np.any(~np.isfinite(values)):
        raise ValueError('Seven strictly positive finite live weights required')
    total = float(np.logaddexp2.reduce(values))
    return values-total, total


def compress_matrix(log_matrix, normalized_log_weights):
    """Return log S with E T <= S E, E=(eZ; (0,w)), sum(w)=1.

    Consequently eZ T1...Tn 1 <= e0 S1...Sn 1, in the original order.
    The terminal vector is exactly ones; there is no repeated normalization
    or extra final live mass. Exact zeros are represented by -infinity.
    """
    matrix = np.asarray(log_matrix, dtype=float)
    weights = np.asarray(normalized_log_weights, dtype=float)
    if (matrix.shape != (8, 8) or np.any(np.isnan(matrix)) or np.any(np.isposinf(matrix))
            or weights.shape != (7,) or np.any(~np.isfinite(weights))
            or abs(float(np.logaddexp2.reduce(weights))) > 1e-12):
        raise ValueError('Eight-state cap and normalized positive live reference required')
    incoming_live = np.logaddexp2.reduce(weights[:, None]+matrix[1:, :], axis=0)
    return np.array([[matrix[0, 0], np.max(matrix[0, 1:]-weights)],
        [incoming_live[0], np.max(incoming_live[1:]-weights)]])


def compress_bank(bank, log_weights=FIXED_LOG_WEIGHTS, *, maximum=18):
    bank = checked_bank(bank, maximum)
    weights, normalization = normalize_live_weights(log_weights)
    return dict(regions={key: compress_matrix(matrix, weights) for key, matrix in bank.items()},
        normalized_log_weights=weights, metadata=dict(
            state_order=['zero', 'weighted-live'], original_state_order=list(STATE_ORDER),
            original_live_order=list(LIVE_ORDER), input_log_weights=list(map(float, log_weights)),
            normalized_log_weights=weights.tolist(), normalization_log2_total_mass=normalization,
            terminal_vector=[1, 1], extra_terminal_factor_log2=0.,
            common_weights_across_regions=True, weights_selected_before_raw_bank=True,
            P_and_G_are_arbitrary_probability_families=True,
            F_is_conditional_flat_second_family=True,
            family_representatives_may_change_between_steps=True,
            probability_rows_assumed=False, raw_regional_banks_built=0,
            parameter_fits=0, directed=False, whole_code_certificate=False))


def entrywise_curves(regions, *, maximum=18):
    """Maximum over f<=1 AFTER compression with the same live reference."""
    regions = checked_bank(regions, maximum, dimension=2)
    return np.stack([np.maximum.reduce([regions[m, f]
        for f in range(min(1, m//2)+1)]).reshape(4) for m in range(maximum+1)], axis=1)
