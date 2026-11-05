"""FLOAT state moments with no outer crossed-pair layer.

The regional common matrix measure and global paired route are unchanged.
Only the outer layer between two RS columns is omitted. For group count r,
an inactive support bit contributes S0=sum_t E(t); an active bit contributes
Sr=sum_t t**r E(t). A fixed mask uses the ordered product of these matrices.
The two column products are then multiplied in chronological order.

All APIs use log2 with -inf as exact zero. No ordinary-domain conversion,
provider, source/root sum, RS normalization, or certificate is included.
"""
import numpy as np

from hd_pair_state_kernel import common_moment_atoms
from hd_post_refresh_quiet_rows import log_matmul


MAX_ROWS = 16
MAX_BATCH = 8192


def _logs(values):
    result = np.asarray(values, dtype=float)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Finite log2 values or exact zero (-inf) required')
    return result


def entry_log_matrices(measure, groups):
    """Return (S0,Sr) in log2; atom masses are not raised to r."""
    if type(groups) is not int or not 1 <= groups <= 18:
        raise ValueError('Group count1..18 required')
    atoms, matrices = common_moment_atoms(measure)
    zero = np.logaddexp2.reduce(matrices, axis=0)
    nonzero = np.logaddexp2.reduce(matrices+groups*atoms[:, None, None], axis=0)
    return zero, nonzero


def mask_table_from_entries(log_zero, log_nonzero, *, rows=16):
    """All masks in integer order, with bit i acting at chronological row i.

    Each entry argument may be one common2x2 matrix or a sequence of rows
    matrices. The latter form permits independent chronology tests.
    """
    if type(rows) is not int or not 1 <= rows <= MAX_ROWS:
        raise ValueError('One to sixteen ordered rows required')
    zero, nonzero = _logs(log_zero), _logs(log_nonzero)
    if zero.shape not in ((2, 2), (rows, 2, 2)) or nonzero.shape not in ((2, 2), (rows, 2, 2)):
        raise ValueError('Common or per-row two-state matrices required')
    zero = np.broadcast_to(zero, (rows, 2, 2))
    nonzero = np.broadcast_to(nonzero, (rows, 2, 2))
    result = np.full((1, 2, 2), -np.inf)
    result[0, 0, 0] = result[0, 1, 1] = 0.
    for row in range(rows):
        result = np.concatenate((log_matmul(result, zero[row]),
                                 log_matmul(result, nonzero[row])), axis=0)
    return result


def mask_log_matrices(measure, groups, *, rows=16):
    """Build the full2**rows table once from the existing common measure."""
    return mask_table_from_entries(*entry_log_matrices(measure, groups), rows=rows)


def branch_log_matrices(first_masks, second_masks, table):
    """Gather and join two natural columns; there is no xor-pair delta."""
    table = _logs(table)
    first, second = np.asarray(first_masks), np.asarray(second_masks)
    size = len(table) if table.ndim else 0
    if (table.ndim != 3 or table.shape[1:] != (2, 2) or not 2 <= size <= 1 << MAX_ROWS
            or size & (size-1) or first.ndim != 1 or second.shape != first.shape
            or len(first) > MAX_BATCH or not np.issubdtype(first.dtype, np.integer)
            or not np.issubdtype(second.dtype, np.integer)
            or np.any(first < 0) or np.any(second < 0)
            or np.any(first >= size) or np.any(second >= size)):
        raise ValueError('Compatible mask table and at most8192 in-range mask pairs required')
    return log_matmul(table[first], table[second])


def metadata():
    return dict(arithmetic='FLOAT-log2-not-directed', outer_crossed_pair_layer=False,
        global_paired_route_unchanged=True, same_common_regional_measure=True,
        empty_entries_from_common_measure=True, atom_masses_powered=False,
        fixed_supports_across_column_rows=True, chronological_columns=True,
        source_count_included=False, rs_normalization_included=False,
        rowset_root_sum_included=False, regional_provider_runs=0,
        native_implementation_changed=False, whole_code_certificate=False)
