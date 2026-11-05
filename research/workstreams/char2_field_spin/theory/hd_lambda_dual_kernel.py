"""FLOAT column bounds from a raw dual support enumerator.

Fix an allowed-input column code C and its dual Cperp, of dimension b.
For atoms 0 <= t_i <= 1, the MacWilliams identity gives

 K_C(t) = q**(-b) sum_T n_T prod_{i not in T}(1+(q-1)t_i)
                                  prod_{i in T}(1-t_i).

Here n_T counts raw dual words, so sum_T n_T = q**b. Categories may
repeat a support mask; they need not be merged. The caller must supply
the actual dual code, with the same column permutations as the primal.

Let E be one common positive matrix measure and r the group count.
For a fixed dual mask T, J_T is the ordered product of the matrices
integral f_{T_i}(t)**r E(dt). Atom masses are NOT raised to r. Entrywise
Minkowski bounds the primal column moment by

 (q**(-b) sum_T n_T J_T**(1/r))**r.

The roots and final powers are entrywise; the products defining J_T
retain chronological matrix order. No normalized-category Jensen step
is used. In particular, the all-one atom gives q**((rows-b)*r), not an
extra dual-mass factor. All arithmetic is FLOAT log2, with -inf for zero.

There is no original-input floor here: a dual weight is not an input
weight. This module includes no RS, source, projective, or row-set sum
factor, no regional provider, and no directed certificate.
"""
from math import log2
from numbers import Integral

import numpy as np

from hd_pair_identity_kernel import mask_table_from_entries
from hd_pair_state_kernel import common_moment_atoms


MAX_ROWS = 16
MAX_GROUPS = 8192


def _logs(values):
    result = np.asarray(values, dtype=float)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Finite log2 values or exact zero (-inf) required')
    return result


def _parameters(groups, q):
    if type(groups) is not int or not 1 <= groups <= MAX_GROUPS:
        raise ValueError('A group count from 1 through 8192 is required')
    if type(q) is not int or not 2 <= q <= 256 or q & (q-1):
        raise ValueError('A characteristic-two field order from 2 through 256 is required')


def entry_logs_from_atoms(log2_atoms, log2_matrices, groups, *, q=256):
    """Return transformed (zero-support, nonzero-support) entry logs.

    At most 32 atoms are accepted. Each positive 2x2 matrix includes its
    measure mass, which is used once. Atoms above one are rejected because
    the dual nonzero factor would be negative.
    """
    _parameters(groups, q)
    atoms, matrices = _logs(log2_atoms), _logs(log2_matrices)
    if (atoms.ndim != 1 or not 1 <= len(atoms) <= 32
            or matrices.shape != (len(atoms), 2, 2) or np.any(atoms > 0)):
        raise ValueError('At most 32 atoms in [0,1] and their 2x2 matrix masses required')
    zero_factors = np.logaddexp2(0., log2(q-1)+atoms)
    # expm1 preserves 1-t near t=1. Endpoints t=0 and t=1 become
    # factor logs 0 and -inf, respectively, without replacing zeros.
    with np.errstate(divide='ignore', invalid='raise'):
        nonzero_factors = np.log2(-np.expm1(np.log(2.)*atoms))
    zero = np.logaddexp2.reduce(matrices+groups*zero_factors[:, None, None], axis=0)
    nonzero = np.logaddexp2.reduce(matrices+groups*nonzero_factors[:, None, None], axis=0)
    return _logs(zero), _logs(nonzero)


def entry_log_matrices(measure, groups, *, q=256):
    """Transform the unchanged common regional matrix measure."""
    return entry_logs_from_atoms(*common_moment_atoms(measure), groups, q=q)


def mask_log_matrices(measure, groups, *, q=256, rows=16):
    """Return log J_T for every mask T, shape (2**rows,2,2).

    Bit i selects the nonzero dual factor at chronological row i. The
    full table is built once; production callers may reuse it for all
    fixed source sets and columns under this same measure and group count.
    """
    return mask_table_from_entries(*entry_log_matrices(measure, groups, q=q), rows=rows)


def _count_logs(raw_counts, dual_arity, q):
    counts = np.asarray(raw_counts, dtype=object)
    if counts.ndim < 1 or not counts.shape[-1]:
        raise ValueError('A nonempty final category axis is required')
    if any(isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral) or n < 0
           for n in counts.flat):
        raise ValueError('Raw category counts must be exact nonnegative integers')
    expected = q**dual_arity
    if any(sum(int(n) for n in row) != expected
           for row in counts.reshape(-1, counts.shape[-1])):
        raise ValueError('Raw dual counts must sum exactly to q**dual_arity')
    return np.fromiter((log2(int(n)) if n else -np.inf for n in counts.flat),
                       dtype=float, count=counts.size).reshape(counts.shape)


def column_log_bound(masks, raw_counts, dual_arity, groups, table, *, q=256):
    """Bound primal column moments, preserving all leading batch axes.

    masks has shape (...,categories). Exact raw_counts must broadcast to
    that shape without adding axes. table is the result of the transformed
    mask construction for the SAME r, q, and common measure. This function
    cannot infer those provenance parameters from an untagged numpy array.

    Return shape is (...,2,2). The implemented expression is
    r*(logsumexp2(log n_T + log J_T/r) - b*log2(q)).
    """
    _parameters(groups, q)
    table, masks = _logs(table), np.asarray(masks)
    size = len(table) if table.ndim else 0
    if (table.ndim != 3 or table.shape[1:] != (2, 2)
            or not 2 <= size <= 1 << MAX_ROWS or size & (size-1)):
        raise ValueError('An ordered mask table of shape (2**rows,2,2) is required')
    rows = size.bit_length()-1
    if type(dual_arity) is not int or not 0 <= dual_arity <= rows:
        raise ValueError('Dual dimension must lie between zero and the row count')
    if (masks.ndim < 1 or not masks.shape[-1]
            or not np.issubdtype(masks.dtype, np.integer)
            or np.any(masks < 0) or np.any(masks >= size)):
        raise ValueError('Integer masks in range with a final category axis required')
    count_logs = _count_logs(raw_counts, dual_arity, q)
    if count_logs.shape[-1] != masks.shape[-1]:
        raise ValueError('Mask and count category axes must agree')
    try:
        count_logs = np.broadcast_to(count_logs, masks.shape)
    except ValueError as error:
        raise ValueError('Raw counts must broadcast to the mask batch') from error
    terms = table[masks]/groups+count_logs[..., None, None]
    return groups*(np.logaddexp2.reduce(terms, axis=-3)-dual_arity*log2(q))


def metadata():
    return dict(arithmetic='FLOAT-log2-not-directed', atom_range='[0,1]',
        dual_support_counts='raw', dual_mass_normalization='q**(-dual_arity)',
        atom_masses_powered=False, entrywise_minkowski=True,
        same_common_regional_measure=True, ordered_matrix_products=True,
        original_input_floor_included=False, rs_normalization_included=False,
        source_count_included=False, projective_factor_included=False,
        rowset_root_sum_included=False, regional_provider_runs=0,
        native_implementation_changed=False, whole_code_certificate=False)
