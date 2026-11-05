"""FLOAT log-domain state kernel for two columns joined by xor-row pairs.

Fix two labelled support masks S,T before the pair layer.  Each pair uses
fresh nonzero alpha and maps (u,v) to (u+alpha*v,u+2*alpha*v).  At scalar
output marks t,u, its support generating functions are1 for type00,
t*u for one active input, and(t+u+(q-3)*t*u)/(q-1) for type11.

A common positive matrix measure assigns each atom t a matrix E(t).
For r groups, the local four-boundary tensor is the sum of
f(t,u)**r E(t)[a,b] E(u)[c,d] over two independent atom indices.  Atom
matrices include their masses; they are NOT raised to r.  The same
measure is used at every position.  It need not have total mass one.

The xor-row tree preserves both natural chronological matrix chains.
At a set bit of delta, only the second chain reverses the two child
intervals.  After the tree, summing T[a,b,b,d] joins the two columns.
No regional matrix commutation, state quotient, source resampling, or
entrywise matrix power is used.

All numeric APIs take and return log2 values, with -inf meaning zero.
These are ordinary FLOAT contractions, not directed certificates.  They
contain no support-count coefficient, Minkowski branch sum, RS factor,
source count, final-state sum, or whole-code claim.
"""
from math import log2

import numpy as np

from hd_post_refresh_quiet_rows import moment_atoms


MAX_ROWS, MAX_BATCH = 16, 8192
EXACT_RATIO = 'fresh-ratio-exact'
SYMMETRIC_MAJORANT = 'symmetric-shear-majorant'


def _logs(values):
    values = np.asarray(values, dtype=float)
    if np.any(np.isnan(values)) or np.any(np.isposinf(values)):
        raise ValueError('Finite log2 values or exact zero(-inf) required')
    return values


def common_moment_atoms(measure):
    """The unchanged five-entry moment measure as atom logs and matrix logs.

    Returns(log2_atoms,log2_matrices).  Each matrix already contains its
    mass.  No normalization by total measure mass is performed.
    """
    masses, _ = moment_atoms(measure)
    atoms = np.array((-np.inf, measure['log2_atom00'], measure['log2_atom01'],
                      measure['log2_atom10'], 0.))
    matrices = np.full((5, 2, 2), -np.inf)
    for index, (row, column) in enumerate(((0, 0), (0, 0), (0, 1), (1, 0), (1, 1))):
        matrices[index, row, column] = masses[index]
    return atoms, matrices


def support_pair_logs(log2_t, log2_u, *, q=256, mode=EXACT_RATIO):
    """Return log f for types00,single,11; the optional mode is an upper law.

    The symmetric-shear expression is(t+u+(q-4)tu)/(q-2).  It is an
    amplitude-independent majorant only for0<=t,u<=1, not an exact law.
    """
    if type(q) is not int or q not in (4, 8, 16, 256):
        raise ValueError('Supported characteristic-two fields are4,8,16,256')
    left, right = np.broadcast_arrays(_logs(log2_t), _logs(log2_u))
    if mode not in (EXACT_RATIO, SYMMETRIC_MAJORANT):
        raise ValueError('Unknown pair support law')
    if mode == SYMMETRIC_MAJORANT and (np.any(left > 0) or np.any(right > 0)):
        raise ValueError('The symmetric-shear majorant requires atoms at most one')
    coefficient, denominator = (q-3, q-1) if mode == EXACT_RATIO else (q-4, q-2)
    single = left+right
    double = np.logaddexp2(left, right)
    if coefficient:
        double = np.logaddexp2(double, log2(coefficient)+single)
    double -= log2(denominator)
    return np.stack((np.zeros(left.shape), single, double))


def leaf_tensors(log2_atoms, log2_matrices, groups, *, q=256, mode=EXACT_RATIO):
    """Return log tensors[00/single/11,a,b,c,d] from a finite common measure."""
    atoms, matrices = _logs(log2_atoms), _logs(log2_matrices)
    if (atoms.ndim != 1 or not 1 <= len(atoms) <= 32
            or matrices.shape != (len(atoms), 2, 2)
            or type(groups) is not int or not 1 <= groups <= 8192):
        raise ValueError('At most32 atoms, their2x2 matrices, and a positive bounded group count required')
    factors = support_pair_logs(atoms[:, None], atoms[None, :], q=q, mode=mode)
    # Axes are(type,first_atom,second_atom,a,b,c,d).  The two atom
    # indices are summed AFTER taking the scalar support factor to r.
    masses = (matrices[:, None, :, :, None, None]
              +matrices[None, :, None, None, :, :])
    terms = masses[None, ...]+groups*factors[:, :, :, None, None, None, None]
    result = np.logaddexp2.reduce(np.logaddexp2.reduce(terms, axis=2), axis=1)
    return _logs(result)


def leaves_from_measure(measure, groups, *, q=256, mode=EXACT_RATIO):
    return leaf_tensors(*common_moment_atoms(measure), groups, q=q, mode=mode)


def _permute_last4(values, order):
    prefix = tuple(range(values.ndim-4))
    return values.transpose(prefix+tuple(values.ndim-4+i for i in order))


def join_tensors(left, right, *, reverse_second=False):
    """Positive contraction with any broadcastable leading batch axes.

    Normal: L[a,b,c,e] R[b,d,e,f] -> result[a,d,c,f].
    Reversed second chain: L[a,b,e,f] R[b,c,d,e] -> result[a,c,d,f].
    In either case only the two shared boundary states are summed.
    """
    left, right = _logs(left), _logs(right)
    if (left.ndim < 4 or right.ndim < 4
            or left.shape[-4:] != (2,)*4 or right.shape[-4:] != (2,)*4
            or type(reverse_second) is not bool):
        raise ValueError('Two four-boundary2x2x2x2 tensors required')
    if reverse_second:
        # Output axes(a,c,d,f), followed by summed axes(b,e).
        a = _permute_last4(left, (0, 3, 1, 2))[..., :, None, None, :, :, :]
        b = _permute_last4(right, (1, 2, 0, 3))[..., None, :, :, None, :, :]
    else:
        # Output axes(a,d,c,f), followed by summed axes(b,e).
        a = _permute_last4(left, (0, 2, 1, 3))[..., :, None, :, None, :, :]
        b = _permute_last4(right, (1, 3, 0, 2))[..., None, :, None, :, :, :]
    result = np.logaddexp2.reduce(np.logaddexp2.reduce(a+b, axis=-1), axis=-1)
    return _logs(result)


def _shape(rows, delta):
    if (type(rows) is not int or not 1 <= rows <= MAX_ROWS or rows & (rows-1)
            or type(delta) is not int or not 0 <= delta < rows):
        raise ValueError('Power-of-two rows through16 and an in-range xor displacement required')


def tensor_tree(tensors, delta):
    """Contract(...,rows,2,2,2,2) leaves, retaining all four boundaries."""
    current = _logs(tensors)
    if current.ndim < 5 or current.shape[-4:] != (2,)*4:
        raise ValueError('A row sequence of four-boundary tensors is required')
    rows = current.shape[-5]
    _shape(rows, delta)
    for bit in range(rows.bit_length()-1):
        current = join_tensors(current[..., 0::2, :, :, :, :],
                               current[..., 1::2, :, :, :, :],
                               reverse_second=bool(delta >> bit & 1))
    return current[..., 0, :, :, :, :]


def join_columns(tensor):
    """Join the end of the first column to the start of the second exactly."""
    tensor = _logs(tensor)
    if tensor.ndim < 4 or tensor.shape[-4:] != (2,)*4:
        raise ValueError('A four-boundary tensor is required')
    return np.logaddexp2(tensor[..., :, 0, 0, :], tensor[..., :, 1, 1, :])


def branch_log_matrices(first_masks, second_masks, delta, leaves, *, rows=16):
    """Batched exact support branches, returning log2(batch,2,2) matrices.

    S_i is paired with T_(i xor delta).  The first and second masks remain
    fixed throughout their columns.  Counts and branch Minkowski weights
    belong to the caller and must not be silently folded into the leaves.
    """
    _shape(rows, delta)
    first, second = np.asarray(first_masks), np.asarray(second_masks)
    leaves = _logs(leaves)
    if (first.ndim != 1 or second.shape != first.shape or len(first) > MAX_BATCH
            or not np.issubdtype(first.dtype, np.integer)
            or not np.issubdtype(second.dtype, np.integer)
            or np.any(first < 0) or np.any(second < 0)
            or np.any(first >= 1 << rows) or np.any(second >= 1 << rows)
            or leaves.shape != (3, 2, 2, 2, 2)):
        raise ValueError('At most8192 matched support-mask pairs and three leaf tensors required')
    positions = np.arange(rows, dtype=np.uint32)[None, :]
    types = (((first.astype(np.uint32)[:, None] >> positions) & 1)
             +((second.astype(np.uint32)[:, None] >> (positions ^ delta)) & 1))
    return join_columns(tensor_tree(leaves[types], delta))


def branch_log_matrix(first_mask, second_mask, delta, leaves, *, rows=16):
    return branch_log_matrices(np.array([first_mask]), np.array([second_mask]),
                               delta, leaves, rows=rows)[0]


def metadata(*, mode=EXACT_RATIO):
    if mode not in (EXACT_RATIO, SYMMETRIC_MAJORANT):
        raise ValueError('Unknown pair support law')
    return dict(arithmetic='FLOAT-log2-not-directed', pair_support_law=mode,
                exact_support_law=mode == EXACT_RATIO,
                requires_atoms_at_most_one=mode == SYMMETRIC_MAJORANT,
                shared_positive_matrix_measure=True, independent_atom_indices_per_region=True,
                scalar_support_factor_powered_before_atom_sum=True,
                atom_matrix_masses_powered=False, fixed_labelled_supports=True,
                both_chronological_state_chains_retained=True,
                columns_joined_once=True, rowset_sum_included=False,
                rs_normalization_included=False, source_count_included=False,
                whole_code_certificate=False)
