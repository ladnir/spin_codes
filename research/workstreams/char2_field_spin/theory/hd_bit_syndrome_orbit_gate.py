"""FLOAT fixed-row contraction with orbit-specific bit-syndrome coefficients.

Only the quiet matroid and outside evaluation-code COUNTS use affine
orbits. Every labelled A keeps its original positions, information cover,
column matrix and whole-word power. No orbit multiplicity replaces that
sum. This module loads no regional bank and supplies no global prefix.

There are two ordered count-DP passes: quiet marks A, nonquiet marks A^c.
Each pass chooses its terminal coefficient vector by the labelled row
set's class. The positive atom measure is the same for both branches and
every A; Young combines branches before the whole-word power. All bounds
here are FLOAT proposals, not directed certificates.
"""
from functools import lru_cache
from math import log2

import numpy as np

from hd_bit_syndrome_gate import (
    ENDPOINTS, LAMBDAS, _logs, checked_lambdas, combine_columns,
    cover_gap_geometry, ordered_cover_products, ordered_gap_products,
    row_geometry, weighted_cover_matrices, weighted_gap_matrices,
)
from hd_bit_syndrome_kernel import affine_eight_orbits, affine_image_counts
from hd_bit_syndrome_outer import CONSTRUCTION, quiet_support_counts
from hd_post_refresh_quiet_rows import (
    combine_row_bounds, log_matmul, log_matrix_power, moment_atoms,
    quiet_group_values, rowset_gaps, transition_geometry,
)

REPRESENTATIVES = (0x00ff, 0x017f, 0x033f, 0x035f)
FULL_QUIET_COUNTS = (4162573845, 0, 16516605, 16386810)
ORBIT_SIZES = (30, 1920, 840, 10080)


@lru_cache(maxsize=1)
def coefficient_tables():
    """Exact maxima per support weight, separately for each A orbit.

    For quiet supports use subsets of the representative A. For nonquiet
    outside supports use im(V_Ac), with its zero word removed. Affine
    relabelling transports these counts, but not chronological matrices.
    """
    if affine_eight_orbits() != tuple(zip(REPRESENTATIVES, ORBIT_SIZES)):
        raise AssertionError('Reviewed four-orbit partition changed')
    quiet_table = quiet_support_counts()
    quiet_caps, outside_caps = [], []
    for representative, full_count in zip(REPRESENTATIVES, FULL_QUIET_COUNTS):
        if quiet_table[representative] != full_count:
            raise AssertionError('Exact full-support orbit invariant changed')
        quiet = [0] * 9
        subset = representative
        while True:
            weight = subset.bit_count()
            quiet[weight] = max(quiet[weight], quiet_table[subset])
            if not subset:
                break
            subset = (subset-1) & representative
        outside = tuple(i for i in range(16) if not representative >> i & 1)
        nonquiet = [0] * 9
        for support, count in affine_image_counts(outside).items():
            if support:
                weight = support.bit_count()
                nonquiet[weight] = max(nonquiet[weight], count)
        quiet_caps.append(tuple(quiet))
        outside_caps.append(tuple(nonquiet))
    return tuple(quiet_caps), tuple(outside_caps)


@lru_cache(maxsize=1)
def labelled_classes():
    """Class of every lexicographically ordered labelled eight-set A."""
    coefficient_tables()
    quiet_table = quiet_support_counts()
    lookup = {value: index for index, value in enumerate(FULL_QUIET_COUNTS)}
    rows = row_geometry()[0]
    result = np.array([lookup[quiet_table[sum(1 << i for i in row)]] for row in rows], dtype=np.int8)
    if tuple(np.bincount(result, minlength=4)) != ORBIT_SIZES:
        raise AssertionError('Class invariant does not partition all labelled A')
    result.setflags(write=False)
    return result


def counted_column_matrices(terminal_log2, class_indices, measure, marked,
                           gap_products, gap_codes, *, batch_size=128,
                           check_time=lambda: None):
    """One ordered DP with class-selected terminal scalar factors.

    terminal_log2[class,count] already includes the group power. Gaps are
    products under the same positive atom measure, in their actual order.
    No per-class DP or row-set multiplicity is needed.
    """
    terminal, products = _logs(terminal_log2), _logs(gap_products)
    classes, codes = np.asarray(class_indices), np.asarray(gap_codes)
    counts, layers = transition_geometry(marked)
    if (terminal.ndim != 2 or terminal.shape[1] != len(counts) or not len(terminal)
            or classes.ndim != 1 or not len(classes) or classes.dtype.kind not in 'iu'
            or np.min(classes) < 0 or np.max(classes) >= len(terminal)
            or products.ndim != 3 or products.shape[1:] != (2, 2)
            or codes.shape != (len(classes), marked+1) or codes.dtype.kind not in 'iu'
            or np.min(codes) < 0 or np.max(codes) >= len(products)
            or type(batch_size) is not int or not 1 <= batch_size <= 512):
        raise ValueError('Complete class/count tables and ordered gap geometry required')
    masses, _ = moment_atoms(measure)
    result = np.full((len(classes), 2, 2), -np.inf)
    for start in range(0, len(classes), batch_size):
        check_time()
        stop = min(len(classes), start+batch_size)
        block = codes[start:stop]
        current = products[block[:, 0]][None, ...]
        for position, (next_size, targets) in enumerate(layers):
            check_time()
            following = np.full((next_size, stop-start, 2, 2), -np.inf)
            after = products[block[:, position+1]]
            for label, ((incoming, outgoing), mass) in enumerate(zip(ENDPOINTS, masses)):
                if mass == -np.inf:
                    continue
                term = current[..., incoming, None] + after[None, :, None, outgoing, :] + mass
                target = targets[label]
                following[target] = np.logaddexp2(following[target], term)
            current = following
        scalar = terminal[classes[start:stop]].T
        result[start:stop] = np.logaddexp2.reduce(current+scalar[:, :, None, None], axis=0)
    return result, dict(rowset_count=len(classes), final_atom_count_classes=len(counts),
        coefficient_classes=len(terminal), dp_passes=1, batch_size=batch_size,
        maximum_dp_layer_bytes=len(counts)*min(batch_size, len(classes))*4*8,
        gap_word_slots=len(products))


def contract(log_quiet_caps, log_outside_caps, row_classes, measure, groups, log2_x, *,
             size=16, columns=16, q=256, lambdas=LAMBDAS, batch_size=128,
             check_time=lambda: None, progress=None, information_cover=False):
    """Generic small-test interface; supplied caps must hold per row class."""
    quiet_caps, outside_caps = _logs(log_quiet_caps), _logs(log_outside_caps)
    if (quiet_caps.ndim != 2 or outside_caps.ndim != 2 or not len(quiet_caps)
            or len(quiet_caps) != len(outside_caps) or not 2 <= quiet_caps.shape[1] <= 9
            or type(groups) is not int or groups < 1 or type(columns) is not int or columns < 1
            or not np.isfinite(log2_x) or log2_x < 0 or type(information_cover) is not bool):
        raise ValueError('Class coefficient tables, positive scope, X>=1 and bool cover required')
    checked_lambdas(lambdas)
    arity = quiet_caps.shape[1]-1
    rows, information, ranks, gap_codes = row_geometry(size, arity)
    marked = size-arity
    classes = np.asarray(row_classes)
    if (outside_caps.shape[1] != marked+1 or np.any(outside_caps[:, 0] != -np.inf)
            or classes.shape != (len(rows),)):
        raise ValueError('All labelled row classes and zero-constant outside caps required')
    check_time()
    quiet_rows, quiet_gaps = rowset_gaps(size, arity)
    if quiet_rows != rows:
        raise AssertionError('Quiet and nonquiet labelled row orders differ')
    _, empty = moment_atoms(measure)
    quiet_products = np.stack([log_matrix_power(empty, n) for n in range(marked+1)])
    quiet_values = np.stack([groups*quiet_group_values(caps, measure, log2_x) for caps in quiet_caps])
    quiet, quiet_cost = counted_column_matrices(quiet_values, classes, measure, arity,
        quiet_products, quiet_gaps, batch_size=batch_size, check_time=check_time)
    if progress:
        progress('quiet_complete')
    covers = None
    if information_cover:
        covers, cover_values, parents, labels, gap_codes = cover_gap_geometry(size, arity)
        bank = weighted_cover_matrices(measure, groups, cover_values, q=q)
        products = ordered_cover_products(bank['cover_matrices'], parents, labels)
    else:
        bank = weighted_gap_matrices(measure, groups, q=q)
        products = ordered_gap_products(bank['information'], bank['constrained'], arity)
    nonquiet_values = np.stack([groups*quiet_group_values(caps, measure, 0.) for caps in outside_caps])
    nonquiet, nonquiet_cost = counted_column_matrices(nonquiet_values, classes, measure, marked,
        products, gap_codes, batch_size=batch_size, check_time=check_time)
    if progress:
        progress('nonquiet_complete')
    combined, shifts = combine_columns(quiet, nonquiet, groups, arity, log2_x, lambdas=lambdas)
    check_time()
    word = log_matrix_power(combined, columns)
    values = np.logaddexp2.reduce(word[:, 0, :], axis=1)
    return dict(rowsets=rows, information_sets=information, outside_ranks=ranks,
        row_classes=classes, information_cover=information_cover, information_cover_exponents=covers,
        per_rowset_word_log2=values, combined_word_log2=combine_row_bounds(values, groups),
        column_log2_matrices=combined, branch_column_log2_matrices=dict(quiet=quiet, nonquiet=nonquiet),
        young_log2_scales=shifts, gap_log2_matrices=bank,
        costs=dict(quiet=quiet_cost, nonquiet=nonquiet_cost, total_dp_passes=2),
        terminal_vector=[1, 1], extra_rowset_factor=False, orbit_multiplicity_factor=False,
        common_positive_measure=True, all_labelled_rowsets=True, branch_sum_before_word_power=True)


def full_row_bounds(measure, groups, log2_x, *, lambdas=LAMBDAS, batch_size=128,
                    check_time=lambda: None, progress=None, information_cover=False):
    """Same raw return interface as the frozen gate; no provider or prefix."""
    check_time()
    quiet, outside = coefficient_tables()
    def logs(table):
        return np.array([[log2(value) if value else -np.inf for value in row] for row in table])
    answer = contract(logs(quiet), logs(outside), labelled_classes(), measure, groups, log2_x,
        size=16, columns=16, q=256, lambdas=lambdas, batch_size=batch_size,
        check_time=check_time, progress=progress, information_cover=information_cover)
    answer['metadata'] = dict(outer_construction=CONSTRUCTION, intended_inner='original-one-state',
        bank_loading='caller-owned; no bank or provider loaded by this engine',
        groups=groups, arity=8, columns=16, field_order=256,
        orbit_representatives=[hex(value) for value in REPRESENTATIVES],
        orbit_sizes=list(ORBIT_SIZES), orbit_full_quiet_counts=list(FULL_QUIET_COUNTS),
        quiet_coefficients=[list(row) for row in quiet], outside_coefficients=[list(row) for row in outside],
        orbit_use='coefficient envelopes only; every labelled row set retains its chronology',
        lexicographic_information_sets=not information_cover, information_cover=information_cover,
        source_rowsets_fixed_across_columns=True, common_positive_atom_measure=True,
        atoms_above_one_supported=True, young_weights=list(map(str, checked_lambdas(lambdas))),
        exact_sparse_restoration=False, global_prefix_included=False, total_dp_passes=2,
        directed=False, whole_code_certificate=False)
    return answer
