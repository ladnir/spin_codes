"""Bounded FLOAT five-edge-count contraction for zero to fifteen marks.

The five labels are 00/atom0, 00/atomB, 01, 10 and 11. Count classes
compress only a symmetric terminal polynomial. Every edge direction and
ordered gap matrix remains in the two-state path contraction. A terminal
table already includes the caller's group power; this engine adds none.

The caller supplies all coefficient caps and one common positive measure.
This module supplies no construction, source count, regional bank, global
prefix or certificate. All numerical operations are ordinary FLOAT.
The reviewed helpers restricted to eight marks are left unchanged.
"""
from functools import lru_cache
from math import comb

import numpy as np

from hd_post_refresh_quiet_rows import moment_atoms
from hd_two_sum_gf16_state_envelope import elementary_logs

ENDPOINTS = ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1))
MAX_MARKED = 15


def _marked(value):
    if type(value) is not int or not 0 <= value <= MAX_MARKED:
        raise ValueError('An exact marked count in0..15 is required')
    return value


def _logs(values):
    result = np.asarray(values, dtype=float)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Finite-or-zero logarithms required')
    return result


def _increment(row, label):
    # Setup only. The hot DP below has five fixed vectorized edge updates.
    values = list(row)
    values[label] += 1
    return tuple(values)


@lru_cache(maxsize=MAX_MARKED+1)
def transition_geometry_general(marked):
    """Lexicographic weak compositions and five injective layer maps.

    Exactly C(marked+4,4) terminal classes occur:715 at9,3876 at15.
    For1..8 the complete tuple and target-index order matches the frozen
    transition_geometry. At0 there is one zero-count class and no layer.
    """
    marked = _marked(marked)
    if marked == 0:
        return ((0, 0, 0, 0, 0),), ()
    previous, layers = transition_geometry_general(marked-1)
    counts = tuple((a, b, c, d, marked-a-b-c-d)
        for a in range(marked+1)
        for b in range(marked-a+1)
        for c in range(marked-a-b+1)
        for d in range(marked-a-b-c+1))
    if len(counts) != comb(marked+4, 4):
        raise AssertionError('Five-label count geometry is incomplete')
    lookup = {row: index for index, row in enumerate(counts)}
    targets = tuple(np.asarray([lookup[_increment(row, label)] for row in previous],
                               dtype=np.int32) for label in range(5))
    for target in targets:
        target.setflags(write=False)
    return counts, layers+((len(counts), targets),)


def quiet_group_values_general(log_counts, measure, log2_x):
    """Log(sum_b caps[b] e_b(atoms) X^b) for every count tuple.

    The historical name does not restrict the constant coefficient or
    the branch: complete image caps with constant1 are accepted. A zero
    constant is represented by -inf. The atom0 labels contribute no
    positive elementary degree, but e0=1 always remains. No coefficient
    multiplicity, group power, edge mass or row-set factor is added here.
    """
    coefficients = _logs(log_counts)
    if coefficients.ndim != 1 or not 1 <= len(coefficients) <= MAX_MARKED+1:
        raise ValueError('One coefficient for each degree through marked0..15 required')
    if type(log2_x) not in (int, float, np.float64, np.float32) or not np.isfinite(log2_x) or log2_x < 0:
        raise ValueError('Finite activity marker X>=1 required')
    moment_atoms(measure)
    counts, _ = transition_geometry_general(len(coefficients)-1)
    result = np.empty(len(counts), dtype=float)
    for index, (_, n00, n01, n10, n11) in enumerate(counts):
        markers = ([measure['log2_atom00']]*n00+[measure['log2_atom01']]*n01
                   +[measure['log2_atom10']]*n10+[0.]*n11)
        elementary = elementary_logs(np.asarray(markers))
        result[index] = np.logaddexp2.reduce(coefficients[:len(elementary)]
            +elementary+log2_x*np.arange(len(elementary)))
    return result


def counted_column_matrices(terminal_log2, class_indices, measure, marked,
                            gap_products, gap_codes, *, batch_size=128,
                            check_time=lambda: None):
    """Class-selected ordered column DP; signature matches the orbit engine.

    terminal_log2[class,count] includes the caller's entire group power.
    Each labelled row has marked+1 ordered gaps. Gap products must use
    the SAME positive measure as the marked edges. No homogeneous-row
    replacement by exact empty operators is made or authorized here.
    """
    marked = _marked(marked)
    terminal, products = _logs(terminal_log2), _logs(gap_products)
    classes, codes = np.asarray(class_indices), np.asarray(gap_codes)
    counts, layers = transition_geometry_general(marked)
    if (terminal.ndim != 2 or terminal.shape[1] != len(counts) or not len(terminal)
            or classes.ndim != 1 or not len(classes) or classes.dtype.kind not in 'iu'
            or np.min(classes) < 0 or np.max(classes) >= len(terminal)
            or products.ndim != 3 or products.shape[1:] != (2, 2) or not len(products)
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
                term = current[..., incoming, None]+after[None, :, None, outgoing, :]+mass
                target = targets[label]
                # Incrementing a fixed coordinate is injective, so this
                # vectorized assignment loses no repeated-index additions.
                following[target] = np.logaddexp2(following[target], term)
            current = following
        scalar = terminal[classes[start:stop]].T
        result[start:stop] = np.logaddexp2.reduce(current+scalar[:, :, None, None], axis=0)
    return result, dict(rowset_count=len(classes), marked_positions=marked,
        final_atom_count_classes=len(counts), coefficient_classes=len(terminal),
        dp_passes=1, batch_size=batch_size,
        maximum_dp_layer_bytes=len(counts)*min(batch_size, len(classes))*4*8,
        gap_word_slots=len(products), numerical_arithmetic='FLOAT log2',
        directed=False, whole_code_certificate=False, global_prefix_included=False)
