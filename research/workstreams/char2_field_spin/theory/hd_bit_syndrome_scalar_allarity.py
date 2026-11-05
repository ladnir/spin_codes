"""FLOAT scalar all-arity bounds for the fixed bit-syndrome outer H*.

Fix a labelled nonempty input set A, its complement O, and C=H*[O,A].
The outside image has weight enumerator C_A.  Since H* is a symmetric
involution, the zero-outside INSIDE fiber is D_A=ker C; its dual is the
outside image for source O.  Exact MacWilliams inversion therefore gives
the weight enumerator Z_A without enumerating inside supports.

For 0<=t<=1 every translated fiber has weighted mass at most Z_A(t).
This follows by positive zero-set expansion: each affine intersection is
empty or has the cardinality of its homogeneous intersection.  Hence
K_A(t)<=C_A(t) Z_A(t).  This is POINTWISE scalar domination, not a
coefficientwise bound on actual output weights.  The dimension-only
alternative replaces Z_A by (1+255t)^dim(D_A).

Both majorants have exact constant1.  Remove that coefficient in integer
arithmetic, then impose >=9 nonzero columns out of16 on each FIXED A.
Only then sum A and pay 255^(-8|A|).  Identical scalar polynomials may be
grouped with their labelled multiplicities; no row set is resampled.

Regional matrices are caller-owned.  A common positive state potential
and a scalar slope t<=1 give R[m,f]v<=eta*t^m*v.  The full-word factor is
eta^256, with one terminal factor 1/min(v) and one projective factor1/255.
All evaluations here are FLOAT proposals, not directed certificates.
No provider, optimizer, native code, or production gate runs at import.
"""
from collections import Counter
from functools import lru_cache
from math import comb, log2

import numpy as np

from hd_bit_syndrome_flat_geometry import (
    FULL, correction_zero_set_counts, flat_hull_mask, rank_cross_masks)
from hd_bit_syndrome_outer import CONSTRUCTION

Q, SIZE, COLUMNS, MINIMUM, INPUT_ROWS = 256, 16, 16, 9, 8
LIMIT = 734003


def _mask(value):
    if type(value) is not int or not 0 <= value <= FULL:
        raise ValueError('A labelled16-coordinate integer mask is required')
    return value


def _affine_dimension(mask):
    # A nonempty affine flat of size2^k has augmented row rank k+1.
    return flat_hull_mask(mask).bit_count().bit_length()


@lru_cache(maxsize=65536)
def outside_weight_counts(source_mask):
    """Exact image WE of H*[A^c,A], for arities0..16.

    Global correction coefficients are counted once by the flat helper.
    Restriction to O has constant fiber size q^(rank(V_A)-rank(C)).
    Dividing each aggregated weight count removes that multiplicity.
    This uses the actual cross rank, including the arity8 exceptions.
    """
    source_mask = _mask(source_mask)
    outside = FULL ^ source_mask
    cross_rank = rank_cross_masks(outside, source_mask)
    divisor = Q ** (_affine_dimension(source_mask)-cross_rank)
    counts = [0] * (outside.bit_count()+1)
    for zeros, count in correction_zero_set_counts(source_mask).items():
        counts[(outside & ~zeros).bit_count()] += count
    if any(count % divisor for count in counts):
        raise AssertionError('Outside restriction has a nonconstant fiber')
    result = tuple(count // divisor for count in counts)
    if result[0] != 1 or sum(result) != Q ** cross_rank:
        raise AssertionError('Outside image mass or unique zero changed')
    return result


def macwilliams_kernel_counts(dual_counts, q=Q):
    """Exact coefficient WE of a kernel from its dual-code WE.

    The input contains one count per weight0..length.  Integer
    Krawtchouk sums are combined before division or FLOAT conversion.
    """
    dual_counts = tuple(dual_counts)
    if (type(q) is not int or q < 2 or not dual_counts
            or any(type(x) is not int or x < 0 for x in dual_counts)
            or dual_counts[0] != 1):
        raise ValueError('A nonnegative exact dual weight enumerator is required')
    length, mass = len(dual_counts)-1, sum(dual_counts)
    dimension_mass = mass
    while dimension_mass > 1 and dimension_mass % q == 0:
        dimension_mass //= q
    if dimension_mass != 1 or q ** length % mass:
        raise ValueError('Dual mass must be a field-power divisor of q^length')
    answer = []
    for weight in range(length+1):
        numerator = 0
        for j, count in enumerate(dual_counts):
            if count:
                krawtchouk = sum((-1) ** ell * (q-1) ** (weight-ell)
                                 * comb(j, ell) * comb(length-j, weight-ell)
                                 for ell in range(max(0, weight-(length-j)),
                                                  min(j, weight)+1))
                numerator += count*krawtchouk
        if numerator < 0 or numerator % mass:
            raise ValueError('Dual data do not yield nonnegative integral coefficients')
        answer.append(numerator // mass)
    if answer[0] != 1 or sum(answer)*mass != q ** length:
        raise AssertionError('MacWilliams zero coefficient or mass failed')
    return tuple(answer)


def positive_nonempty_coefficients(outside_counts, fiber_counts):
    """Multiply positive integer WEs, then delete their exact constant1."""
    left, right = tuple(outside_counts), tuple(fiber_counts)
    if (not left or not right or left[0] != 1 or right[0] != 1
            or any(type(x) is not int or x < 0 for x in left+right)):
        raise ValueError('Positive integer polynomials with constant1 required')
    product = [0] * (len(left)+len(right)-1)
    for i, value in enumerate(left):
        if value:
            for j, other in enumerate(right):
                product[i+j] += value*other
    if product[0] != 1:
        raise AssertionError('The zero column was not counted exactly once')
    product[0] = 0
    return tuple(product)


def _log_coefficients(values):
    return tuple(log2(value) if value else -np.inf for value in values)


def prepare_geometry(*, source_masks=None, max_masks=FULL, check_time=lambda: None):
    """Prepare compressed exact classes; None requests every nonempty A.

    No numeric bound is evaluated.  A caller may pass a small explicit
    subset for tests; metadata then forbids a complete all-arity claim.
    check_time is called periodically, including between geometry classes.
    """
    if type(max_masks) is not int or not 1 <= max_masks <= FULL:
        raise ValueError('An explicit positive mask budget up to65535 is required')
    masks = tuple(range(1, FULL+1)) if source_masks is None else tuple(source_masks)
    if (not masks or len(masks) > max_masks or len(set(masks)) != len(masks)
            or any(type(a) is not int or not 1 <= a <= FULL for a in masks)):
        raise ValueError('Unique nonempty labelled masks within the budget required')
    compressed = Counter()
    rank_histogram = Counter()
    for index, a in enumerate(masks):
        if index % 64 == 0:
            check_time()
        arity = a.bit_count()
        outside_counts = outside_weight_counts(a)
        dual_counts = outside_weight_counts(FULL ^ a)
        cross_rank = rank_cross_masks(FULL ^ a, a)
        if sum(dual_counts) != Q ** cross_rank:
            raise AssertionError('Primal and dual cross ranks disagree')
        compressed[arity, outside_counts, dual_counts, cross_rank] += 1
        rank_histogram[arity, cross_rank] += 1
    classes = []
    for (arity, outside, dual, cross_rank), multiplicity in sorted(compressed.items()):
        check_time()
        dimension = arity-cross_rank
        exact_fiber = macwilliams_kernel_counts(dual)
        if sum(exact_fiber) != Q ** dimension:
            raise AssertionError('Actual zero-outside fiber dimension changed')
        coarse_fiber = tuple(comb(dimension, j)*(Q-1)**j for j in range(dimension+1))
        exact_nonempty = positive_nonempty_coefficients(outside, exact_fiber)
        coarse_nonempty = positive_nonempty_coefficients(outside, coarse_fiber)
        if sum(exact_nonempty) != Q ** arity-1 or sum(coarse_nonempty) != Q ** arity-1:
            raise AssertionError('Column mass at t=1 changed')
        classes.append(dict(arity=arity, counts_by_outputweight=outside,
                            dual_counts_by_weight=dual, d=dimension, cross_rank=cross_rank,
                            multiplicity=multiplicity, zero_fiber_counts=exact_fiber,
                            nonempty_coefficients=exact_nonempty,
                            dimension_nonempty_coefficients=coarse_nonempty,
                            nonempty_log2_coefficients=_log_coefficients(exact_nonempty),
                            dimension_nonempty_log2_coefficients=_log_coefficients(coarse_nonempty)))
    complete = len(masks) == FULL
    return dict(classes=classes, source_mask_count=len(masks),
                counts_by_arity={a: sum(c['multiplicity'] for c in classes if c['arity'] == a)
                                 for a in range(1, SIZE+1)},
                rank_histogram={f'{a}:{rank}': count for (a, rank), count in sorted(rank_histogram.items())},
                metadata=dict(construction=CONSTRUCTION, field_order=Q, float_only=True,
                              all_labelled_nonempty_rowsets=complete, labels_retained=True,
                              fixed_rowset_across_columns=True, rowset_resampling=False,
                              exact_local_enumerators=True, scalar_domination_only=True,
                              valid_scalar_interval=[0, 1], whole_code_certificate=False))


def log_polynomial(log_coefficients, log2_t):
    """Stable positive polynomial evaluation; permits t=0 via log2_t=-inf."""
    if np.isnan(log2_t) or log2_t > 0:
        raise ValueError('Scalar atom must lie in[0,1]')
    coefficients = np.asarray(log_coefficients, dtype=float)
    if (coefficients.ndim != 1 or not coefficients.size
            or np.any(np.isnan(coefficients)) or np.any(np.isposinf(coefficients))):
        raise ValueError('Finite-or-negative-infinity coefficient logs required')
    if np.isneginf(log2_t):
        return float(coefficients[0])
    return float(np.logaddexp2.reduce(coefficients + np.arange(coefficients.size)*log2_t))


def nonzero_column_tail_log(log_nonempty, *, columns=COLUMNS, minimum=MINIMUM):
    """Positive binomial sum; the zero-column coefficient is exactly1."""
    if (type(columns) is not int or type(minimum) is not int
            or not 1 <= minimum <= columns or np.isnan(log_nonempty) or np.isposinf(log_nonempty)):
        raise ValueError('A positive finite column range and valid log mass required')
    return float(np.logaddexp2.reduce([log2(comb(columns, j))+j*log_nonempty
                                      for j in range(minimum, columns+1)]))


def group_logs(log2_t, prepared, *, fiber_bound='exact'):
    """Per-arity group log bounds, including all A and RS normalization.

    The multiplicity multiplies each fixed-A WHOLE-WORD binomial tail.
    It is never inserted inside the16-column power.
    """
    if fiber_bound not in ('exact', 'dimension'):
        raise ValueError('fiber_bound must be exact or dimension')
    if np.isnan(log2_t) or log2_t > 0:
        raise ValueError('Scalar atom must lie in[0,1]')
    key = 'nonempty_log2_coefficients' if fiber_bound == 'exact' else 'dimension_nonempty_log2_coefficients'
    result = {a: -np.inf for a in range(1, SIZE+1)}
    for item in prepared['classes']:
        a = item['arity']
        value = (nonzero_column_tail_log(log_polynomial(item[key], log2_t))
                 - INPUT_ROWS*a*log2(Q-1)+log2(item['multiplicity']))
        result[a] = float(np.logaddexp2(result[a], value))
    return result


def scalar_envelope(log_curves, log2_potential, log2_t):
    """Caller-owned all-f entry maxima -> common R_m v<=eta*t^m*v.

    Rows are00,01,10,11 and columns are every m=0..maximum.  The
    caller must use the same declared source-conditioned f range for all
    curves.  No regional provider or slope/potential fit is performed.
    """
    curves = np.asarray(log_curves, dtype=float)
    if (curves.ndim != 2 or curves.shape[0] != 4 or curves.shape[1] < 2
            or np.any(np.isnan(curves)) or np.any(np.isposinf(curves))
            or not np.isfinite(log2_potential) or log2_potential > 0
            or not np.isfinite(log2_t) or log2_t > 0):
        raise ValueError('4x(maximum+1) log curves, 0<p<=1 and 0<t<=1 required')
    rows = np.vstack((np.logaddexp2(curves[0], curves[1]+log2_potential),
                      np.logaddexp2(curves[2]-log2_potential, curves[3])))
    residual = rows - np.arange(curves.shape[1])[None, :]*log2_t
    eta = float(np.max(residual))
    if not np.isfinite(eta):
        raise ValueError('A finite positive regional scalar bound is required')
    return dict(log2_eta=eta, log2_t=float(log2_t), log2_potential=float(log2_potential),
                endpoint_log2=-float(log2_potential), maximum=curves.shape[1]-1,
                row_contraction_log2=rows.tolist(),
                contacts=[dict(state=int(i), occupancy=int(m))
                          for i, m in np.argwhere(residual >= eta-1e-10)])


def allarity_bound(prepared, envelope, *, groups, source_count, source_scope, mark):
    """FLOAT all-cohort bound at one weight, within caller-declared sources."""
    if (not prepared['metadata']['all_labelled_nonempty_rowsets']
            or type(groups) is not int or not 1 <= groups <= envelope['maximum']
            or type(source_count) is not int or source_count < 1
            or not isinstance(source_scope, str) or not source_scope
            or not 0 < mark < 1):
        raise ValueError('Complete geometry and explicit valid source/weight scope required')
    logs = group_logs(envelope['log2_t'], prepared)
    group_total = float(np.logaddexp2.reduce(tuple(logs.values())))
    prefix = (log2(source_count)-LIMIT*log2(mark)-log2(Q-1)
              +SIZE*COLUMNS*envelope['log2_eta']+envelope['endpoint_log2'])
    return dict(bound_log2=prefix+groups*group_total, prefix_log2=prefix,
                group_log2=group_total, group_logs_by_arity=logs, groups=groups,
                source_count=str(source_count), source_scope=source_scope,
                all_arities=True, all_cohorts=True, float_only=True,
                requires_fresh_directed_replay=True, whole_code_certificate=False)
