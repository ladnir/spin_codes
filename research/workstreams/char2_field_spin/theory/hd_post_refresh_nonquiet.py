"""Positive nonquiet remainder for the saved post-refresh moment interface.

For arity eight among sixteen coordinates, masks of weight <=8 are quiet
and masks of weight >=9 are nonquiet. Keep that event as an OR flag through
the group product and the ordered column-matrix product. No subtraction
of full and quiet bounds occurs. Raw fixed-rowset Holder weights remain
unchanged; this module only separates a positive part of their majorant.

No provider, parameter fit, artifact load, or production gate runs here.
The quiet flag can be replaced by a separately proved coherent-rowset
bound, then added once to this nonquiet remainder. Neither component alone
is a whole-class certificate. All arithmetic here is FLOAT screening.
"""
from math import comb, expm1, log, log2

import numpy as np

from hd_two_sum_gf16_mean_column import raw_holder_weights
from hd_two_sum_gf16_state_envelope import elementary_logs, path_classes
from hd_two_sum_late_alpha_quiet import log_convolve


QUIET, NONQUIET = 0, 1
OR_PRODUCTS = ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 1))


def _logs(values):
    values = np.asarray(values, dtype=float)
    if np.any(np.isnan(values)) or np.any(np.isposinf(values)):
        raise ValueError('Finite-or-zero log coefficients required')
    return values


def flag_product(left, right, check_time=lambda: None):
    """Positive OR product of two scalar polynomials, shape(flag,degree)."""
    left, right = _logs(left), _logs(right)
    if (left.ndim != 2 or right.ndim != 2 or left.shape[0] != 2
            or right.shape[0] != 2 or not left.shape[1] or not right.shape[1]):
        raise ValueError('Two complete flagged scalar polynomials required')
    result = np.full((2, left.shape[1]+right.shape[1]-1), -np.inf)
    for first, second, target in OR_PRODUCTS:
        if np.any(np.isfinite(left[first])) and np.any(np.isfinite(right[second])):
            term = log_convolve(left[first], right[second], check_time)
            np.logaddexp2(result[target], term, out=result[target])
    return result


def flag_power(values, exponent, check_time=lambda: None):
    values = _logs(values)
    if (values.ndim != 2 or values.shape[0] != 2 or not values.shape[1]
            or type(exponent) is not int or exponent < 0):
        raise ValueError('Flagged scalar polynomial and nonnegative exponent required')
    result = np.array([[0.], [-np.inf]])
    base = values
    while exponent:
        check_time()
        if exponent & 1:
            result = flag_product(result, base, check_time)
        exponent >>= 1
        if exponent:
            base = flag_product(base, base, check_time)
    return result


def group_polynomials(weights, counts, active00, measure, arity):
    """Flag by true output-mask weight BEFORE merging surrogate degrees."""
    weights = _logs(weights)
    if (len(counts) != 4 or any(type(n) is not int or n < 0 for n in counts)
            or weights.shape != (sum(counts)+1,) or type(active00) is not int
            or not 0 <= active00 <= counts[0] or type(arity) is not int
            or not 1 <= arity < len(weights)):
        raise ValueError('Complete mask weights and feasible atom counts required')
    _, n01, n10, n11 = counts
    markers = np.array([measure['log2_atom00']]*active00
        +[measure['log2_atom01']]*n01+[measure['log2_atom10']]*n10+[0.]*n11)
    if np.any(~np.isfinite(markers)):
        raise ValueError('Finite nonzero atom logs required')
    elementary = elementary_logs(markers)
    result = np.full((2, arity+1), -np.inf)
    for weight, value in enumerate(elementary):
        flag, degree = int(weight > arity), min(weight, arity)
        result[flag, degree] = np.logaddexp2(result[flag, degree], weights[weight]+value)
    return result


def column_matrices(weights, groups, measure, *, arity, check_time=lambda: None):
    """Ordered column matrix with a positive mask-event flag.

    The state-path multiplicities and atom masses stay outside the group
    polynomial. Its two flag coefficients are raised to groups in the OR
    semiring. Each cache key keeps directed 01/10 counts distinct.
    """
    weights = _logs(weights)
    size = len(weights)-1
    if (weights.ndim != 1 or not 1 <= size <= 16 or type(groups) is not int
            or groups < 1 or type(arity) is not int or not 1 <= arity <= size):
        raise ValueError('Complete weights and positive group/arity scope required')
    logp = measure['log2_probability']
    if not np.isfinite(logp) or logp > 0:
        raise ValueError('00 nonzero-atom probability must belong to (0,1]')
    if any(not np.isfinite(measure[key]) for key in
           ('log2_scale01', 'log2_scale10', 'log2_scale11')):
        raise ValueError('Finite positive entry scales required')
    logother = -np.inf if logp == 0 else log2(-expm1(logp*log(2)))
    result = np.full((2, 2, 2, groups*arity+1), -np.inf)
    cached, terms = {}, 0
    classes = path_classes(size)
    for initial, final, counts, multiplicity in classes:
        check_time()
        n00, n01, n10, n11 = counts
        state_mass = (log2(multiplicity)+n01*measure['log2_scale01']
            +n10*measure['log2_scale10']+n11*measure['log2_scale11'])
        for active00 in range(n00+1):
            if logp == 0 and active00 != n00:
                continue
            atom_mass = (log2(comb(n00, active00))+active00*logp
                +((n00-active00)*logother if active00 < n00 else 0.))
            key = active00, n01, n10, n11
            if key not in cached:
                polynomial = group_polynomials(weights, counts, active00, measure, arity)
                cached[key] = flag_power(polynomial, groups, check_time)
            target = result[:, initial, final]
            np.logaddexp2(target, cached[key]+atom_mass+state_mass, out=target)
            terms += 1
    return result, dict(path_classes=len(classes), positive_path_atom_terms=terms,
        distinct_group_polynomials=len(cached), separate_cross_counts_in_cache=True)


def matrix_product(left, right, check_time=lambda: None):
    """Chronological 2x2 matrix product over positive flagged polynomials."""
    left, right = _logs(left), _logs(right)
    if (left.ndim != 4 or right.ndim != 4 or left.shape[:3] != (2, 2, 2)
            or right.shape[:3] != (2, 2, 2) or not left.shape[3] or not right.shape[3]):
        raise ValueError('Complete flag/row/column/degree arrays required')
    result = np.full((2, 2, 2, left.shape[3]+right.shape[3]-1), -np.inf)
    for row in range(2):
        for middle in range(2):
            for column in range(2):
                check_time()
                term = flag_product(left[:, row, middle], right[:, middle, column], check_time)
                target = result[:, row, column]
                np.logaddexp2(target, term, out=target)
    return result


def matrix_power(matrix, exponent, check_time=lambda: None):
    matrix = _logs(matrix)
    if (matrix.ndim != 4 or matrix.shape[:3] != (2, 2, 2) or not matrix.shape[3]
            or type(exponent) is not int or exponent < 0):
        raise ValueError('Flagged 2x2 matrix polynomial and nonnegative exponent required')
    result = np.full((2, 2, 2, 1), -np.inf)
    result[0, 0, 0, 0] = result[0, 1, 1, 0] = 0.
    base = matrix
    while exponent:
        check_time()
        if exponent & 1:
            result = matrix_product(result, base, check_time)
        exponent >>= 1
        if exponent:
            base = matrix_product(base, base, check_time)
    return result


def tail_logs(word, floor):
    word = _logs(word)
    if (word.ndim != 4 or word.shape[:3] != (2, 2, 2) or type(floor) is not int
            or not 0 <= floor < word.shape[3]):
        raise ValueError('Complete flagged word and an in-range nonnegative floor required')
    return np.array([np.logaddexp2.reduce(word[flag, 0, :, floor:].ravel())
                     for flag in range(2)])


def post_refresh_remainder(measure, prefix, *, check_time=lambda: None):
    """Fixed r18/a8 component under a caller-validated raw post-refresh law.

    Prefix includes the exact Fmax1 source count, z^-LIMIT, projective255^-1
    and RS255^-1152 exactly once. This routine neither authenticates that
    witness nor evaluates a provider. Its nonquiet result must be added
    to an independently valid quiet bound before claiming this class.
    """
    if not np.isfinite(prefix):
        raise ValueError('One finite once-only prefix required')
    column, costs = column_matrices(raw_holder_weights(), 18, measure,
        arity=8, check_time=check_time)
    word = matrix_power(column, 16, check_time)
    tails = tail_logs(word, 1296)
    return dict(quiet_holder_log2=float(prefix+tails[0]),
        nonquiet_holder_log2=float(prefix+tails[1]),
        full_holder_log2=float(prefix+np.logaddexp2(tails[0], tails[1])),
        component='at least one nonquiet output mask among all groups and columns',
        flag_order=['all-quiet', 'at-least-one-nonquiet'],
        active_groups=18, arity=8, columns=16, aggregate_surrogate_floor=1296,
        fixed_rowsets_across_columns=True, raw_holder_exponent=16,
        extra_rowset_factor=False, initial_state=0, terminal_vector=[1, 1],
        positive_OR_arithmetic=True, subtraction_used=False,
        complete_homogeneous_class=False, quiet_complement_required=True,
        raw_regional_banks_built=0, parameter_fits=0, directed=False,
        whole_code_certificate=False, **costs)
