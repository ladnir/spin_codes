"""Positive full-a8 column contraction under one scalar moment envelope.

No provider, fit, artifact loading, or production gate runs here. Inputs
are reviewed raw alpha1 regional log matrices and a common two-atom law.
Fixed row-set choices enter one16-column Holder weight, never a fresh
row-set sample. Raw weights deliberately EXCLUDE the RS denominator.
"""
from math import comb, expm1, log, log2

import numpy as np

from hd_two_sum_closed_kernel import arity_power_sum, zero_two_moments
from hd_two_sum_gf16_single_dense_matrix import log_matmul
from hd_two_sum_late_alpha_quiet import log_power


def _logs(values):
    result = np.asarray(values, dtype=float)
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Finite-or-zero log entries required')
    return result


def raw_holder_weights(*, quiet=False, group=16, arity=8, columns=16):
    """log2(sum_A N_A(T)^columns)/columns; NO RS denominator included.

    For half arity the quiet masks have w<=arity, and every nonquiet
    mask has w>=arity+1. Hence kappa(w)=min(w,arity) is independent ofA.
    """
    if (group not in (8, 16) or type(arity) is not int or arity*2 != group
            or type(columns) is not int or columns < 1 or type(quiet) is not bool):
        raise ValueError('Half-arity g8/g16 and positive Holder exponent required')
    result = []
    for weight in range(group+1):
        if quiet:
            moment = (comb(group-weight, arity-weight)*zero_two_moments(weight)**columns
                      if weight <= arity else 0)
        else:
            moment = arity_power_sum(arity, weight, group, columns)
        result.append(log2(moment)/columns if moment else -np.inf)
    return np.asarray(result)


def two_atom_moments(maximum, log_probability, log_atom):
    """L=atom with probabilityp, else0; moment0 is exactly1."""
    if (type(maximum) is not int or maximum < 1
            or not np.isfinite(log_probability) or log_probability > 0
            or not np.isfinite(log_atom)):
        raise ValueError('Positive finite two-atom parameters with p<=1 required')
    result = log_probability+log_atom*np.arange(maximum+1)
    result[0] = 0.
    return result


def common_matrix(regions, log_probability, log_atom):
    """Smallest entrywise C for the GIVEN M_m, using all declared m.

    This is a finite envelope calculation, not a slope fit. Zero entries
    remain zero when they vanish at every occupancy. Nonzero C01/C10
    can create empty transitions: this deliberate relaxation is visible.
    """
    regions = _logs(regions)
    if regions.ndim != 3 or regions.shape[1:] != (2, 2) or len(regions) < 2:
        raise ValueError('Complete occupancy0..r log2x2 matrices required')
    moments = two_atom_moments(len(regions)-1, log_probability, log_atom)
    ratios = regions-moments[:, None, None]
    matrix = np.max(ratios, axis=0)
    return matrix, moments


def state_factor(matrix, steps):
    """log2(e0 C^steps ones); state chronology is ordinary multiplication."""
    matrix = _logs(matrix)
    if matrix.shape != (2, 2) or type(steps) is not int or steps < 0:
        raise ValueError('One2x2 matrix and nonnegative integer steps required')
    answer = np.array([[0., -np.inf], [-np.inf, 0.]])
    base = matrix
    while steps:
        if steps & 1:
            answer = log_matmul(answer, base)
        steps >>= 1
        if steps:
            base = log_matmul(base, base)
    return float(np.logaddexp2(answer[0, 0], answer[0, 1]))


def column_polynomial(weights, groups, log_probability, log_atom, *,
                      arity=8, check_time=lambda: None):
    """Exact positive latent contraction in FLOAT, raw RS normalization.

    Q(X)=sum_d Bin(n,p)[d]*(sum_w<=d C(d,w)W_w B^w
                           X^min(w,arity))^groups.
    Binomial factors count actual labelled atom subsets and mask subsets;
    no extra row-set factor or fractional count power is introduced.
    """
    weights = _logs(weights)
    if (weights.ndim != 1 or len(weights) < 2 or type(groups) is not int
            or groups < 1 or type(arity) is not int or not 1 <= arity < len(weights)):
        raise ValueError('One complete mask-weight table and positive groups/arity required')
    two_atom_moments(groups, log_probability, log_atom)
    size = len(weights)-1
    log_other = (-np.inf if log_probability == 0
                 else log2(-expm1(log_probability*log(2))))
    answer = np.full(groups*arity+1, -np.inf)
    for active in range(size+1):
        check_time()
        if log_probability == 0 and active != size:
            continue
        log_mass = (log2(comb(size, active))+active*log_probability
                    +(size-active)*log_other if active < size
                    else size*log_probability)
        factor = np.full(arity+1, -np.inf)
        for weight in range(active+1):
            degree = min(weight, arity)
            term = log2(comb(active, weight))+weights[weight]+weight*log_atom
            factor[degree] = np.logaddexp2(factor[degree], term)
        value = log_power(factor, groups, check_time)
        np.logaddexp2(answer, value+log_mass, out=answer)
    return answer


def whole_word_screen(weights, matrix, groups, log_probability, log_atom, *,
                      log_source_count, mark, limit, columns=16, arity=8,
                      dimension=8, minimum=9, check_time=lambda: None):
    """Complete declared homogeneous class under supplied raw envelopes.

    The caller must establish the source stratum and R1<=M*C. This
    function cannot authenticate a bank. It is FLOAT, never a certificate.
    RS normalization255^(-dimension*arity*groups) and projective255^-1
    each occur once. There is NO external C(group,arity)^groups factor.
    """
    if (type(columns) is not int or columns < 1 or type(dimension) is not int
            or not 1 <= dimension < columns or type(minimum) is not int
            or not 1 <= minimum <= columns or not 0 < mark < 1
            or type(limit) is not int or limit < 0 or not np.isfinite(log_source_count)):
        raise ValueError('Valid constituent parameters and one finite witness required')
    weights = _logs(weights)
    raw = column_polynomial(weights, groups, log_probability, log_atom,
                            arity=arity, check_time=check_time)
    word = log_power(raw, columns, check_time)
    floor = minimum*arity*groups
    tail = float(np.logaddexp2.reduce(word[floor:]))
    states = state_factor(matrix, (len(weights)-1)*columns)
    denominator = dimension*arity*groups*log2(255)
    prefix = log_source_count-limit*log2(mark)-log2(255)-denominator
    return dict(log2_screen=prefix+tail+states, raw_word_tail_log2=tail,
        state_factor_log2=states, prefix_log2=prefix, rs_normalizer_log2=denominator,
        aggregate_surrogate_floor=floor, surrogate='min(output_weight,arity)',
        columns=columns, regions_per_column=len(weights)-1, groups=groups,
        row_sets_fixed_across_columns=True, rowset_choices_inside_holder=True,
        extra_rowset_factor=False, raw_weights_exclude_rs_denominator=True,
        alpha=1, initial_state=0, final_vector=[1, 1],
        exact_offdiagonal_empty_zeros_retained=False,
        individual_row_floors=False, directed=False, whole_code_certificate=False)
