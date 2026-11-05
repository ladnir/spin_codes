"""FLOAT diagnostics for the zero-empty-pair boundary of frozen all-pair caps.

Let Q(x,y) be the existing paired matrix polynomial, of ordinary total
degree at most64 (weighted occupancy degree at most128). For B(t) equal to
its degree64 homogeneous part evaluated at(1,t), the regional boundary cap
coefficient is [t**f] B(t)**64 / C(4096,f), f=r-4096. Matrix products remain
chronological; the raw occupancy product was taken before entrywise alpha
inside Q. No extra singleton-orientation factor occurs.

A positive finite marker t bounds the coefficient by B(t)**64/t**f.
The convolution result below is only a FLOAT estimate, not a directed
upper bound. Underflow and computed-zero ambiguities are explicit.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import comb, isfinite, log2
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

import hd_field128_gamma_all_pair_screen as provider


NEG = float('-inf')


def source_pins():
    pins = provider.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def top_logs(composed, *, slots=32):
    """Extract ordinary degree2*slots, not maximum weighted occupancy."""
    side = 2*slots+1
    values = np.asarray(composed, dtype=float)
    if (type(slots) is not int or not 1 <= slots <= 32 or values.shape != (4, side, side)
            or np.any(np.isnan(values)) or np.any(np.isposinf(values))):
        raise ValueError('Invalid complete paired log coefficients')
    degree = np.arange(side)
    return values[:, 2*slots-degree, degree].copy()


def positive_marker_logs(log_coefficients, exponent, degree, log_marker):
    """FLOAT evaluation of the analytic entrywise positive-marker bound."""
    logs = np.asarray(log_coefficients, dtype=float)
    if (logs.ndim != 2 or logs.shape[0] != 4 or logs.shape[1] < 1
            or np.any(np.isnan(logs)) or np.any(np.isposinf(logs))
            or type(exponent) is not int or not 1 <= exponent <= 128
            or type(degree) is not int or not 0 <= degree <= (logs.shape[1]-1)*exponent
            or not isfinite(log_marker) or abs(log_marker) > 64):
        raise ValueError('Invalid bounded univariate marker request')
    value = np.logaddexp2.reduce(logs+np.arange(logs.shape[1])*log_marker, axis=1).reshape(2, 2)
    return np.asarray(provider.alias.power(value, exponent)).ravel()-degree*log_marker


def coefficient_logs(log_coefficients, exponent, degree, *, log_marker=0., reverse=True):
    """Scaled positive np.convolve estimate, truncated at the requested degree.

    Each matrix entry has its own positive coefficient vector and log scale.
    Scaling keeps products finite but does not make binary64 directed. Exact
    reversal of this same polynomial reduces degree; it does not complement
    the encoder, local caps, source weight, or outside prefix.
    """
    logs = np.asarray(log_coefficients, dtype=float).copy()
    # Reuse the interface checks without assuming a positive target coefficient.
    positive_marker_logs(logs, exponent, degree, log_marker)
    maximum = (logs.shape[1]-1)*exponent
    reversed_polynomial = bool(reverse and degree > maximum-degree)
    if reversed_polynomial:
        logs = logs[:, ::-1]
        degree = maximum-degree
        log_marker = -log_marker
    diagnostics = dict(input_underflow_terms=0, scale_underflow_events=0,
                       normalization_underflow_terms=0, addition_underflow_terms=0,
                       convolution_underflow_possible=False,
                       convolution_zero_with_positive_factors=False,
                       nonfinite=False, zero_coefficient_ambiguity=False)
    def normalized(values, scale):
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            diagnostics['nonfinite'] = True
            raise ArithmeticError('Nonfinite or negative positive-convolution result')
        maximum_value = float(np.max(values))
        if maximum_value == 0:
            return np.zeros(1), NEG
        with np.errstate(under='ignore'):
            scaled = values/maximum_value
        diagnostics['normalization_underflow_terms'] += int(np.count_nonzero((values > 0) & (scaled == 0)))
        return scaled, scale+log2(maximum_value)
    def make(row):
        row = row[:degree+1]
        finite = np.isfinite(row)
        if not np.any(finite):
            return np.zeros(1), NEG
        scale = float(np.max(row[finite]))
        with np.errstate(under='ignore'):
            values = np.exp2(row-scale)
        diagnostics['input_underflow_terms'] += int(np.count_nonzero(finite & (values == 0)))
        return normalized(values, scale)
    def product(a, b):
        if a[1] == NEG or b[1] == NEG:
            return np.zeros(1), NEG
        minimum_log = log2(float(np.min(a[0][a[0] > 0])))+log2(float(np.min(b[0][b[0] > 0])))
        diagnostics['convolution_underflow_possible'] |= minimum_log < -1074
        with np.errstate(under='ignore'):
            value = np.convolve(a[0], b[0])[:degree+1]
        if not np.any(value) and np.any(a[0]) and np.any(b[0]):
            # May also be an exact zero after truncation; never assert which.
            diagnostics['convolution_zero_with_positive_factors'] = True
        return normalized(value, a[1]+b[1])
    def add(a, b):
        if a[1] == NEG:
            return b
        if b[1] == NEG:
            return a
        scale = max(a[1], b[1])
        value = np.zeros(max(len(a[0]), len(b[0])))
        for row, old_scale in (a, b):
            with np.errstate(under='ignore'):
                factor = float(np.exp2(old_scale-scale))
            diagnostics['scale_underflow_events'] += int(factor == 0)
            with np.errstate(under='ignore'):
                scaled = factor*row
            diagnostics['addition_underflow_terms'] += int(np.count_nonzero((row > 0) & (scaled == 0)))
            value[:len(row)] += scaled
        return normalized(value, scale)
    def matrix_product(a, b):
        return tuple(add(product(a[2*i], b[j]), product(a[2*i+1], b[2+j]))
                     for i in range(2) for j in range(2))
    tilt = np.arange(logs.shape[1])*log_marker
    matrix = tuple(make(row+tilt) for row in logs)
    zero, one = (np.zeros(1), NEG), (np.ones(1), 0.)
    result, power = (one, zero, zero, one), exponent
    while power:
        if power & 1:
            result = matrix_product(result, matrix)
        power //= 2
        if power:
            matrix = matrix_product(matrix, matrix)
    answer = []
    for values, scale in result:
        value = float(values[degree]) if degree < len(values) else 0.
        if value == 0:
            diagnostics['zero_coefficient_ambiguity'] = True
            answer.append(NEG)
        else:
            answer.append(log2(value)+scale-degree*log_marker)
    diagnostics.update(polynomial_reversed=reversed_polynomial,
                       computed_degree=degree, maximum_retained_vector_length=degree+1,
                       intermediate_convolution_length_bound=2*degree+1,
                       omitted_underflow_correction=True, directed=False,
                       convolution_underflow_test='conservative minimum-positive-factor range; not a lost-term count',
                       rounding_direction_unknown=True,
                       estimated_zeros_are_not_proved=True)
    return np.array(answer), diagnostics


def evaluate_boundary(z, alpha, weight, *, state=30, max_seconds=20., max_calls=96):
    """One boundary regional matrix only; never builds a regional alias bank."""
    if (type(weight) is not int or not 4096 <= weight <= 8192
            or not isfinite(max_seconds) or not 0 < max_seconds <= 60
            or type(max_calls) is not int or not 8 <= max_calls <= 256):
        raise ValueError('Expected one boundary weight and bounded resources')
    z, alpha = F(z), F(alpha)
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError('Invalid exact mark proposals')
    pins, directory = source_pins(), Path(__file__).resolve().parent
    budget = provider.old.Budget(max_seconds)
    coefficients = provider.AliasCoefficients(float(z), float(alpha), state=state, _budget=budget)
    logs = top_logs(coefficients.composed)
    full, single = weight-4096, 8192-weight
    normalizer = log2(comb(4096, full))
    cache, rows = {}, []
    def values(eta):
        eta = float(eta)
        budget.check()
        if eta not in cache:
            cache[eta] = positive_marker_logs(logs, 64, full, eta)
        return cache[eta]
    for entry in range(4):
        best_eta, best_value = 0., float(values(0.)[entry])
        calls = 0
        def objective(eta):
            nonlocal best_eta, best_value, calls
            calls += 1
            value = float(values(eta)[entry])
            if value < best_value:
                best_eta, best_value = float(eta), value
            return value
        fit = minimize_scalar(objective, bounds=(-64., 64.), method='bounded',
                              options=dict(maxiter=max_calls, xatol=1e-10))
        for endpoint in (-64., 64.):
            value = float(values(endpoint)[entry])
            if value < best_value:
                best_eta, best_value = endpoint, value
        budget.check()
        estimate, diagnostic = coefficient_logs(logs, 64, full, log_marker=best_eta)
        rows.append(dict(entry=entry, marker_log2_exact=str(F(best_eta)),
                         positive_marker_exact=str(F(2.**best_eta)),
                         positive_marker_role='rounded exp2 proposal; displayed FLOAT value uses fitted log marker',
                         marker_log2_float=best_eta, marker_bound_log2_FLOAT=best_value-normalizer,
                         coefficient_estimate_log2_FLOAT=float(estimate[entry]-normalizer),
                         marker_calls=calls, optimizer_success=bool(fit.success),
                         estimate_diagnostics=diagnostic))
    budget.check()
    if source_pins() != pins or any(sha256((directory/name).read_bytes()).hexdigest() != digest
                                   for name, digest in pins.items()):
        raise RuntimeError('Boundary source changed during the diagnostic')
    return dict(schema='gamma-empty-boundary-univariate-FLOAT-v1',
                construction=provider.gamma.construction(state), state_dimension=state,
                active_groups=weight, single_pairs=single, full_pairs=full, empty_pairs=0,
                paired_polynomial_power=64, regional_slots=4096,
                normalization='binomial(4096,full_pairs); no singleton orientation factor',
                witness=dict(z=str(z), alpha=str(alpha)), entries=rows,
                top_log_coefficients_sha256=sha256(logs.tobytes()).hexdigest(),
                source_sha256=pins, maximum_marker_calls_per_entry=max_calls+3,
                distinct_marker_evaluations=len(cache), max_seconds=max_seconds,
                full_regional_bank_built=False, chronological_product_before_alpha=True,
                original_weight_prefix_required=True, coefficient_estimate_is_upper_bound=False,
                directed=False, selected_weight_certificate=False, whole_code_certificate=False,
                scope='Only the zero-empty regional matrix; no all-profile sum or distance certificate')
