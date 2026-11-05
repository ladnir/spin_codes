"""FLOAT empty-count extraction at a fixed positive full-pair marker.

For Q(x,y) of ordinary degree D, C(u,t)=u**D Q(1/u,t/u)
marks empty pairs by u and full pairs by t. For P=D*h slots and
histogram (e,s,f), the regional numerator is [u**e t**f] C(u,t)**h.
At fixed t>0 it is bounded by [u**e] C(u,t)**h / t**f.

This diagnostic estimates that coefficient using positive binary64
convolution. It retains underflow diagnostics and is NOT a directed bound.
The old polynomial, setup, original input weight, and prefix are unchanged.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import comb, isfinite, log2
from pathlib import Path

import numpy as np

import hd_field128_gamma_empty_boundary_screen as boundary
import hd_field128_gamma_midpoint_dense_calibrate as dense


def source_pins():
    pins = {**boundary.source_pins(), **dense.source_pins()}
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def empty_layer_logs(composed, log_full_marker, *, slots=32):
    """Positive log coefficients of C(u,t) at fixed t; no coefficient deletion."""
    side = 2*slots+1
    values = np.asarray(composed, dtype=float)
    if (type(slots) is not int or not 1 <= slots <= 32
            or values.shape != (4, side, side) or np.any(np.isnan(values))
            or np.any(np.isposinf(values)) or not isfinite(log_full_marker)
            or abs(log_full_marker) > 64):
        raise ValueError('Invalid polynomial or finite full marker')
    single, full = np.indices((side, side))
    if np.any(np.isfinite(values[:, single+full > 2*slots])):
        raise ValueError('Polynomial exceeds the declared ordinary total degree')
    layers = []
    for empty in range(side):
        degree = np.arange(side-empty)
        terms = values[:, 2*slots-empty-degree, degree]+degree*log_full_marker
        layers.append(np.logaddexp2.reduce(terms, axis=1))
    return np.stack(layers, axis=1)


def regional_estimate(composed, weight, full, log_full_marker, log_empty_marker=0.,
                      *, slots=32, exponent=64):
    pairs = 2*slots*exponent
    empty, single, full = dense.histogram(pairs, weight, full)
    logs = empty_layer_logs(composed, log_full_marker, slots=slots)
    coefficient, diagnostics = boundary.coefficient_logs(
        logs, exponent, empty, log_marker=log_empty_marker)
    normalizer = log2(comb(pairs, empty)*comb(pairs-empty, full))
    estimate = coefficient-full*log_full_marker-normalizer
    positive = boundary.positive_marker_logs(logs, exponent, empty, log_empty_marker)
    positive -= full*log_full_marker+normalizer
    return estimate, positive, diagnostics


def original_weight_prefix(z, alpha, weight):
    if type(weight) is not int or not 1 <= weight <= 8192:
        raise ValueError('Expected original weight1..8192')
    z, alpha = float(z), float(alpha)
    if not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError('Invalid finite marks')
    old = dense.provider.old
    forbidden = int(F(2, 5)*(1 << 21))
    density = old.outer_log_density(1 << 128, 'invertible-toeplitz')
    return alpha*(density*weight-128-old.log1mexp2(-128)-forbidden*log2(z)), forbidden


def evaluate_profile(z, alpha, weight, full, *, max_seconds=30., max_calls=128):
    """One symmetric profile, fresh local coefficients, no regional alias bank."""
    empty, single, full = dense.histogram(4096, weight, full)
    if (not 4096 < weight < 8192 or empty < 1 or single < 1
            or not isfinite(max_seconds) or not 0 < max_seconds <= 60
            or type(max_calls) is not int or not 8 <= max_calls <= 128):
        raise ValueError('Expected a bounded dense interior profile')
    z, alpha = F(z), F(alpha)
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError('Invalid exact marks')
    provider = dense.provider
    pins = source_pins()
    budget = provider.old.Budget(max_seconds)
    coeff = provider.AliasCoefficients(float(z), float(alpha), state=30, _budget=budget)
    def numerator(point):
        budget.check()
        return np.asarray(provider.alias.power(coeff.pair_logs(point, _budget=budget), 64)).ravel()
    rows = []
    for entry in range(4):
        positive, fit = dense.fit_entry(numerator, 4096, weight, full, entry, [], budget,
                                        max_calls=max_calls)
        x, y = map(float, positive)
        eta, tilt = log2(y/x), -log2(x)
        if abs(eta) > 64 or abs(tilt) > 64:
            raise ArithmeticError('Fitted chart exceeds the diagnostic convolution bounds')
        estimate, point, diagnostic = regional_estimate(coeff.composed, weight, full, eta, tilt)
        if not isfinite(float(estimate[entry])):
            raise ArithmeticError('Coefficient estimate is zero or nonfinite; not a usable estimate')
        rows.append(dict(entry=entry, fit=fit, full_marker_log2_FLOAT=eta,
                         empty_tilt_log2_FLOAT=tilt,
                         coefficient_estimate_log2_FLOAT=float(estimate[entry]),
                         positive_marker_log2_FLOAT=float(point[entry]),
                         estimate_diagnostics=diagnostic))
    counts = dense.mode.MidpointCounts(4096, weight, max_singletons=64, budget=budget)
    count = counts.count(full, full)
    transfer = provider.profiles.transfer
    estimated = np.array([r['coefficient_estimate_log2_FLOAT'] for r in rows]).reshape(2, 2)
    positive = np.array([r['positive_marker_log2_FLOAT'] for r in rows]).reshape(2, 2)
    prefix, forbidden = original_weight_prefix(z, alpha, weight)
    def profile_score(matrix):
        product = transfer.power(matrix, 256)
        return prefix+log2(count)+float(np.logaddexp2(product[0, 0], product[0, 1]))
    directory = Path(__file__).resolve().parent
    if source_pins() != pins or any(sha256((directory/name).read_bytes()).hexdigest() != digest
                                   for name, digest in pins.items()):
        raise RuntimeError('Diagnostic source changed')
    return dict(schema='gamma-dense-empty-extraction-FLOAT-v1',
                construction=provider.gamma.construction(30), state_dimension=30,
                target_distance='2/5', active_groups=weight, even_full=full, odd_full=full,
                histogram=dict(empty=empty, single=single, full=full),
                witness=dict(z=str(z), alpha=str(alpha)), source_sha256=pins, entries=rows,
                log2_outer_prefix=prefix, forbidden_output_weight=forbidden,
                log2_profile_count_upper=log2(count),
                positive_marker_profile_score_FLOAT=profile_score(positive),
                coefficient_estimate_profile_score_FLOAT=profile_score(estimated),
                full_regional_bank_built=False, original_weight_prefix_used=True,
                count_complement_only=True, chronological_product_before_alpha=True,
                coefficient_estimate_is_upper_bound=False, directed=False,
                selected_weight_certificate=False, whole_code_certificate=False,
                scope='One symmetric interior profile only; positive convolution estimate, not a certificate')
