"""Support-checked scaled floating placement, for witness proposals only.

Ordinary floating addition can discard a summand negligible relative to its
destination. Log-sum-exp has the same limitation. This helper permits that
rounding, but tracks the exact Boolean support of every completed matrix
coefficient. If a mathematically positive entry becomes zero, it refuses the
result. The caller must then use full log-domain arithmetic or an Arb replay.
There is no magnitude cutoff, pruning of states, or outward-rounding claim.

Positive hypergeometric coefficients give an exact support recurrence:
the support of the next coefficient is the Boolean union of the ordered
matrix products that contribute. Matrix products remain in chronological
order. All positive final entries are logged before subsequent mixtures.
"""
from math import comb, log

import numpy as np
from scipy.special import gammaln

from packet_regional_log import log_placement, regional_uniform_log, log_power_matrix
from packet_outer_geometry_proposal import upper_arrays
from packet_outer_cost_search import geometry, _occupancies, _tilt


def placement_logs(operators, degree, *, epochs, windows=32, rowwise=False):
    values = np.asarray(operators, dtype=float)
    if (values.ndim != 3 or len(values) != windows+1 or not values.shape[1]
            or values.shape[1] != values.shape[2] or not np.isfinite(values).all()
            or np.any(values < 0) or type(degree) is not int or degree < 0
            or type(epochs) is not int or epochs < 1 or type(windows) is not int
            or windows < 1 or degree > epochs*windows):
        raise ValueError('nonnegative complete local family and feasible placement required')
    size = values.shape[1]
    local_scale = values.max(axis=(1, 2))
    if np.any(local_scale <= 0):
        raise FloatingPointError('positive local scaling factors required')
    local_support = values > 0
    with np.errstate(over='raise', invalid='raise', divide='raise', under='ignore'):
        local = values / local_scale[:, None, None]
        if not np.array_equal(local > 0, local_support):
            raise FloatingPointError('local normalization lost a positive entry')
        local_logs = np.log(local_scale)
        current = np.eye(size)[None, :, :]
        scales = np.zeros((1, size)) if rowwise else np.zeros(1)
        for epoch in range(1, epochs+1):
            limit = min(degree, epoch*windows)
            total, previous = epoch*windows, (epoch-1)*windows
            terms = []
            base = np.full((limit+1, size) if rowwise else limit+1, -np.inf)
            for k in range(min(windows, limit)+1):
                js = np.arange(k, min(limit, k+len(current)-1)+1)
                old = js-k
                weights = (log(comb(windows, k))+gammaln(previous+1)
                    -gammaln(old+1)-gammaln(previous-old+1)-gammaln(total+1)
                    +gammaln(js+1)+gammaln(total-js+1))
                exponents = (weights[:, None] if rowwise else weights)+scales[old]+local_logs[k]
                base[js] = np.maximum(base[js], exponents)
                terms.append((k, js, old, exponents))
            if not np.isfinite(base).all():
                raise FloatingPointError('nonfinite placement normalization')
            nxt = np.zeros((limit+1, size, size))
            support = np.zeros(nxt.shape, dtype=bool)
            for k, js, old, exponents in terms:
                # Every hypergeometric coefficient here is mathematically
                # positive, even if this individual scaled summand rounds off.
                support[js] |= np.matmul(current[old] > 0, local_support[k])
                weights = np.exp(exponents-base[js])
                nxt[js] += (weights[:, :, None] if rowwise else weights[:, None, None])*(current[old] @ local[k])
            if not np.array_equal(nxt > 0, support):
                raise FloatingPointError('a completed positive coefficient entry became zero')
            scale = nxt.max(axis=2 if rowwise else (1, 2))
            if not np.isfinite(scale).all() or np.any(scale <= 0):
                raise FloatingPointError('nonfinite or zero placement coefficient')
            current = nxt/(scale[:, :, None] if rowwise else scale[:, None, None])
            scales = base+np.log(scale)
            if not np.array_equal(current > 0, support):
                raise FloatingPointError('renormalization lost a positive entry')
        result = np.full(current.shape, -np.inf)
        positive = current > 0
        result[positive] = np.log(current[positive])
        result += scales[:, :, None] if rowwise else scales[:, None, None]
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise FloatingPointError('invalid final placement logarithm')
    return result


def estimate(local, *, K, envelope, occupancies, tilt, windows=32):
    """Support-checked proposal, automatically falling back to full log DP."""
    from packet_regional_log import upper_logs
    meta = geometry(K, envelope, windows=windows)
    qs, tilt = _occupancies(occupancies, meta['groups']), _tilt(tilt)
    backend = 'support-checked-row-scaled'
    try:
        regional = placement_logs(upper_arrays(local), qs[-1],
            epochs=meta['epochs_per_region'], windows=windows, rowwise=True)
    except FloatingPointError:
        backend = 'log'
        regional = log_placement(upper_logs(local), qs[-1],
            epochs=meta['epochs_per_region'], windows=windows)
    log_beta = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    witnesses = {}
    for q in qs:
        moment = log_power_matrix(regional_uniform_log(regional, q), meta['regions'])
        value = log(comb(meta['groups'], q))+q*log_beta+float(tilt)*meta['cutoff']+moment
        witnesses[str(q)] = dict(estimated_margin_bits=-value/log(2), log_moment=moment)
    return dict(meta, schema='larger-outer-support-checked-proposal-1', tilt=str(tilt),
        backend=backend, witnesses=witnesses, proposal_only=True, whole_code_certificate=False)
