"""Independent floating diagnosis of all-zero-state trajectories.

The calculation isolates a positive summand of the current uniform-input
outer envelope. It is not a certificate or a codeword counterexample.
Selected group slots carry uniform bytes INCLUDING zero, so Fourier factors
are nonnegative and no signed cancellation is used in this diagnostic.
"""
from __future__ import annotations

import argparse
from math import comb, log

import numpy as np
from scipy.special import logsumexp

import independent_checks as algebra


def character_weights(points):
    states = np.arange(65536, dtype=np.uint32)
    low, high = states & 255, states >> 8
    return np.stack([np.bitwise_count(low ^ np.asarray(
        [algebra.adjoint_multiply(alpha, byte) for byte in range(256)],
        dtype=np.uint32)[high]) for alpha in points], axis=1)


def marked_zero_polynomial(weights, tilt):
    """Coefficients C(W,r) E[z^wtX 1{CX=0}] for r uniformly placed slots."""
    z = np.exp(-np.longdouble(tilt))
    packet = (1+z)**(8-weights)*(1-z)**weights/256
    states, width = weights.shape
    polynomial = np.zeros((states, width+1), dtype=np.longdouble)
    polynomial[:, 0] = 1
    for position in range(width):
        # Descending coefficients retain the preceding physical-packet stage.
        for degree in range(position+1, 0, -1):
            polynomial[:, degree] += packet[:, position]*polynomial[:, degree-1]
    result = polynomial.mean(axis=0)
    if not np.isfinite(result).all() or np.any(result <= 0):
        raise ArithmeticError('positive finite zero-sector coefficients required')
    return result


def region_logs(polynomial, max_q, *, groups=512):
    width = len(polynomial)-1
    if groups % width:
        raise ValueError('whole physical steps per region required')
    local = np.log(np.asarray(polynomial, dtype=float))
    result = np.array([0.])
    for epoch in range(groups//width):
        count = min(max_q, width*(epoch+1))+1
        terms = np.full((width+1, count), -np.inf)
        for k in range(min(width, count-1)+1):
            length = min(len(result), count-k)
            terms[k, k:k+length] = result[:length]+local[k]
        result = logsumexp(terms, axis=0)
    return result-np.asarray([log(comb(groups, q)) for q in range(len(result))])


def screen(points=tuple(range(8)), *, occupancies=(119,), tilts=(.2, .3, .4, .5, .6, .8, 1.)):
    weights = character_weights(points)
    beta_log = 256*log(2)-8*log(65535)
    result = []
    for tilt in tilts:
        polynomial = marked_zero_polynomial(weights, tilt)
        regional = region_logs(polynomial, max(occupancies))
        margins = {str(q): -(log(comb(512, q))+q*beta_log+tilt*13107+32*regional[q])/log(2)
                   for q in occupancies}
        result.append(dict(tilt=tilt, zero_sector_margin_bits=margins,
                           regional_log_moments={str(q): float(regional[q]) for q in occupancies},
                           marked_local_coefficients=list(map(float, polynomial))))
    return dict(proposal_only=True, whole_code_certificate=False,
                interpretation='one positive all-zero-state summand of the outer comparison envelope',
                K=65536, N=131072, groups=512, regions=32, packet_bits=8,
                physical_step_bits=8*len(points), state_bits=16, points=list(points),
                occupancies=list(occupancies), trials=result)


def envelope_ablations(*, q=119, tilt=.4, name='byte_native'):
    """Selected-point operator diagnostics, explicitly not certificates."""
    import maps
    import screen as proof
    data, _ = maps.prepare(name)
    weighted, emission, _ = proof.moments(data, np.exp(-tilt))
    local = proof.operators(weighted, emission)
    variations = {'full': local}
    short = local.copy()
    short[:, 1, :] = 0
    short[:, :, 1] = 0
    variations['restrict_no_uniform_state_paths'] = short
    no_returns = local.copy()
    no_returns[:, 1:, 0] = 0
    variations['INVALID_suppress_returns'] = no_returns
    short_small_returns = short.copy()
    short_small_returns[3:, 1:, 0] = 0
    variations['restrict_short_returns_only_j1_j2'] = short_small_returns
    beta_log = 256*log(2)-8*log(65535)
    values = {}
    for label, family in variations.items():
        regional, _ = proof.placement(family, q, epochs=64, windows=8)
        moment = proof.logarithmic.log_power_matrix(proof.uniform_mixture(regional, q), 32)
        values[label] = -(log(comb(512, q))+q*beta_log+tilt*13107+moment)/log(2)
    return dict(proposal_only=True, q=q, tilt=tilt, map=name,
                interpretation='restricted path contributions except explicitly INVALID optimistic ablation',
                margin_bits=values)


def closure_slack_bound(tilt, *, physical_steps=2048, packet_bits=8):
    """Analytic comparison factor, evaluated in float for a diagnostic.

    One injective active packet bounds each tilted syndrome atom by eta.
    Two packets of full joint rank improve the zero-return correction to
    eta squared. No asymptotic claim or numerical certificate is produced.
    """
    from math import exp, log2
    normalizer = (1+exp(-tilt))**packet_bits-1
    eta = 1/normalizer
    if eta >= 1:
        raise ValueError('this simple sandwich is nontrivial only when eta < 1')
    return dict(tilt=tilt, physical_steps=physical_steps,
                tilted_syndrome_atom_upper=eta,
                whole_closure_slack_bits=-physical_steps*log2(1-eta),
                exact_return_subtraction_slack_bits=-physical_steps*log2(1-eta*eta))


if __name__ == '__main__':
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', type=int, default=list(range(8)))
    parser.add_argument('--q', nargs='+', type=int, default=[113, 119])
    parser.add_argument('--tilts', nargs='+', type=float, default=[.1, .2, .3, .4, .5, .6, .8, 1.])
    args = parser.parse_args()
    print(json.dumps(screen(tuple(args.points), occupancies=tuple(args.q), tilts=tuple(args.tilts)), indent=2))
