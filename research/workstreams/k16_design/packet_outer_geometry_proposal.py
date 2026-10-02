"""Floating, noncertificate comparisons of RS outer block geometries.

The caller supplies a complete macro operator family and an independently
justified UniformInputEnvelope.  State persists across all routing regions.
Nothing here transfers a certificate, resets the state, or takes a scalar
power of a zero-start regional moment.  Selected witnesses need Arb replay.
"""
from fractions import Fraction as Q
from math import comb, isfinite, log

import numpy as np
from packet_rs_length_proposal import scaled_placement, regional_uniform
from packet_rs_s20_proposal import log_power
from rs_uniform_envelope import UniformInputEnvelope


def upper_arrays(family):
    size = family[0].nrows()
    arrays = np.empty((len(family), size, size))
    for k, matrix in enumerate(family):
        if matrix.nrows() != size or matrix.ncols() != size:
            raise ValueError('consistent square operators required')
        for i in range(size):
            for j in range(size):
                endpoint = matrix[i, j].upper()
                value = float(endpoint)
                if not endpoint.is_finite() or endpoint < 0 or not isfinite(value):
                    raise FloatingPointError('nonfinite or negative local endpoint')
                if endpoint > 0 and value == 0:
                    raise FloatingPointError('local conversion lost positive mass')
                arrays[k, i, j] = value
    return arrays


def estimate(local, *, K, envelope, occupancies, tilt, windows=32, distance=Q(1, 10)):
    """Evaluate one rational tilt for an explicitly chosen outer geometry."""
    if not isinstance(envelope, UniformInputEnvelope) or envelope.packet_bits != 4:
        raise ValueError('an explicit four-bit-packet uniform envelope is required')
    if type(K) is not int or K <= 0 or K % envelope.message_bits:
        raise ValueError('K must be a positive multiple of the outer dimension')
    L, R = K // envelope.message_bits, envelope.regions
    if type(windows) is not int or windows < 1 or L % windows:
        raise ValueError('each region must contain complete macros')
    qs = tuple(occupancies)
    if not qs or tuple(sorted(set(qs))) != qs or any(type(q) is not int or not 1 <= q <= L for q in qs):
        raise ValueError('ordered distinct feasible occupancies required')
    tilt, distance = Q(tilt), Q(distance)
    if tilt <= 0 or not 0 < distance < 1:
        raise ValueError('positive tilt and relative distance in (0,1) required')
    N = L * envelope.output_bits
    cutoff = (N * distance.numerator) // distance.denominator
    regional = scaled_placement(local, qs[-1], epochs=L // windows, windows=windows)
    beta = envelope.beta
    log_beta = log(beta.numerator) - log(beta.denominator)
    witnesses = {}
    for q in qs:
        mixed, scale = regional_uniform(regional, q)
        with np.errstate(over='raise', invalid='raise', divide='raise', under='raise'):
            moment = log_power(mixed, R) + R * scale
        value = log(comb(L, q)) + q * log_beta + float(tilt) * cutoff + moment
        if not isfinite(value):
            raise FloatingPointError('nonfinite proposal objective')
        witnesses[str(q)] = dict(estimated_margin_bits=-value / log(2), log_moment=moment)
    return dict(schema='rs-outer-geometry-proposal-1', proposal_only=True,
        whole_code_certificate=False, K=K, N=N, groups=L, regions=R,
        outer_dimension=envelope.message_bits, outer_length=envelope.output_bits,
        outer_symbols=[envelope.n, envelope.k], symbol_bits=envelope.symbol_bits,
        windows=windows, epochs_per_region=L // windows, cutoff=cutoff,
        distance=str(distance), tilt=str(tilt), zero_initial_state=True,
        terminal='sum of all coordinates', continuous_state_across_regions=True,
        witnesses=witnesses,
        scope='Floating proposal only; no certificate or transfer from another geometry.')
