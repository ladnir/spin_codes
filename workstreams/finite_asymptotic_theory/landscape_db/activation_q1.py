"""Binary64 Q1 screen with zero/arbitrary-live/density-bounded-live classes.

The transfer is the activation-aware envelope documented in the BCH bridge's
README (ba80 worktree), generalized to arbitrary complete region lengths.
This module has no dependency on that worktree. It is not an outward verifier.
"""
from __future__ import annotations

import math
import numpy as np


def epoch_logs(t: int, s: int, spectrum: dict[int, int], lambdas: np.ndarray):
    """Return a batch of row-to-column transfers, indexed by positive tilt."""
    if s < 2 or sum(spectrum.values()) != (1 << s) - 1:
        raise ValueError("expected the complete nonzero A spectrum, with s >= 2")
    if any(w <= 0 or w > t or c <= 0 for w, c in spectrum.items()):
        raise ValueError("invalid nonzero A spectrum")
    lam = np.asarray(lambdas, dtype=float)
    if lam.ndim != 1 or not np.all(np.isfinite(lam) & (lam > 0)):
        raise ValueError("positive finite tilt vector required")
    den = (1 << s) - 1
    log_den = math.log(den)
    kappa = math.log1p(1 / (den - 1))
    m0 = np.full(lam.shape, -np.inf)
    m1 = np.full(lam.shape, -np.inf)
    for w, count in spectrum.items():
        mass = math.log(count) - log_den
        np.logaddexp(m0, mass - lam * w, out=m0)
        np.logaddexp(m1, mass + math.log(w / t) - lam * (w - 1), out=m1)
        if w < t:
            np.logaddexp(m1, mass + math.log((t - w) / t) - lam * (w + 1), out=m1)
    d = min(spectrum)
    zero = np.full((len(lam), 3, 3), -np.inf)
    one = np.full_like(zero, -np.inf)
    zero[:, 0, 0] = 0
    zero[:, 1, 2] = -lam * d
    zero[:, 2, 2] = kappa + m0
    one[:, 0, 1] = -lam
    one[:, 1, 0] = -lam * (d - 1) - log_den
    one[:, 1, 2] = -lam * (d - 1)
    one[:, 2, 0] = kappa + m1 - log_den
    one[:, 2, 2] = kappa + m1
    return zero, one


def matrix_product(a, b):
    # Fixed three-state contraction, batched across tilts.
    return np.logaddexp(
        np.logaddexp(a[:, :, 0, None] + b[:, None, 0, :],
                     a[:, :, 1, None] + b[:, None, 1, :]),
        a[:, :, 2, None] + b[:, None, 2, :])


def region_logs(zero, one, epochs: int):
    """Degree-one matrix-polynomial powering in O(log epochs)."""
    if epochs < 1:
        raise ValueError("a region must contain at least one complete epoch")
    rz = np.full_like(zero, -np.inf)
    rz[:, 0, 0] = rz[:, 1, 1] = rz[:, 2, 2] = 0
    ra = np.full_like(one, -np.inf)
    bz, ba = zero, one
    remaining = epochs
    while remaining:
        if remaining & 1:
            ra = np.logaddexp(matrix_product(ra, bz), matrix_product(rz, ba))
            rz = matrix_product(rz, bz)
        remaining >>= 1
        if remaining:
            ba = np.logaddexp(matrix_product(ba, bz), matrix_product(bz, ba))
            bz = matrix_product(bz, bz)
    return rz, ra - math.log(epochs)


def vector_product(v, m, out):
    # v is [tilt, weight, state]. Keep the state contraction explicit.
    np.logaddexp(v[:, :, 0, None] + m[:, None, 0, :],
                 v[:, :, 1, None] + m[:, None, 1, :], out=out)
    np.logaddexp(out, v[:, :, 2, None] + m[:, None, 2, :], out=out)


def coefficient_logs(rz, ra, length: int):
    """Uniform outer-weight support average; state starts at zero."""
    if length < 1:
        raise ValueError("positive outer length required")
    shape = (len(rz), length + 1, 3)
    current = np.full(shape, -np.inf)
    updated = np.empty(shape)
    marked = np.empty(shape)
    current[:, 0, 0] = 0
    for n in range(length):
        vector_product(current[:, :n + 1], rz, updated[:, :n + 1])
        vector_product(current[:, :n + 1], ra, marked[:, :n + 1])
        updated[:, n + 1] = -np.inf
        np.logaddexp(updated[:, 1:n + 2], marked[:, :n + 1], out=updated[:, 1:n + 2])
        current, updated = updated, current
    normalizers = np.array([math.log(math.comb(length, w)) for w in range(length + 1)])
    return np.logaddexp.reduce(current, axis=2) - normalizers


def screen(t, s, a_spectrum, length, dimension, message_exponent, log_counts, tilts):
    message_bits = 1 << message_exponent
    if message_bits % dimension or (message_bits // dimension) % t:
        raise ValueError("non-native length: dimension | k and t | k/dimension required")
    rows = message_bits // dimension
    cutoff = (length * rows) // 10
    lam = np.exp(np.asarray(tilts))
    moments = coefficient_logs(*region_logs(*epoch_logs(t, s, a_spectrum, lam), rows // t), length)
    values = np.minimum(0., moments + cutoff * lam[:, None])
    witnesses = np.argmin(values, axis=0)
    best = values[witnesses, np.arange(length + 1)]
    terms = np.array([math.log(rows) + count + best[w] for w, count in sorted(log_counts.items())])
    weights = sorted(log_counts)
    dominant = weights[int(np.argmax(terms))]
    margin = -float(np.logaddexp.reduce(terms)) / math.log(2)
    if not math.isfinite(margin):
        raise ArithmeticError("nonfinite margin")
    return dict(margin_bits=margin, dominant_weight=dominant,
                dominant_log_surprisal=float(tilts[witnesses[dominant]]),
                dominant_witness_at_grid_edge=bool(witnesses[dominant] in (0, len(tilts) - 1)),
                message_bits=message_bits, output_bits=length * rows, outer_rows=rows,
                epochs_per_region=rows // t, bad_weight=cutoff)
