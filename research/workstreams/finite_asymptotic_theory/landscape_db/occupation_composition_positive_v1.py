"""Positive composition contractions with a log-domain fallback.

The envelope is a maximum of positive linear forms. Scaling each matrix
entry commutes with this operation. Rare coefficients or lost support send
that entry through the original log evaluator; no FFT is used.
"""
import math

import numpy as np

import occupation_composition_lazy_v1 as lazy
import occupation_refresh_v1 as transfer


def shifted_mixture(regions, distribution, remaining):
    count = len(distribution)
    if len(regions) < count+remaining:
        raise ValueError('missing region coefficients')
    shape = regions.shape[1:]
    flat = regions[:count+remaining].reshape(count+remaining, -1)
    output = np.full((remaining+1, flat.shape[1]), -np.inf)
    dscale = float(np.max(distribution))
    if not math.isfinite(dscale):
        return output.reshape((remaining+1, *shape))
    weights = np.exp(distribution-dscale)
    for entry in range(flat.shape[1]):
        column = flat[:, entry]
        scale = float(np.max(column))
        if not math.isfinite(scale):
            continue
        raw = np.correlate(np.exp(column-scale), weights, mode='valid')
        with np.errstate(divide='ignore'):
            output[:, entry] = np.log(raw)+scale+dscale
        # An underflowed summand cannot be material relative to an ordinary
        # output above this threshold. Recompute all tiny outputs directly.
        for j in np.flatnonzero(raw < 1e-130):
            output[j, entry] = np.logaddexp.reduce(column[j:j+count]+distribution)
    return output.reshape((remaining+1, *shape))


def log_fold(logs, envelope):
    current = logs.copy()
    while len(current) > 1:
        current = envelope.apply(current[:-1], current[1:])
    return current[0]


def envelope_fold(logs, envelope):
    if len(logs) == 1:
        return logs[0].copy()
    max_a, max_b = float(np.max(envelope.log_a)), float(np.max(envelope.log_b))
    # Over eight steps these guards bound both growth and the decay of a
    # surviving maximum, preventing normalization from magnifying underflow.
    if max(max_a, max_b) > 10 or min(max_a, max_b) < -10:
        return log_fold(logs, envelope)
    shape = logs.shape[1:]
    original = logs.reshape(len(logs), -1)
    maxima = np.max(original, axis=0)
    scales = np.where(np.isfinite(maxima), maxima, 0.)
    centered = original-scales
    support = np.isfinite(original)
    unsafe = np.any(support & (centered < -250), axis=0)
    if np.all(unsafe):
        return log_fold(logs, envelope)
    current = np.exp(centered)
    output = np.empty_like(current)
    left_term = np.empty_like(current)
    right_term = np.empty_like(current)
    next_support = np.empty_like(support)
    a, b = np.exp(envelope.log_a), np.exp(envelope.log_b)
    # The support recurrence below requires at least one positive coefficient
    # on each side; endpoint-only envelopes use the original evaluator.
    if not np.any(a > 0) or not np.any(b > 0):
        return log_fold(logs, envelope)
    size = len(current)
    steps = 0
    while size > 1:
        size -= 1; steps += 1
        left, right = current[:size], current[1:size+1]
        target, tmp, tmp2 = output[:size], left_term[:size], right_term[:size]
        np.multiply(left, a[0], out=target)
        np.multiply(right, b[0], out=tmp)
        np.add(target, tmp, out=target)
        for index in range(1, len(a)):
            np.multiply(left, a[index], out=tmp)
            np.multiply(right, b[index], out=tmp2)
            np.add(tmp, tmp2, out=tmp)
            np.maximum(target, tmp, out=target)
        np.logical_or(support[:size], support[1:size+1], out=next_support[:size])
        if steps % 8 == 0 or size == 1:
            scale = np.max(target, axis=0)
            scale = np.where(scale > 0, scale, 1.)
            target /= scale
            scales += np.log(scale)
            unsafe |= np.any(next_support[:size] & (target < math.exp(-250)), axis=0)
        current, output = output, current
        support, next_support = next_support, support
    with np.errstate(divide='ignore'):
        result = np.log(current[0])+scales
    if np.any(unsafe):
        result[unsafe] = log_fold(original[:, unsafe], envelope)
    return result.reshape(shape)


class CompositionBoxes(lazy.CompositionBoxes):
    def matrix(self, regions, lower, total):
        lower = np.asarray(lower, dtype=np.int64)
        if (len(lower) != len(self.bands) or np.any(lower < 0)
                or lower.sum() > total or total > self.maximum):
            raise ValueError('invalid lower type counts')
        remaining = total-int(lower.sum()); key = tuple(map(int, lower))
        distribution = self.distribution_cache.get(key)
        if distribution is None:
            distribution = np.array([0.])
            for group, count in enumerate(lower):
                distribution = transfer.boxes.positive_log_convolve(
                    distribution, self.binomials[group][int(count)])
            if len(self.distribution_cache) >= 2048:
                self.distribution_cache.clear()
            self.distribution_cache[key] = distribution
        current = shifted_mixture(regions, distribution, remaining)
        current += float(lower@self.log_gamma)/self.block
        return envelope_fold(current, self.envelope)
