"""Floating proposals for RS16/GL20 output tilts and marker probabilities.

The input is freshly prepared t64/s20 data, not a saved numerical endpoint.
For each positive output tilt, cache the 33 local occupancy operators as
floating matrices. For marker probability p, their Binomial(32,15*p/16)
mixture is the iid macro operator B_p. The objective includes the regional
conditioning loss for exactly q of L markers in each of R regions:

    log C(L,q) + q log beta + lambda floor(N/10)
      - R log[C(L,q) p^q (1-p)^(L-q)]
      + log(e_zero B_p^(N/128) 1).

The matrix power is rescaled during every multiplication. We retain its
zero starting state and all terminal mass instead of substituting an
asymptotic eigenvalue. Optimization in logit(p) avoids coarse marker grids.
Optional midpoint refinement searches between supplied output tilts.

All scores are floating proposals, not upper bounds or certificates. The
returned rational tilt and p lists can be passed to the existing outward
fugacity evaluator. That evaluator must recompute every numerical endpoint.
This module neither prepares a full state census nor writes result files.
"""
from __future__ import annotations

from dataclasses import asdict
from fractions import Fraction as Q
from math import comb, exp, isfinite, log, log1p, sqrt

import numpy as np
from scipy.optimize import minimize_scalar
from flint import ctx
import packet_q1 as q1
import packet_rs_state_sparse as sparse
from rs_uniform_envelope import UniformInputEnvelope


def log_power(matrix, exponent):
    """Floating log(e_zero M^exponent 1), rescaled before the first square."""
    matrix = np.asarray(matrix, dtype=float)
    if (matrix.ndim != 2 or not len(matrix) or matrix.shape[0] != matrix.shape[1]
            or not np.all(np.isfinite(matrix)) or np.any(matrix < 0)
            or type(exponent) is not int or exponent < 0):
        raise ValueError('finite nonnegative square matrix and nonnegative integer exponent required')
    if exponent == 0:
        return 0.0
    scale = float(matrix.max())
    if scale <= 0:
        raise FloatingPointError('zero matrix has no positive finite log moment')
    power, power_log = matrix / scale, log(scale)
    row = np.zeros(len(matrix)); row[0] = 1.0
    row_log = 0.0
    while exponent:
        if exponent & 1:
            row = row @ power
            scale = float(row.max())
            if not isfinite(scale) or scale <= 0:
                raise FloatingPointError('zero or nonfinite floating path mass; replay with Arb')
            row /= scale
            row_log += power_log + log(scale)
        exponent >>= 1
        if exponent:
            power = power @ power
            scale = float(power.max())
            if not isfinite(scale) or scale <= 0:
                raise FloatingPointError('zero or nonfinite floating matrix power; replay with Arb')
            power /= scale
            power_log = 2 * power_log + log(scale)
    result = row_log + log(float(row.sum()))
    if not isfinite(result):
        raise FloatingPointError('nonfinite floating log moment')
    return result


def iid_macro(local, marker_probability):
    """Floating mixture of the complete 0..32 nonzero-packet macro family."""
    local = np.asarray(local, dtype=float)
    p = float(marker_probability)
    if (local.ndim != 3 or len(local) != 33 or local.shape[1] != local.shape[2]
            or not local.shape[1] or not np.all(np.isfinite(local)) or np.any(local < 0)
            or not isfinite(p) or not 0 < p <= 1):
        raise ValueError('complete finite macro family and marker probability in (0,1] required')
    theta = 15 * p / 16
    weights = np.array([comb(32, j) * theta**j * (1-theta)**(32-j) for j in range(33)])
    if abs(float(weights.sum()) - 1) > 1e-12:
        raise FloatingPointError('binomial mixture failed floating normalization')
    return np.tensordot(weights, local, axes=1)


def _probability(value):
    """Round a proposal to a reproducible rational, preserving open endpoints."""
    if not isfinite(value) or not 0 < value <= 1:
        raise ValueError('positive finite probability at most one required')
    result = Q(format(value, '.12g'))
    if not 0 < result <= 1:
        raise ValueError('rounded marker probability left its domain')
    return result


def _sigmoid(value):
    if value >= 0:
        return 1 / (1 + exp(-value))
    e = exp(value)
    return e / (1 + e)


class ProposalModel:
    """Reuse one prepared map and cached local families across occupancies."""
    def __init__(self, data, map_record, *, geometry=sparse.DEFAULT_GEOMETRY, precision=192):
        if sparse.validated(data, map_record) != 20:
            raise ValueError('this proposal model requires the declared twenty-bit inner')
        if (not isinstance(geometry, q1.Geometry) or geometry.regions != 64
                or geometry.group_dimension != 128 or type(precision) is not int or precision < 128):
            raise ValueError('RS16 packet geometry and integer precision >=128 required')
        self.data, self.map_record, self.geometry = data, map_record, geometry
        self.precision = precision
        beta = UniformInputEnvelope(16, 8, 4, 4).beta
        self.log_beta = log(beta.numerator) - log(beta.denominator)
        self.sources = sparse.source_snapshot(map_record)
        self.locals = {}
        self.moments = {}
        self.log_counts = {}

    def local(self, tilt):
        tilt = Q(tilt)
        if tilt <= 0:
            raise ValueError('positive output tilt required')
        if tilt not in self.locals:
            prior_precision = ctx.prec
            try:
                ctx.prec = self.precision
                family = q1.kernel_t64.local_operators(self.data, tilt, activity=Q(1, 2))
                size = family[0].nrows()
                if len(family) != 33 or any(m.nrows() != size or m.ncols() != size for m in family):
                    raise ValueError('complete square macro-occupancy family required')
                arrays = np.array([[[float(m[i, j]) for j in range(size)] for i in range(size)]
                                   for m in family])
                if (not np.all(np.isfinite(arrays)) or np.any(arrays < 0)
                        or any(arrays[k, i, j] == 0 and family[k][i, j] > 0
                               for k in range(33) for i in range(size) for j in range(size))):
                    raise FloatingPointError('local float conversion lost finite positive coefficients')
                self.locals[tilt] = arrays
            finally:
                ctx.prec = prior_precision
        return self.locals[tilt]

    def score(self, occupancy, tilt, marker_probability):
        """Return a floating log first-moment proposal, never an outward bound."""
        groups = self.geometry.group_count
        if type(occupancy) is not int or not 1 <= occupancy <= groups:
            raise ValueError('occupancy must lie in 1..L')
        tilt, probability = Q(tilt), Q(marker_probability)
        p = float(probability)
        if tilt <= 0 or not 0 < probability <= 1 or (probability == 1 and occupancy != groups):
            raise ValueError('positive tilt and marker probability admitting exactly q markers required')
        if not 0 < p <= 1 or p == 1 and probability != 1:
            raise FloatingPointError('marker probability cannot be represented faithfully as float')
        key = tilt, probability
        if key not in self.moments:
            mixed = iid_macro(self.local(tilt), p)
            self.moments[key] = log_power(mixed, self.geometry.N // 128)
        if occupancy not in self.log_counts:
            self.log_counts[occupancy] = log(comb(groups, occupancy))
        lc = self.log_counts[occupancy]
        log_condition = 0.0 if probability == 1 else (
            lc + occupancy * log(p) + (groups-occupancy) * log1p(-p))
        result = (lc + occupancy * self.log_beta + float(tilt) * (self.geometry.N // 10)
                  - self.geometry.regions * log_condition + self.moments[key])
        if not isfinite(result):
            raise FloatingPointError('nonfinite objective; no witness proposed')
        return result

    def _marker_fit(self, occupancy, tilt, markers, max_iterations):
        groups = self.geometry.group_count
        # The logit grid spans sparse and nearly complete marker populations.
        # Include the unconstrained binomial mode q/L and every caller proposal.
        candidates = {_probability(_sigmoid(float(x))) for x in np.linspace(-18, 18, 25)}
        candidates.update(markers)
        candidates.add(Q(occupancy, groups))
        candidates = sorted(p for p in candidates if p < 1 or occupancy == groups)
        ranked = [(self.score(occupancy, tilt, p), p) for p in candidates]
        best = min(ranked)
        # Refine both strongest brackets; row-choice envelopes need not be smooth.
        interior = [p for p in candidates if p < 1]
        ranks = sorted(range(len(interior)), key=lambda i: self.score(occupancy, tilt, interior[i]))[:2]
        for index in ranks:
            lo = max(0, index-1); hi = min(len(interior)-1, index+1)
            if lo == hi:
                continue
            bounds = tuple(log(float(interior[i])) - log1p(-float(interior[i])) for i in (lo, hi))
            def objective(logit):
                p = _probability(_sigmoid(logit))
                return self.score(occupancy, tilt, p)
            fit = minimize_scalar(objective, bounds=bounds, method='bounded',
                                  options=dict(xatol=1e-7, maxiter=max_iterations))
            # Re-evaluate the rounded rational; never report the unrounded fit.fun.
            p = _probability(_sigmoid(float(fit.x)))
            candidate = self.score(occupancy, tilt, p), p
            best = min(best, candidate)
        return best

    def propose(self, occupancies, tilts, *, marker_probabilities=(), tilt_refinements=0,
                max_marker_iterations=80):
        """Optimize markers per tilt, optionally bisecting neighboring tilt intervals.

        Each refinement round adds at most two local-operator families per q.
        The default only searches the supplied tilt grid. There is no full census.
        """
        qs = tuple(occupancies)
        tilts = tuple(map(Q, tilts))
        markers = tuple(map(Q, marker_probabilities))
        if (not qs or any(type(q) is not int or not 1 <= q <= self.geometry.group_count for q in qs)
                or tuple(sorted(set(qs))) != qs or not tilts or min(tilts) <= 0
                or len(set(tilts)) != len(tilts) or any(not 0 < p <= 1 for p in markers)
                or type(tilt_refinements) is not int or not 0 <= tilt_refinements <= 8
                or type(max_marker_iterations) is not int or max_marker_iterations < 1):
            raise ValueError('explicit ordered occupancies, distinct positive tilts, and bounded search settings required')
        witnesses = {}
        for q in qs:
            fitted = {t: self._marker_fit(q, t, markers, max_marker_iterations) for t in sorted(tilts)}
            for _ in range(tilt_refinements):
                grid = sorted(fitted)
                winner = min(grid, key=lambda t: fitted[t][0])
                index = grid.index(winner)
                neighbors = [grid[j] for j in (index-1, index+1) if 0 <= j < len(grid)]
                for neighbor in neighbors:
                    middle = Q(format(sqrt(float(winner) * float(neighbor)), '.12g'))
                    if middle > 0 and middle not in fitted:
                        fitted[middle] = self._marker_fit(q, middle, markers, max_marker_iterations)
            tilt = min(fitted, key=lambda t: fitted[t][0])
            value, probability = fitted[tilt]
            witnesses[str(q)] = dict(tilt=str(tilt), marker_probability=str(probability),
                estimated_margin_bits=-value/log(2),
                tilt_at_search_boundary=tilt in (min(tilts), max(tilts)),
                requires_outward_replay=True)
        if self.sources != sparse.source_snapshot(self.map_record):
            raise RuntimeError('loaded source or map base changed during proposal search')
        used_tilts = sorted({Q(w['tilt']) for w in witnesses.values()})
        used_markers = sorted({Q(w['marker_probability']) for w in witnesses.values()})
        return dict(schema='rs16-s20-fugacity-proposals-1', proposal_only=True,
            whole_code_certificate=False, has_numerical_upper_endpoints=False,
            occupancy_values=list(qs), geometry=asdict(self.geometry), state_bits=20,
            map_sha256=self.data['map_sha256'], source_sha256=self.sources,
            precision_for_local_envelopes=self.precision, witnesses=witnesses,
            tilts=list(map(str, used_tilts)), marker_probabilities=list(map(str, used_markers)),
            local_tilts_evaluated=len(self.locals), cached_moments=len(self.moments),
            replay='Pass returned tilts and marker_probabilities to the existing outward '
                   'run_tail(method="fugacity") for the requested occupancies; recompute all endpoints.',
            scope='Floating search proposals only. Scores may be optimistic; no distance '
                  'or failure-probability claim follows before fresh outward replay.')
