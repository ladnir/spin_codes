"""Floating tilt proposals for RS16 packet codes at several message lengths.

The caller supplies freshly prepared t64/GL_s maps, with s at most 22.
For K divisible by 4096, L=K/128 groups occupy 64 routing regions.
Each region contains L/32 ordered macros, with 32 packet slots per macro.

At positive tilt lambda, local Arb operators describe each macro occupancy.
We convert their upper endpoints to floats, then average ordered products
over exact without-replacement placements. This is not an iid approximation.
For q active group labels, the uniform input majorant has J~Bin(q,15/16)
nonzero packets per region. The floating objective is

    log C(L,q) + q log beta + lambda floor(2K/10)
      + log(e_zero (E[R_J(lambda)])^64 1),

where beta=2^256/(2^16-1)^8. State persists across regions; all terminal
coordinates count mass. The local row-selection activity is 1/2, not the
packet activity 15/16. Per-occupancy scaling protects long regional products.
Binomial mixing takes place in the log domain, including rare occupancies.

Scores are proposals, not bounds, despite using local Arb upper endpoints.
Every selected rational tilt requires fresh outward replay with the same
maps, length, occupancy, and ideal independent-uniform setup distribution.
This module does not prepare a state census, write receipts, or certify a code.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from fractions import Fraction as Q
from math import comb, isfinite, log

import numpy as np
from scipy.special import gammaln, logsumexp
from flint import ctx
import packet_q1 as q1
import packet_rs_state_sparse as sparse
from packet_rs_s20_proposal import log_power
from rs_uniform_envelope import UniformInputEnvelope


def geometry(K):
    """Return the fixed RS16 geometry without preparing any maps."""
    if type(K) is not int or K <= 0 or K % 4096:
        raise ValueError('K must be a positive integer multiple of 4096')
    return q1.Geometry(K // 128, 64, 128)


def _family(operators):
    operators = np.asarray(operators, dtype=float)
    if (operators.ndim != 3 or not len(operators) or not operators.shape[1]
            or operators.shape[1] != operators.shape[2]
            or not np.isfinite(operators).all() or np.any(operators < 0)):
        raise ValueError('finite nonnegative square operator family required')
    if np.any(operators.max(axis=(1, 2)) <= 0):
        raise FloatingPointError('zero local operator has no positive scaling factor')
    return operators


@dataclass(frozen=True)
class ScaledPlacement:
    """R_j = exp(log_scales[j])*matrices[j], in floating arithmetic only."""
    matrices: np.ndarray
    log_scales: np.ndarray


def scaled_placement(operators, degree, *, epochs, windows=32):
    """Conditional ordered placement, truncated to degree, with rowwise scaling.

    Appending a macro to e-1 macros gives weight
    C(w,k)*C((e-1)w,j-k)/C(ew,j) to R_(j-k)*T_k.
    This is the normalized hypergeometric recurrence also used by
    occupancy_sensitivity.float_placement. Scaling is separate for each j.
    The recurrence keeps the multiplication order, including noncommuting T_k.
    """
    operators = _family(operators)
    if (type(degree) is not int or degree < 0 or type(epochs) is not int or epochs < 1
            or type(windows) is not int or windows < 1 or len(operators) != windows + 1
            or degree > epochs * windows):
        raise ValueError('complete local family and feasible integer placement geometry required')
    size = operators.shape[1]
    local_scale = operators.max(axis=(1, 2))
    with np.errstate(over='raise', invalid='raise', divide='raise', under='raise'):
        local = operators / local_scale[:, None, None]
        local_logs = np.log(local_scale)
        current = np.eye(size)[None, :, :]
        current_logs = np.zeros(1)
        for epoch in range(1, epochs + 1):
            limit = min(degree, epoch * windows)
            total, previous = epoch * windows, (epoch - 1) * windows
            terms = []
            base = np.full(limit + 1, -np.inf)
            for k in range(min(windows, limit) + 1):
                js = np.arange(k, min(limit, k + len(current) - 1) + 1)
                old = js - k
                log_weight = (log(comb(windows, k)) + gammaln(previous + 1)
                    - gammaln(old + 1) - gammaln(previous - old + 1)
                    - gammaln(total + 1) + gammaln(js + 1) + gammaln(total - js + 1))
                exponent = log_weight + current_logs[old] + local_logs[k]
                base[js] = np.maximum(base[js], exponent)
                terms.append((k, js, old, exponent))
            if not np.isfinite(base).all():
                raise FloatingPointError('missing finite positive placement term')
            nxt = np.zeros((limit + 1, size, size))
            for k, js, old, exponent in terms:
                weights = np.exp(exponent - base[js])
                nxt[js] += weights[:, None, None] * (current[old] @ local[k])
            scale = nxt.max(axis=(1, 2))
            if not np.isfinite(scale).all() or np.any(scale <= 0):
                raise FloatingPointError('zero or nonfinite regional placement; use Arb replay')
            current = nxt / scale[:, None, None]
            current_logs = base + np.log(scale)
        if not np.isfinite(current_logs).all() or not np.isfinite(current).all():
            raise FloatingPointError('nonfinite scaled regional family')
    current.setflags(write=False)
    current_logs.setflags(write=False)
    return ScaledPlacement(current, current_logs)


def regional_uniform(regional, occupancy):
    """Return normalized E[R_J] and its log scale, J~Bin(q,15/16).

    Log-domain summation retains rare J even when their probabilities alone
    underflow. Positive entries that cannot survive final float normalization
    cause an error; they are not silently discarded to improve a score.
    """
    if not isinstance(regional, ScaledPlacement):
        raise ValueError('scaled placement family required')
    matrices, scales = regional.matrices, regional.log_scales
    if (type(occupancy) is not int or not 0 <= occupancy < len(matrices)
            or matrices.ndim != 3 or not matrices.shape[1]
            or matrices.shape[1] != matrices.shape[2]
            or scales.shape != (len(matrices),) or not np.isfinite(scales).all()
            or not np.isfinite(matrices).all() or np.any(matrices < 0)):
        raise ValueError('finite scaled family and covered occupancy required')
    q = occupancy
    js = np.arange(q + 1)
    log_weights = (gammaln(q + 1) - gammaln(js + 1) - gammaln(q - js + 1)
                   + js * log(15) - q * log(16))
    normalization = float(logsumexp(log_weights))
    if not isfinite(normalization) or abs(normalization) > 1e-8:
        raise FloatingPointError('binomial weights failed floating normalization')
    log_weights -= normalization
    # The maximum is chosen per matrix entry, not just at the binomial mode.
    # A rare occupancy can dominate an entry of this tilted expectation.
    logs = np.full(matrices[:q + 1].shape, -np.inf)
    positive = matrices[:q + 1] > 0
    logs[positive] = np.log(matrices[:q + 1][positive])
    terms = logs + (scales[:q + 1] + log_weights)[:, None, None]
    entries = logsumexp(terms, axis=0)
    finite = np.isfinite(entries)
    if not np.any(finite) or np.any(np.isnan(entries)) or np.any(np.isposinf(entries)):
        raise FloatingPointError('zero or nonfinite binomial matrix mixture')
    scale = float(entries[finite].max())
    result = np.zeros(entries.shape)
    if np.any(entries[finite] - scale < log(np.finfo(float).tiny)):
        raise FloatingPointError('positive mixed entry would become subnormal or underflow')
    with np.errstate(over='raise', invalid='raise', divide='raise', under='raise'):
        result[finite] = np.exp(entries[finite] - scale)
    return result, scale


def _tilts(values):
    if isinstance(values, (str, bytes)):
        raise ValueError('sequence of distinct positive rational tilts required')
    try:
        values = tuple(Q(value) for value in values)
    except (ValueError, TypeError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('finite rational tilts required') from error
    if not values or min(values) <= 0 or len(set(values)) != len(values):
        raise ValueError('distinct positive rational tilts required')
    for value in values:
        if not isfinite(float(value)) or float(value) <= 0:
            raise FloatingPointError('tilt cannot be represented as positive finite float')
    return tuple(sorted(values))


class ProposalModel:
    """Reuse fresh maps and local operators across lengths and q grids."""
    def __init__(self, data, map_record, *, precision=192):
        bits = sparse.validated(data, map_record)
        if not 1 <= bits <= 22:
            raise ValueError('this search supports physical state dimensions 1..22')
        if type(precision) is not int or precision < 192:
            raise ValueError('local Arb precision must be an integer at least 192')
        self.data, self.map_record, self.state_bits = data, map_record, bits
        self.precision = precision
        beta = UniformInputEnvelope(16, 8, 4, 4).beta
        self.log_beta = log(beta.numerator) - log(beta.denominator)
        self.sources = sparse.source_snapshot(map_record)
        self.locals, self.regionals, self.moments, self.log_counts = {}, {}, {}, {}
        self.placement_evaluations = 0

    def local(self, tilt):
        """Convert finite nonnegative local Arb upper endpoints, not midpoints."""
        tilt = _tilts([tilt])[0]
        if tilt not in self.locals:
            previous = ctx.prec
            try:
                ctx.prec = self.precision
                family = q1.kernel_t64.local_operators(self.data, tilt, activity=Q(1, 2))
                if not family or len(family) != 33 or not family[0].nrows():
                    raise ValueError('complete nonempty 0..32 local family required')
                size = family[0].nrows()
                if any(m.nrows() != size or m.ncols() != size for m in family):
                    raise ValueError('consistent square local operators required')
                arrays = np.empty((33, size, size))
                for k, matrix in enumerate(family):
                    for i in range(size):
                        for j in range(size):
                            endpoint = matrix[i, j].upper()
                            if not endpoint.is_finite() or endpoint < 0:
                                raise FloatingPointError('nonfinite or negative local upper endpoint')
                            value = float(endpoint)
                            if not isfinite(value) or (endpoint > 0 and value == 0):
                                raise FloatingPointError('local conversion lost finite positive mass')
                            arrays[k, i, j] = value
                _family(arrays)
                arrays.setflags(write=False)
                self.locals[tilt] = arrays
            finally:
                ctx.prec = previous
        return self.locals[tilt]

    def regional(self, K, tilt, maximum_occupancy):
        """Cache the largest requested degree per length/tilt; reuse its prefix."""
        shape = geometry(K)
        tilt = _tilts([tilt])[0]
        if type(maximum_occupancy) is not int or not 0 <= maximum_occupancy <= shape.group_count:
            raise ValueError('maximum occupancy must lie in 0..L')
        key = K, tilt
        cached = self.regionals.get(key)
        if cached is None or len(cached.matrices) <= maximum_occupancy:
            cached = scaled_placement(self.local(tilt), maximum_occupancy,
                epochs=shape.macros_per_region, windows=shape.macro_windows)
            self.regionals[key] = cached
            self.placement_evaluations += 1
        return cached

    def score(self, K, occupancy, tilt):
        """Return the floating log first-moment objective, never an upper bound."""
        shape = geometry(K)
        if type(occupancy) is not int or not 1 <= occupancy <= shape.group_count:
            raise ValueError('occupancy must lie in 1..L')
        tilt = _tilts([tilt])[0]
        key = K, occupancy, tilt
        if key not in self.moments:
            regional = self.regional(K, tilt, occupancy)
            mixed, scale = regional_uniform(regional, occupancy)
            with np.errstate(over='raise', invalid='raise', divide='raise', under='raise'):
                self.moments[key] = log_power(mixed, shape.regions) + shape.regions * scale
        count_key = K, occupancy
        if count_key not in self.log_counts:
            self.log_counts[count_key] = log(comb(shape.group_count, occupancy))
        result = (self.log_counts[count_key] + occupancy * self.log_beta
                  + float(tilt) * (shape.N // 10) + self.moments[key])
        if not isfinite(result):
            raise FloatingPointError('nonfinite objective; no witness proposed')
        return result

    def propose(self, lengths, occupancies, tilts):
        """Search a rational tilt grid for each requested (K,q).

        Occupancies are either one ordered grid for every K, or a mapping
        from each requested K to its own ordered grid. All q lie in 3..L.
        Local operators are reused across K. For each (K,tilt), a single
        regional recurrence covers every q through the largest requested q.
        """
        lengths = tuple(lengths)
        if not lengths or len(set(lengths)) != len(lengths):
            raise ValueError('nonempty distinct message lengths required')
        shapes = {K: geometry(K) for K in lengths}
        if isinstance(occupancies, Mapping):
            if set(occupancies) != set(lengths):
                raise ValueError('occupancy mapping must cover exactly the requested lengths')
            grids = {K: tuple(occupancies[K]) for K in lengths}
        else:
            shared = tuple(occupancies)
            grids = {K: shared for K in lengths}
        for K, qs in grids.items():
            if (not qs or any(type(q) is not int or not 3 <= q <= shapes[K].group_count for q in qs)
                    or tuple(sorted(set(qs))) != qs):
                raise ValueError('ordered distinct occupancies in 3..L required for every K')
        tilts = _tilts(tilts)
        best = {K: {} for K in lengths}
        for tilt in tilts:
            for K in lengths:
                self.regional(K, tilt, grids[K][-1])
                for q in grids[K]:
                    candidate = self.score(K, q, tilt), tilt
                    if q not in best[K] or candidate < best[K][q]:
                        best[K][q] = candidate
        if self.sources != sparse.source_snapshot(self.map_record):
            raise RuntimeError('loaded source or map source changed during proposal search')
        records = {}
        for K in lengths:
            choices = {str(q): dict(tilt=str(best[K][q][1]),
                estimated_margin_bits=-best[K][q][0] / log(2),
                tilt_at_search_boundary=best[K][q][1] in (tilts[0], tilts[-1]),
                requires_outward_replay=True) for q in grids[K]}
            records[str(K)] = dict(K=K, N=shapes[K].N, geometry=asdict(shapes[K]),
                threshold=shapes[K].N // 10, macros_per_region=shapes[K].macros_per_region,
                occupancy_values=list(grids[K]), witnesses=choices,
                tilts=list(map(str, sorted({best[K][q][1] for q in grids[K]}))))
        return dict(schema='rs16-length-regional-proposals-1', proposal_only=True,
            whole_code_certificate=False, has_numerical_upper_endpoints=False,
            state_bits=self.state_bits, physical_t=64, inner_distribution='uniform_gl',
            distance='1/10', zero_initial_state=True, final_flush=False,
            terminal='sum of all envelope coordinates', local_selection_activity='1/2',
            packet_activity='15/16', map_record=self.map_record,
            map_sha256=self.data['map_sha256'], source_sha256=self.sources,
            precision_for_local_envelopes=self.precision, input_tilts=list(map(str, tilts)),
            lengths=records, cached_local_tilts=len(self.locals),
            cached_regional_families=len(self.regionals),
            placement_evaluations=self.placement_evaluations,
            replay='Freshly recompute the selected (K,q,tilt) values with the exact '
                   'regional Arb evaluator and identical maps; these scores are not endpoints.',
            scope='Floating witness search only. Neither a probability bound nor a distance '
                  'certificate follows until fresh outward replay covers the desired occupancies.')
