"""Compare explicit larger MDS outers without changing the t64 inner.

This isolated helper preserves the retained proof sources.  ``estimate``
uses logarithmic floating arithmetic and makes no certificate claim.
``replay`` freshly bounds only its listed occupancies with Arb arithmetic.
Both apply the same uniform-input measure majorant, including artificial
zero words, at every occupancy.  In particular, q=1 and q=2 are permitted;
no exact-shell endpoint from a different outer is imported.

An active outer group has expected multiplicity measure at most beta times
uniform binary input.  For q active groups, each region therefore has a
Bin(q,15/16) count of nonzero four-bit packets.  The exact placement average
retains distinct regional slots.  Inner state starts at zero, persists
through every physical step and region, and contributes all final mass.

The caller must justify the MDS envelope for its concrete outer.  Eight
parallel GF16 RS[16,8] rows and four parallel GF256 RS[16,8] rows give the
same expected output measure after independent uniform GL32 symbol maps.
For any zero-coordinate set of size j, both have 2^(32*max(8-j,0)) words.
Inclusion-exclusion determines the same exact symbol-support counts; GL32
then makes every nonzero symbol value uniform.  This observation concerns
the averaged first moment, not equality of sampled code distributions.
"""
from __future__ import annotations

from fractions import Fraction as Q
from math import comb, isfinite, log

from flint import arb, ctx
import numpy as np
import packet_q1 as q1
import packet_rs_state_sparse as sparse
from packet_regional_log import upper_logs, log_placement, regional_uniform_log, log_power_matrix
from packet_regional_power import placement_power
from packet_uniform_tail import regional_uniform
from packet_outer_geometry_proposal import upper_arrays
from packet_rs_length_proposal import scaled_placement, regional_uniform as scaled_uniform
from packet_rs_s20_proposal import log_power
from rs_uniform_envelope import UniformInputEnvelope


def geometry(K, envelope, *, windows=32, distance=Q(1, 10)):
    """Validate an explicit rate-one-half four-bit-packet configuration."""
    if (not isinstance(envelope, UniformInputEnvelope) or envelope.packet_bits != 4
            or envelope.n != 2*envelope.k):
        raise ValueError('explicit rate-one-half four-bit-packet MDS envelope required')
    if type(K) is not int or K <= 0 or K % envelope.message_bits:
        raise ValueError('K must be a positive multiple of the outer dimension')
    groups = K // envelope.message_bits
    if type(windows) is not int or windows <= 0 or groups % windows:
        raise ValueError('complete macro steps in each region required')
    distance = Q(distance)
    if not 0 < distance < Q(1, 2):
        raise ValueError('target distance must lie strictly between zero and one half')
    N = groups*envelope.output_bits
    return dict(K=K, N=N, groups=groups, regions=envelope.regions,
        outer_dimension=envelope.message_bits, outer_length=envelope.output_bits,
        outer_symbols=[envelope.n, envelope.k], symbol_bits=envelope.symbol_bits,
        windows=windows, epochs_per_region=groups//windows, distance=str(distance),
        cutoff=(N*distance.numerator)//distance.denominator,
        beta=str(envelope.beta), packet_activity='15/16', zero_initial_state=True,
        continuous_state_across_regions=True, final_flush=False,
        terminal='sum of all coordinates')


def _occupancies(values, groups):
    values = tuple(values)
    if (not values or tuple(sorted(set(values))) != values
            or any(type(q) is not int or not 1 <= q <= groups for q in values)):
        raise ValueError('ordered distinct feasible positive occupancies required')
    return values


def _tilt(value):
    value = Q(value)
    if value <= 0 or not isfinite(float(value)):
        raise ValueError('finite positive rational tilt required')
    return value


def estimate(local, *, K, envelope, occupancies, tilt, windows=32, distance=Q(1, 10), backend='log'):
    """Log-domain proposal from complete local Arb upper-envelope matrices."""
    meta = geometry(K, envelope, windows=windows, distance=distance)
    qs, tilt = _occupancies(occupancies, meta['groups']), _tilt(tilt)
    if backend not in ('log', 'scaled-or-log'):
        raise ValueError('backend must be log or scaled-or-log')
    moments = {}
    actual_backend = 'log'
    if backend == 'scaled-or-log':
        try:
            regional = scaled_placement(upper_arrays(local), qs[-1],
                epochs=meta['epochs_per_region'], windows=windows)
            for q in qs:
                matrix, scale = scaled_uniform(regional, q)
                with np.errstate(over='raise', invalid='raise', divide='raise', under='raise'):
                    moments[q] = log_power(matrix, envelope.regions) + envelope.regions*scale
            actual_backend = 'scaled'
        except FloatingPointError:
            moments = {}
    if not moments:
        regional = log_placement(upper_logs(local), qs[-1],
            epochs=meta['epochs_per_region'], windows=windows)
        moments = {q: log_power_matrix(regional_uniform_log(regional, q), envelope.regions) for q in qs}
    beta = envelope.beta
    log_beta = log(beta.numerator)-log(beta.denominator)
    result = {}
    for q in qs:
        moment = moments[q]
        value = log(comb(meta['groups'], q))+q*log_beta+float(tilt)*meta['cutoff']+moment
        if not isfinite(value):
            raise FloatingPointError('nonfinite floating objective')
        result[str(q)] = dict(estimated_margin_bits=-value/log(2), log_moment=moment)
    return dict(meta, schema='larger-outer-log-cost-proposal-1', tilt=str(tilt), backend=actual_backend,
        witnesses=result, proposal_only=True, whole_code_certificate=False,
        scope='Floating witness search only; selected witnesses require outward replay.')


def replay(local, *, K, envelope, occupancies, tilt, windows=32, distance=Q(1, 10)):
    """Fresh positive Arb bound from complete authenticated local operators.

    This low-level evaluator does not authenticate maps or pin sources;
    callers must use Model.replay or provide those checks themselves.
    """
    meta = geometry(K, envelope, windows=windows, distance=distance)
    qs, tilt = _occupancies(occupancies, meta['groups']), _tilt(tilt)
    if ctx.prec < 192:
        raise ValueError('outward replay requires at least 192-bit interval precision')
    family = placement_power(local, epochs=meta['epochs_per_region'],
        windows=windows, maximum_groups=qs[-1])
    aq, up = q1.kernel_t64.aq, q1.kernel_t64.up
    values = {}
    for q in qs:
        matrix = regional_uniform(family, q)**meta['regions']
        moment = sum((matrix[0, j] for j in range(matrix.ncols())), arb(0))
        value = up(comb(meta['groups'], q)*aq(envelope.beta)**q*
            (aq(tilt)*meta['cutoff']).exp()*moment)
        if not value.is_finite() or not value > 0:
            raise ArithmeticError('positive finite interval endpoint required')
        values[str(q)] = dict(upper=q1.endpoint(value), margin_bits=str(-value.log()/arb(2).log()))
    return dict(meta, schema='larger-outer-uniform-occupancy-component-1',
        tilt=str(tilt), precision=ctx.prec, occupancy_covered=list(qs),
        values=values, fresh_computation=True, whole_code_certificate=False,
        scope='Listed-occupancy uniform-majorant component only; not a whole certificate.')


class Model:
    """Cache fresh local operators across explicit outer and length choices."""
    def __init__(self, data, map_record, *, precision=192):
        if type(precision) is not int or precision < 192:
            raise ValueError('at least 192-bit local interval precision required')
        self.bits = sparse.validated(data, map_record)
        if not 16 <= self.bits <= 22:
            raise ValueError('selected t64/s16..22 inner required')
        self.data, self.map_record, self.precision = data, map_record, precision
        self.sources = sparse.source_snapshot(map_record)
        self.locals = {}

    def authenticate(self):
        sparse.validated(self.data, self.map_record)
        if self.sources != sparse.source_snapshot(self.map_record):
            raise RuntimeError('loaded mathematical sources changed during comparison')

    def local(self, tilt):
        self.authenticate()
        tilt = _tilt(tilt)
        if tilt not in self.locals:
            previous = ctx.prec
            try:
                ctx.prec = self.precision
                self.locals[tilt] = q1.kernel_t64.local_operators(self.data, tilt, activity=Q(1, 2))
            finally:
                ctx.prec = previous
            self.authenticate()
        return self.locals[tilt]

    def estimate(self, *, tilt, **kwargs):
        if kwargs.get('windows', 32) != 32:
            raise ValueError('the authenticated t64 two-step model has exactly 32 packet slots')
        result = estimate(self.local(tilt), tilt=tilt, **kwargs)
        self.authenticate()
        return dict(result, state_bits=self.bits, physical_t=64,
            map_record=self.map_record, source_sha256=self.sources)

    def replay(self, *, tilt, **kwargs):
        if kwargs.get('windows', 32) != 32:
            raise ValueError('the authenticated t64 two-step model has exactly 32 packet slots')
        local = self.local(tilt)
        previous = ctx.prec
        try:
            ctx.prec = self.precision
            result = replay(local, tilt=tilt, **kwargs)
        finally:
            ctx.prec = previous
        self.authenticate()
        return dict(result, state_bits=self.bits, physical_t=64,
            map_record=self.map_record, source_sha256=self.sources)


def prepare(bits):
    """Freshly prepare one selected inner; no retained numerical table is read."""
    if type(bits) is not int or not 16 <= bits <= 22:
        raise ValueError('integer state dimension in16..22 required')
    if bits == 16:
        return q1.kernel_t64.prepare(birth_density='capped')
    if bits <= 18:
        import packet_inner_small_extension as adapter
    else:
        import packet_inner_quadratic_extension as adapter
    return adapter.prepare(bits, birth_density='capped')
