"""Dense-occupancy bounds for the ideal RS16 packet ensemble at explicit K.

Each outer group maps 128 message bits to 64 four-bit packets. The map uses
four parallel GF16 RS[16,8] words and independent uniform GL16 symbol maps.
Group setups, group packet shuffles, and regional shuffles are independent.
The selected t64/s16 inner starts at zero, retains state between regions,
and uses a fresh independent uniform GL16 update at each physical step.
There is no final flush. The default geometry has L=8192 groups and R=64
regions, hence K=2^20, N=2^21, and 256 macro steps per region.

Fix q active message groups. Their expected output counting measure is
dominated by beta^q times independent uniform 256-bit group outputs, where
beta=2^256/(2^16-1)^8. The comparison includes artificial zero outputs without
removing active-group labels. Each region therefore places q uniform four-bit
packets in a uniform q-subset of its L slots. These regional laws are
independent, but the inner state is not reset between regions.

To bound this fixed-q law, fix any 0<p<1, mark each slot independently with
probability p, and give each marked slot an independent uniform four-bit packet.
Conditioning on exactly q markers gives the required regional input law.
The event has probability c_q=C(L,q)p^q(1-p)^(L-q). For every nonnegative
path functional F, E[F | exactly q markers] <= E[F]/c_q. The conditioning
events are independent between regions, so the loss is D_q^R, D_q=1/c_q.
The endpoint p=1 covers only q=L, with D_L=1. The original run() uses p=q/L;
run_fugacity() compares a shared grid of p values and retains the best bound.

Let T_j(lambda) be the fixed nonnegative local upper-envelope operator for
j nonzero packets in a uniform j-subset of the 32 macro slots. The iid
marker experiment gives packet activity theta=15p/16. Its macro operator
is B_p=sum_j Binomial(32,theta)(j) T_j(lambda). The same T_j family bounds
all paths before and after conditioning. Chronological products retain state.
For cutoff floor(N/10), the expected low-weight count at occupancy q is at
most

    C(L,q) beta^q exp(lambda*cutoff) D_q^R
        * e_zero B_p^(R*L/32) 1.

For a fixed p, the same matrix moment serves every q. Hypergeometric
composition of the two physical halves gives B_p=M_p^2, where M_p is the
Binomial(16,theta) mixture of their common physical occupancy operators.
run_fugacity() uses this identity and selects complete valid birth rows with
a floating continuation potential at activity theta. Row selection changes
the envelope representation, not the input or setup law. Cached row choices
are never combined through entrywise minima.

The comparison is between exact polynomials in the fixed T_j matrices,
not independently rounded regional bounds. Arb rounds the final numerical
calculation outward. Every terminal envelope coordinate represents mass,
including the uniform-density coordinate whose state-count factor is already
included. Thus the terminal vector is all ones, without an extra factor.

This entry point records only its requested q interval. It never asserts a
whole-code certificate or a guarantee for a seeded implementation. The
uniform outer envelope and local transition envelopes are upper bounds,
not exact failure probabilities. No K16 receipt or numerical bound is reused.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
from time import monotonic

from flint import arb, arb_mat, ctx
import numpy as np
import packet_q1 as q1
from rs_outer import expected_group_support_counts
from rs_uniform_envelope import UniformInputEnvelope


SCHEMA = 'rs16-conditioned-iid-dense-occupancy-1'
FUGACITY_SCHEMA = 'rs16-fugacity-iid-dense-occupancy-1'
DEFAULT_TILTS = ('.0032', '.0064', '.0128', '.0256', '.0512', '.1024',
                 '.2048', '.4096', '.8192', '1.6384', '2.1972246')
DEFAULT_MARKER_PROBABILITIES = (
    '1/256', '3/512', '1/128', '3/256', '1/64', '3/128', '1/32', '3/64',
    '1/16', '3/32', '1/8', '3/16', '1/4', '3/8', '1/2', '5/8', '3/4',
    '7/8', '15/16', '31/32', '1')


def geometry(group_count=8192):
    """Use complete macro steps and the fixed 128-to-256 RS16 outer group."""
    return q1.Geometry(group_count, 64, 128)


def _occupancy(groups, q):
    if type(groups) is not int or groups < 1 or type(q) is not int or not 0 <= q <= groups:
        raise ValueError('integer L>=1 and 0<=q<=L required')


def _marker(groups, q, marker_probability):
    _occupancy(groups, q)
    p = Q(q, groups) if marker_probability is None else Q(marker_probability)
    if not 0 <= p <= 1 or (p == 0 and q != 0) or (p == 1 and q != groups):
        raise ValueError('marker probability must give the conditioning event positive probability')
    return p


def conditioning_loss_exact(groups, q, marker_probability=None):
    """Exact D_q for algebra tests; large production values use its Arb log."""
    p = _marker(groups, q, marker_probability)
    return 1 / (comb(groups, q) * p**q * (1 - p)**(groups - q))


def log_conditioning_loss(groups, q, *, log_binomial=None, marker_probability=None):
    """Enclose log D_q, including both certain-event endpoints without 0 log 0."""
    p_exact = _marker(groups, q, marker_probability)
    if p_exact in (0, 1):
        return arb(0)
    if log_binomial is None:
        log_binomial = arb(comb(groups, q)).log()
    p = q1.kernel_t64.aq(p_exact)
    complement = q1.kernel_t64.aq(1 - p_exact)
    # Retain the complete interval: negating an upper endpoint of log C is unsafe.
    return -log_binomial - q * p.log() - (groups - q) * complement.log()


def binomial_weights(degree, probability):
    """Return exact rational weights, including probability zero and one."""
    if type(degree) is not int or degree < 0:
        raise ValueError('nonnegative integer binomial degree required')
    probability = Q(probability)
    if not 0 <= probability <= 1:
        raise ValueError('binomial probability must lie in [0,1]')
    a, b = probability.numerator, probability.denominator
    denominator = b**degree
    weights = tuple(Q(comb(degree, j) * a**j * (b - a)**(degree - j), denominator)
                    for j in range(degree + 1))
    if sum(weights, Q(0)) != 1:
        raise ArithmeticError('exact binomial weights do not normalize')
    return weights


def iid_macro(local, probability, *, matrix=arb_mat,
              rational=q1.kernel_t64.aq, rounding=q1.rounded):
    """Average a complete local occupancy family with exact binomial weights."""
    if not local or local[0].nrows() < 1:
        raise ValueError('nonempty local matrix family required')
    size = local[0].nrows()
    if any(value.nrows() != size or value.ncols() != size for value in local):
        raise ValueError('matching square local matrices required')
    result = matrix(size, size)
    for weight, value in zip(binomial_weights(len(local) - 1, probability), local):
        if weight:
            result += value * rational(weight)
    return rounding(result)


def physical_candidates(data, tilt):
    """Cache complete outward physical birth-row candidates once per tilt."""
    physical = q1.kernel_t64.authenticate(data)
    if physical['windows'] != 16 or Q(tilt) <= 0 or physical.get('birth_density') != 'capped':
        raise ValueError('capped t64 physical data and positive tilt required')
    z = (-q1.kernel_t64.aq(Q(tilt))).exp()
    local = q1.kernel_t64.sparse_kernel.outward_at_z(physical, z)
    density = q1.kernel_t64.kernel_birth_density
    numerators, denominators = density.conditional_laws(physical, z)
    size = local[0].nrows()
    rows = tuple(tuple(tuple(row) for row in density.outward_rows(physical,
        [matrix[0, k] for k in range(size)], numerators[:, j], denominators[j]))
        for j, matrix in enumerate(local))
    arrays = np.array([[[float(matrix[i, j]) for j in range(size)]
                        for i in range(size)] for matrix in local])
    row_arrays = tuple(np.array([[float(value) for value in row] for row in options])
                       for options in rows)
    return dict(local=local, rows=rows, arrays=arrays, row_arrays=row_arrays,
                precision=ctx.prec, tilt=str(Q(tilt)), map_sha256=data['map_sha256'])


def select_physical(candidates, activity, *, selection_activity=None):
    """Choose whole birth rows, then mix at the actual input activity.

    The optional selection activity affects only the floating proposal.
    Every selected row is a complete outward representation of the same law.
    """
    activity = Q(activity)
    selection_activity = activity if selection_activity is None else Q(selection_activity)
    if (not 0 <= activity <= 1 or not 0 <= selection_activity <= 1
            or candidates['precision'] != ctx.prec):
        raise ValueError('valid activities and unchanged candidate precision required')
    proposal_weights = np.array([float(weight) for weight in
                                 binomial_weights(16, selection_activity)])
    potential = q1.kernel_t64.kernel_birth_density.continuation(
        np.tensordot(proposal_weights, candidates['arrays'], axes=1))
    indices = tuple(int(np.argmin(options @ potential)) for options in candidates['row_arrays'])
    local = candidates['local']
    mixed = iid_macro(local, activity, rounding=lambda value: value)
    actual_weights = tuple(q1.kernel_t64.aq(weight) for weight in binomial_weights(16, activity))
    for target in range(mixed.ncols()):
        mixed[0, target] = sum((weight * candidates['rows'][j][indices[j]][target]
                               for j, weight in enumerate(actual_weights)), arb(0))
    mixed = q1.rounded(mixed)
    if any(not mixed[i, j].is_finite() or not mixed[i, j] >= 0
           for i in range(mixed.nrows()) for j in range(mixed.ncols())):
        raise ArithmeticError('finite nonnegative physical mixture required')
    return mixed, indices


def coefficient_log_bound(*, finite_geometry, q, tilt, beta, marker_probability,
                          log_moment, log_binomial=None):
    """Apply the marker coefficient bound to an enclosed common matrix moment."""
    groups, regions = finite_geometry.group_count, finite_geometry.regions
    p = _marker(groups, q, marker_probability)
    if q == 0 or p == 0 or Q(tilt) <= 0 or Q(beta) <= 0:
        raise ValueError('positive occupancy, marker probability, tilt and beta required')
    if log_binomial is None:
        log_binomial = arb(comb(groups, q)).log()
    base = log_moment + q1.kernel_t64.aq(Q(tilt)) * (finite_geometry.N // 10)
    log_beta = q1.kernel_t64.aq(Q(beta)).log()
    if p == 1:
        return base + groups * log_beta
    log_p, log_complement = q1.kernel_t64.aq(p).log(), q1.kernel_t64.aq(1 - p).log()
    return (base - regions * groups * log_complement + (1 - regions) * log_binomial
            + q * (log_beta - regions * (log_p - log_complement)))


def occupancy_upper(local, *, finite_geometry, q, tilt, beta, log_binomial=None):
    """Evaluate one valid occupancy bound; all local operators use the same tilt."""
    groups = finite_geometry.group_count
    _occupancy(groups, q)
    tilt, beta = Q(tilt), Q(beta)
    if q == 0 or tilt <= 0 or beta <= 0 or len(local) != finite_geometry.macro_windows + 1:
        raise ValueError('positive occupancy/tilt/beta and complete macro family required')
    if log_binomial is None:
        log_binomial = arb(comb(groups, q)).log()
    log_loss = log_conditioning_loss(groups, q, log_binomial=log_binomial)
    activity = Q(15 * q, 16 * groups)
    mixed = iid_macro(local, activity)
    total_macros = finite_geometry.regions * finite_geometry.macros_per_region
    power = mixed**total_macros
    moment = sum((power[0, j] for j in range(power.ncols())), arb(0))
    if not moment.is_finite() or not moment > 0:
        raise ArithmeticError('positive finite moment enclosure required; increase precision if needed')
    logarithm = (log_binomial + q * q1.kernel_t64.aq(beta).log()
                 + q1.kernel_t64.aq(tilt) * (finite_geometry.N // 10)
                 + finite_geometry.regions * log_loss + moment.log())
    upper = q1.kernel_t64.up(logarithm.exp())
    if not upper.is_finite() or not upper > 0:
        raise ArithmeticError('positive finite occupancy endpoint required')
    return upper


def _options(group_count, q_min, q_max, precision, tilts, output):
    finite_geometry = geometry(group_count)
    if q_max is None:
        q_max = group_count
    if (type(q_min) is not int or type(q_max) is not int or not 1 <= q_min <= q_max <= group_count
            or type(precision) is not int or precision < 128):
        raise ValueError('valid nonzero occupancy interval and integer precision >=128 required')
    if isinstance(tilts, (str, bytes)):
        raise ValueError('a sequence of distinct positive rational tilts is required')
    try:
        values = tuple(Q(t) for t in tilts)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('finite rational tilts required') from error
    if not values or any(t <= 0 for t in values) or len(set(values)) != len(values):
        raise ValueError('distinct positive rational tilts required')
    path = Path(output) if output is not None else None
    if path is not None and path.exists():
        raise ValueError('output must be fresh')
    return finite_geometry, q_max, tuple(map(str, values)), path


def run(*, group_count=8192, q_min=33, q_max=None, precision=256,
        tilts=DEFAULT_TILTS, output=None, data=None, map_record=None):
    """Evaluate the interval; shared data/map_record must come from fresh prepare()."""
    finite_geometry, q_max, tilts, output = _options(
        group_count, q_min, q_max, precision, tilts, output)
    if (data is None) != (map_record is None):
        raise ValueError('shared prepared data and map record must be supplied together')
    start = monotonic()
    ctx.prec = precision
    outer = UniformInputEnvelope(n=16, k=8, packet_bits=4, packets_per_symbol=4)
    outer.verify_shell_domination()
    counts = expected_group_support_counts(n=16, k=8, packet_bits=4, packets_per_symbol=4)
    if (outer.regions != finite_geometry.regions or outer.message_bits != finite_geometry.group_dimension
            or sum(counts) != (1 << finite_geometry.group_dimension) - 1 or counts[0] != 0):
        raise ArithmeticError('RS16 exact outer counts do not match finite geometry')
    if data is None:
        data, map_record = q1.kernel_t64.prepare(birth_density='capped')
    if not map_record:
        raise ArithmeticError('fresh selected-map record required')
    q1.kernel_t64.authenticate(data)
    if (data['physical_step_bits'] != 64 or data['bits'] != 16
            or data['macro_windows'] != 32 or data['birth_density'] != 'capped'
            or map_record.get('map_sha256') != data['map_sha256']):
        raise ValueError('matching freshly prepared capped t64/s16 map and record required')
    sources = q1.source_snapshot()
    best = {q: None for q in range(q_min, q_max + 1)}
    choices = dict.fromkeys(best)
    # These interval logs are independent of the tilt; do not round them upward
    # before passing them to the conditioning expression, where their sign flips.
    log_counts = {}
    count = comb(group_count, q_min)
    for q in best:
        log_counts[q] = arb(count).log()
        if q < q_max:
            count = count * (group_count - q) // (q + 1)
    record = dict(schema=SCHEMA, outer='rs16', K=finite_geometry.K, N=finite_geometry.N,
        geometry=asdict(finite_geometry), group_output_bits=256,
        groups=group_count, regions=finite_geometry.regions,
        group_dimension=finite_geometry.group_dimension,
        physical_t=64, state_bits=16, physical_steps=finite_geometry.N // 64,
        macro_t=128, macro_steps=finite_geometry.N // 128,
        macros_per_region=finite_geometry.macros_per_region, physical_steps_per_macro=2,
        threshold=finite_geometry.N // 10, distance='1/10',
        q_min=q_min, q_max=q_max, occupancy_covered=[q_min, q_max],
        evaluated_every_integer_occupancy=True, whole_code_certificate=False,
        zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        terminal='sum of all mass-envelope coordinates; no flush',
        outer_comparison=outer.metadata(), beta=str(outer.beta), every_shell_checked=True,
        count_kind='shells',
        count_sha256=hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest(),
        conditioning='per-region iid markers p=q/L conditioned on exactly q markers',
        conditioning_loss='[1/(C(L,q)*(q/L)^q*(1-q/L)^(L-q))]^regions; D_L=1',
        iid_packet_activity='15*q/(16*L)', local_selection_activity='1/2',
        local_selection_activity_role='valid birth-row proposal only, not the input law',
        regional_placement='domination by conditioned iid marker law, not exact equality',
        arithmetic='Arb intervals with outward final dyadic upper endpoints',
        precision=precision, tilts=list(tilts), source_sha256=sources,
        map_record=map_record, trials=[],
        scope='Expected low-weight count for the explicit occupancy interval in the ideal '
              'RS16 ensemble, with independent uniform GL16 symbol maps, group packet '
              'shuffles, regional shuffles and physical inner updates. Not a bound '
              'for every setup, a seeded implementation, or the complete code.')
    print(f'DENSE RS16 K={finite_geometry.K} N={finite_geometry.N} L={group_count} '
          f'q={q_min}..{q_max} macros_per_region={finite_geometry.macros_per_region}', flush=True)
    for tilt in tilts:
        trial_start = monotonic()
        local = q1.kernel_t64.local_operators(data, Q(tilt), activity=Q(1, 2))
        if (len(local) != 33 or any(value.nrows() != local[0].nrows()
                or value.ncols() != local[0].nrows() for value in local)
                or any(not value[i, j].is_finite() or not value[i, j] >= 0
                       for value in local for i in range(value.nrows()) for j in range(value.ncols()))):
            raise ArithmeticError('complete finite nonnegative local upper envelopes required')
        for q in best:
            upper = occupancy_upper(local, finite_geometry=finite_geometry, q=q,
                                    tilt=tilt, beta=outer.beta, log_binomial=log_counts[q])
            if best[q] is None or upper < best[q]:
                best[q], choices[q] = upper, tilt
            if (q - q_min + 1) % 512 == 0:
                print(f'tilt={tilt} evaluated through q={q} '
                      f'elapsed={monotonic() - trial_start:.2f}s', flush=True)
        total = q1.kernel_t64.up(sum(best.values(), arb(0)))
        worst = sorted(best, key=lambda q: float(best[q].log()), reverse=True)[:8]
        margin = str(-total.log() / arb(2).log())
        record['trials'].append(dict(tilt=tilt, union_upper=q1.endpoint(total),
            margin_bits=margin, worst_q=worst, elapsed_seconds=monotonic() - trial_start))
        if sources != q1.source_snapshot():
            raise RuntimeError('loaded local mathematical source changed during evaluation')
        print(f'tilt={tilt} complete q={q_min}..{q_max} union_margin={margin} '
              f'worst_q={worst} elapsed={monotonic() - trial_start:.2f}s', flush=True)
    record.update(union_upper=q1.endpoint(total), margin_bits=margin,
        occupancy_uppers={str(q): q1.endpoint(value) for q, value in best.items()},
        occupancy_choices={str(q): tilt for q, tilt in choices.items()},
        elapsed_seconds=monotonic() - start)
    if sources != q1.source_snapshot():
        raise RuntimeError('loaded local mathematical source changed before receipt completion')
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(record, indent=2) + '\n')
    return record


def run_fugacity(*, group_count=8192, q_min=33, q_max=None, precision=256,
        tilts=DEFAULT_TILTS, marker_probabilities=DEFAULT_MARKER_PROBABILITIES,
        selection_activity=None, distance=Q(1, 10), output=None, data=None, map_record=None):
    """Evaluate a shared (tilt,p) grid with one physical matrix power per pair.

    All q values in the interval are evaluated for each p<1. The endpoint
    p=1 applies only to q=L. Optional shared data must come from fresh prepare().
    Optional selection_activity fixes a proposal control, not the input law.
    """
    finite_geometry, q_max, tilts, output = _options(
        group_count, q_min, q_max, precision, tilts, output)
    distance = Q(distance)
    if not 0 < distance < Q(1, 2):
        raise ValueError('exact distance in (0,1/2) required')
    threshold = finite_geometry.N * distance.numerator // distance.denominator
    if (data is None) != (map_record is None):
        raise ValueError('shared prepared data and map record must be supplied together')
    if isinstance(marker_probabilities, (str, bytes)):
        raise ValueError('a sequence of distinct marker probabilities is required')
    try:
        markers = tuple(Q(p) for p in marker_probabilities)
        selection_activity = None if selection_activity is None else Q(selection_activity)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('finite rational marker/proposal probabilities required') from error
    if (not markers or len(set(markers)) != len(markers) or any(not 0 < p <= 1 for p in markers)
            or (q_min < group_count and not any(p < 1 for p in markers))
            or (selection_activity is not None and not 0 <= selection_activity <= 1)):
        raise ValueError('valid marker grid covering the interval and selection activity required')
    start = monotonic()
    ctx.prec = precision
    outer = UniformInputEnvelope(n=16, k=8, packet_bits=4, packets_per_symbol=4)
    outer.verify_shell_domination()
    counts = expected_group_support_counts(n=16, k=8, packet_bits=4, packets_per_symbol=4)
    if (outer.regions != finite_geometry.regions or outer.message_bits != finite_geometry.group_dimension
            or sum(counts) != (1 << finite_geometry.group_dimension) - 1 or counts[0] != 0):
        raise ArithmeticError('RS16 counts do not match finite geometry')
    if data is None:
        data, map_record = q1.kernel_t64.prepare(birth_density='capped')
    q1.kernel_t64.authenticate(data)
    if (not map_record or data['physical_step_bits'] != 64 or data['bits'] != 16
            or data['macro_windows'] != 32 or data['birth_density'] != 'capped'
            or map_record.get('map_sha256') != data['map_sha256']):
        raise ValueError('matching freshly prepared capped t64/s16 map and record required')
    sources = q1.source_snapshot()
    best = {q: None for q in range(q_min, q_max + 1)}
    choices = dict.fromkeys(best)
    log_counts = {}
    count = comb(group_count, q_min)
    for q in best:
        log_counts[q] = arb(count).log()
        if q < q_max:
            count = count * (group_count - q) // (q + 1)
    record = dict(schema=FUGACITY_SCHEMA, outer='rs16', K=finite_geometry.K, N=finite_geometry.N,
        geometry=asdict(finite_geometry), group_output_bits=256, group_dimension=128,
        groups=group_count, regions=finite_geometry.regions,
        physical_t=64, state_bits=16, physical_steps=finite_geometry.N // 64,
        macro_t=128, macro_steps=finite_geometry.N // 128,
        physical_steps_per_macro=2, macros_per_region=finite_geometry.macros_per_region,
        threshold=threshold, distance=str(distance), q_min=q_min, q_max=q_max,
        occupancy_covered=[q_min, q_max], evaluated_every_integer_occupancy=True,
        whole_code_certificate=False, fresh_computation=True,
        zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        terminal='sum of all mass-envelope coordinates; no flush',
        outer_comparison=outer.metadata(), beta=str(outer.beta), every_shell_checked=True,
        count_kind='shells',
        count_sha256=hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest(),
        conditioning='per-region iid markers with freely chosen p, conditioned on exactly q markers',
        conditioning_loss='[1/(C(L,q)*p^q*(1-p)^(L-q))]^regions; p=1 only for q=L',
        iid_packet_activity='15*p/16',
        local_selection_activity='15*p/16' if selection_activity is None else str(selection_activity),
        local_selection_activity_role='valid complete birth-row proposal only, not the input law',
        regional_placement='iid coefficient majorant with one inverse conditioning probability per region',
        matrix_identity='Binomial32 mixture of two-half convolution equals the square of Binomial16 physical mixture',
        arithmetic='Arb intervals with outward final dyadic upper endpoints', precision=precision,
        tilts=list(tilts), marker_probabilities=list(map(str, markers)),
        source_sha256=sources, map_record=map_record, trials=[],
        scope='Expected low-weight count over the explicit occupancy interval in the ideal '
              'RS16 ensemble with independent uniform GL16 symbol maps, group packet shuffles, '
              'regional shuffles, and physical inner updates. Not a seed-specific or whole-code certificate.')
    log_beta = q1.kernel_t64.aq(outer.beta).log()
    total_physical = finite_geometry.N // 64
    print(f'FUGACITY RS16 K={finite_geometry.K} N={finite_geometry.N} L={group_count} '
          f'q={q_min}..{q_max} p_grid={len(markers)}', flush=True)
    for tilt in tilts:
        candidates = physical_candidates(data, Q(tilt))
        for p in markers:
            if p == 1 and q_max != group_count:
                continue
            trial_start = monotonic()
            activity = Q(15, 16) * p
            mixed, indices = select_physical(candidates, activity, selection_activity=selection_activity)
            power = mixed**total_physical
            moment = sum((power[0, j] for j in range(power.ncols())), arb(0))
            if not moment.is_finite() or not moment > 0:
                raise ArithmeticError('positive finite moment required; increase precision if needed')
            log_moment = moment.log()
            base = log_moment + q1.kernel_t64.aq(Q(tilt)) * threshold
            selected_qs = (group_count,) if p == 1 else best
            if p < 1:
                log_p = q1.kernel_t64.aq(p).log()
                log_complement = q1.kernel_t64.aq(1 - p).log()
                base -= finite_geometry.regions * group_count * log_complement
                slope = log_beta - finite_geometry.regions * (log_p - log_complement)
            for q in selected_qs:
                logarithm = (base + group_count * log_beta if p == 1 else
                    base + q * slope + (1 - finite_geometry.regions) * log_counts[q])
                upper = q1.kernel_t64.up(logarithm.exp())
                if not upper.is_finite() or not upper > 0:
                    raise ArithmeticError('positive finite occupancy endpoint required')
                if best[q] is None or upper < best[q]:
                    best[q] = upper
                    choices[q] = dict(tilt=tilt, marker_probability=str(p),
                        selection_activity=str(activity if selection_activity is None else selection_activity))
            complete = all(value is not None for value in best.values())
            trial = dict(tilt=tilt, marker_probability=str(p), iid_packet_activity=str(activity),
                selection_activity=str(activity if selection_activity is None else selection_activity),
                selected_birth_rows=list(indices), moment_upper=q1.endpoint(moment),
                elapsed_seconds=monotonic() - trial_start, current_interval_complete=complete)
            if complete:
                total = q1.kernel_t64.up(sum(best.values(), arb(0)))
                worst = sorted(best, key=lambda q: float(best[q].log()), reverse=True)[:8]
                margin = str(-total.log() / arb(2).log())
                trial.update(union_upper=q1.endpoint(total), margin_bits=margin, worst_q=worst)
                print(f'tilt={tilt} p={p} union_margin={margin} worst_q={worst} '
                      f'elapsed={monotonic() - trial_start:.2f}s', flush=True)
            record['trials'].append(trial)
        if sources != q1.source_snapshot():
            raise RuntimeError('loaded local mathematical source changed during fugacity evaluation')
    if any(value is None for value in best.values()):
        raise ArithmeticError('marker grid did not cover every requested occupancy')
    record.update(union_upper=q1.endpoint(total), margin_bits=margin,
        occupancy_uppers={str(q): q1.endpoint(value) for q, value in best.items()},
        occupancy_choices= {str(q): choice for q, choice in choices.items()},
        elapsed_seconds=monotonic() - start)
    if sources != q1.source_snapshot():
        raise RuntimeError('loaded local mathematical source changed before receipt completion')
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(record, indent=2) + '\n')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', type=int, default=8192)
    parser.add_argument('--q-min', type=int, default=33)
    parser.add_argument('--q-max', type=int)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--tilts', nargs='+', default=DEFAULT_TILTS)
    parser.add_argument('--fugacity', action='store_true', help='use a shared tunable marker probability grid')
    parser.add_argument('--marker-probabilities', nargs='+', default=DEFAULT_MARKER_PROBABILITIES)
    parser.add_argument('--selection-activity', help='optional fixed birth-row proposal activity; default15p/16')
    parser.add_argument('--distance', default='1/10', help='exact distance for fugacity experiments; default1/10')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.fugacity:
        run_fugacity(group_count=args.groups, q_min=args.q_min, q_max=args.q_max,
            precision=args.precision, tilts=args.tilts, output=args.output,
            marker_probabilities=args.marker_probabilities, selection_activity=args.selection_activity,
            distance=args.distance)
    else:
        if args.selection_activity is not None or Q(args.distance) != Q(1, 10):
            parser.error('--selection-activity and nondefault --distance require --fugacity')
        run(group_count=args.groups, q_min=args.q_min, q_max=args.q_max,
            precision=args.precision, tilts=args.tilts, output=args.output)


if __name__ == '__main__':
    main()
