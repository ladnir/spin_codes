"""Fresh RS16 proof replay with an explicitly selected state-17..22 adapter.

The outer law is four GF16 RS[16,8] rows per 128-bit group, independent
uniform GL16 symbol maps, and independent uniform group/regional shuffles.
The selected quadratic-prefix adapter defines the actual t64/s inner maps.
Each physical update is independent uniform GL(s,2). State starts at zero,
continues across all steps and regions, and is not flushed.

Only fresh replay covers every occupancy and may certify the 10%/40-bit
target. Component runs cover their explicit occupancy lists only. Numeric
helpers are reused with their actual state-s data; no module globals are
rebound, and no saved result or state-size label replaces fresh computation.
Exact placement uses the recorded linear or binary backend. Both compute
upper envelopes for the same chronological regional transition operator.
"""
from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
from time import monotonic

from flint import arb, ctx
import packet_inner_small_extension as small_inner
import packet_inner_quadratic_extension as quadratic_inner
import packet_regional_power as regional_power
import packet_rs_small_state as shared

q1, lengths, dense, low = shared.q1, shared.lengths, shared.dense, shared.low
arithmetic = shared.arithmetic
TAIL_SCHEMA = 'rs16-state-ladder-tail-1'
WHOLE_SCHEMA = 'rs16-state-ladder-complete-replay-1'


def adapter_for(bits):
    """Select a frozen constructor explicitly; never change another module."""
    if type(bits) is not int:
        raise ValueError('integer state dimension from 17 through 22 required')
    if bits in small_inner.SUPPORTED_BITS:
        return small_inner
    if bits in quadratic_inner.SUPPORTED_BITS:
        return quadratic_inner
    raise ValueError('state dimension must lie in 17..22')


def _prepare(adapter, bits, precision, data=None, map_record=None):
    if adapter is not adapter_for(bits) or (data is None) != (map_record is None):
        raise ValueError('matching explicit adapter and paired fresh data/record required')
    before = q1.source_snapshot()
    ctx.prec = precision
    if data is None:
        data, map_record = adapter.prepare(bits, birth_density='capped')
    adapter.authenticate(data, map_record)
    if (data['bits'] != bits or map_record['s'] != bits or data['birth_density'] != 'capped' or
            map_record['return_denominator'] != (1 << bits)-1):
        raise ValueError('actual state dimension and uniform-GL return denominator required')
    sources = q1.source_snapshot()
    if any(sources.get(path) != digest for path, digest in before.items()):
        raise RuntimeError('preexisting proof source changed during preparation')
    return data, map_record, sources


def _placement_backend(value):
    if type(value) is not str or value not in ('linear', 'binary'):
        raise ValueError('placement backend must be linear or binary')
    return value


def _options(K, bits, occupancies, tilts, precision, output, method, markers, placement_backend='linear'):
    _placement_backend(placement_backend)
    geometry, adapter = lengths.geometry(K), adapter_for(bits)
    tilts, output = lengths.checked_options(tilts, precision, output)
    if isinstance(occupancies, (str, bytes)) or isinstance(markers, (str, bytes)):
        raise ValueError('explicit occupancy and marker sequences required')
    qs = tuple(occupancies)
    if (not qs or any(type(q) is not int or not 3 <= q <= geometry.group_count for q in qs) or
            qs != tuple(sorted(set(qs))) or method not in ('exact', 'fugacity')):
        raise ValueError('sorted distinct tail occupancies and exact/fugacity method required')
    try:
        markers = tuple(map(Q, markers))
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('finite rational marker probabilities required') from error
    if (len(set(markers)) != len(markers) or any(not 0 < p <= 1 for p in markers) or
            (method == 'fugacity' and (not markers or (qs[0] < geometry.group_count and not any(p < 1 for p in markers))))):
        raise ValueError('distinct valid marker probabilities covering every requested q required')
    return geometry, adapter, qs, tuple(map(Q, tilts)), markers, output


def run_tail(*, K, bits, occupancies, tilts, precision=192, output=None,
             method='exact', marker_probabilities=(), placement_backend='linear', data=None, map_record=None):
    """Bound listed occupancies; placement_backend selects only exact placement."""
    geometry, adapter, qs, tilts, markers, output = _options(
        K, bits, occupancies, tilts, precision, output, method, marker_probabilities, placement_backend)
    start = monotonic()
    beta, counts = lengths.exact_outer()
    data, map_record, sources = _prepare(adapter, bits, precision, data, map_record)
    best, choices = dict.fromkeys(qs), dict.fromkeys(qs)
    logs = {q: arb(comb(geometry.group_count, q)).log() for q in qs} if method == 'fugacity' else {}
    log_beta = q1.kernel_t64.aq(beta).log() if method == 'fugacity' else None
    record = dict(shared._common(geometry, bits, map_record, precision, sources), schema=TAIL_SCHEMA,
        method=method, return_denominator=(1 << bits)-1, birth_density='capped',
        occupancy_covered=list(qs), every_requested_occupancy_checked=True,
        beta=str(beta), every_shell_checked=True, count_kind='exact_expected_shells',
        count_sha256=lengths.count_hash(counts), tilts=list(map(str, tilts)),
        marker_probabilities=list(map(str, markers)), trials=[],
        scope=f'Listed-occupancy first moment for ideal independent outer GL16 and inner GL({bits},2) '
              'setup with uniform routing, continuous zero-initial state, and no flush. No whole-code certificate.')
    if method == 'exact':
        record['placement_backend'] = placement_backend

    def consider(q, bound, witness):
        endpoint = arithmetic.dyadic(lengths.positive_endpoint(bound))
        if best[q] is None or arithmetic.dyadic_less(endpoint, best[q]):
            best[q], choices[q] = endpoint, witness

    for tilt in tilts:
        trial_start = monotonic()
        if method == 'exact':
            local = q1.kernel_t64.local_operators(data, tilt, activity=Q(1, 2))
            placement = q1.placement if placement_backend == 'linear' else regional_power.placement_power
            regional = placement(local, epochs=geometry.macros_per_region,
                windows=32, rounding=q1.rounded, maximum_groups=qs[-1])
            for q in qs:
                consider(q, lengths.occupancy_upper(regional, K=K, occupancy=q, beta=beta, tilt=tilt), dict(tilt=str(tilt)))
        else:
            candidates = dense.physical_candidates(data, tilt)
            for p in markers:
                if p == 1 and geometry.group_count not in best:
                    continue
                mixed, indices = dense.select_physical(candidates, Q(15, 16)*p)
                power = mixed ** (geometry.N//64)
                moment = sum((power[0, j] for j in range(power.ncols())), arb(0))
                if not moment.is_finite() or not moment > 0:
                    raise ArithmeticError('positive finite all-terminal-state moment required')
                base, slope = shared._fugacity_affine(geometry, tilt, p, moment.log(), log_beta)
                for q in ((geometry.group_count,) if p == 1 else qs):
                    logarithm = base if p == 1 else base + (1-geometry.regions)*logs[q] + q*slope
                    consider(q, q1.kernel_t64.up(logarithm.exp()),
                        dict(tilt=str(tilt), marker_probability=str(p), selected_birth_rows=list(indices)))
        if any(value is None for value in best.values()):
            raise ArithmeticError('requested occupancy lacks a valid endpoint')
        upper, margin = shared._summary(best.values(), precision)
        record['trials'].append(dict(tilt=str(tilt), union_upper=upper, margin_bits=margin,
                                    elapsed_seconds=monotonic()-trial_start))
        shared._stable(sources, precision)
        print(f'RS16/S{bits} K={K} {method} q={qs[0]}..{qs[-1]} tilt={tilt} margin_bits={margin}', flush=True)
    adapter.authenticate(data, map_record)
    shared._stable(sources, precision)
    record.update(union_upper=upper, margin_bits=margin,
        occupancy_uppers={str(q): list(best[q]) for q in qs}, occupancy_choices={str(q): choices[q] for q in qs},
        elapsed_seconds=monotonic()-start)
    if output is not None:
        arithmetic._write_new(output, record)
    return record


def recipe(K, bits, *, placement_backend='linear', **options):
    adapter_for(bits)
    placement_backend = _placement_backend(placement_backend)
    options.setdefault('q1_tilts', shared.DEFAULT_Q1)
    return dict(arithmetic.recipe(K, **options), state_bits=bits, placement_backend=placement_backend)


def _endpoints(record, geometry, bits, map_record, sources, precision, count_hash, beta, placement_backend=None):
    if placement_backend is not None:
        _placement_backend(placement_backend)
    if record.get('schema') in ('finite-packet-rs-state-q1-1', 'finite-packet-rs-state-q2-1'):
        return shared._component_endpoints(record, geometry, bits, map_record, sources, precision, count_hash)
    expected = shared._common(geometry, bits, map_record, precision, sources)
    if (json.dumps({k: record.get(k) for k in expected}, sort_keys=True) != json.dumps(expected, sort_keys=True) or
            record.get('schema') != TAIL_SCHEMA or record.get('method') not in ('exact', 'fugacity') or
            record.get('return_denominator') != (1 << bits)-1 or record.get('birth_density') != 'capped' or
            record.get('every_requested_occupancy_checked') is not True or record.get('every_shell_checked') is not True or
            record.get('count_kind') != 'exact_expected_shells' or record.get('count_sha256') != count_hash or
            record.get('beta') != str(beta)):
        raise ValueError('fresh actual-state tail metadata, sources, and outer counts must match')
    if record['method'] == 'exact':
        recorded_backend = _placement_backend(record.get('placement_backend'))
        if placement_backend is not None and recorded_backend != placement_backend:
            raise ValueError('exact placement backend must match the replay recipe')
    elif 'placement_backend' in record:
        raise ValueError('fugacity components do not use a placement backend')
    qs = record.get('occupancy_covered')
    if (not isinstance(qs, list) or not qs or any(type(q) is not int or not 3 <= q <= geometry.group_count for q in qs) or
            qs != sorted(set(qs)) or not isinstance(record.get('occupancy_uppers'), dict) or
            set(record['occupancy_uppers']) != {str(q) for q in qs}):
        raise ValueError('exactly one endpoint for every explicit tail occupancy required')
    return {q: arithmetic.dyadic(record['occupancy_uppers'][str(q)]) for q in qs}


def replay(output, *, K, bits, **options):
    """Prepare once and recompute all q; only a passing full union certifies."""
    start, plan = monotonic(), recipe(K, bits, **options)
    output, geometry, adapter = Path(output), lengths.geometry(K), adapter_for(bits)
    names = ['q1', 'q2'] + (['exact'] if plan['exact_occupancies'] else []) + ['dense']
    paths = {name: output.with_name(output.stem+'-'+name+'.json') for name in names}
    if output.exists() or any(path.exists() for path in paths.values()):
        raise ValueError('whole replay and every component output must be fresh')
    precision = plan['precision']
    data, map_record, sources = _prepare(adapter, bits, precision)
    fresh = dict(data=data, map_record=map_record, precision=precision)
    records = []
    for q in (1, 2):
        adapter.authenticate(data, map_record)
        records.append(low.run_low(q, geometry=geometry, tilts=plan['tilts'][f'q{q}'], output=paths[f'q{q}'], **fresh))
        shared._stable(sources, precision)
    if plan['exact_occupancies']:
        records.append(run_tail(K=K, bits=bits, occupancies=plan['exact_occupancies'],
            tilts=plan['tilts']['exact'], method='exact', placement_backend=plan['placement_backend'],
            output=paths['exact'], **fresh))
    records.append(run_tail(K=K, bits=bits, occupancies=range(plan['dense_interval'][0], geometry.group_count+1),
        tilts=plan['tilts']['dense'], method='fugacity', marker_probabilities=plan['marker_probabilities'],
        output=paths['dense'], **fresh))
    beta, counts = lengths.exact_outer()
    count_hash, values, selected = lengths.count_hash(counts), {}, {}
    for index, record in enumerate(records):
        for q, bound in _endpoints(record, geometry, bits, map_record, sources, precision,
                                  count_hash, beta, plan['placement_backend']).items():
            if q < 3 and q in values:
                raise ValueError('duplicate low-occupancy component')
            if q not in values or arithmetic.dyadic_less(bound, values[q]):
                values[q], selected[q] = bound, index
    if set(values) != set(range(1, geometry.group_count+1)):
        raise ValueError('full replay must cover every occupancy 1..L')
    upper, margin = shared._summary(values.values(), precision)
    tail, tail_margin = shared._summary((values[q] for q in range(3, geometry.group_count+1)), precision)
    adapter.authenticate(data, map_record)
    shared._stable(sources, precision)
    passes = arithmetic.dyadic_less(arithmetic.dyadic(upper), (1, -40))
    result = dict(shared._common(geometry, bits, map_record, precision, sources), schema=WHOLE_SCHEMA,
        target_margin_bits=40, target_minimum_distance=geometry.N//10+1, return_denominator=(1 << bits)-1,
        birth_density='capped', count_sha256=count_hash, occupancy_covered=[1, geometry.group_count],
        all_occupancies_covered=True, fresh_replay=True, target_met=passes, recipe=plan,
        occupancy_uppers={str(q): list(values[q]) for q in sorted(values)}, selected_components={str(q): selected[q] for q in sorted(values)},
        q1_upper=list(values[1]), q2_upper=list(values[2]), tail_upper=tail, union_upper=upper,
        margin_bits=margin, tail_margin_bits=tail_margin, elapsed_seconds=monotonic()-start,
        arithmetic='Exact sum of minimum dyadic endpoints; one final upward rounding.',
        components=[dict(path=str(paths[name].resolve()), schema=record['schema'],
            sha256=hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
            hash_scope='canonical record JSON', occupancy_covered=record['occupancy_covered']) for name, record in zip(names, records)],
        scope=f'Full first moment under ideal independent outer GL16 symbol maps, uniform group/regional routing, '
              f'and inner GL({bits},2) updates. Continuous zero-initial state, no flush, no seeded guarantee.')
    result['whole_code_certificate'] = passes
    arithmetic._write_new(output, result)
    print(f'RS16/S{bits} K={K} COMPLETE UNION margin_bits={margin} certificate={passes}', flush=True)
    return result
