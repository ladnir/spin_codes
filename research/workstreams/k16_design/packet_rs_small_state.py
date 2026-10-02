"""Fresh occupancy bounds and full replay for RS16 with state 17 or 18.

Each 128-bit group uses four GF16 RS[16,8] rows and independent uniform
GL(16,2) symbol maps. Independent uniform group-packet and regional shuffles
route four-bit packets through 64 regions. The selected16 quadratic-prefix
extension defines the t64/s inner, for s=17 or 18. Each physical update is
independent uniform GL(s,2). State starts at zero, persists across steps
and regions, and is discarded without a flush.

Q1/Q2 use exact expected outer shells and all ordered support placements.
The tail compares each active group's expected measure to beta times
uniform 256-bit output, retaining its active label even at artificial zero.
Exact regional placement and conditioned iid-marker bounds both preserve
chronological state transitions and count every terminal mass coordinate.

run_tail() covers only its explicit occupancy list. replay() freshly prepares
all 2^s states once, recomputes every component, and includes all q=1..K/128.
Only a source-stable replay with union below 2^-40 marks a certificate for
distance exceeding floor(2K/10), over the ideal setup law above. No saved
component, seeded guarantee, or S16 numerical result is accepted as a replay.
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

from flint import arb, ctx
import packet_inner_small_extension as inner
import packet_rs_length_whole as arithmetic
import packet_rs_state_sparse as low

q1, lengths, dense = arithmetic.q1, arithmetic.lengths, arithmetic.dense
TAIL_SCHEMA = 'rs16-small-state-tail-1'
WHOLE_SCHEMA = 'rs16-small-state-complete-replay-1'
CONTINUITY = arithmetic.CONTINUITY
DEFAULT_Q1 = (*arithmetic.retained.DEFAULT_Q1, '.0024')


def _prepare(bits, precision, data=None, map_record=None):
    inner._bits(bits)
    if (data is None) != (map_record is None):
        raise ValueError('shared prepared data and map record must be supplied together')
    before = q1.source_snapshot()
    ctx.prec = precision
    if data is None:
        data, map_record = inner.prepare(bits, birth_density='capped')
    inner.authenticate(data, map_record)
    if data['bits'] != bits or map_record['s'] != bits or data['birth_density'] != 'capped':
        raise ValueError('matching fresh capped state17/18 data required')
    sources = q1.source_snapshot()
    if any(sources.get(path) != digest for path, digest in before.items()):
        raise RuntimeError('preexisting mathematical source changed during preparation')
    return data, map_record, sources


def _stable(sources, precision):
    if ctx.prec != precision or q1.source_snapshot() != sources:
        raise RuntimeError('precision or mathematical source changed during computation')


def _common(geometry, bits, map_record, precision, sources):
    return dict(outer='rs16', K=geometry.K, N=geometry.N, geometry=asdict(geometry),
        groups=geometry.group_count, regions=64, group_dimension=128, group_output_bits=256,
        distance='1/10', threshold=geometry.N//10, physical_t=64, state_bits=bits,
        physical_steps=geometry.N//64, macro_t=128, macro_steps=geometry.N//128,
        macros_per_region=geometry.macros_per_region, physical_steps_per_macro=2,
        inner_distribution='uniform_gl', inner_group=f'GL({bits},2)',
        physical_updates_independent=True, outer_symbol_group='GL(16,2)',
        zero_initial_state=True, final_flush=False, state_continuity=CONTINUITY,
        terminal='sum of all mass-envelope coordinates; no flush',
        map_record=map_record, source_sha256=sources, precision=precision,
        fresh_computation=True, whole_code_certificate=False)


def _options(K, bits, occupancies, tilts, precision, output, method, markers):
    geometry = lengths.geometry(K)
    inner._bits(bits)
    tilts, output = lengths.checked_options(tilts, precision, output)
    if isinstance(occupancies, (str, bytes)):
        raise ValueError('explicit sorted occupancy sequence required')
    qs = tuple(occupancies)
    if (not qs or any(type(q) is not int or not 3 <= q <= geometry.group_count for q in qs)
            or qs != tuple(sorted(set(qs)))):
        raise ValueError('explicit increasing distinct tail occupancies in3..L required')
    if method not in ('exact', 'fugacity') or isinstance(markers, (str, bytes)):
        raise ValueError('exact/fugacity method and a marker sequence required')
    try:
        markers = tuple(map(Q, markers))
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('finite rational markers required') from error
    if (len(set(markers)) != len(markers) or any(not 0 < p <= 1 for p in markers) or
            (method == 'fugacity' and (not markers or
             (qs[0] < geometry.group_count and not any(p < 1 for p in markers))))):
        raise ValueError('distinct valid markers covering every requested occupancy required')
    return geometry, qs, tuple(map(Q, tilts)), markers, output


def _summary(values, precision):
    endpoint = arithmetic.rounded_endpoint(arithmetic.dyadic_sum(values), precision)
    return endpoint, arithmetic._margin(endpoint)


def _fugacity_affine(geometry, tilt, p, log_moment, log_beta):
    """Cache the q-independent terms of dense.coefficient_log_bound exactly."""
    base = log_moment + q1.kernel_t64.aq(tilt)*(geometry.N//10)
    if p == 1:
        return base + geometry.group_count*log_beta, None
    log_p, log_complement = q1.kernel_t64.aq(p).log(), q1.kernel_t64.aq(1-p).log()
    return (base - geometry.regions*geometry.group_count*log_complement,
            log_beta - geometry.regions*(log_p-log_complement))


def run_tail(*, K, bits, occupancies, tilts, precision=192, output=None,
             method='exact', marker_probabilities=(), data=None, map_record=None):
    """Recompute every listed tail occupancy with actual authenticated state-s data."""
    geometry, qs, tilts, markers, output = _options(
        K, bits, occupancies, tilts, precision, output, method, marker_probabilities)
    start = monotonic()
    beta, counts = lengths.exact_outer()
    data, map_record, sources = _prepare(bits, precision, data, map_record)
    best, choices = dict.fromkeys(qs), dict.fromkeys(qs)
    log_counts = {q: arb(comb(geometry.group_count, q)).log() for q in qs} if method == 'fugacity' else {}
    log_beta = q1.kernel_t64.aq(beta).log() if method == 'fugacity' else None
    record = dict(_common(geometry, bits, map_record, precision, sources), schema=TAIL_SCHEMA,
        method=method, return_denominator=(1 << bits)-1, birth_density='capped',
        occupancy_covered=list(qs), every_requested_occupancy_checked=True,
        beta=str(beta), every_shell_checked=True, count_kind='exact_expected_shells',
        count_sha256=lengths.count_hash(counts), tilts=list(map(str, tilts)),
        marker_probabilities=list(map(str, markers)), trials=[],
        scope=f'Expected low-weight message count for the listed occupancies under independent '
              f'ideal outer GL16, uniform routing, and inner GL({bits},2) setup. '
              'Continuous zero-initial state, no flush; not a whole-code certificate.')

    def consider(q, value, witness):
        endpoint = arithmetic.dyadic(lengths.positive_endpoint(value))
        if best[q] is None or arithmetic.dyadic_less(endpoint, best[q]):
            best[q], choices[q] = endpoint, witness

    for tilt in tilts:
        trial_start = monotonic()
        if method == 'exact':
            local = q1.kernel_t64.local_operators(data, tilt, activity=Q(1, 2))
            regional = q1.placement(local, epochs=geometry.macros_per_region,
                windows=32, rounding=q1.rounded, maximum_groups=qs[-1])
            for q in qs:
                consider(q, lengths.occupancy_upper(regional, K=K, occupancy=q, beta=beta, tilt=tilt),
                         dict(tilt=str(tilt)))
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
                log_moment = moment.log()
                base, slope = _fugacity_affine(geometry, tilt, p, log_moment, log_beta)
                for q in ((geometry.group_count,) if p == 1 else qs):
                    logarithm = base if p == 1 else base + (1-geometry.regions)*log_counts[q] + q*slope
                    consider(q, q1.kernel_t64.up(logarithm.exp()),
                             dict(tilt=str(tilt), marker_probability=str(p), selected_birth_rows=list(indices)))
        if any(value is None for value in best.values()):
            raise ArithmeticError('every requested occupancy must receive a bound')
        upper, margin = _summary(best.values(), precision)
        record['trials'].append(dict(tilt=str(tilt), union_upper=upper, margin_bits=margin,
            elapsed_seconds=monotonic()-trial_start))
        _stable(sources, precision)
        print(f'RS16/S{bits} K={K} {method} q={qs[0]}..{qs[-1]} '
              f'count={len(qs)} tilt={tilt} margin_bits={margin}', flush=True)
    inner.authenticate(data, map_record)
    _stable(sources, precision)
    record.update(union_upper=upper, margin_bits=margin,
        occupancy_uppers={str(q): list(best[q]) for q in qs},
        occupancy_choices={str(q): choices[q] for q in qs},
        elapsed_seconds=monotonic()-start)
    if output is not None:
        arithmetic._write_new(output, record)
    return record


def recipe(K, bits, **options):
    """Validate a complete coverage plan and its witnesses before any census."""
    inner._bits(bits)
    options.setdefault('q1_tilts', DEFAULT_Q1)
    return dict(arithmetic.recipe(K, **options), state_bits=bits)


def _component_endpoints(record, geometry, bits, map_record, sources, precision, count_hash):
    expected = _common(geometry, bits, map_record, precision, sources)
    if json.dumps({key: record.get(key) for key in expected}, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise ValueError('fresh component scope, geometry, map, sources, and precision must match')
    if (record.get('count_kind') != 'exact_expected_shells' or
            record.get('count_sha256') != count_hash or
            record.get('whole_code_certificate') is not False or
            record.get('fresh_computation') is not True):
        raise ValueError('exact outer shell hash and fresh noncertificate component required')
    schema = record.get('schema')
    if schema in ('finite-packet-rs-state-q1-1', 'finite-packet-rs-state-q2-1'):
        q = 1 if schema.endswith('q1-1') else 2
        if (record.get('occupancy_covered') != [q] or record.get('all_supports_covered') is not True or
                record.get('all_two_group_support_pairs_covered') is not (q == 2)):
            raise ValueError('complete q1/q2 support coverage required')
        return {q: arithmetic.dyadic(record.get(f'q{q}_upper'))}
    if (schema != TAIL_SCHEMA or record.get('method') not in ('exact', 'fugacity') or
            record.get('return_denominator') != (1 << bits)-1 or record.get('birth_density') != 'capped' or
            record.get('every_requested_occupancy_checked') is not True or
            record.get('every_shell_checked') is not True or record.get('beta') != str(lengths.exact_outer()[0])):
        raise ValueError('matching explicit exact/fugacity tail component required')
    qs = record.get('occupancy_covered')
    if (not isinstance(qs, list) or not qs or any(type(q) is not int or not 3 <= q <= geometry.group_count for q in qs)
            or qs != sorted(set(qs)) or not isinstance(record.get('occupancy_uppers'), dict) or
            set(record['occupancy_uppers']) != {str(q) for q in qs}):
        raise ValueError('exactly one endpoint per explicitly listed tail occupancy required')
    return {q: arithmetic.dyadic(record['occupancy_uppers'][str(q)]) for q in qs}


def replay(output, *, K, bits, **options):
    """Freshly recompute a complete union; there is no saved-component input."""
    start, plan = monotonic(), recipe(K, bits, **options)
    output, geometry = Path(output), lengths.geometry(K)
    names = ['q1', 'q2'] + (['exact'] if plan['exact_occupancies'] else []) + ['dense']
    paths = {name: output.with_name(output.stem+'-'+name+'.json') for name in names}
    if output.exists() or any(path.exists() for path in paths.values()):
        raise ValueError('whole replay and all component paths must be fresh')
    precision = plan['precision']
    data, map_record, sources = _prepare(bits, precision)
    shared = dict(data=data, map_record=map_record, precision=precision)
    records = []
    for q in (1, 2):
        inner.authenticate(data, map_record)
        records.append(low.run_low(q, geometry=geometry, tilts=plan['tilts'][f'q{q}'],
                                   output=paths[f'q{q}'], **shared))
        _stable(sources, precision)
    if plan['exact_occupancies']:
        records.append(run_tail(K=K, bits=bits, occupancies=plan['exact_occupancies'],
            tilts=plan['tilts']['exact'], output=paths['exact'], method='exact', **shared))
    records.append(run_tail(K=K, bits=bits,
        occupancies=range(plan['dense_interval'][0], geometry.group_count+1),
        tilts=plan['tilts']['dense'], marker_probabilities=plan['marker_probabilities'],
        output=paths['dense'], method='fugacity', **shared))
    _, counts = lengths.exact_outer()
    count_hash = lengths.count_hash(counts)
    values, selected = {}, {}
    for index, record in enumerate(records):
        endpoints = _component_endpoints(record, geometry, bits, map_record, sources, precision, count_hash)
        for q, value in endpoints.items():
            if q < 3 and q in values:
                raise ValueError('duplicate q1/q2 component')
            if q not in values or arithmetic.dyadic_less(value, values[q]):
                values[q], selected[q] = value, index
    if set(values) != set(range(1, geometry.group_count+1)):
        raise ValueError('fresh replay must cover every occupancy1..L')
    upper, margin = _summary(values.values(), precision)
    tail, tail_margin = _summary((values[q] for q in range(3, geometry.group_count+1)), precision)
    inner.authenticate(data, map_record)
    _stable(sources, precision)
    passes = arithmetic.dyadic_less(arithmetic.dyadic(upper), (1, -40))
    result = dict(_common(geometry, bits, map_record, precision, sources), schema=WHOLE_SCHEMA,
        target_margin_bits=40, target_minimum_distance=geometry.N//10+1,
        return_denominator=(1 << bits)-1, birth_density='capped', count_sha256=count_hash,
        occupancy_covered=[1, geometry.group_count], all_occupancies_covered=True,
        fresh_replay=True, target_met=passes, recipe=plan,
        occupancy_uppers={str(q): list(values[q]) for q in sorted(values)},
        selected_components={str(q): selected[q] for q in sorted(values)},
        q1_upper=list(values[1]), q2_upper=list(values[2]), tail_upper=tail, union_upper=upper,
        margin_bits=margin, tail_margin_bits=tail_margin,
        arithmetic='Exact sum of minimum dyadic endpoints; one final upward rounding.',
        components=[dict(path=str(paths[name].resolve()), schema=record['schema'],
            sha256=hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
            hash_scope='canonical record JSON', occupancy_covered=record['occupancy_covered'])
            for name, record in zip(names, records)], elapsed_seconds=monotonic()-start,
        scope=f'Full first-moment union under independent ideal outer GL16 symbol maps, '
              f'uniform group/regional routing, and inner GL({bits},2) updates; '
              'zero initial state, continuous state, no final flush. No seeded guarantee.')
    result['whole_code_certificate'] = passes
    arithmetic._write_new(output, result)
    print(f'RS16/S{bits} K={K} COMPLETE UNION margin_bits={margin} certificate={passes}', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--K', type=int, required=True)
    parser.add_argument('--bits', type=int, choices=(17, 18), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--exact-occupancies', type=int, nargs='*', default=list(arithmetic.DEFAULT_EXACT))
    parser.add_argument('--dense-min', type=int, default=3)
    for name, default in (('q1', DEFAULT_Q1), ('q2', arithmetic.retained.DEFAULT_Q2),
                          ('exact', arithmetic.retained.DEFAULT_SPARSE), ('dense', arithmetic.retained.DEFAULT_DENSE)):
        parser.add_argument(f'--{name}-tilts', nargs='+', default=list(default))
    parser.add_argument('--marker-probabilities', nargs='+', default=list(dense.DEFAULT_MARKER_PROBABILITIES))
    parser.add_argument('--dry-run', action='store_true')
    args = vars(parser.parse_args())
    output, dry = args.pop('output'), args.pop('dry_run')
    if dry:
        print(json.dumps(recipe(**args), indent=2))
    else:
        replay(output, **args)


if __name__ == '__main__':
    main()
