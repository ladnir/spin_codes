"""Fresh occupancy bounds for RS16 with an isolated 20-bit recursive state.

The outer remains four parallel GF16 RS[16,8] words per128-bit group, with
independent uniform GL16 symbol mixers. Group shuffles and64 regional
shuffles remain independent and uniform. Only the inner changes: its fixed
t64 maps use20 state coordinates and each physical update is uniform GL20.
State starts at zero, persists across physical steps and regions, and is not
flushed. All probabilities concern this ideal independently sampled setup.

For q active groups the outer expected measure is dominated by beta^q times
uniform256-bit group outputs. The exact method averages all distinct regional
packet placements. The fugacity method removes conditioning of iid markers
on exactlyq markers perregion, paying the inverse conditioning probability
in each of64 regions. Both keep chronological matrix products and all final
state mass. Numerical endpoints cover only the explicitly listed occupancies;
no result from this module alone is a whole-code certificate.
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
import packet_q1 as q1
import packet_uniform_tail as uniform
import packet_rs_k20_dense as density
import packet_inner_s20 as inner


GEOMETRY = q1.Geometry(8192, 64, 128)


def authenticate(data, map_record):
    physical = inner.authenticate(data)
    if (data['bits'] != 20 or data['physical_step_bits'] != 64 or
            data['macro_windows'] != 32 or physical['distribution'] != 'uniform_gl' or
            map_record.get('s') != 20 or map_record.get('t') != 64 or
            map_record.get('map_sha256') != data['map_sha256'] or
            map_record.get('distribution') != 'uniform_gl'):
        raise ValueError('authenticated t64/s20 map record and uniform GL20 wrapper required')
    source = map_record.get('source')
    if (not isinstance(source, dict) or not source.get('path') or
            hashlib.sha256(Path(source['path']).read_bytes()).hexdigest() != source.get('sha256')):
        raise ValueError('current base-map source must match the construction record')
    inner._check_pins(map_record)
    if (map_record.get('full_state_census') is not True or
            map_record.get('state_count') != 1 << 20):
        raise ValueError('fresh complete twenty-bit state census required')
    return physical


def options(occupancies, tilts, precision, output, method, markers):
    if type(precision) is not int or precision < 128:
        raise ValueError('integer precision >=128 required')
    qs = tuple(occupancies)
    if (not qs or any(type(q) is not int or not 3 <= q <= GEOMETRY.group_count for q in qs)
            or tuple(sorted(set(qs))) != qs):
        raise ValueError('explicit increasing distinct occupancies in3..8192 required')
    if isinstance(tilts, (str, bytes)):
        raise ValueError('distinct positive rational tilt sequence required')
    tilts = tuple(map(Q, tilts))
    if not tilts or any(t <= 0 for t in tilts) or len(set(tilts)) != len(tilts):
        raise ValueError('distinct positive rational tilts required')
    if method not in ('exact', 'fugacity'):
        raise ValueError('explicit supported proof method required')
    if isinstance(markers, (str, bytes)):
        raise ValueError('marker probabilities must be a sequence')
    markers = tuple(map(Q, markers))
    if (method == 'fugacity' and (not markers or len(set(markers)) != len(markers)
            or any(not 0 < p <= 1 for p in markers)
            or (qs[0] < GEOMETRY.group_count and not any(p < 1 for p in markers)))):
        raise ValueError('valid marker grid covering every requested occupancy required')
    path = None if output is None else Path(output)
    if path is not None and path.exists():
        raise ValueError('fresh output required')
    return qs, tilts, markers, path


def run_tail(*, occupancies, tilts, precision=192, output=None, method='exact',
             marker_probabilities=(), data=None, map_record=None, birth_density='classes'):
    """Compute new bounds, never substitute stored endpoints or state-size labels."""
    qs, tilts, markers, output = options(occupancies, tilts, precision, output,
                                       method, marker_probabilities)
    if (data is None) != (map_record is None):
        raise ValueError('prepared data and map record must be supplied together')
    start = monotonic()
    ctx.prec = precision
    if data is None:
        data, map_record = inner.prepare(birth_density=birth_density)
    authenticate(data, map_record)
    beta, counts = uniform.uniform_envelope(16, 8, 4)
    sources = q1.source_snapshot()
    threshold = GEOMETRY.N // 10
    best, choices = dict.fromkeys(qs), dict.fromkeys(qs)
    log_counts = {q: arb(comb(GEOMETRY.group_count, q)).log() for q in qs}
    contiguous = qs == tuple(range(qs[0], qs[-1] + 1))
    record = dict(schema='rs16-s20-occupancy-component-1', outer='rs16',
        K=GEOMETRY.K, N=GEOMETRY.N, geometry=asdict(GEOMETRY), state_bits=20,
        physical_t=64, macro_t=128, distance='1/10', threshold=threshold,
        outer_symbol_update='independent uniform GL16',
        inner_state_update='independent uniform GL20', method=method,
        return_denominator=(1 << 20) - 1, birth_density=map_record['birth_density'],
        occupancy_values=list(qs), contiguous_occupancy_interval=contiguous,
        occupancy_covered=[qs[0], qs[-1]] if contiguous else None,
        q_min=qs[0], q_max=qs[-1], every_requested_occupancy_checked=True,
        zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        beta=str(beta), count_kind='exact_expected_shells',
        count_sha256=hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest(),
        map_record=map_record, source_sha256=sources, precision=precision,
        tilts=list(map(str, tilts)), marker_probabilities=list(map(str, markers)),
        whole_code_certificate=False, fresh_computation=True, trials=[],
        scope='Expected low-weight message count for only the listed occupancies, '
              'under the independent uniform RS16 outer/routing and GL20 inner setup. '
              'Continuous state, no flush. Not a whole-code certificate.')
    for tilt in tilts:
        trial_start = monotonic()
        local = q1.kernel_t64.local_operators(data, tilt, activity=Q(1, 2))
        if method == 'exact':
            regional = q1.placement(local, epochs=GEOMETRY.macros_per_region,
                windows=32, rounding=q1.rounded, maximum_groups=qs[-1])
            candidates = []
            for q in qs:
                power = uniform.regional_uniform(regional, q) ** GEOMETRY.regions
                moment = sum((power[0, j] for j in range(power.ncols())), arb(0))
                if not moment.is_finite() or not moment > 0:
                    raise ArithmeticError('positive finite regional moment required')
                log_bound = (log_counts[q] + q * q1.kernel_t64.aq(beta).log() +
                             q1.kernel_t64.aq(tilt) * threshold + moment.log())
                candidates.append((q, log_bound, dict(tilt=str(tilt))))
        else:
            candidates = []
            for p in markers:
                if p == 1 and GEOMETRY.group_count not in best:
                    continue
                mixed = density.iid_macro(local, Q(15, 16) * p)
                power = mixed ** (GEOMETRY.N // 128)
                moment = sum((power[0, j] for j in range(power.ncols())), arb(0))
                if not moment.is_finite() or not moment > 0:
                    raise ArithmeticError('positive finite iid moment required')
                log_moment = moment.log()
                for q in ((GEOMETRY.group_count,) if p == 1 else qs):
                    log_bound = density.coefficient_log_bound(finite_geometry=GEOMETRY,
                        q=q, tilt=tilt, beta=beta, marker_probability=p,
                        log_moment=log_moment, log_binomial=log_counts[q])
                    candidates.append((q, log_bound, dict(tilt=str(tilt), marker_probability=str(p))))
        for q, log_bound, witness in candidates:
            upper = q1.kernel_t64.up(log_bound.exp())
            if not upper.is_finite() or not upper > 0:
                raise ArithmeticError('positive finite outward occupancy endpoint required')
            if best[q] is None or upper < best[q]:
                best[q], choices[q] = upper, witness
        if any(v is None for v in best.values()):
            raise ArithmeticError('missing requested occupancy')
        total = q1.kernel_t64.up(sum(best.values(), arb(0)))
        worst = sorted(best, key=lambda q: float(best[q].log()), reverse=True)[:8]
        margin = str(-total.log() / arb(2).log())
        record['trials'].append(dict(tilt=str(tilt), union_upper=q1.endpoint(total),
            margin_bits=margin, worst_q=worst, elapsed_seconds=monotonic() - trial_start))
        if ctx.prec != precision or sources != q1.source_snapshot():
            raise RuntimeError('precision or loaded mathematical source changed during computation')
        print(f'RS16/S20 {method} q={qs[0]}..{qs[-1]} count={len(qs)} '
              f'tilt={tilt} margin_bits={margin} worst_q={worst}', flush=True)
    record.update(union_upper=q1.endpoint(total), margin_bits=margin,
        occupancy_uppers={str(q): q1.endpoint(v) for q, v in best.items()},
        occupancy_choices={str(q): v for q, v in choices.items()},
        elapsed_seconds=monotonic() - start)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(record, indent=2) + '\n')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q-min', type=int, default=3)
    parser.add_argument('--q-max', type=int, default=32)
    parser.add_argument('--points', type=int, nargs='+')
    parser.add_argument('--method', choices=('exact', 'fugacity'), default='exact')
    parser.add_argument('--tilts', nargs='+', default=['.0072', '.008', '.0096'])
    parser.add_argument('--marker-probabilities', nargs='+', default=list(density.DEFAULT_MARKER_PROBABILITIES))
    parser.add_argument('--birth-density', choices=('classes', 'capped'), default='classes')
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    occupancies = args.points if args.points is not None else range(args.q_min, args.q_max + 1)
    run_tail(occupancies=occupancies, tilts=args.tilts, precision=args.precision,
        output=args.output, method=args.method, marker_probabilities=args.marker_probabilities,
        birth_density=args.birth_density)


if __name__ == '__main__':
    main()
