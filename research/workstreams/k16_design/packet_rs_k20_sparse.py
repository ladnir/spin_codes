"""Fresh sparse-occupancy bounds for the K=2^20 RS16 packet ensemble.

There are L=8192 independent outer groups. Each group encodes 128 message
bits into 64 four-bit packets using four parallel GF16 RS[16,8] words.
Independent uniform GL(16,2) maps randomize its sixteen symbol positions.
An independent uniform permutation shuffles each group's packets. In each
of 64 regions, another independent uniform permutation routes the L packets.
The selected t64/s16 inner starts at zero and retains state across regions.
Its physical updates use independent uniform GL(16,2) maps; no flush follows.

This module bounds the expected number of low-weight outputs with q active
outer groups. Q1 and Q2 reuse the unchanged exact support-placement routines.
Their group counts are expected shell counts, not cumulative caps or counts
for every realized outer code. Independence between different messages is
unnecessary; group setups must be independent.

For q>=3, let beta=2^256/(2^16-1)^8. Each active group's expected output
measure is at most beta times the uniform measure on all 256-bit vectors.
The comparison includes an artificial zero vector without removing its
active-group label. Its packet activities are independent Bernoulli(15/16).
Thus a region has J~Bin(q,15/16) active packets in distinct uniform slots.
Let R_j average the local transition upper envelopes over those placements,
and let R(q)=E[R_J]. A region contains 256 ordered 128-bit macro steps.
For tilt lambda>0, the occupancy contribution is at most

    C(8192,q) beta^q exp(lambda*209715) e_zero R(q)^64 1.

The matrix power retains the inner state; every terminal coordinate is mass.
Minimizing valid upper bounds over tilts is valid separately for each q.
The returned component covers only its explicit occupancy interval. It is
not a whole-code certificate, a claim about every setup, or a seed guarantee.
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
import packet_q2 as q2
import packet_uniform_tail as uniform


GEOMETRY = q1.Geometry(8192, 64, 128)
THRESHOLD = GEOMETRY.N // 10
DEFAULT_TILTS = ('.00016', '.00032', '.00064', '.00128', '.00256', '.00512')
TAIL_SCHEMA = 'k20-rs16-uniform-regional-tail-1'


def exact_outer():
    """Return exact beta and the 65 expected RS16 shells, excluding zero."""
    return uniform.uniform_envelope(n=16, k=8, packets_per_symbol=4)


def count_hash(counts, *, cumulative=False):
    if cumulative:
        total, values = Q(0), []
        for count in counts:
            total += count
            values.append(total)
        counts = values
    return hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest()


def checked_options(tilts, precision, output):
    """Validate exact positive tilts and reserve no output or numerical state."""
    if type(precision) is not int or precision < 128:
        raise ValueError('integer precision >=128 required')
    if isinstance(tilts, (str, bytes)):
        raise ValueError('a nonempty sequence of distinct positive rational tilts is required')
    try:
        values = tuple(Q(t) for t in tilts)
    except (TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('finite rational tilts required') from error
    if not values or any(t <= 0 for t in values) or len(set(values)) != len(values):
        raise ValueError('distinct positive rational tilts required')
    path = Path(output) if output is not None else None
    if path is not None and path.exists():
        raise ValueError('output must be fresh')
    return tuple(map(str, values)), path


def positive_endpoint(value):
    """Serialize a positive finite outward dyadic upper endpoint."""
    if not value.is_finite() or not value > 0:
        raise ArithmeticError('positive finite component endpoint required')
    upper = q1.kernel_t64.up(value)
    if not upper.is_finite() or not upper > 0:
        raise ArithmeticError('positive finite rounded endpoint required')
    return q1.endpoint(upper)


def endpoint_value(pair):
    if (not isinstance(pair, list) or len(pair) != 2 or
            any(type(v) is not int for v in pair) or pair[0] <= 0 or abs(pair[1]) > 1000000):
        raise ValueError('positive bounded-exponent dyadic endpoint required')
    return Q(pair[0]) * Q(2) ** pair[1]


def _premises():
    return dict(outer='four parallel GF16 RS[16,8] words',
        label_mixing='independent uniform GL16 per aligned four-nibble symbol',
        independent_setups_between_groups=True, exact_expected_shell_counts=True,
        independent_uniform_group_packet_shuffles=True,
        independent_uniform_regional_shuffles=True)


def _prepare(precision, data=None, map_record=None):
    ctx.prec = precision
    if (data is None) != (map_record is None):
        raise ValueError('shared prepared data and map record must be supplied together')
    if data is None:
        data, map_record = q1.kernel_t64.prepare(birth_density='capped')
    if not map_record:
        raise ArithmeticError('freshly enumerated physical-map record required')
    return data, map_record, q1.source_snapshot()


def _finish(record, sources, output):
    if sources != q1.source_snapshot():
        raise RuntimeError('loaded mathematical source changed during K20 computation')
    if record.get('whole_code_certificate') is not False:
        raise ArithmeticError('a sparse component cannot claim a whole-code certificate')
    record['source_sha256'] = sources
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        # Refuse replacement even if another process created the path meanwhile.
        with output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(record, indent=2) + '\n')
    return record


def _run_low(occupancy, *, tilts, precision, output, data=None, map_record=None):
    tilts, output = checked_options(tilts, precision, output)
    _, counts = exact_outer()
    data, map_record, sources = _prepare(precision, data, map_record)
    kwargs = dict(tilts=tilts, precision=precision, data=data,
                  map_record=map_record, metadata=_premises(), output=None)
    if occupancy == 1:
        record = q1.evaluate_q1(counts, group_count=GEOMETRY.group_count,
            regions=GEOMETRY.regions, epochs_per_region=GEOMETRY.macros_per_region,
            group_dimension=GEOMETRY.group_dimension, count_kind='shells', **kwargs)
    else:
        record = q2.evaluate_q2(counts, geometry=GEOMETRY, **kwargs)
    if (record.get('geometry') != asdict(GEOMETRY) or record.get('K') != GEOMETRY.K
            or record.get('N') != GEOMETRY.N or record.get('threshold') != THRESHOLD
            or record.get('distance') != '1/10' or record.get('map_record') != map_record
            or record.get('zero_initial_state') is not True or record.get('final_flush') is not False
            or record.get('occupancy_covered') != [occupancy]
            or record.get('count_sha256') != count_hash(counts, cumulative=occupancy == 1)
            or ctx.prec != precision):
        raise ArithmeticError('inconsistent K20 sparse evaluator receipt')
    endpoint_value(record[f'q{occupancy}_upper'])
    if record.get('source_sha256') != sources:
        raise ArithmeticError('underlying evaluator source pins must remain unchanged')
    record.update(outer='rs16', k20_sparse_component=f'q{occupancy}',
                  fresh_computation=True, macros_per_region=GEOMETRY.macros_per_region)
    return _finish(record, sources, output)


def run_q1(*, tilts=DEFAULT_TILTS, precision=192, output=None, data=None, map_record=None):
    """Freshly bound q1; optional data/map_record must come from current prepare()."""
    return _run_low(1, tilts=tilts, precision=precision, output=output,
                    data=data, map_record=map_record)


def run_q2(*, tilts=DEFAULT_TILTS, precision=192, output=None, data=None, map_record=None):
    """Freshly bound q2; optional data/map_record must come from current prepare()."""
    return _run_low(2, tilts=tilts, precision=precision, output=output,
                    data=data, map_record=map_record)


def occupancy_upper(regional, *, occupancy, beta, tilt):
    """Evaluate the stated K20 uniform-envelope expression with outward arithmetic."""
    if type(occupancy) is not int or not 3 <= occupancy <= GEOMETRY.group_count:
        raise ValueError('occupancy must lie in 3..8192')
    if Q(beta) <= 0 or Q(tilt) <= 0:
        raise ValueError('positive beta and tilt required')
    matrix = uniform.regional_uniform(regional, occupancy) ** GEOMETRY.regions
    moment = sum((matrix[0, j] for j in range(matrix.ncols())), arb(0))
    factor = (q1.kernel_t64.aq(Q(tilt)) * THRESHOLD).exp()
    upper = q1.kernel_t64.up(comb(GEOMETRY.group_count, occupancy) *
                            q1.kernel_t64.aq(Q(beta)) ** occupancy * factor * moment)
    positive_endpoint(upper)
    return upper


def run_tail(*, q_min=3, q_max=32, tilts=DEFAULT_TILTS, precision=192, output=None,
             data=None, map_record=None):
    """Freshly bound q_min..q_max; shared data/map_record must come from prepare()."""
    tilts, output = checked_options(tilts, precision, output)
    if (type(q_min) is not int or type(q_max) is not int or
            not 3 <= q_min <= q_max <= GEOMETRY.group_count):
        raise ValueError('explicit occupancy subinterval of 3..8192 required')
    start = monotonic()
    beta, counts = exact_outer()
    data, map_record, sources = _prepare(precision, data, map_record)
    best = dict.fromkeys(range(q_min, q_max + 1))
    choices = dict.fromkeys(best)
    record = dict(schema=TAIL_SCHEMA, outer='rs16', K=GEOMETRY.K, N=GEOMETRY.N,
        geometry=asdict(GEOMETRY), threshold=THRESHOLD, distance='1/10',
        q_min=q_min, q_max=q_max, occupancy_covered=[q_min, q_max],
        whole_code_certificate=False, fresh_computation=True,
        physical_t=64, state_bits=16, macro_t=128, physical_steps_per_macro=2,
        macros_per_region=GEOMETRY.macros_per_region,
        zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        regional_placement='exact without-replacement average of local upper envelopes',
        outer_comparison='beta times uniform full 256-bit group output; artificial zero retains active label',
        beta=str(beta), every_shell_checked=True, count_sha256=count_hash(counts),
        count_premises=_premises(), precision=precision, tilts=list(tilts),
        source_sha256=sources, map_record=map_record, trials=[],
        scope='Expected-message first moment over ideal independent uniform GL16 '
              'outer symbol maps, group packet shuffles, regional shuffles, and '
              'physical GL16 updates. Only the explicit occupancy interval is '
              'covered; no whole-code or deterministic-seed certificate.')
    for tilt in tilts:
        trial_start = monotonic()
        # This 1/2 selects the valid local density envelope. It does not replace
        # the comparison's Bernoulli(15/16) packet activity or nonzero label law.
        local = q1.kernel_t64.local_operators(data, Q(tilt), activity=Q(1, 2))
        regional = q1.placement(local, epochs=GEOMETRY.macros_per_region,
            windows=GEOMETRY.macro_windows, rounding=q1.rounded, maximum_groups=q_max)
        for q in best:
            upper = occupancy_upper(regional, occupancy=q, beta=beta, tilt=tilt)
            if best[q] is None or upper < best[q]:
                best[q], choices[q] = upper, tilt
        total = q1.kernel_t64.up(sum(best.values(), arb(0)))
        if ctx.prec != precision:
            raise ArithmeticError('precision changed during K20 computation')
        margin = str(-total.log() / arb(2).log())
        record['trials'].append(dict(tilt=tilt, union_upper=positive_endpoint(total),
            margin_bits=margin, elapsed_seconds=monotonic() - trial_start))
        record.update(union_upper=positive_endpoint(total), margin_bits=margin,
            occupancy_uppers={str(q): positive_endpoint(v) for q, v in best.items()},
            occupancy_choices={str(q): v for q, v in choices.items()},
            elapsed_seconds=monotonic() - start)
        if sources != q1.source_snapshot():
            raise RuntimeError('loaded mathematical source changed during K20 tail computation')
        print(f'K20 RS16 q={q_min}..{q_max} tilt={tilt} margin_bits={margin}', flush=True)
    return _finish(record, sources, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('component', choices=('q1', 'q2', 'tail'))
    parser.add_argument('--q-min', type=int, default=3)
    parser.add_argument('--q-max', type=int, default=32)
    parser.add_argument('--tilts', nargs='+', default=list(DEFAULT_TILTS))
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    kwargs = dict(tilts=args.tilts, precision=args.precision, output=args.output)
    if args.component == 'tail':
        run_tail(q_min=args.q_min, q_max=args.q_max, **kwargs)
    elif args.component == 'q1':
        run_q1(**kwargs)
    else:
        run_q2(**kwargs)


if __name__ == '__main__':
    main()
