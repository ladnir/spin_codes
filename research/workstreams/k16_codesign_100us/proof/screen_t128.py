"""Fresh bounded screens for two explicit t128/s16 candidates at K=65536.

The outer is the retained four-row GF16 RS[16,8] construction. Each symbol
has its own independent uniform-nonzero-image randomizer. Routing uses
independent within-group packet shuffles and independent regional shuffles.
Each physical inner step emits x+A*a and updates a to M*a+A^T*x. The M's
are independent uniform GL16 maps (or any independent exact transitive
linear family). Initial state is zero, state persists, and there is no flush.

prefix16 uses the first sixteen rows of the pinned t128/s20 map. repeat64
duplicates each selected t64/s16 expansion row into both 64-bit halves.
Both are new constructions, evaluated as ONE physical step per 128 bits.
No prior numerical endpoint is reused. Outputs are partial bound receipts,
never whole-code certificates, even if all requested occupancies pass.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
DESIGN = HERE.parents[1] / 'k16_design'
sys.path.insert(0, str(DESIGN))

from flint import arb, ctx
import packet_q1 as q1
import packet_q2 as q2
import packet_uniform_tail as tail
import packet_inner_t128_extension as t128

aq, up = q1.kernel_t64.aq, q1.kernel_t64.up


def prepare(candidate):
    if candidate == 'prefix16':
        source = t128.BASE_MAP
        if hashlib.sha256(source.read_bytes()).hexdigest() != t128.BASE_SHA256:
            raise ArithmeticError('pinned t128/s20 map changed')
        rows = tuple(int(x, 16) for x in json.loads(source.read_bytes())['generator_rows_hex'][:16])
    elif candidate == 'repeat64':
        source = q1.kernel_t64.SELECTED_MAP.resolve()
        _, _, base = q1.kernel_t64._source_maps(source)
        rows = tuple(int(x, 16) | (int(x, 16) << 64) for x in base['expansion_rows_hex'])
    else:
        raise ValueError('explicit prefix16 or repeat64 candidate required')
    columns = t128.columns_from_rows(rows)
    rank = q1.kernel_t64.s16_maps.binary_rank
    if (len(rows) != 16 or rank(rows) != 16 or rank(columns) != 16
            or any((a & b).bit_count() & 1 for a in rows for b in rows)
            or any(rank(columns[i:i+4]) != 4 for i in range(0, 128, 4))):
        raise ArithmeticError('rank, CA=0, or full-rank-packet check failed')
    images = tuple(q1.kernel_t64.s16_maps.images_from_rows(rows))
    spectrum = dict(sorted(Counter(x.bit_count() for x in images).items()))
    if len(set(images)) != 1 << 16 or spectrum.get(0) != 1:
        raise ArithmeticError('full injective state census failed')
    data = q1.kernel_t64.kernel_maps.prepare_maps(images, columns, bits=16,
        distribution='uniform_gl', birth_density='capped')
    record = dict(candidate=candidate, physical_t=128, state_bits=16,
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=columns,
        feedback_definition='C=A^T', expansion_rank=16, feedback_rank=16,
        feedback_times_expansion_zero=True, full_state_census=True,
        minimum_expansion_weight=min(w for w in spectrum if w),
        expansion_spectrum={str(w): n for w, n in spectrum.items()},
        map_sha256=data['map_sha256'], independent_transitive_update_per_step=True,
        source=dict(path=str(source), sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    return data, record


def source_pins(map_record):
    result = q1.source_snapshot()
    result[str(HERE / 'screen_t128.py')] = hashlib.sha256((HERE / 'screen_t128.py').read_bytes()).hexdigest()
    result[map_record['source']['path']] = map_record['source']['sha256']
    return result


def local_operators(data, tilt):
    q1.kernel_t64.kernel_maps.authenticate(data)
    if data['bits'] != 16 or data['windows'] != 32:
        raise ValueError('one physical t128/s16 step required')
    z = (-aq(tilt)).exp()
    local = q1.kernel_t64.sparse_kernel.outward_at_z(data, z)
    return q1.kernel_t64.kernel_birth_density.refine_local(data, local, z, Q(1, 2))


def run(candidate, tilts, max_q, precision, output):
    tilts = tuple(map(Q, tilts))
    if (not tilts or min(tilts) <= 0 or len(set(tilts)) != len(tilts)
            or not 1 <= max_q <= 512 or precision < 128 or output.exists()):
        raise ValueError('positive distinct tilts, q range, precision, and fresh output required')
    ctx.prec = precision
    start = monotonic()
    data, maps = prepare(candidate)
    beta, counts = tail.uniform_envelope(16, 8, 4)
    pins = source_pins(maps)
    best_one = [arb(1) for _ in range(65)]
    best_two = [[arb(1) for _ in range(v + 1)] for v in range(65)] if max_q >= 2 else None
    best_tail = {q: None for q in range(3, max_q + 1)}
    choices = {}
    trials = []
    print(f'Preparing {candidate}: minimum expansion weight={maps["minimum_expansion_weight"]}', flush=True)
    for tilt in tilts:
        local = local_operators(data, tilt)
        regional = q1.placement(local, epochs=16, windows=32,
            rounding=q1.rounded, maximum_groups=max_q)
        factor = (aq(tilt) * 13107).exp()
        one_moments = q1.support_moments(regional[0], regional[1], 64)
        for v, moment in enumerate(one_moments):
            best_one[v] = min(best_one[v], up(factor * moment))
        one = up(512 * sum((aq(c) * m for c, m in zip(counts, best_one)), arb(0)))
        values = {1: one}
        if max_q >= 2:
            moments = q2.pair_support_moments(regional[:3], 64)
            for v, column in enumerate(moments):
                for u, moment in enumerate(column):
                    best_two[v][u] = min(best_two[v][u], up(factor * moment))
            values[2] = q2.fold_shell_pairs(counts, best_two, 512)
        for q in best_tail:
            matrix = tail.regional_uniform(regional, q) ** 64
            moment = sum((matrix[0, j] for j in range(matrix.ncols())), arb(0))
            value = up(comb(512, q) * aq(beta) ** q * factor * moment)
            if best_tail[q] is None or value < best_tail[q]:
                best_tail[q], choices[q] = value, str(tilt)
        values.update(best_tail)
        if any(not v.is_finite() or not v > 0 for v in values.values()):
            raise ArithmeticError('finite positive outward endpoints required')
        total = up(sum(values.values(), arb(0)))
        margins = {str(q): str(-v.log() / arb(2).log()) for q, v in values.items()}
        trials.append(dict(tilt=str(tilt), occupancy_margin_bits=margins))
        if source_pins(maps) != pins:
            raise ArithmeticError('mathematical source changed during screen')
        result = dict(schema='k16-t128-s16-bounded-screen-1', K=65536, N=131072,
            threshold=13107, distance='1/10', groups=512, regions=64,
            group_dimension=128, physical_t=128, state_bits=16, physical_steps=1024,
            physical_steps_per_region=16, zero_initial_state=True, final_flush=False,
            state_continuity='retained_between_every_step_and_region',
            outer='four parallel GF16 RS[16,8] rows with independent transitive GL16 symbol maps',
            map_record=maps, source_sha256=pins, precision=precision,
            tilts=list(map(str, tilts)), trials=trials, fresh_computation=True,
            occupancy_covered=[1, max_q], whole_code_certificate=False,
            occupancy_uppers={str(q): q1.endpoint(v) for q,v in values.items()},
            occupancy_margin_bits=margins, tail_choices=choices,
            union_upper=q1.endpoint(total), margin_bits=str(-total.log() / arb(2).log()),
            elapsed_seconds=monotonic() - start,
            scope='Fresh partial first-moment bound for the explicit ideal ensemble; no whole-code certificate or timing claim.')
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2) + '\n')
        print(f'tilt={tilt} q1={margins["1"]} union={result["margin_bits"]}', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', choices=('prefix16', 'repeat64'), required=True)
    parser.add_argument('--tilts', nargs='+', default=['.00256', '.00512', '.01024', '.02048'])
    parser.add_argument('--max-q', type=int, default=1)
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.candidate, args.tilts, args.max_q, args.precision, args.output)


if __name__ == '__main__':
    main()
