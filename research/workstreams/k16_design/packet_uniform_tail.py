"""RS occupancy screen using a pointwise uniform outer comparison.

For each active group, the expected measure on complete 256-bit outputs is
dominated by beta times the uniform measure. This is a comparison measure,
not the setup distribution. Its artificial zero output keeps the active label.
The exact regional average retains state across all regions. q=1 and q=2 are
excluded even when q=3..512 is covered; no whole-code certificate is claimed.
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
import packet_q1 as q1
import rs_outer


def uniform_envelope(n=8, k=4, packets_per_symbol=8):
    """Check every shell of beta times uniform bits against exact RS counts.

For h active symbols, the shortened MDS code has dimension h-(n-k).
An information set injects its exact-support words into nonzero tuples, so
there are at most (Q-1)^(h-(n-k)) such words per specified support. Independent
uniform GL maps therefore give each nonzero n-symbol vector mass at most
(Q-1)^-(n-k). The bound also holds at the zero vector after excluding the
zero message. Shuffling preserves the uniform majorant.
"""
    if n * packets_per_symbol != 64 or n != 2 * k:
        raise ValueError('64 four-bit packets and half-rate RS geometry required')
    alphabet = 1 << (4 * packets_per_symbol)
    beta = Q(1 << 256, (alphabet - 1)**(n - k))
    counts = rs_outer.expected_group_support_counts(n=n, k=k, packet_bits=4,
                                                    packets_per_symbol=packets_per_symbol)
    for u, count in enumerate(counts):
        if count > beta * Q(comb(64, u) * 15**u, 16**64):
            raise ArithmeticError(f'uniform shell domination failed at {u}')
    return beta, counts


def regional_uniform(regional, active_groups, *, matrix=arb_mat,
                     rational=q1.kernel_t64.aq, rounding=q1.rounded):
    """Average R_j with J~Bin(active_groups,15/16), preserving matrix order."""
    if type(active_groups) is not int or not 0 <= active_groups < len(regional):
        raise ValueError('complete regional occupancy operators required')
    size = regional[0].nrows()
    value = matrix(size, size)
    total = Q(0)
    for j in range(active_groups + 1):
        probability = Q(comb(active_groups, j) * 15**j, 16**active_groups)
        total += probability
        value += regional[j] * rational(probability)
    if total != 1:
        raise ArithmeticError('exact binomial occupancy normalization failed')
    return rounding(value)


def run(*, outer='rs8', q_min=3, q_max=32, precision=192,
        tilts=('.00512', '.01024', '.0256', '.0512', '.1024', '.2048'), output=None):
    if (outer not in ('rs8', 'rs16') or type(q_min) is not int or type(q_max) is not int
            or not 3 <= q_min <= q_max <= 512 or precision < 128 or not tilts
            or any(Q(t) <= 0 for t in tilts) or len(set(map(Q, tilts))) != len(tilts)
            or (output is not None and output.exists())):
        raise ValueError('valid outer, q=3..512 subrange, positive tilts, precision, and fresh output required')
    start = monotonic()
    ctx.prec = precision
    n, k, packets = (8, 4, 8) if outer == 'rs8' else (16, 8, 4)
    beta, counts = uniform_envelope(n, k, packets)
    geometry = q1.Geometry(512, 64, 128)
    threshold = geometry.N // 10
    data, map_record = q1.kernel_t64.prepare(birth_density='capped')
    sources = q1.source_snapshot()
    best = {q: None for q in range(q_min, q_max + 1)}
    choices = dict.fromkeys(best)
    record = dict(schema='rs-uniform-regional-tail-diagnostic-1', outer=outer,
        K=geometry.K, N=geometry.N, geometry=asdict(geometry), threshold=threshold,
        distance='1/10', q_min=q_min, q_max=q_max, occupancy_covered=[q_min, q_max],
        whole_code_certificate=False, tail_cover_complete=q_min == 3 and q_max == 512,
        zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        regional_placement='exact without-replacement average of local upper envelopes',
        outer_comparison='beta times uniform full 256-bit group output, artificial zero keeps active label',
        beta=str(beta), every_shell_checked=True,
        count_sha256=hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest(),
        precision=precision, tilts=list(tilts), source_sha256=sources, map_record=map_record, trials=[],
        scope='RS expected outer measure over independent GL mixing, independent uniform '
              'column and regional shuffles, selected t64/s16 maps with fresh ideal uniform '
              'GL16 updates. Exact conditional placement of transition upper envelopes. '
              'No seed-specific or whole-code certification; q=1 and q=2 excluded.')
    print(f'UNIFORM RS TAIL {outer} q={q_min}..{q_max} K={geometry.K} N={geometry.N}', flush=True)
    for tilt in tilts:
        trial_start = monotonic()
        local = q1.kernel_t64.local_operators(data, Q(tilt), activity=Q(1, 2))
        print(f'tilt={tilt} building regional operators through q={q_max}', flush=True)
        regional = q1.placement(local, epochs=geometry.macros_per_region, windows=32,
                                 rounding=q1.rounded, maximum_groups=q_max)
        factor = (q1.kernel_t64.aq(Q(tilt)) * threshold).exp()
        for q in best:
            matrix = regional_uniform(regional, q)**geometry.regions
            moment = sum((matrix[0, j] for j in range(matrix.ncols())), arb(0))
            upper = q1.kernel_t64.up(comb(geometry.group_count, q) *
                                    q1.kernel_t64.aq(beta)**q * factor * moment)
            if not upper.is_finite() or not upper > 0:
                raise ArithmeticError('positive finite occupancy endpoint required')
            if best[q] is None or upper < best[q]:
                best[q], choices[q] = upper, tilt
        total = q1.kernel_t64.up(sum(best.values(), arb(0)))
        worst = sorted(best, key=lambda q: float(best[q].log()), reverse=True)[:8]
        margin = str(-total.log() / arb(2).log())
        record['trials'].append(dict(tilt=tilt, union_upper=q1.endpoint(total), margin_bits=margin,
            worst_q=worst, elapsed_seconds=monotonic() - trial_start))
        record.update(union_upper=q1.endpoint(total), margin_bits=margin,
            occupancy_uppers={str(q): q1.endpoint(v) for q, v in best.items()},
            occupancy_choices={str(q): value for q, value in choices.items()},
            elapsed_seconds=monotonic() - start)
        if sources != q1.source_snapshot():
            raise RuntimeError('loaded local mathematical source changed during screen')
        if output is not None:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(record, indent=2) + '\n')
        print(f'tilt={tilt} complete_q={q_min}..{q_max} union_margin={margin} '
              f'worst_q={worst} elapsed={monotonic() - trial_start:.2f}s', flush=True)
    print('Whole-code certificate: false; q=1 and q=2 are excluded.', flush=True)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', choices=('rs8', 'rs16'), default='rs8')
    parser.add_argument('--q-min', type=int, default=3)
    parser.add_argument('--q-max', type=int, default=32)
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--tilts', nargs='+', default=['.00512', '.01024', '.0256', '.0512', '.1024', '.2048'])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    run(outer=args.outer, q_min=args.q_min, q_max=args.q_max, precision=args.precision,
        tilts=args.tilts, output=args.output)


if __name__ == '__main__':
    main()
