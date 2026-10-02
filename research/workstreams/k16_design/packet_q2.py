"""Occupancy-two RS diagnostic with ordered, continuous-state placement.

The regional placement averages are exact averages of local upper-envelope
operators. They are not exact transition probabilities. This file covers all
two-active-group support pairs, not other occupancies or a whole-code union.
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

import packet_q1 as q1
import rs_outer
from flint import arb, arb_mat, ctx


def pair_support_moments(regional, regions, *, matrix=arb_mat,
                         scalar=arb, rounding=q1.kernel_t64.up,
                         progress=None):
    """Return triangular moments[v][u], 0 <= u <= v <= regions.

For R_j the j-active-packet regional operator, retain the coefficients of
e_zero (R_0 + (x+y) R_1 + xy R_2)^regions as row vectors. Divide coefficient
(u,v) by C(regions,u) C(regions,v), then sum terminal mass coordinates.
The recurrence appends each region on the right: it never resets state or
commutes the R_j. Symmetry only exchanges the two labeled active groups.

All coordinates of the selected t64/s16 class envelope are terminal mass.
Use exact rational matrices and identity rounding for small algebraic tests.
"""
    if type(regions) is not int or regions < 0 or len(regional) != 3:
        raise ValueError('three regional operators and a nonnegative region count required')
    size = regional[0].nrows()
    if (size < 1 or any(m.nrows() != size or m.ncols() != size for m in regional)
            or any(not m[i, j] >= 0 or
                   (hasattr(m[i, j], 'is_finite') and not m[i, j].is_finite()) for m in regional
                   for i in range(size) for j in range(size))):
        raise ValueError('matching finite nonnegative square operators required')
    r0, r1, r2 = regional
    previous = [[matrix([[1] + [0] * (size - 1)])]]
    for r in range(1, regions + 1):
        following = []
        for v in range(r + 1):
            column = []
            for u in range(v + 1):
                value = matrix(1, size)
                if v < r:
                    value += previous[v][u] * r0
                # Combine both single-group predecessors before multiplying.
                singles = None
                if u and v < r:
                    singles = previous[v][u - 1]
                if v and u < r:
                    other = previous[v - 1][u] if u < v else previous[u][v - 1]
                    singles = other if singles is None else singles + other
                if singles is not None:
                    value += singles * r1
                if u:
                    value += previous[v - 1][u - 1] * r2
                column.append(matrix([[rounding(value[0, j]) for j in range(size)]]))
            following.append(column)
        previous = following
        if progress is not None:
            progress(r)
    binomials = [comb(regions, u) for u in range(regions + 1)]
    return [[rounding(sum((row[0, j] for j in range(size)), scalar(0)) /
                      (binomials[u] * binomials[v]))
             for u, row in enumerate(column)]
            for v, column in enumerate(previous)]


def fold_shell_pairs(counts, weights, group_count):
    """Exact expected shell-count products; these inputs must not be CDF caps."""
    total = arb(0)
    for v, column in enumerate(weights):
        if not counts[v]:
            continue
        for u, probability in enumerate(column):
            if counts[u]:
                # The physical group pair is unordered, its two messages labeled.
                count = counts[u] * counts[v] * (1 if u == v else 2)
                total += q1.kernel_t64.aq(count) * probability
    return q1.kernel_t64.up(comb(group_count, 2) * total)


def evaluate_q2(counts, *, geometry, tilts, precision=192, data=None,
                map_record=None, metadata=None, output=None):
    """Bound all q=2 messages under exact independent expected shell counts.

Each group's support must be uniform conditional on its size, with uniform
independent nonzero packet labels. The two groups' outer setups and column
shuffles are independent. Regional routing samples distinct slots uniformly.
The initial state is zero; every subsequent region retains the previous state.
"""
    if not isinstance(geometry, q1.Geometry) or geometry.group_count < 2:
        raise ValueError('checked half-rate packet geometry with at least two groups required')
    if (type(precision) is not int or precision < 128 or not tilts
            or any(Q(t) <= 0 for t in tilts) or len(set(map(Q, tilts))) != len(tilts)):
        raise ValueError('precision >=128 and distinct positive rational tilts required')
    counts = tuple(map(Q, counts))
    if (len(counts) != geometry.regions + 1 or counts[0] != 0
            or any(c < 0 for c in counts)
            or sum(counts) != (1 << geometry.group_dimension) - 1):
        raise ValueError('complete exact nonzero-message shell counts required, not a CDF')
    if output is not None and output.exists():
        raise ValueError('output must be fresh')
    start = monotonic()
    regions = geometry.regions
    threshold = geometry.N // 10
    ctx.prec = precision
    print(f'Q2 DIAGNOSTIC K={geometry.K} N={geometry.N} L={geometry.group_count} '
          f'regions={regions} macros_per_region={geometry.macros_per_region} '
          f'cutoff={threshold}', flush=True)
    if data is None:
        data, map_record = q1.kernel_t64.prepare(birth_density='capped')
    sources = q1.source_snapshot()
    best = [[arb(1) for _ in range(v + 1)] for v in range(regions + 1)]
    choices = [[None for _ in range(v + 1)] for v in range(regions + 1)]
    record = dict(schema='finite-packet-q2-diagnostic-1', K=geometry.K, N=geometry.N,
        geometry=asdict(geometry), groups=geometry.group_count, regions=regions,
        group_output_bits=4 * regions, group_dimension=geometry.group_dimension,
        distance='1/10', threshold=threshold, physical_t=64, state_bits=16,
        physical_steps=geometry.N // 64, macro_t=128, macro_steps=geometry.N // 128,
        macros_per_region=geometry.macros_per_region, physical_steps_per_macro=2,
        zero_initial_state=True, state_continuity='retained_between_all_steps_and_regions',
        final_flush=False, occupancy_covered=[2], occupancies_not_covered=[[1, 1], [3, geometry.group_count]],
        all_two_group_support_pairs_covered=True, whole_code_certificate=False,
        precision=precision, tilts=list(tilts), count_kind='exact_expected_shells',
        count_sha256=hashlib.sha256(json.dumps(list(map(str, counts))).encode()).hexdigest(),
        count_premises=metadata, map_record=map_record, source_sha256=sources,
        support_pair_storage='row v contains u=0..v; off-diagonal count products doubled',
        regional_placement='exact without-replacement average of upper-envelope operators',
        support_placement='ordered bivariate coefficient / (binom(R,u)*binom(R,v))',
        scope='Independent outer setups and uniform column shuffles for the two labeled '
              'active groups; independent uniform regional shuffles, uniform nonzero '
              'four-bit labels, selected t64/s16 maps, independent ideal uniform GL16 '
              'updates. Upper-envelope diagnostic, not seed-specific or whole-code certification.',
        trials=[])
    for tilt in tilts:
        trial_start = monotonic()
        local = q1.kernel_t64.local_operators(data, Q(tilt), activity=Q(1, 2))
        regional = q1.placement(local, epochs=geometry.macros_per_region,
            windows=geometry.macro_windows, rounding=q1.rounded, maximum_groups=2)
        print(f'tilt={tilt} regional R0/R1/R2 ready, state_coordinates={regional[0].nrows()}', flush=True)

        def progress(r):
            if r % 16 == 0:
                print(f'tilt={tilt} ordered regions={r}/{regions} '
                      f'elapsed={monotonic() - trial_start:.2f}s', flush=True)

        moments = pair_support_moments(regional, regions, progress=progress)
        factor = (q1.kernel_t64.aq(Q(tilt)) * threshold).exp()
        for v, column in enumerate(moments):
            for u, moment in enumerate(column):
                candidate = q1.kernel_t64.up(factor * moment)
                if candidate < best[v][u]:
                    best[v][u], choices[v][u] = candidate, tilt
        upper = fold_shell_pairs(counts, best, geometry.group_count)
        if ctx.prec != precision or not upper.is_finite() or not upper > 0:
            raise ArithmeticError('positive finite diagnostic endpoint required')
        margin = str(-upper.log() / arb(2).log())
        terms = [(q1.kernel_t64.aq(counts[u] * counts[v] * (1 if u == v else 2)) *
                  best[v][u], u, v)
                 for v in range(regions + 1) for u in range(v + 1) if counts[u] and counts[v]]
        dominant = [[u, v] for _, u, v in sorted(terms, key=lambda item: float(item[0]), reverse=True)[:8]]
        record['trials'].append(dict(tilt=tilt, q2_upper=q1.endpoint(upper),
            margin_bits=margin, dominant_support_pairs=dominant,
            state_coordinates=regional[0].nrows(), elapsed_seconds=monotonic() - trial_start))
        record.update(q2_upper=q1.endpoint(upper), margin_bits=margin,
            q2_below_2_minus_40=bool(upper < arb(2) ** -40),
            support_pair_choices=choices,
            support_pair_probability_uppers=[[q1.endpoint(p) for p in column] for column in best],
            elapsed_seconds=monotonic() - start)
        if sources != q1.source_snapshot():
            raise RuntimeError('loaded local mathematical source changed during diagnostic')
        if output is not None:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(record, indent=2) + '\n')
        print(f'tilt={tilt} q2_margin_bits={margin} dominant_support_pairs={dominant}', flush=True)
    print('Only q=2 is covered; q=1, q>=3, and the whole-code union are not certified.', flush=True)
    return record


def run(tilts, precision=192, output=None):
    geometry = q1.Geometry(rs_outer.GROUP_COUNT, rs_outer.REGION_COUNT, rs_outer.GROUP_DIMENSION)
    if geometry.K != 1 << 16:
        raise ValueError('this entry point is specifically the proposed K16 RS geometry')
    return evaluate_q2(rs_outer.expected_group_support_counts(), geometry=geometry,
        tilts=tilts, precision=precision, output=output,
        metadata=dict(outer='four parallel GF256 RS[8,4] words',
            label_mixing='independent uniform GL32 per aligned four-byte symbol',
            independent_setups_between_groups=True, exact_expected_shell_counts=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts', nargs='+', default=['.00512', '.01024', '.0256', '.0512'])
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    run(args.tilts, args.precision, args.output)


if __name__ == '__main__':
    main()
