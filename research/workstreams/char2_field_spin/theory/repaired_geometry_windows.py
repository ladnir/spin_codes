"""Bounded threshold-aware window screens for two repaired outer geometries.

Anchors and relaxed-parity ladders are initialization heuristics only. Every
reported dimension uses the geometry's TRUE outer redundancy. A no-find is
not a distance certificate or a probability bound. Frozen old fixtures are
imported neither for modification nor as repaired-geometry witnesses.
"""

from hashlib import sha256
import argparse
import json

import numpy as np

from window_dimension_search import (generate_route, count_hits, dimension, coordinate_ascent,
                                     exchange_lengths, strict_window_geometry)


GEOMETRIES = {
    '33/22': {'rows': 47680, 'regions': 33, 'parity': 11},
    '64/44': {'rows': 23872, 'regions': 64, 'parity': 20},
}


def inverse_positions(route):
    result = np.empty_like(route)
    for inverse, region in zip(result, route):
        inverse[region] = np.arange(len(region))
    return result


def anchored_starts(inverse, anchors, width, alignment=16):
    """Aligned windows maximizing a lexicographic weighted anchor coverage.

    Weights are powers of two, so covering the first anchor beats every
    combination excluding it. Therefore that row is present in ALL regions.
    Additional anchors choose among windows containing the first. Ties prefer
    the first anchor's centered window; this is not a true-dimension score.
    """
    regions, rows = inverse.shape
    anchors = tuple(int(a) for a in anchors)
    if (not 1 <= len(anchors) <= 4 or len(set(anchors)) != len(anchors)
            or any(a not in range(rows) for a in anchors)
            or not 1 <= width <= rows or alignment < 1 or width % alignment or rows % alignment):
        raise ValueError('require one to four distinct anchors and a valid aligned width')
    choices = np.arange(0, rows-width+1, alignment, dtype=np.int64)
    weights = tuple(1 << (len(anchors)-1-j) for j in range(len(anchors)))
    starts = np.empty(regions, dtype=np.int64)
    for j in range(regions):
        positions = inverse[j, list(anchors)]
        scores = np.zeros(len(choices), dtype=np.int32)
        for position, weight in zip(positions, weights):
            scores += weight*((choices <= position) & (position < choices+width))
        centered = max(0, min(int(positions[0])-width//2, rows-width))//alignment*alignment
        best = int(scores.argmax())
        starts[j] = centered if scores[centered//alignment] == scores[best] else choices[best]
        assert starts[j] <= positions[0] < starts[j]+width
    return starts


def parity_ladder(route, width, initial, true_parity, ladder, *, alignment=16,
                  max_sweeps=12, seed=0):
    """Continue surrogate maxima, but retain candidates by true dimension.

    The returned AscentResult always has its dimension/history measured at
    true_parity. Intermediate surrogate values are separately named records.
    """
    if not ladder or any(not 0 <= p <= true_parity for p in ladder):
        raise ValueError('nonempty relaxed-parity ladder required')
    generator = np.random.Generator(np.random.PCG64(seed))
    starts = np.asarray(initial, dtype=np.int64).copy()
    best_starts = starts.copy()
    best_true = dimension(count_hits(route, starts, width), true_parity)
    records = []
    for relaxed in ladder:
        result = coordinate_ascent(route, width, relaxed, starts, alignment=alignment,
                                    max_sweeps=max_sweeps, generator=generator)
        true_value = dimension(result.counts, true_parity)
        records.append({'surrogate_parity': relaxed, 'surrogate_dimension': result.dimension,
                        'true_parity': true_parity, 'true_dimension': true_value,
                        'sweeps': len(result.history)-1})
        if true_value > best_true:
            best_true, best_starts = true_value, result.starts.copy()
        starts = result.starts
    polished = coordinate_ascent(route, width, true_parity, best_starts, alignment=alignment,
                                 max_sweeps=max_sweeps, generator=generator)
    assert polished.dimension >= best_true
    return polished, records


def initialize(route, width, seed, alignment=16):
    """Eight finite starts: two random, three centered, and two/three/four anchors."""
    generator = np.random.Generator(np.random.PCG64(seed))
    inverse = inverse_positions(route)
    regions, rows = route.shape
    choices = (rows-width)//alignment+1
    result = []
    for j in range(2):
        starts = generator.integers(0, choices, size=regions, dtype=np.int64)*alignment
        result.append((f'random-{j}', starts, ()))
    for j in range(3):
        anchor = int(generator.integers(rows))
        result.append((f'centered-{j}', anchored_starts(inverse, (anchor,), width, alignment), (anchor,)))
    for number in (2, 3, 4):
        anchors = tuple(int(a) for a in generator.choice(rows, size=number, replace=False))
        result.append((f'weighted-{number}', anchored_starts(inverse, anchors, width, alignment), anchors))
    return result


def screen(geometry, route_seed, *, search_seed=193711, exchange_seed=9173,
           exchanges=512, max_sweeps=12, alignment=16, state=4, requested_source=2**20):
    rows, regions, parity = (geometry[key] for key in ('rows', 'regions', 'parity'))
    width, e = strict_window_geometry(rows, regions, alignment)
    route = generate_route(rows, regions, route_seed)
    starts_list = initialize(route, width, search_seed, alignment)
    best = None
    records = []
    for index, (kind, starts, anchors) in enumerate(starts_list):
        initial_true = dimension(count_hits(route, starts, width), parity)
        if anchors:
            assert initial_true >= regions-parity
        generator = np.random.Generator(np.random.PCG64(search_seed+101+index))
        direct = coordinate_ascent(route, width, parity, starts, alignment=alignment,
                                   max_sweeps=max_sweeps, generator=generator)
        records.append({'initialization': kind, 'anchors': list(anchors), 'initial_true_dimension': initial_true,
                        'direct_true_dimension': direct.dimension})
        if best is None or direct.dimension > best.dimension:
            best, best_kind = direct, kind+'/direct'
        # Two bounded continuation paths: one random and the four-anchor path.
        if index in (0, len(starts_list)-1):
            ladder = tuple(range(min(10, parity), parity+1))
            continued, stages = parity_ladder(route, width, starts, parity, ladder,
                                              alignment=alignment, max_sweeps=max_sweeps,
                                              seed=search_seed+1001+index)
            records[-1]['ladder'] = stages
            records[-1]['ladder_polished_true_dimension'] = continued.dimension
            if continued.dimension > best.dimension:
                best, best_kind = continued, kind+'/ladder'
    before_exchange = best.dimension
    best, exchange_log = exchange_lengths(route, best, parity, alignment=alignment,
                                          seed=exchange_seed, proposals=exchanges,
                                          max_sweeps=max_sweeps)
    # Reconstruct by a different counting path, not incremental search state.
    selected = np.concatenate([route[j, first:first+length]
                               for j, (first, length) in enumerate(zip(best.starts, best.widths))])
    counts = np.bincount(selected, minlength=rows)
    assert dimension(counts, parity) == best.dimension
    assert int(best.widths.sum()) <= e
    outer_k = regions-parity
    padding = max(0, rows*outer_k-requested_source)
    excluded = (padding+outer_k-1)//outer_k
    retained_dimension = dimension(counts[:rows-excluded], parity)
    report = {'geometry': f'{regions}/{outer_k}', 'route_seed': route_seed,
              'numpy_version': np.__version__, 'rows': rows, 'regions': regions, 'true_parity': parity,
              'state': state, 'exit_constraints_upper': regions*state,
              'initial_window_width': width, 'total_window_symbols': int(best.widths.sum()),
              'bad_output_limit': e, 'true_dimension_before_exchange': before_exchange,
              'max_true_dimension_found': best.dimension, 'best_initialization': best_kind,
              'excluded_tail_rows': excluded, 'tail_excluded_true_dimension': retained_dimension,
              'realized_bad_route': best.dimension > regions*state,
              'retained_source_bad_route': retained_dimension > regions*state,
              'width_range': [int(best.widths.min()), int(best.widths.max())],
              'positive_dimension_rows': int((counts > parity).sum()),
              'max_row_hits': int(counts.max()), 'initialization_records': records,
              'length_exchanges': exchange_log,
              'route_sha256': sha256(route.astype('<i4').tobytes()).hexdigest(),
              'starts_sha256': sha256(best.starts.astype('<i8').tobytes()).hexdigest(),
              'widths_sha256': sha256(best.widths.astype('<i8').tobytes()).hexdigest(),
              'no_find_is_not_certificate': True, 'probability_claimed': False}
    if report['realized_bad_route']:
        report['compact_witness'] = {'starts': best.starts.tolist(), 'widths': best.widths.tolist(),
                                     'positive_dimension_row_hits': [[int(i), int(counts[i])]
                                                                     for i in np.flatnonzero(counts > parity)]}
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometries', nargs='+', choices=tuple(GEOMETRIES), default=tuple(GEOMETRIES))
    parser.add_argument('--route-seeds', nargs='+', type=int, default=(0, 1, 2))
    parser.add_argument('--exchanges', type=int, default=512)
    parser.add_argument('--max-sweeps', type=int, default=12)
    args = parser.parse_args()
    for name in args.geometries:
        for seed in args.route_seeds:
            print(json.dumps(screen(GEOMETRIES[name], seed, exchanges=args.exchanges,
                                     max_sweeps=args.max_sweeps), sort_keys=True), flush=True)
