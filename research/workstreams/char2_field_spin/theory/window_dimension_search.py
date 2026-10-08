"""Bounded search for a realized-route shortened-dimension obstruction.

Choose one aligned window in each coordinate region. A route with
sum_i max(c_i-parity,0) > regions*state has a nonzero source supported in
those windows whose state vanishes after every window. This rejects that
realized route for every linear inner coefficient assignment. Search success
alone is not a distribution-level lower bound on setup failure.
"""

from dataclasses import dataclass
from hashlib import sha256
import argparse
import json

import numpy as np


@dataclass
class AscentResult:
    starts: np.ndarray
    widths: np.ndarray
    counts: np.ndarray
    dimension: int
    history: tuple
    coordinate_updates: int
    converged: bool


def generate_route(rows, regions, seed):
    """Deterministic PCG64 replay; arrays map positions to outer-row labels."""
    if rows < 1 or regions < 1 or regions > 32767:
        raise ValueError('positive geometry with at most 32767 regions required')
    generator = np.random.Generator(np.random.PCG64(seed))
    route = np.empty((regions, rows), dtype=np.int32)
    for region in range(regions):
        route[region] = generator.permutation(rows)
    return route


def validate_route(route):
    if route.ndim != 2 or not all(route.shape):
        raise ValueError('route must be a nonempty matrix')
    rows = route.shape[1]
    expected = np.arange(rows)
    if any(not np.array_equal(np.sort(region), expected) for region in route):
        raise ValueError('every region must be a permutation of the row labels')


def dimension(counts, parity):
    if parity < 0:
        raise ValueError('nonnegative outer redundancy required')
    return int(np.maximum(np.asarray(counts, dtype=np.int64)-parity, 0).sum())


def strict_window_geometry(rows, regions, alignment=16):
    """Largest equal aligned windows fitting ceil(N/10)-1, not floor(N/10)."""
    if rows < 1 or regions < 1 or alignment < 1:
        raise ValueError('positive geometry required')
    bad_output_limit = (rows*regions+9)//10-1
    width = alignment*(bad_output_limit//(regions*alignment))
    return width, bad_output_limit


def count_hits(route, starts, width):
    regions, rows = route.shape
    widths = np.full(regions, width, dtype=np.int64) if np.ndim(width) == 0 else np.asarray(width, dtype=np.int64)
    if len(starts) != regions or len(widths) != regions or np.any((widths < 1) | (widths > rows)):
        raise ValueError('invalid starts or window width')
    counts = np.zeros(rows, dtype=np.int32)
    for region, first, length in zip(route, starts, widths):
        if not 0 <= first <= rows-length:
            raise ValueError('window outside its region')
        counts[region[int(first):int(first+length)]] += 1
    return counts


def aligned_scores(region, counts_without_region, parity, width, alignment):
    """Exact dimension increments for every aligned start in this region."""
    if not 1 <= width <= len(region) or alignment < 1:
        raise ValueError('invalid width/alignment')
    hot = (counts_without_region >= parity)[region]
    prefix = np.empty(len(region)+1, dtype=np.int32)
    prefix[0] = 0
    np.cumsum(hot, dtype=np.int32, out=prefix[1:])
    starts = np.arange(0, len(region)-width+1, alignment, dtype=np.int64)
    return prefix[starts+width]-prefix[starts]


def coordinate_ascent(route, width, parity, initial_starts, *, alignment=16,
                      max_sweeps=24, generator=None):
    """Strict coordinate improvement; keep the current start on every tie.

    The prefix score is exact because adding a row coordinate increases
    max(c-parity,0) precisely when its previous count is at least parity.
    Every sweep is independently checked by a fresh exact dimension sum.
    """
    regions, rows = route.shape
    widths = np.full(regions, width, dtype=np.int64) if np.ndim(width) == 0 else np.asarray(width, dtype=np.int64).copy()
    if (len(widths) != regions or np.any((widths < 1) | (widths > rows))
            or alignment < 1 or np.any(widths % alignment)
            or max_sweeps < 1 or parity < 0):
        raise ValueError('invalid geometry or search limit')
    starts = np.asarray(initial_starts, dtype=np.int64).copy()
    if np.any(starts % alignment):
        raise ValueError('initial starts must be aligned')
    counts = count_hits(route, starts, widths)
    candidates_by_width = {int(length): np.arange(0, rows-length+1, alignment, dtype=np.int64)
                           for length in np.unique(widths)}
    hot = np.empty(rows, dtype=np.bool_)
    ordered = np.empty(rows, dtype=np.bool_)
    prefix = np.empty(rows+1, dtype=np.int32)
    prefix[0] = 0
    value = dimension(counts, parity)
    history = [value]
    updates = 0
    converged = False
    for _ in range(max_sweeps):
        changed = False
        order = range(regions) if generator is None else generator.permutation(regions)
        for region_index in order:
            region = route[region_index]
            length = int(widths[region_index])
            candidates = candidates_by_width[length]
            previous = int(starts[region_index])
            counts[region[previous:previous+length]] -= 1
            np.greater_equal(counts, parity, out=hot)
            np.take(hot, region, out=ordered)
            np.cumsum(ordered, dtype=np.int32, out=prefix[1:])
            scores = prefix[candidates+length]-prefix[candidates]
            best_index = int(scores.argmax())
            old_score = int(scores[previous//alignment])
            gain = int(scores[best_index])-old_score
            if gain > 0:
                previous = int(candidates[best_index])
                starts[region_index] = previous
                value += gain
                changed = True
                updates += 1
            counts[region[previous:previous+length]] += 1
        exact = dimension(counts, parity)
        if exact != value or exact < history[-1]:
            raise AssertionError('coordinate score failed exact monotonicity check')
        history.append(exact)
        if not changed:
            converged = True
            break
    return AscentResult(starts, widths, counts, value, tuple(history), updates, converged)


def search_route(route, width, parity, *, alignment=16, search_seed=0,
                 restarts=8, perturbations=2, perturb_regions=4, max_sweeps=24):
    """Finite deterministic-seed multistart search, never a global optimizer."""
    regions, rows = route.shape
    if restarts < 1 or perturbations < 0 or not 1 <= perturb_regions <= regions:
        raise ValueError('invalid finite restart/perturbation limits')
    if not 1 <= width <= rows or alignment < 1 or width % alignment:
        raise ValueError('invalid aligned window geometry')
    validate_route(route)
    generator = np.random.Generator(np.random.PCG64(search_seed))
    choices = (rows-width)//alignment+1
    best = None
    records = []
    for restart in range(restarts):
        starts = generator.integers(0, choices, size=regions, dtype=np.int64)*alignment
        for stage in range(perturbations+1):
            if stage:
                starts = current.starts.copy()
                changed_regions = generator.choice(regions, size=perturb_regions, replace=False)
                starts[changed_regions] = generator.integers(0, choices, size=perturb_regions)*alignment
            current = coordinate_ascent(route, width, parity, starts, alignment=alignment,
                                        max_sweeps=max_sweeps, generator=generator)
            record = {'restart': restart, 'stage': stage, 'initial_dimension': current.history[0],
                      'dimension': current.dimension, 'sweeps': len(current.history)-1,
                      'coordinate_updates': current.coordinate_updates, 'converged': current.converged}
            records.append(record)
            if best is None or current.dimension > best.dimension:
                best = current
                best_record = record
    # Reconstruct the witness from its starts; do not trust incremental counts.
    replay = count_hits(route, best.starts, width)
    if not np.array_equal(replay, best.counts) or dimension(replay, parity) != best.dimension:
        raise AssertionError('best witness failed independent count reconstruction')
    return best, best_record, records


def exchange_lengths(route, initial, parity, *, alignment=16, seed=9173,
                     proposals=256, amounts=(16, 64, 256, 1024), max_sweeps=12):
    """Bounded pair-length exchanges with exact coordinate polishing.

    Each proposal preserves the total length and keeps every window nonempty.
    Accept only a nondecreasing exact dimension; ties may change the layout.
    """
    regions, rows = route.shape
    if proposals < 0 or any(amount < alignment or amount % alignment for amount in amounts) or not amounts:
        raise ValueError('invalid bounded aligned exchanges')
    generator = np.random.Generator(np.random.PCG64(seed))
    current = initial
    total = int(initial.widths.sum())
    accepted = 0
    improvements = []
    for proposal in range(proposals):
        grow, shrink = generator.choice(regions, size=2, replace=False)
        amount = int(amounts[int(generator.integers(len(amounts)))])
        lengths = current.widths.copy()
        if lengths[grow]+amount > rows or lengths[shrink]-amount < alignment:
            continue
        lengths[grow] += amount
        lengths[shrink] -= amount
        starts = np.minimum(current.starts, (rows-lengths)//alignment*alignment)
        candidate = coordinate_ascent(route, lengths, parity, starts, alignment=alignment,
                                      max_sweeps=max_sweeps, generator=generator)
        if candidate.dimension >= current.dimension:
            if candidate.dimension > current.dimension:
                improvements.append({'proposal': proposal, 'dimension': candidate.dimension})
            current = candidate
            accepted += 1
    reconstructed = count_hits(route, current.starts, current.widths)
    if int(current.widths.sum()) != total or dimension(reconstructed, parity) != current.dimension:
        raise AssertionError('length exchange failed budget or dimension replay')
    return current, {'seed': seed, 'proposals': proposals, 'accepted': accepted,
                     'amounts': list(amounts), 'improvements': improvements}


def summarize(route_seed, search_seed, route, width, parity, state, alignment,
              best, record, records, bad_output_limit, exclude_tail_rows=0,
              include_witness=False):
    regions, rows = route.shape
    histogram = np.bincount(best.counts, minlength=regions+1)
    total_length = int(best.widths.sum())
    if not 0 <= exclude_tail_rows <= rows:
        raise ValueError('invalid excluded tail row count')
    tail_free_dimension = dimension(best.counts[:rows-exclude_tail_rows], parity)
    summary = {'route_seed': route_seed, 'search_seed': search_seed, 'numpy_version': np.__version__,
            'regions': regions, 'rows': rows, 'outer_redundancy': parity, 'state': state,
            'initial_window_width': width, 'width_range': [int(best.widths.min()), int(best.widths.max())],
            'alignment': alignment, 'total_window_symbols': total_length,
            'bad_output_limit': bad_output_limit, 'max_dimension_found': best.dimension,
            'full_code_source_symbols': rows*(regions-parity),
            'excluded_tail_rows': exclude_tail_rows, 'tail_excluded_dimension': tail_free_dimension,
            'exit_constraints_upper': regions*state, 'remaining_dimension_lower': best.dimension-regions*state,
            'realized_bad_route_witness': best.dimension > regions*state and total_length <= bad_output_limit,
            'equal_window_search_location': record, 'local_search_runs': len(records),
            'positive_dimension_rows': int(histogram[parity+1:].sum()),
            'high_hit_histogram': {str(c): int(histogram[c]) for c in range(parity+1, regions+1) if histogram[c]},
            'starts_sha256': sha256(best.starts.astype('<i8').tobytes()).hexdigest(),
            'widths_sha256': sha256(best.widths.astype('<i8').tobytes()).hexdigest(),
            'route_sha256': sha256(route.astype('<i4').tobytes()).hexdigest(),
            'distribution_failure_lower_bound_claimed': False}
    if include_witness and summary['realized_bad_route_witness']:
        summary['compact_witness'] = {'starts': best.starts.tolist(), 'widths': best.widths.tolist(),
                                      'positive_dimension_row_hits': [[int(i), int(best.counts[i])]
                                                                      for i in np.flatnonzero(best.counts > parity)]}
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--route-seeds', type=int, nargs='+', default=(0, 1, 2))
    parser.add_argument('--rows', type=int, default=47680)
    parser.add_argument('--regions', type=int, default=32)
    parser.add_argument('--parity', type=int, default=10)
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--alignment', type=int, default=16)
    parser.add_argument('--width', type=int)
    parser.add_argument('--restarts', type=int, default=8)
    parser.add_argument('--perturbations', type=int, default=2)
    parser.add_argument('--perturb-regions', type=int, default=4)
    parser.add_argument('--max-sweeps', type=int, default=24)
    parser.add_argument('--search-seed', type=int, default=193711)
    parser.add_argument('--length-exchanges', type=int, default=0)
    parser.add_argument('--exchange-seed', type=int, default=9173)
    parser.add_argument('--exclude-tail-rows', type=int, default=0)
    parser.add_argument('--include-witness', action='store_true')
    args = parser.parse_args()
    default_width, e = strict_window_geometry(args.rows, args.regions, args.alignment)
    width = args.width if args.width is not None else default_width
    for seed in args.route_seeds:
        route = generate_route(args.rows, args.regions, seed)
        best, record, records = search_route(route, width, args.parity, alignment=args.alignment,
                                             search_seed=args.search_seed, restarts=args.restarts,
                                             perturbations=args.perturbations, perturb_regions=args.perturb_regions,
                                             max_sweeps=args.max_sweeps)
        exchanges = None
        if args.length_exchanges:
            best, exchanges = exchange_lengths(route, best, args.parity, alignment=args.alignment,
                                                seed=args.exchange_seed, proposals=args.length_exchanges)
        summary = summarize(seed, args.search_seed, route, width, args.parity, args.state,
                            args.alignment, best, record, records, e,
                            args.exclude_tail_rows, args.include_witness)
        if exchanges is not None:
            summary['length_exchange_search'] = exchanges
        print(json.dumps(summary, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
