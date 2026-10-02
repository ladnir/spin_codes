"""Exact t64/s16 maps spanned by affine functions and nine monomials.

Coordinate x is the six-bit vector (x_0,...,x_5), least significant first.
Rows are 1,x_0,...,x_5 followed by x_i*x_j for a sorted nine-edge graph.
The edge (0,1) gives rank four on each consecutive four-coordinate packet.
All products of two rows have degree at most four in six variables, hence
their sums over all coordinates vanish: C A=0 for C=A^T.
"""
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path

import screen_t128 as screen

EDGES = tuple(combinations(range(6), 2))
HERE = Path(__file__).resolve()


def quadratic_rank(mask):
    rows = [0] * 6
    for k, (i, j) in enumerate(EDGES):
        if mask >> k & 1:
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    return screen.q1.kernel_t64.s16_maps.binary_rank(rows)


def graph_census():
    """Rank spectra for every nine-edge graph containing (0,1), exactly."""
    counts = {rank: [int(quadratic_rank(mask) == rank) for mask in range(1 << 15)]
              for rank in (0, 2, 4, 6)}
    for values in counts.values():
        for bit in range(15):
            for mask in range(1 << 15):
                if mask >> bit & 1:
                    values[mask] += values[mask ^ (1 << bit)]
    result = []
    for mask in range(1 << 15):
        if mask & 1 and mask.bit_count() == 9:
            ranks = tuple(counts[rank][mask] for rank in (0, 2, 4, 6))
            result.append((ranks[1], ranks[2], mask, ranks))
    return sorted(result)


def edge_set(name):
    if name == 'minrank':
        mask = graph_census()[0][2]
        return tuple(e for k, e in enumerate(EDGES) if mask >> k & 1)
    if name == 'bipartite':
        return tuple(e for e in EDGES if (e[0] % 2) != (e[1] % 2))
    if name == 'prism':
        return ((0, 1), (0, 2), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (3, 5), (4, 5))
    if name == 'lex':
        return EDGES[:9]
    raise ValueError('named monomial graph required')


def prepare(name):
    edges = edge_set(name)
    rows = ((1 << 64) - 1,) + tuple(sum(1 << x for x in range(64) if x >> i & 1) for i in range(6))
    rows += tuple(sum(1 << x for x in range(64) if (x >> i & 1) and (x >> j & 1)) for i,j in edges)
    maps = screen.q1.kernel_t64.s16_maps
    columns = tuple(sum(((r >> x) & 1) << j for j,r in enumerate(rows)) for x in range(64))
    if (len(rows) != 16 or maps.binary_rank(rows) != 16 or maps.binary_rank(columns) != 16
            or any((a & b).bit_count() & 1 for a in rows for b in rows)
            or any(maps.binary_rank(columns[i:i+4]) != 4 for i in range(0,64,4))):
        raise ArithmeticError('dimension, CA=0, or packet rank failed')
    images = tuple(maps.images_from_rows(rows))
    spectrum = dict(sorted(Counter(x.bit_count() for x in images).items()))
    physical = screen.q1.kernel_t64.kernel_maps.prepare_maps(images, columns, bits=16,
        distribution='uniform_gl', birth_density='capped')
    data = screen.q1.kernel_t64.wrap(physical)
    record = dict(name=name, monomial_pairs=list(map(list, edges)),
        monomial_masks=[(1 << i) | (1 << j) for i,j in edges],
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=list(columns),
        expansion_spectrum={str(w): n for w,n in spectrum.items()},
        full_state_census=True, state_count=65536, physical_t=64, state_bits=16,
        feedback_definition='C=A^T', feedback_times_expansion_zero=True,
        independent_transitive_update_per_step=True, macro_physical_steps=2,
        minimum_expansion_weight=min(w for w in spectrum if w), map_sha256=physical['map_sha256'],
        source_sha256={str(HERE): hashlib.sha256(HERE.read_bytes()).hexdigest()})
    return data, record


if __name__ == '__main__':
    census = graph_census()
    minimum = census[0][:2]
    best = [entry for entry in census if entry[:2] == minimum]
    named_masks = {name: sum(1 << EDGES.index(edge) for edge in edge_set(name))
                   for name in ('bipartite', 'prism', 'lex', 'minrank')}
    print(json.dumps(dict(graphs=len(census), best_rank2=minimum[0], best_rank4=minimum[1],
        minimizers=len(best), first_edges=edge_set('minrank'),
        named_rank_counts={name: Counter(quadratic_rank(mask) for mask in range(1 << 15)
            if not mask & ~allowed) for name,allowed in named_masks.items()}), indent=2))
