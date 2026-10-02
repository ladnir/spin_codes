"""Pure-rational search scheduling; no cell is accepted by these helpers.

Broad unresolved intervals can be subdivided before invoking an expensive
proposal routine whose strongest refinement is disabled at those widths.
Existing accepted witnesses stay unchanged and still require fresh replay.
"""
import copy
from fractions import Fraction as Q


def path_cell(root, path):
    root = tuple(map(Q, root))
    if len(root) != 2 or not 0 <= root[0] < root[1] <= 1:
        raise ValueError('nonempty rational mean domain in [0,1] required')
    if not isinstance(path, str) or len(path) > 64 or set(path)-{'0', '1'}:
        raise ValueError('binary subdivision path of length at most64 required')
    lo, hi = root
    for bit in path:
        mid = (lo+hi)/2
        if bit == '0':
            hi = mid
        else:
            lo = mid
    return lo, hi


def partition(root, leaves, unresolved):
    """Reconstruct exact intervals and reject gaps, overlaps, or stale cells."""
    root = tuple(map(Q, root))
    path_cell(root, '')
    if not isinstance(leaves, dict) or not isinstance(unresolved, dict) or set(leaves) & set(unresolved):
        raise ValueError('disjoint accepted and unresolved maps required')
    rows = dict(leaves, **unresolved)
    if not rows:
        raise ValueError('complete nonempty partition required')
    cells = {}
    for path, row in rows.items():
        if not isinstance(row, dict):
            raise ValueError('partition rows must be dictionaries')
        cells[path] = path_cell(root, path)
        if 'cell' in row and tuple(map(Q, row['cell'])) != cells[path]:
            raise ValueError('saved coordinates differ from the exact path')
    previous = root[0]
    for lo, hi in sorted(cells.values()):
        if lo != previous:
            raise ValueError('partition has a gap or overlap')
        previous = hi
    if previous != root[1]:
        raise ValueError('partition misses the right endpoint')
    return cells


def presplit(root, cover, maximum_width=Q(1, 1024), maximum_depth=22):
    """Refine only unresolved cells, retaining the exact complete partition.

    Return (new_cover, statistics), without mutating cover. Newly split cells
    contain geometry only: neither a parent's score nor its variance-specific
    witness is valid as a child's certificate. The usual search/replay must
    evaluate every descendant before accepting it. Work counters are left
    unchanged because no numerical model evaluation occurred.
    """
    maximum_width = Q(maximum_width)
    if (not isinstance(cover, dict) or not 0 < maximum_width <= 1
            or type(maximum_depth) is not int or not 0 <= maximum_depth <= 64):
        raise ValueError('cover, positive rational width, and depth0..64 required')
    cells = partition(root, cover['leaves'], cover['unresolved'])
    result = copy.deepcopy(cover)
    result['unresolved'] = {}
    skipped = 0
    for original, row in cover['unresolved'].items():
        pending = [(original, cells[original])]
        while pending:
            path, (lo, hi) = pending.pop()
            if hi-lo <= maximum_width:
                result['unresolved'][path] = (copy.deepcopy(row) if path == original else
                    dict(cell=[str(lo), str(hi)], presplit_from=original))
                continue
            if len(path) >= maximum_depth:
                raise ValueError('maximum depth cannot reach the requested width gate')
            mid = (lo+hi)/2
            pending.extend(((path+'1', (mid, hi)), (path+'0', (lo, mid))))
            skipped += 1
    partition(root, result['leaves'], result['unresolved'])
    statistics = dict(method='geometry-only unresolved presplit', maximum_width=str(maximum_width),
        maximum_depth=maximum_depth, accepted_retained=len(result['leaves']),
        pending_before=len(cover['unresolved']), pending_after=len(result['unresolved']),
        skipped_internal_cells=skipped, numerical_evaluations=0)
    return result, statistics


def retile(root, cover, maximum_width=Q(1, 1024), maximum_depth=22):
    """Coalesce unresolved siblings, then apply an exact width-limited tiling.

    Accepted cells are barriers and remain untouched. Any coalesced parent
    receives geometry only; presplit likewise discards parent witnesses on
    newly split descendants. Thus even a merge followed by a split back to
    an earlier path cannot accidentally restore that path's old numeric data.
    Return (new_cover, statistics), with no numerical model evaluations and
    no mutation of the caller's cover or its work counters.
    """
    root = tuple(map(Q, root))
    maximum_width = Q(maximum_width)
    if (not isinstance(cover, dict) or not 0 < maximum_width <= 1
            or type(maximum_depth) is not int or not 0 <= maximum_depth <= 64):
        raise ValueError('cover, positive rational width, and depth0..64 required')
    partition(root, cover['leaves'], cover['unresolved'])
    coalesced = copy.deepcopy(cover)
    pending = coalesced['unresolved']
    merges = 0
    # Descending depth permits each newly created parent to continue upward.
    # Only members of the unresolved frontier participate in these merges.
    for start in sorted(tuple(pending), key=lambda path: (-len(path), path)):
        current = start
        while current and current in pending:
            sibling = current[:-1]+('1' if current[-1] == '0' else '0')
            if sibling not in pending:
                break
            parent = current[:-1]
            if parent in coalesced['leaves']:
                raise ValueError('cannot coalesce across an accepted cell')
            del pending[current]
            del pending[sibling]
            pending[parent] = dict(cell=list(map(str, path_cell(root, parent))))
            merges += 1
            current = parent
    partition(root, coalesced['leaves'], pending)
    result, split_stats = presplit(root, coalesced, maximum_width, maximum_depth)
    partition(root, result['leaves'], result['unresolved'])
    statistics = dict(method='geometry-only unresolved retile', maximum_width=str(maximum_width),
        maximum_depth=maximum_depth, accepted_retained=len(result['leaves']),
        pending_before=len(cover['unresolved']), pending_after_coalesce=len(pending),
        coalesced_internal_cells=merges, presplit_internal_cells=split_stats['skipped_internal_cells'],
        pending_after=len(result['unresolved']), numerical_evaluations=0)
    return result, statistics
