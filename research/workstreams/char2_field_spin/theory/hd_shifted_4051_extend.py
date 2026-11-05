"""Bounded FLOAT-only extension of the frozen4096 boundary partition.

Recover cached affine intercepts from the original accepted-box screen
values, checking repeated values for consistency. No operator is built.
The output contains ordinary exact-replay proposals, never certificates.
"""
from hashlib import sha256
import heapq
import json
from math import log2
from pathlib import Path

import numpy as np

from hd_shifted_bulk_search import LOG_DENSITY, box_geometry, compact_plan, split_box
from hd_shifted_cover_backend import profile_count
from hd_shifted_swap_bulk import restore_identity
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION


def recovered_bank(plan):
    slopes = []
    for w in plan['witnesses']:
        (s0, f0), (s1, f1) = w['marker_log_odds']
        first = (w['alpha']*LOG_DENSITY+1-128*(s0+s1), -2+256*s0-128*f0, 256*s1-128*f1)
        slopes.append((first, np.asarray(first)+(0, 2, -2)))
    slopes = np.asarray(slopes)
    intercepts = np.full(len(slopes), np.nan)
    disagreement = 0.
    for leaf in plan['accepted_boxes']:
        index = leaf['witness_id']
        box = tuple(map(tuple, leaf['ranges']))
        corners, entropy = box_geometry(box)
        offset = np.min(np.max(slopes[index]@corners.T+entropy, axis=1))
        value = leaf['union_log2_upper_screen']-log2(profile_count(box))-offset
        if np.isfinite(intercepts[index]):
            disagreement = max(disagreement, abs(value-intercepts[index]))
        else:
            intercepts[index] = value
    if not np.all(np.isfinite(intercepts)) or disagreement > 1e-7:
        raise ArithmeticError('Inconsistent or missing recovered FLOAT intercept')
    # Cushion only guards reconstruction roundoff; this remains a FLOAT screen.
    return slopes, intercepts+1e-7, float(disagreement)


def extend(plan, *, first=4051, last=4095, threshold=-65., max_leaves=8192, progress=None):
    if (plan.get('construction') != SWAP_CONSTRUCTION or plan['requested_weight_range'][0] != 4096
            or not 1 <= first <= last < 4096 or max_leaves < 1):
        raise ValueError('Expected the complete frozen4096 boundary and a lower-weight extension')
    slopes, intercepts, disagreement = recovered_bank(plan)
    def evaluate(box):
        corners, entropy = box_geometry(box)
        values = intercepts+np.min(np.max(slopes@corners.T+entropy, axis=2), axis=1)
        index = int(np.argmin(values))
        return float(values[index]+log2(profile_count(box))), index
    pending, accepted, failed = [], [], []
    serial = 0
    def enqueue(box):
        nonlocal serial
        count = profile_count(box)
        if not count:
            return
        value, index = evaluate(box)
        row = dict(ranges=box, witness_id=index, union_log2_upper_screen=value,
                   integer_triples=count, inherited=False)
        if value <= threshold:
            accepted.append(row)
        else:
            heapq.heappush(pending, (-value, serial, row))
            serial += 1
    initial = []
    for row in plan['accepted_boxes']:
        if row['ranges'][0][0] <= 4096 <= row['ranges'][0][1]:
            box = ((first, last), tuple(row['ranges'][1]), tuple(row['ranges'][2]))
            if profile_count(box):
                initial.append(box)
                enqueue(box)
    expected = profile_count(((first, last), (0, 4096), (0, 4096)))
    if sum(profile_count(box) for box in initial) != expected:
        raise ArithmeticError('Boundary extension does not preserve all feasible profiles')
    leaf_count, splits = len(initial), 0
    while pending and leaf_count < max_leaves:
        _, _, row = heapq.heappop(pending)
        box = row['ranges']
        if all(lo == hi for lo, hi in box):
            failed.append(row)
            continue
        children = [child for child in split_box(box) if profile_count(child)]
        leaf_count += len(children)-1
        splits += 1
        for child in children:
            enqueue(child)
        if progress and splits % 256 == 0:
            progress(dict(splits=splits, leaves=leaf_count, accepted=len(accepted),
                          pending=len(pending), fixed_failures=len(failed)))
    failed.extend(row for _, _, row in pending)
    maximum = max((row['union_log2_upper_screen'] for row in accepted), default=-float('inf'))
    union = maximum+log2(sum(2**(row['union_log2_upper_screen']-maximum) for row in accepted)) if accepted else None
    result = dict(profile_coordinates=['r', 'a', 'b'], all_even_odd_splits=True,
        requested_weight_range=[first, last], accepted_boxes=accepted,
        unresolved_boxes=[dict(ranges=row['ranges'], reason='frozen-bank-or-leaf-budget',
                               best_log2_upper_screen=row['union_log2_upper_screen']) for row in failed],
        witnesses=plan['witnesses'], accepted_integer_triples=sum(row['integer_triples'] for row in accepted),
        feasible_domain_triples=expected, accepted_domain_fraction=sum(row['integer_triples'] for row in accepted)/expected,
        union_log2_upper_screen=union, threshold=threshold, directed=False, whole_code_certificate=False,
        reconstruction=dict(max_intercept_disagreement=disagreement, intercept_cushion=1e-7,
                            initial_boxes=len(initial), max_leaves=max_leaves, splits=splits,
                            final_leaves=leaf_count, operator_builds=0, witness_fits=0,
                            requires_fresh_directed_replay=True))
    # compact_plan expects these historical bookkeeping fields.
    result.update(proposal_budget=0, proposals_used=0, slab_width=45, full_width=256, minimum_width=1)
    return restore_identity(compact_plan(result))


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=Path('hd_shifted_swap_bulk_plan3.json'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--max-leaves', type=int, default=8192)
    p.add_argument('--first', type=int, default=4051)
    p.add_argument('--last', type=int, default=4095)
    args = p.parse_args()
    if args.output.exists():
        p.error('Frozen output already exists')
    raw = args.source.read_bytes()
    result = extend(json.loads(raw), first=args.first, last=args.last, max_leaves=args.max_leaves,
                    progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(raw).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, separators=(',', ':'), allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('accepted_boxes', 'unresolved_boxes', 'witnesses')}), flush=True)
    print(json.dumps(dict(path=str(args.output), accepted_boxes=len(result['accepted_boxes']),
                          unresolved_boxes=len(result['unresolved_boxes']), witnesses=len(result['witnesses']))))
