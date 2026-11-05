"""Bounded FLOAT positive-operator box proposals for the swap ensemble.

The old partition engine is reused without identifying the constructions.
A separate subdivision continuation loads only explicitly labeled swap
plans. Every output retains a complete partition, including unresolved boxes.
"""
from pathlib import Path
from math import log2
import heapq

import numpy as np

from hd_shifted_bulk_search import (AdaptiveAlphaProposer, best_witness, box_geometry,
    compact_plan, initial_partition, search, split_box, volume)
from hd_shifted_cover_backend import profile_count
from hd_shifted_swap_screen import CONSTRUCTION, SwapPositiveWitness


class SwapProposer(AdaptiveAlphaProposer):
    witness_type = SwapPositiveWitness

    def __init__(self, frontier=(), *, seed_boxes=(), threshold=-96., **kwargs):
        super().__init__(**kwargs)
        self.frontier = tuple(frontier)
        self.seed_boxes = tuple(seed_boxes)
        self.threshold = threshold
        self.seed_evaluations = 0

    def initial_bank(self):
        result = []
        for row in self.frontier:
            if self.exhausted:
                break
            witness = row['positive']
            result.append(self.create(witness['mark'], witness['alpha'],
                                      witness['marker_log_odds'], row['profile']))
        result.extend(super().initial_bank())
        # The original alpha proposer fixes a mark tuned to r>=5120.
        # Supply density-scaled witnesses before its local refinements.
        # Each candidate still checks the entire clipped coarse box.
        for box in self.seed_boxes:
            if self.exhausted:
                break
            _, score, _ = best_witness(result, box)
            if score <= self.threshold:
                continue
            centre = tuple(np.mean(box_geometry(box)[0], axis=0))
            r, a, b = centre
            mark = 2**(-r/1000)
            for alpha in (.1, .25, .5):
                if self.exhausted:
                    break
                odds = []
                for full in (a, b):
                    empty, singleton = 4096-r+full, r-2*full
                    odds.append((log2(max(1, singleton)/max(1, empty))+.5*alpha/.1,
                                 log2(max(1, full)/max(1, empty))+.9*alpha/.1))
                candidate = self.create(mark, alpha, odds, centre)
                self.seed_evaluations += 1
                value, _ = candidate.evaluate(box)
                if value < score:
                    result.append(candidate)
                    score = value
                if score <= self.threshold:
                    break
        return result

    def statistics(self):
        return dict(super().statistics(), density_seed_boxes=len(self.seed_boxes),
                    density_seed_evaluations=self.seed_evaluations)


def restore_identity(plan):
    """compact_plan intentionally keeps a minimal old-format witness schema."""
    plan = dict(plan, construction=CONSTRUCTION, pair_swaps=True,
                orientation_aggregation='independent-fair-swaps')
    plan['witnesses'] = [dict(row, construction=CONSTRUCTION, pair_swaps=True,
                             orientation_aggregation='independent-fair-swaps')
                         for row in plan['witnesses']]
    return plan


def bulk(frontier=(), *, first=3072, last=5119, proposals=300,
         max_evaluations=6000, marker_maxfev=80, weight_maxfev=16,
         threshold=-96., progress=None):
    initial, boundaries = initial_partition(first, last, slab_width=256, full_width=256)
    seed_boxes = [box for box, _ in initial]+boundaries
    adaptive = SwapProposer(frontier, seed_boxes=seed_boxes, threshold=threshold,
                            max_evaluations=max_evaluations,
                            marker_maxfev=marker_maxfev, weight_maxfev=weight_maxfev)
    result = search(first=first, last=last, proposals=proposals, threshold=threshold,
                    slab_width=256, full_width=256, minimum_width=16,
                    adaptive=adaptive, progress=progress)
    return restore_identity(compact_plan(result))


def subdivide(saved, frontier=(), *, bank_plans=(), first=4096, last=4351, threshold=-96.,
              max_splits=4096, max_candidates=2000, centre_boxes=False, progress=None):
    """Refine a clipped band using a frozen bank, without local optimization.

    With centre_boxes, translate the best witness to each failing box and
    try a bounded analytic density family before splitting. The default
    retains the earlier centre-only fallback for frozen reproduction.
    All accepted leaves are checked on the full clipped histogram domain.
    """
    if (saved.get('construction') != CONSTRUCTION or max_splits < 0 or max_candidates < 0
            or not saved['requested_weight_range'][0] <= first <= last <= saved['requested_weight_range'][1]):
        raise ValueError('Invalid construction, bounded work, or continuation band')
    proposals = list(saved['witnesses'])
    for extra in bank_plans:
        if (extra.get('construction') != CONSTRUCTION or any(
                w.get('construction') != CONSTRUCTION or w.get('pair_swaps') is not True
                or w.get('orientation_aggregation') != 'independent-fair-swaps'
                for w in extra.get('witnesses', []))):
            raise ValueError('An extra bank must explicitly select independent fair swaps')
        proposals.extend(extra['witnesses'])
    # An extra bank supplies witnesses only, never its partition or claims.
    # Deduplicate before building positive matrices; translated-marker seeds
    # retain their own centres even when the other witness fields agree.
    unique = {}
    for row in proposals:
        key = (row['mark'], row['alpha'], tuple(map(tuple, row['marker_log_odds'])),
               tuple(row['proposal_centre']))
        unique.setdefault(key, row)
    bank = [SwapPositiveWitness.create(w['mark'], w['alpha'], w['marker_log_odds'],
                                        w['proposal_centre']) for w in unique.values()]
    bank.extend(SwapPositiveWitness.create(row['positive']['mark'], row['positive']['alpha'],
                                           row['positive']['marker_log_odds'], row['profile'])
                for row in frontier)
    initial_bank_size = len(bank)
    pending, serial = [], 0
    def enqueue(box):
        nonlocal serial
        if profile_count(box):
            heapq.heappush(pending, (-volume(box), serial, box))
            serial += 1
    for row in saved['accepted_boxes']+saved['unresolved_boxes']:
        (low, high), a, b = row['ranges']
        box = ((max(first, low), min(last, high)), tuple(a), tuple(b))
        if box[0][0] <= box[0][1]:
            enqueue(box)
    accepted, unresolved = [], []
    splits = candidates = inspected = 0
    while pending:
        _, _, box = heapq.heappop(pending)
        witness_id, score, worst = best_witness(bank, box)
        if score > threshold and candidates < max_candidates:
            centre = tuple(np.mean(box_geometry(box)[0], axis=0))
            r = round(centre[0])
            a, b = (min(r//2, max(0, r-4096, round(value))) for value in centre[1:])
            if a+b == r and r < 8192:
                b = max(0, b-1)
            if not centre_boxes:
                point_box = tuple((value, value) for value in (r, a, b))
                _, point_score, _ = best_witness(bank, point_box)
            if centre_boxes:
                source = bank[witness_id]
                sr, sa, sb = source.centre
                translated = []
                for index, (full, old_full) in enumerate(zip((a, b), (sa, sb))):
                    empty, singleton = max(1, 4096-r+full), r-2*full
                    old_empty, old_singleton = max(1, 4096-sr+old_full), sr-2*old_full
                    target_counts, old_counts = (singleton, full), (old_singleton, old_full)
                    row = []
                    for axis, (target_count, old_count) in enumerate(zip(target_counts, old_counts)):
                        offset = (source.odds[index][axis]-log2(old_count/old_empty)
                                  if old_count > 0 else (.5, .9)[axis]*source.alpha/.1)
                        row.append(log2(max(1, target_count)/empty)+offset)
                    translated.append(row)
                trials = [(source.z, source.alpha, translated)]
                for alpha in dict.fromkeys((source.alpha, .1, .25, .5)):
                    odds = []
                    for full in (a, b):
                        empty, singleton = max(1, 4096-r+full), r-2*full
                        odds.append((log2(max(1, singleton)/empty)+.5*alpha/.1,
                                     log2(max(1, full)/empty)+.9*alpha/.1))
                    trials.append((2**(-r/1000), alpha, odds))
                for mark, alpha, odds in trials:
                    if candidates >= max_candidates:
                        break
                    witness = SwapPositiveWitness.create(mark, alpha, odds, (r, a, b))
                    candidates += 1
                    value, vertex = witness.evaluate(box)
                    if value < score:
                        bank.append(witness)
                        witness_id, score, worst = len(bank)-1, value, vertex
                    if score <= threshold:
                        break
            elif point_score > threshold:
                for alpha in (.1, .25, .5):
                    if candidates >= max_candidates:
                        break
                    odds = []
                    for full in (a, b):
                        empty, singleton = 4096-r+full, r-2*full
                        odds.append((log2(max(1, singleton)/max(1, empty))+.5*alpha/.1,
                                     log2(max(1, full)/max(1, empty))+.9*alpha/.1))
                    witness = SwapPositiveWitness.create(2**(-r/1000), alpha, odds, (r, a, b))
                    candidates += 1
                    value, _ = witness.evaluate(point_box)
                    if value < point_score:
                        bank.append(witness)
                        point_score = value
                    if point_score <= threshold:
                        break
                witness_id, score, worst = best_witness(bank, box)
        if score <= threshold:
            accepted.append(dict(ranges=box, witness_id=witness_id,
                                 union_log2_upper_screen=score, worst_corner=worst,
                                 integer_triples=profile_count(box), inherited=False))
        elif splits >= max_splits or all(lo == hi for lo, hi in box):
            unresolved.append(dict(ranges=box, reason='split-budget' if splits >= max_splits else 'fixed-point',
                                   best_log2_upper_screen=score))
        else:
            for child in split_box(box):
                enqueue(child)
            splits += 1
        inspected += 1
        if progress and inspected % 256 == 0:
            progress(dict(stage='swap-subdivide', inspected=inspected, splits=splits,
                          candidates=candidates, accepted=len(accepted), pending=len(pending)))
    covered = sum(row['integer_triples'] for row in accepted)
    domain = profile_count(((first, last), (0, 4096), (0, 4096)))
    maximum = max((row['union_log2_upper_screen'] for row in accepted), default=-float('inf'))
    union = maximum+log2(sum(2**(row['union_log2_upper_screen']-maximum) for row in accepted)) if accepted else None
    result = dict(profile_coordinates=['r', 'a', 'b'], all_even_odd_splits=True,
                  requested_weight_range=[first, last], proposal_budget=max_candidates,
                  proposals_used=candidates, threshold=threshold, slab_width=256,
                  full_width=256, minimum_width=1, accepted_boxes=accepted, unresolved_boxes=unresolved,
                  witnesses=[w.serializable() for w in bank], accepted_integer_triples=covered,
                  feasible_domain_triples=domain, accepted_domain_fraction=covered/domain,
                  union_log2_upper_screen=union, directed=False, whole_code_certificate=False,
                  subdivision=dict(max_splits=max_splits, splits=splits, inspected=inspected,
                                   max_candidates=max_candidates, candidates=candidates,
                                   initial_bank_size=initial_bank_size,
                                   whole_box_centering=centre_boxes, local_optimizations=0))
    return restore_identity(compact_plan(result))


if __name__ == '__main__':
    import argparse, json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--frontier', type=Path, default=Path(__file__).with_name('hd_shifted_swap_frontier.json'))
    parser.add_argument('--first', type=int, default=3072)
    parser.add_argument('--last', type=int, default=5119)
    parser.add_argument('--proposals', type=int, default=300)
    parser.add_argument('--max-evaluations', type=int, default=6000)
    parser.add_argument('--marker-maxfev', type=int, default=80)
    parser.add_argument('--weight-maxfev', type=int, default=16)
    parser.add_argument('--subdivide-plan', type=Path, help='Refine a saved swap partition with a frozen bank')
    parser.add_argument('--bank-plan', type=Path, action='append', default=[],
                        help='Add swap witnesses only from this frozen plan; keep the source partition')
    parser.add_argument('--max-splits', type=int, default=4096)
    parser.add_argument('--max-candidates', type=int, default=2000)
    parser.add_argument('--output', type=Path, help='Create a new JSON result artifact; refuse to overwrite it')
    parser.add_argument('--centre-boxes', action='store_true', help='Try whole-box centered markers even when its centre passes')
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output already exists; frozen plans are never overwritten')
    if args.bank_plan and not args.subdivide_plan:
        parser.error('--bank-plan requires --subdivide-plan')
    with args.frontier.open(encoding='utf-8') as source:
        saved = json.load(source)
    if saved.get('construction') != CONSTRUCTION:
        parser.error('The frontier must explicitly name the swap construction')
    progress = lambda row: print(json.dumps(row), flush=True)
    if args.subdivide_plan:
        with args.subdivide_plan.open(encoding='utf-8') as source:
            previous = json.load(source)
        bank_plans = []
        for path in args.bank_plan:
            with path.open(encoding='utf-8') as source:
                bank_plans.append(json.load(source))
        result = subdivide(previous, saved['results'], bank_plans=bank_plans,
                           first=args.first, last=args.last,
                           max_splits=args.max_splits, max_candidates=args.max_candidates,
                           centre_boxes=args.centre_boxes, progress=progress)
    else:
        result = bulk(saved['results'], first=args.first, last=args.last, proposals=args.proposals,
                      max_evaluations=args.max_evaluations, marker_maxfev=args.marker_maxfev,
                      weight_maxfev=args.weight_maxfev, progress=progress)
    if args.output:
        # This is generated numerical output, not a source-code edit. An
        # exclusive create protects every frozen plan against replacement.
        with args.output.open('x', encoding='utf-8') as destination:
            json.dump(result, destination, separators=(',', ':'), allow_nan=False)
            destination.write('\n')
        print(json.dumps(dict(stage='swap-bulk-saved', path=str(args.output.resolve()),
                              accepted_boxes=len(result['accepted_boxes']),
                              unresolved_boxes=len(result['unresolved_boxes']),
                              accepted_integer_triples=result['accepted_integer_triples'],
                              feasible_domain_triples=result['feasible_domain_triples'],
                              union_log2_upper_screen=result['union_log2_upper_screen'])), flush=True)
    else:
        print(json.dumps(dict(stage='swap-bulk-result', result=result), allow_nan=False), flush=True)
