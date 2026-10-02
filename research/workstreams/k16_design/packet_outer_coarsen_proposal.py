"""Reduce an exploratory tilt grid to a cheaper fresh-replay schedule.

Only floating proposal scores are used here. The result remains proposal-only:
the whole replay accepts its tilts, discards all supplied numeric scores, and
recomputes every occupancy with outward arithmetic. The interval cover minimizes
a simple cost proxy for local preparation plus quadratic coefficient powering.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import isfinite
from pathlib import Path


def select(trials, occupancies, minimum_margin=90.0):
    qs = tuple(occupancies)
    if not qs or qs != tuple(range(qs[0], qs[-1]+1)) or qs[0] < 1:
        raise ValueError('a nonempty consecutive positive occupancy interval is required')
    if not isfinite(minimum_margin):
        raise ValueError('finite proposal margin required')
    available = {}
    for trial in trials:
        tilt = Q(trial['tilt'])
        if tilt <= 0:
            raise ValueError('positive rational tilt required')
        scores = available.setdefault(tilt, {})
        for q in qs:
            witness = trial['witnesses'].get(str(q))
            if witness is not None:
                score = float(witness['estimated_margin_bits'])
                if isfinite(score) and score >= minimum_margin:
                    scores[q] = max(scores.get(q, float('-inf')), score)
    # Every edge covers a consecutive interval with a single qualifying tilt.
    # Dynamic programming chooses an exact cover; it never infers values at
    # an occupancy absent from that trial.
    cost, previous = [0.0]+[float('inf')]*len(qs), [None]*(len(qs)+1)
    for i, q in enumerate(qs):
        if not isfinite(cost[i]):
            continue
        for tilt, scores in available.items():
            for j in range(i, len(qs)):
                if qs[j] not in scores:
                    break
                candidate = cost[i] + 10000 + qs[j]*qs[j]
                if candidate < cost[j+1]:
                    cost[j+1], previous[j+1] = candidate, (i, tilt)
    if previous[-1] is None:
        missing = [q for q in qs if not any(q in scores for scores in available.values())]
        raise ValueError(f'no qualifying proposal cover; individually missing occupancies: {missing}')
    chosen, cursor = set(), len(qs)
    while cursor:
        cursor, tilt = previous[cursor]
        chosen.add(tilt)
    best = {}
    for q in qs:
        tilt = max((t for t in chosen if q in available[t]), key=lambda t: available[t][q])
        best[str(q)] = dict(margin=available[tilt][q], tilt=str(tilt))
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('proposal')
    parser.add_argument('--output', required=True)
    parser.add_argument('--minimum-margin', type=float, default=90.0)
    args = parser.parse_args()
    raw = Path(args.proposal).read_bytes()
    data = json.loads(raw)
    if data.get('proposal_only') is not True or data.get('whole_code_certificate') is not False:
        raise ValueError('an explicitly noncertifying proposal record is required')
    occupancies = sorted(map(int, data['best']))
    best = select(data['trials'], occupancies, args.minimum_margin)
    result = {key: data[key] for key in ('K', 'state_bits', 'envelope')}
    result.update(proposal_only=True, whole_code_certificate=False, best=best,
        minimum_proposal_margin=args.minimum_margin,
        source_proposal=dict(path=str(Path(args.proposal).resolve()),
            sha256=hashlib.sha256(raw).hexdigest()),
        scope='Rational tilt choices only; all scores must be discarded by the outward replay.')
    with Path(args.output).open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(occupancies=len(best),
        tilts=len({entry['tilt'] for entry in best.values()}),
        weakest=min(entry['margin'] for entry in best.values()))))


if __name__ == '__main__':
    main()
