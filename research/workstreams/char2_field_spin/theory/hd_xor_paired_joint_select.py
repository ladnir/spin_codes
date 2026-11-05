"""Finite FLOAT witness selection for a saved, fixed r17 geometry partition.

Only cells weaker than the requested cutoff are revisited. The fixed finite
mark bank changes proof witnesses, not the construction. Output is a screen;
hd_xor_paired_joint_gate must independently replay all saved literal witnesses.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_paired_joint_gate import (GROUPS, CONSTRUCTION, LIMIT, Q, REGION_LENGTH,
    literal_schedules, ordered_tables, matching_joint_counts, geometry_stratum_counts,
    pair_regions_screen, MultiplicityCapsFloat, joint_log_count)

MARKS = tuple(F(n, 10**9) for n in
              (996248217, 996398217, 996448217, 996498217, 996548217, 996573217, 996598217))


def select(seed, cutoff=-60., progress=None, marks=MARKS, potentials=None):
    if (seed.get('construction') != CONSTRUCTION or seed.get('active_groups') != GROUPS
            or seed.get('active_rows_per_group') != 2 or seed.get('directed') is not False):
        raise ValueError('Expected the matching FLOAT r17 homogeneous screen')
    counts = matching_joint_counts(REGION_LENGTH, GROUPS)
    pairs = geometry_stratum_counts(GROUPS)
    schedules, complete = literal_schedules(set(counts), witness_rows=seed['entries'])
    if not complete:
        raise ValueError('The saved partition is incomplete')
    marks = tuple(F(mark) for mark in marks)
    if not marks or any(not 0 < mark < 1 for mark in marks) or not np.isfinite(cutoff):
        raise ValueError('A nonempty finite mark bank and finite cutoff are required')
    potentials = None if potentials is None else tuple(F(value) for value in potentials)
    if potentials is not None and (not potentials or any(value > 0 for value in potentials)):
        raise ValueError('A nonempty bank of nonpositive log potentials is required')
    started = perf_counter()
    raw, entries, evaluations = {}, [], 0
    by_source = {}
    for row in seed['entries']:
        by_source.setdefault(tuple(row['source_matching_counts']), []).append(row)
    for (a, b), rows in sorted(by_source.items()):
        evaluators = {}
        for previous in rows:
            low, high = previous['mu_range']
            old_logp = F(previous['log2_potential'])
            geometries = sum(pairs[mu] for mu in range(low, high+1))
            best = dict(previous)
            if previous['failure_log2_screen'] > cutoff:
                trials = []
                witness_bank = {(z, p) for z in marks for p in
                                ((old_logp,) if potentials is None else potentials)}
                witness_bank.add((F(previous['mark']), old_logp))
                for z, logp in sorted(witness_bank):
                    if z not in raw:
                        raw[z] = pair_regions_screen(float(z), GROUPS)
                    key = z, logp
                    if key not in evaluators:
                        evaluators[key] = MultiplicityCapsFloat(
                            ordered_tables(raw[z], a, b), float(logp), GROUPS)
                    value = (log2(counts[a, b])+log2(geometries)-log2(Q-1)
                             -LIMIT*log2(float(z))-float(logp)
                             +joint_log_count(evaluators[key](high), GROUPS))
                    trials.append(dict(mark=str(z), log2_potential=str(logp), failure_log2_screen=value))
                    evaluations += 1
                winner = min(trials, key=lambda x: x['failure_log2_screen'])
                best.update(winner)
                best['finite_trials'] = trials
                if progress:
                    progress(dict(source=[a, b], mu=[low, high], before=previous['failure_log2_screen'],
                                  after=winner['failure_log2_screen'], mark=winner['mark'],
                                  log2_potential=winner['log2_potential']))
            entries.append(best)
    result = dict(seed)
    result.update(entries=entries, directed=False, whole_code_certificate=False,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])),
        elapsed_seconds=perf_counter()-started, selector_cutoff=cutoff,
        selector_marks=list(map(str, marks)), selector_new_evaluations=evaluations)
    result['selector_log2_potentials'] = None if potentials is None else list(map(str, potentials))
    result.pop('source_sha256', None)
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=Path, required=True)
    parser.add_argument('--cutoff', type=float, default=-60.)
    parser.add_argument('--marks', nargs='+', type=F, default=MARKS)
    parser.add_argument('--potentials', nargs='+', type=F)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen artifacts are never overwritten')
    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    seed_bytes = args.seed.read_bytes()
    result = select(json.loads(seed_bytes), args.cutoff, marks=args.marks, potentials=args.potentials,
                    progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result.update(selector_sha256=source_hash, seed_sha256=sha256(seed_bytes).hexdigest())
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, allow_nan=False, indent=2)
        target.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()), bound=result['failure_log2_screen'],
                          evaluations=result['selector_new_evaluations'],
                          elapsed_seconds=result['elapsed_seconds'])), flush=True)
