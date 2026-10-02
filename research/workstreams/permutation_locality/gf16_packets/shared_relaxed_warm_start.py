"""Seed a stronger-count search with old rational witnesses, not old bounds.

Require identical activity geometry and componentwise smaller masses. Every
old leaf becomes unresolved. The parallel search must freshly check it with
--reuse-parent-witness, and a full certificate still needs fresh assembly.
The optional append-points mode instead preserves a partial partition and
adds only rational proposals; continuation still rechecks its retained leaves.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import shared_relaxed_strategy as proof


def seed(source, comparison):
    old = proof.validate_record(source, 0, complete=True)
    new = proof.validate_record(comparison, 0)
    for key in ('updates', 'K', 'N', 'minimum_groups', 'base_tilt', 'variance_bins'):
        if source[key] != comparison[key]:
            raise ValueError('identical construction and activity geometry required')
    before = {Q(v['activity']): Q(v['mass']) for v in source['mixture']}
    after = {Q(v['activity']): Q(v['mass']) for v in comparison['mixture']}
    if (len(before) != len(source['mixture']) or len(after) != len(comparison['mixture'])
            or set(before) != set(after) or any(not 0 < after[p] <= before[p] for p in before)
            or tuple(map(Q, old['root'])) != tuple(map(Q, new['root']))):
        raise ValueError('unchanged activities/root and positive nonincreasing masses required')
    cells = proof.sc.partition(SimpleNamespace(root=tuple(map(Q, old['root']))), old['cover']['leaves'], {})
    result = copy.deepcopy(comparison)
    result['results'][0]['cover'] = dict(leaves={}, visited=0,
        unresolved={path: dict(cell=list(map(str, cell)),
            witness=copy.deepcopy(old['cover']['leaves'][path]['witness'])) for path, cell in cells.items()})
    result['results'][0].pop('rechecked_precision', None)
    result['proof_status'] = ('Warm start only: all cells unresolved. Old rational witnesses, '
        'not numerical bounds, were retained. Recheck every cell and regenerate full assembly.')
    return result


def append_points(source, comparison):
    """Add rational point proposals to an unchanged partial/full partition."""
    old = proof.validate_record(source)
    new = proof.validate_record(comparison)
    for key in ('updates', 'K', 'N', 'minimum_groups', 'base_tilt', 'variance_bins',
                'mixture', 'cap_sha256', 'count_witnesses', 'count_refinement_iterations'):
        if source.get(key) != comparison.get(key):
            raise ValueError('point proposals require the identical comparison')
    if old['threshold'] != new['threshold'] or old['root'] != new['root']:
        raise ValueError('point proposals require the identical cutoff and root')
    proof.sc.partition(SimpleNamespace(root=tuple(map(Q, old['root']))),
                       old['cover']['leaves'], old['cover']['unresolved'])
    result = copy.deepcopy(source)
    points = result['results'][0].setdefault('probes', [])
    means = {Q(p['mean']) for p in points}
    for point in new.get('probes', []):
        if not isinstance(point, dict) or not isinstance(point.get('witness'), dict):
            raise ValueError('a rational witness is required for every appended point')
        mean = Q(point['mean'])
        if not 0 < mean < 1 or mean in means:
            raise ValueError('distinct interior point means required')
        means.add(mean)
        points.append(dict(mean=str(mean), witness=copy.deepcopy(point['witness'])))
    result['proof_status'] = ('Point proposals added without importing numerical bounds. '
        'The partition is unchanged; continuation and complete assembly recheck every leaf.')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('comparison', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--append-points', action='store_true',
        help='Preserve a partition and add point proposals from an identical-comparison record')
    args = parser.parse_args()
    if args.output.exists(): parser.error('preserve existing evidence; choose a new output path')
    source, comparison = args.source.read_bytes(), args.comparison.read_bytes()
    result = (append_points if args.append_points else seed)(json.loads(source), json.loads(comparison))
    result['warm_start_sources'] = dict(source_sha256=hashlib.sha256(source).hexdigest(),
        comparison_sha256=hashlib.sha256(comparison).hexdigest())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('Prepared', len(result['results'][0]['cover']['unresolved']), 'unresolved cells;',
          len(result['results'][0].get('probes', [])), 'point proposals; no bounds transferred.')


if __name__ == '__main__': main()
