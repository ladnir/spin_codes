"""Merge search fragments against one immutable canonical GL32 snapshot.

This checks provenance and exact coverage, not numerical bounds. Every
original pending root must be assigned exactly once. A fragment may leave
children unresolved. The ordinary dense output retains witnesses only as
search candidates; fresh dense and whole-code replay remain mandatory.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import dense_cover
import search_strategy as geometry


PARTIAL_SCHEMA = 'packed-canonical-gl32-partial-search-1'
SCOPE_FIELDS = (
    'ensemble', 'K', 'N', 'updates', 'block_width', 'distance', 'threshold',
    'minimum_groups', 'maximum_groups', 'root', 'comparison',
    'expected_cdf_sha256', 'comparison_caps_sha256', 'last_lp', 'refined',
    'base_tilt', 'variance_bins', 'regional_count', 'outer_premises', 'mixture',
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _exact(value):
    require(type(value) in (str, int), 'exact rational string or integer required')
    return Q(value)


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def _hash(value):
    return (type(value) is str and len(value) == 64
            and not set(value)-set('0123456789abcdef'))


def _source(source):
    require(isinstance(source, dict) and set(source) == {'path', 'sha256'}
            and type(source['path']) is str and Path(source['path']).is_absolute()
            and _hash(source['sha256']), 'absolute source path and SHA256 required')


def read_source(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    record = json.loads(raw)
    require(isinstance(record, dict), 'JSON object receipt required')
    return record, dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest())


def snapshot_scope(snapshot):
    """Validate the source interface, without authenticating its old bounds."""
    require(isinstance(snapshot, dict), 'dense snapshot required')
    dense_cover.validate_record(snapshot)
    require(all(key in snapshot for key in SCOPE_FIELDS), 'complete canonical scope required')
    require(snapshot['minimum_groups'] == 33 and snapshot['maximum_groups'] == 2048,
            'dense occupancy scope must be exactly33..2048')
    require(isinstance(snapshot['root'], list) and len(snapshot['root']) == 2,
            'two exact global root endpoints required')
    tuple(map(_exact, snapshot['root']))
    require(_exact(snapshot['distance']) > 0 and _exact(snapshot['base_tilt']) == Q(3, 16),
            'exact canonical claim and base tilt required')
    for key in ('expected_cdf_sha256', 'comparison_caps_sha256'):
        require(_hash(snapshot[key]), 'comparison SHA256 required')
    require(isinstance(snapshot['outer_premises'], dict) and bool(snapshot['outer_premises']),
            'outer count premises required')
    mixture = snapshot['mixture']
    require(isinstance(mixture, list) and bool(mixture), 'positive rational mixture required')
    centers = []
    for row in mixture:
        require(isinstance(row, dict) and set(row) == {'mass', 'activity'}, 'exact mixture row required')
        mass, activity = _exact(row['mass']), _exact(row['activity'])
        require(mass > 0 and 0 < activity <= 1, 'positive mixture mass and activity required')
        centers.append(activity)
    require(len(set(centers)) == len(centers), 'distinct mixture centers required')
    require(type(snapshot.get('precision')) is int and snapshot['precision'] >= 128
            and type(snapshot.get('target_bits')) is int and snapshot['target_bits'] >= 1,
            'snapshot precision and search target required')
    cover = snapshot['cover']
    geometry.partition(snapshot['root'], cover['leaves'], cover['unresolved'])
    require(all(isinstance(row.get('witness'), dict) for row in cover['leaves'].values()),
            'retained accepted rows require witnesses')
    return {key: copy.deepcopy(snapshot[key]) for key in SCOPE_FIELDS}


def _visits(cover):
    value = cover.get('visited', 0)
    require(type(value) is int and value >= 0, 'nonnegative integer work count required')
    return value


def fragment_cover(snapshot, source, partial):
    """Check each assigned subtree and return sanitized global-path rows."""
    require(isinstance(partial, dict) and partial.get('schema') == PARTIAL_SCHEMA,
            'canonical partial-search schema required')
    require(partial.get('source') == source, 'partial source path or hash differs')
    require(all(key in partial for key in SCOPE_FIELDS)
            and _canonical({key: partial[key] for key in SCOPE_FIELDS})
                == _canonical({key: snapshot[key] for key in SCOPE_FIELDS}),
            'partial construction, claim, premises, or mixture differs')
    require(type(partial.get('precision')) is int and partial['precision'] == snapshot['precision']
            and type(partial.get('target_bits')) is int and partial['target_bits'] == snapshot['target_bits'],
            'partial precision or search target differs')
    assigned = partial.get('assigned_roots')
    require(isinstance(assigned, list) and bool(assigned)
            and all(type(path) is str for path in assigned)
            and len(set(assigned)) == len(assigned)
            and set(assigned) <= set(snapshot['cover']['unresolved']),
            'distinct original pending roots must be assigned')
    cover = partial.get('cover')
    require(isinstance(cover, dict) and isinstance(cover.get('leaves'), dict)
            and isinstance(cover.get('unresolved'), dict)
            and not set(cover['leaves']) & set(cover['unresolved']),
            'disjoint accepted and unresolved fragment maps required')
    local = {root: ({}, {}) for root in assigned}
    cleaned = ({}, {})
    for index, name in enumerate(('leaves', 'unresolved')):
        for path, row in cover[name].items():
            require(type(path) is str and isinstance(row, dict), 'global binary path and row required')
            cell = geometry.path_cell(snapshot['root'], path)
            owners = [root for root in assigned if path.startswith(root)]
            require(len(owners) == 1, 'fragment path escapes or overlaps its assignments')
            require(isinstance(row.get('cell'), list) and len(row['cell']) == 2
                    and tuple(map(_exact, row['cell'])) == cell, 'fragment cell geometry differs')
            clean = dict(cell=list(map(str, cell)))
            if index == 0:
                require(isinstance(row.get('witness'), dict), 'accepted fragment row requires a witness')
                clean['witness'] = copy.deepcopy(row['witness'])
            root = owners[0]
            local[root][index][path[len(root):]] = clean
            cleaned[index][path] = clean
    for root, (leaves, unresolved) in local.items():
        geometry.partition(geometry.path_cell(snapshot['root'], root), leaves, unresolved)
    return assigned, cleaned[0], cleaned[1], _visits(cover)


def merge_records(snapshot, source, partials):
    """Merge (partial_record, parsed-byte provenance) pairs; never sum bounds."""
    scope = snapshot_scope(snapshot)
    _source(source)
    leaves, unresolved = copy.deepcopy(snapshot['cover']['leaves']), {}
    assigned, provenance, seen_files = set(), [], {source['path']}
    visits = _visits(snapshot['cover'])
    for partial, metadata in partials:
        _source(metadata)
        require(metadata['path'] not in seen_files, 'duplicate partial source file')
        seen_files.add(metadata['path'])
        roots, new_leaves, new_unresolved, work = fragment_cover(snapshot, source, partial)
        require(not assigned.intersection(roots), 'pending root assigned more than once')
        require(not (set(leaves) | set(unresolved)) & (set(new_leaves) | set(new_unresolved)),
                'duplicate global frontier path')
        assigned.update(roots)
        leaves.update(new_leaves)
        unresolved.update(new_unresolved)
        visits += work
        provenance.append(dict(metadata, assigned_roots=list(roots)))
    require(assigned == set(snapshot['cover']['unresolved']), 'not all original pending roots were assigned')
    geometry.partition(snapshot['root'], leaves, unresolved)
    result = dict(scope, schema=dense_cover.SCHEMA, precision=snapshot['precision'],
        target_bits=snapshot['target_bits'], cover=dict(leaves=leaves, unresolved=unresolved, visited=visits),
        source=copy.deepcopy(source), partial_sources=provenance,
        merge=dict(method='exact assigned-subtree replacement', assigned_roots=sorted(assigned),
            accepted_retained=len(snapshot['cover']['leaves']), numerical_evaluations=0),
        final_replay_required=True,
        proof_status='Merged search witnesses only; no bounds authenticated by merging. '
            'Fresh uncached dense replay and sparse/whole-code assembly are mandatory.')
    dense_cover.validate_record(result)
    return result


def run(source_path, partial_paths, output_path):
    output_path = Path(output_path).resolve()
    require(not output_path.exists(), 'a new output file is required')
    snapshot, source = read_source(source_path)
    partials = [read_source(path) for path in partial_paths]
    result = merge_records(snapshot, source, partials)
    for metadata in [source, *(metadata for _, metadata in partials)]:
        require(hashlib.sha256(Path(metadata['path']).read_bytes()).hexdigest() == metadata['sha256'],
                'an input receipt changed during merging')
    encoded = json.dumps(result, indent=2, allow_nan=False)+'\n'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('x', encoding='utf-8') as stream:
        stream.write(encoded)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--partials', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.source, args.partials, args.output)
    print('MERGED SEARCH ONLY:', len(result['cover']['leaves']), 'witness leaves;',
        len(result['cover']['unresolved']), 'unresolved; final fresh replay required.')


if __name__ == '__main__':
    main()
