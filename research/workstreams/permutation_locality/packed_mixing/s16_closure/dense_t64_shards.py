"""Search and merge disjoint t64 dense subtrees; never certify saved bounds.

Workers freshly authenticate their model and bind one immutable search
snapshot. Merging checks exact global coverage and retains witnesses only.
Every merged leaf still requires the ordinary fresh whole-code replay.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import dense_t64 as dense

search = dense.prior.cell_search
geometry = dense.geometry
SCHEMA = 'packed-gl32-t64-s16-dense-shard-1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON object key: '+key)
        result[key] = value
    return result


def read(path):
    path = Path(path).resolve(); raw = path.read_bytes()
    record = json.loads(raw, object_pairs_hook=_object)
    require(isinstance(record, dict), 'JSON object required')
    return record, dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest())


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def metadata(value):
    require(isinstance(value, dict) and set(value) == {'path', 'sha256'}
        and isinstance(value['path'], str) and Path(value['path']).is_absolute()
        and isinstance(value['sha256'], str) and len(value['sha256']) == 64
        and not set(value['sha256'])-set('0123456789abcdef'),
        'absolute source path and SHA256 required')
    return value


def unchanged(value):
    metadata(value)
    dense.unchanged(value)


def snapshot(record):
    cells = dense.validate_cover_record(record)
    require(type(record.get('precision')) is int and record['precision'] >= 128
        and type(record.get('cell_target_bits')) is int and 20 <= record['cell_target_bits'] <= 256,
        'snapshot precision and cell target required')
    metadata(record.get('source'))
    _visits(record['cover'])
    return cells


def _visits(cover):
    value = cover.get('visited', 0)
    require(type(value) is int and value >= 0, 'nonnegative integer visit count required')
    return value


def assigned_paths(record, shard_index, shard_count):
    snapshot(record)
    pending = sorted(record['cover']['unresolved'])
    require(type(shard_count) is int and 1 <= shard_count <= len(pending)
        and type(shard_index) is int and 0 <= shard_index < shard_count,
        'nonempty shard count and zero-based index required')
    return pending[shard_index::shard_count]


def _fresh_output(path):
    path = Path(path).resolve()
    temporary = path.with_name(path.name+'.tmp')
    require(not path.exists() and not temporary.exists(), 'new output and temporary paths required')
    return path, temporary


def _write(path, temporary, record):
    encoded = json.dumps(record, indent=2, allow_nan=False)+'\n'
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_text(encoded, encoding='utf-8')
    temporary.replace(path)


def run_shard(source_path, shard_index, shard_count, output, *,
              max_cells=200, max_depth=24, checkpoint_every=5):
    output, temporary = _fresh_output(output)
    old, source = read(source_path)
    assigned = assigned_paths(old, shard_index, shard_count)
    scope = old['scope']; proposal = old['source']
    unchanged(proposal)
    model, fresh_scope, fresh_source = dense.fresh_model(proposal['path'],
        precision=old['precision'], variance_bins=scope['variance_bins'],
        distance=Q(scope['distance']), birth_density=scope['birth_density'])
    require(canonical(fresh_scope) == canonical(scope) and fresh_source == proposal,
        'fresh construction, maps, count premises or mixture source differ from snapshot')
    require(tuple(map(Q, model.root)) == tuple(map(Q, scope['root'])), 'fresh global root differs')
    model.proposal_stop_bits = old['cell_target_bits']+2
    record = dict(schema=SCHEMA, snapshot_source=source, source=copy.deepcopy(proposal),
        scope=copy.deepcopy(scope), precision=old['precision'], cell_target_bits=old['cell_target_bits'],
        assigned_roots=assigned, shard_index=shard_index, shard_count=shard_count,
        max_cells=max_cells, max_depth=max_depth, whole_code_certificate=False,
        final_replay_required=True, proof_status='Partial search witnesses only; no certificate.')
    def checkpoint(cover):
        search.validate_cover(scope['root'], assigned, cover)
        require(all(isinstance(row.get('witness'), dict) for row in cover['leaves'].values()),
            'every accepted cell needs a witness')
        unchanged(source); unchanged(proposal)
        record.update(cover=cover, complete_assigned=not cover['unresolved'])
        _write(output, temporary, record)
        print('T64 SHARD', shard_index, '/', shard_count, 'visited', cover['visited'],
              'accepted', len(cover['leaves']), 'pending', len(cover['unresolved']), flush=True)
    search.search(model, assigned, precision=old['precision'], target_bits=old['cell_target_bits'],
        max_cells=max_cells, max_depth=max_depth, checkpoint_every=checkpoint_every, checkpoint=checkpoint)
    unchanged(source); unchanged(proposal)
    return record


def _clean(root, rows, accepted):
    result = {}
    for path, row in rows.items():
        cell = geometry.path_cell(root, path)
        cleaned = dict(cell=list(map(str, cell)))
        if accepted:
            require(isinstance(row.get('witness'), dict), 'accepted search witness required')
            cleaned['witness'] = copy.deepcopy(row['witness'])
        result[path] = cleaned
    return result


def merge_records(old, source, partials):
    """Merge parsed fragments with provenance; never evaluate or sum bounds."""
    snapshot(old); metadata(source)
    root, pending = old['scope']['root'], set(old['cover']['unresolved'])
    leaves = _clean(root, old['cover']['leaves'], True)
    unresolved, assigned, provenance, seen_files = {}, set(), [], {source['path']}
    visits = _visits(old['cover'])
    for partial, origin in partials:
        metadata(origin)
        require(origin['path'] not in seen_files, 'duplicate shard source file')
        seen_files.add(origin['path'])
        require(partial.get('schema') == SCHEMA and partial.get('snapshot_source') == source
            and partial.get('source') == old['source']
            and canonical(partial.get('scope')) == canonical(old['scope'])
            and type(partial.get('precision')) is int and partial['precision'] == old['precision']
            and type(partial.get('cell_target_bits')) is int and partial['cell_target_bits'] == old['cell_target_bits'],
            'shard snapshot, construction, maps, precision or target differ')
        expected = assigned_paths(old, partial.get('shard_index'), partial.get('shard_count'))
        roots = partial.get('assigned_roots')
        require(roots == expected and not assigned.intersection(roots),
            'shard roots differ from deterministic assignment or overlap')
        cover = partial.get('cover')
        search.validate_cover(root, roots, cover)
        new_leaves = _clean(root, cover['leaves'], True)
        new_unresolved = _clean(root, cover['unresolved'], False)
        require(not (set(leaves) | set(unresolved)) & (set(new_leaves) | set(new_unresolved)),
            'duplicate global frontier path')
        assigned.update(roots); leaves.update(new_leaves); unresolved.update(new_unresolved)
        visits += _visits(cover)
        provenance.append(dict(origin, assigned_roots=list(roots)))
    require(assigned == pending, 'every original unresolved root must be assigned exactly once')
    geometry.partition(root, leaves, unresolved)
    result = dict(schema=dense.SCHEMA, scope=copy.deepcopy(old['scope']),
        source=copy.deepcopy(old['source']), precision=old['precision'],
        cell_target_bits=old['cell_target_bits'], whole_code_certificate=False,
        cover=dict(leaves=leaves, unresolved=unresolved, visited=visits),
        complete_search_partition=not unresolved, final_replay_required=True,
        snapshot_source=copy.deepcopy(source), partial_sources=provenance,
        merge=dict(method='exact disjoint pending-subtree replacement', assigned_roots=sorted(assigned),
            accepted_retained=len(old['cover']['leaves']), numerical_evaluations=0),
        proof_status='Merged search witnesses only; all leaves require fresh whole-code replay.')
    dense.validate_cover_record(result)
    return result


def run_merge(source_path, partial_paths, output):
    output, temporary = _fresh_output(output)
    old, source = read(source_path)
    unchanged(old['source'])
    partials = [read(path) for path in partial_paths]
    result = merge_records(old, source, partials)
    for origin in [source, old['source'], *(origin for _, origin in partials)]:
        unchanged(origin)
    _write(output, temporary, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    worker = sub.add_parser('search')
    worker.add_argument('--input', type=Path, required=True)
    worker.add_argument('--shard-index', type=int, required=True)
    worker.add_argument('--shard-count', type=int, required=True)
    worker.add_argument('--max-cells', type=int, default=200)
    worker.add_argument('--max-depth', type=int, default=24)
    worker.add_argument('--checkpoint-every', type=int, default=5)
    worker.add_argument('--output', type=Path, required=True)
    merger = sub.add_parser('merge')
    merger.add_argument('--input', type=Path, required=True)
    merger.add_argument('--partials', type=Path, nargs='+', required=True)
    merger.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.mode == 'search':
        run_shard(args.input, args.shard_index, args.shard_count, args.output,
            max_cells=args.max_cells, max_depth=args.max_depth, checkpoint_every=args.checkpoint_every)
    else:
        run_merge(args.input, args.partials, args.output)


if __name__ == '__main__':
    main()
