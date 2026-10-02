"""Independent search of assigned pending subtrees, not a full certificate.

Fresh outer authentication and the ordinary global-root model are retained.
Only the assigned unresolved frontier paths are searched; unrelated accepted
leaves are neither imported as numerical bounds nor replayed. A separate
merger and final independent whole-code replay remain necessary.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import heapq
import importlib.util
import json
import math
from pathlib import Path
import sys

from flint import arb, ctx
import search_strategy as geometry

# The BCH authenticator imports a DIFFERENT legacy dense_cover.Model. Never
# register the packed CLI interface under that shared global module name.
_import_path = sys.path[:]
try:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    _dense_spec = importlib.util.spec_from_file_location(
        'packed_gl32_dense_interface', Path(__file__).resolve().with_name('dense_cover.py'))
    dense = importlib.util.module_from_spec(_dense_spec)
    sys.modules[_dense_spec.name] = dense
    _dense_spec.loader.exec_module(dense)
finally:
    sys.path[:] = _import_path
    del _import_path

SCHEMA = 'packed-canonical-gl32-partial-search-1'
SCOPE_FIELDS = ('ensemble', 'K', 'N', 'updates', 'block_width', 'distance', 'threshold',
    'minimum_groups', 'maximum_groups', 'root', 'comparison', 'expected_cdf_sha256',
    'comparison_caps_sha256', 'last_lp', 'refined', 'base_tilt', 'variance_bins',
    'regional_count', 'outer_premises', 'mixture')


def assigned_paths(source, paths=None, shard_index=None, shard_count=None):
    """Choose exact pending frontier roots; sharding is contiguous by path."""
    cover = dense.validate_record(source)
    geometry.partition(source['root'], cover['leaves'], cover['unresolved'])
    pending = sorted(cover['unresolved'])
    if paths is not None:
        if shard_index is not None or shard_count is not None:
            raise ValueError('choose explicit paths or one deterministic shard')
        selected = list(paths)
    else:
        if (type(shard_count) is not int or not 1 <= shard_count <= len(pending)
                or type(shard_index) is not int or not 0 <= shard_index < shard_count):
            raise ValueError('nonempty shard count and zero-based shard index required')
        selected = pending[len(pending)*shard_index//shard_count:
            len(pending)*(shard_index+1)//shard_count]
    if not selected or len(set(selected)) != len(selected) or any(p not in cover['unresolved'] for p in selected):
        raise ValueError('distinct original unresolved frontier paths required')
    return sorted(selected)


def validate_cover(root, assigned, cover):
    """Check exact coverage of each assigned subtree in global coordinates."""
    if (not isinstance(cover, dict) or not isinstance(cover.get('leaves'), dict)
            or not isinstance(cover.get('unresolved'), dict)
            or set(cover['leaves']) & set(cover['unresolved'])
            or not assigned or len(set(assigned)) != len(assigned)):
        raise ValueError('disjoint partial maps and distinct assigned roots required')
    roots = {p: geometry.path_cell(root, p) for p in assigned}
    for a in assigned:
        if any(a != b and a.startswith(b) for b in assigned):
            raise ValueError('assigned roots must be disjoint')
    all_rows = dict(cover['leaves'], **cover['unresolved'])
    for path, row in all_rows.items():
        if sum(path.startswith(p) for p in assigned) != 1:
            raise ValueError('partial path lies outside its unique assigned subtree')
        if (not isinstance(row, dict) or 'cell' not in row
                or tuple(map(Q, row['cell'])) != geometry.path_cell(root, path)):
            raise ValueError('every partial row must retain exact global coordinates')
    for prefix, cell in roots.items():
        leaves = {p[len(prefix):]: row for p, row in cover['leaves'].items() if p.startswith(prefix)}
        unresolved = {p[len(prefix):]: row for p, row in cover['unresolved'].items() if p.startswith(prefix)}
        geometry.partition(cell, leaves, unresolved)


def search(model, assigned, *, precision=256, target_bits=36, max_cells=200,
           max_depth=22, checkpoint_every=10, checkpoint=None):
    """Search assigned global paths without replaying other parts of the tree."""
    if (type(precision) is not int or precision < 128
            or type(target_bits) is not int or target_bits < 20
            or type(max_cells) is not int or max_cells < 0
            or type(max_depth) is not int or not 1 <= max_depth <= 64
            or type(checkpoint_every) is not int or checkpoint_every < 1):
        raise ValueError('valid integer precision, target, and work limits required')
    root = tuple(map(Q, model.root))
    pending = [(len(p), p, geometry.path_cell(root, p)) for p in assigned]
    if any(depth > max_depth for depth, _, _ in pending):
        raise ValueError('assigned root exceeds the depth limit')
    heapq.heapify(pending)
    leaves, unresolved, visited = {}, {}, 0

    def state():
        result = dict(leaves=copy.deepcopy(leaves), unresolved={**copy.deepcopy(unresolved),
            **{p: dict(cell=list(map(str, cell))) for _, p, cell in pending}}, visited=visited)
        validate_cover(root, assigned, result)
        return result

    initial = state()
    if checkpoint:
        checkpoint(initial)
    while pending and visited < max_cells:
        depth, path, cell = heapq.heappop(pending)
        ctx.prec = precision
        score, witness = model.proposal(cell)
        if not math.isfinite(score) or not isinstance(witness, dict):
            raise ArithmeticError('finite proposal and rational witness required')
        accepted = False
        if score < -(target_bits+2):
            ctx.prec = precision
            upper = model.outward(cell, witness)
            if ctx.prec != precision or not upper.is_finite() or not upper > 0:
                raise ArithmeticError('finite positive outward bound at requested precision required')
            if upper < arb(2)**-target_bits:
                leaves[path] = dict(cell=list(map(str, cell)), witness=copy.deepcopy(witness),
                    proposal=float(score), upper=[int(v) for v in upper.upper().man_exp()])
                accepted = True
        if not accepted:
            if depth == max_depth:
                unresolved[path] = dict(cell=list(map(str, cell)), witness=copy.deepcopy(witness),
                    proposal=float(score))
            else:
                mid = sum(cell)/2
                heapq.heappush(pending, (depth+1, path+'0', (cell[0], mid)))
                heapq.heappush(pending, (depth+1, path+'1', (mid, cell[1])))
        visited += 1
        if tuple(map(Q, model.root)) != root:
            raise ArithmeticError('global model domain changed during partial search')
        if checkpoint and visited % checkpoint_every == 0:
            checkpoint(state())
    result = state()
    if checkpoint:
        checkpoint(result)
    return result


def fresh_model(source, precision, target_bits):
    """Authenticate the source mixture, then prepare ordinary search wrappers."""
    dense.validate_record(source)
    parent = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(parent))
    sys.path.insert(0, str(parent/'gf16_packets'))
    from canonical_counts import authenticated_bch_cdf
    caps, premises = authenticated_bch_cdf()
    if dense.fingerprint(caps) != source['expected_cdf_sha256'] or premises != source['outer_premises']:
        raise ValueError('fresh BCH/GL32 premises differ from the source')
    comparison = caps
    if source['comparison'] == 'direct-expected-shell-majorant':
        from monotone import transport_shells
        from local_models import full_block
        comparison = transport_shells(premises['canonical_cdf'], full_block(8))
    if dense.fingerprint(comparison) != source['comparison_caps_sha256']:
        raise ValueError('fresh shell comparison differs from the source')
    mixture, _ = dense.exact_mixture(comparison, source['mixture'])
    import birth_classes
    import scalar_cover as sc
    import shared_mixture
    data = birth_classes.actual(2)
    ctx.prec = precision
    model = sc.Model(shared_mixture.as_components(mixture), data, source['threshold'],
        source['minimum_groups'], Q(3, 16), inner=birth_classes,
        variance_shuffle=True, variance_bins=16, regional_count=True)
    if model.root != tuple(map(Q, source['root'])):
        raise ValueError('fresh global root differs from the source')
    model.proposal_stop_bits = target_bits+2
    if source.get('point_hint_sources'):
        from point_reuse import load_points, PointReuseModel
        from region_cache import RegionalCache
        hints = source['point_hint_sources']
        probes, actual = load_points([row['path'] for row in hints])
        if [row['sha256'] for row in actual] != [row['sha256'] for row in hints]:
            raise ValueError('point hint source changed')
        cache = RegionalCache()
        base = model
        model = PointReuseModel(base, probes, target_bits,
            outward=lambda cell, witness: cache.bound(base, cell, witness),
            regional_warm=lambda witness: cache.is_warm(base, witness),
            maximum_optimization_width=Q(1, 1024))
    return model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--paths', nargs='+')
    parser.add_argument('--shard-index', type=int)
    parser.add_argument('--shard-count', type=int)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--target-bits', type=int, default=36)
    parser.add_argument('--max-cells', type=int, default=200)
    parser.add_argument('--max-depth', type=int, default=22)
    parser.add_argument('--checkpoint-every', type=int, default=10)
    args = parser.parse_args()
    temporary = args.output.with_name(args.output.name+'.tmp')
    if (args.output.exists() or temporary.exists() or args.precision < 128 or args.target_bits < 20
            or args.max_cells < 0 or not 1 <= args.max_depth <= 64 or args.checkpoint_every < 1):
        parser.error('new output path and valid proof/work limits required')
    raw = args.source.read_bytes()
    source = json.loads(raw)
    if args.precision != source.get('precision') or args.target_bits != source.get('target_bits'):
        parser.error('partial precision and target must match the immutable source snapshot')
    assigned = assigned_paths(source, args.paths, args.shard_index, args.shard_count)
    record = {key: copy.deepcopy(source[key]) for key in SCOPE_FIELDS}
    record.update(schema=SCHEMA, source=dict(path=str(args.source.resolve()),
        sha256=hashlib.sha256(raw).hexdigest()), assigned_roots=assigned,
        precision=args.precision, target_bits=args.target_bits,
        max_depth=args.max_depth, proof_status='Partial subtree search only; no full certificate.')
    model = fresh_model(source, args.precision, args.target_bits)

    def save(cover):
        if args.source.read_bytes() != raw:
            raise ValueError('source snapshot changed during partial search')
        record['cover'] = cover
        record['complete_assigned'] = not cover['unresolved']
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(json.dumps(record, indent=2)+'\n')
        temporary.replace(args.output)
        print('PACKED PARTIAL visited', cover['visited'], 'accepted', len(cover['leaves']),
            'pending', len(cover['unresolved']), flush=True)

    search(model, assigned, precision=args.precision, target_bits=args.target_bits,
        max_cells=args.max_cells, max_depth=args.max_depth,
        checkpoint_every=args.checkpoint_every, checkpoint=save)


if __name__ == '__main__':
    main()
