"""Scoped point/interval witness hints for hill search, never final replay.

The caller first authenticates the current unwrapped model. Only rational
witness parameters and exact source geometry cross into this wrapper. All
hint bounds are evaluated anew on each target cell; no saved score or
endpoint is read. The regional cache lives only in this search process.
"""
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

_import_path = sys.path[:]
try:
    import hill_cover
    from point_reuse import PointReuseModel, WITNESS_FIELDS, rebase_witness
    from region_cache import RegionalCache, LOCAL_SELECTORS, COUNT_SELECTORS
finally:
    sys.path[:] = _import_path
    del _import_path


def canonical_scope(scope):
    hill_cover.validate_scope(scope)
    return {key: copy.deepcopy(scope[key])
            for key in ('schema', *hill_cover.cell_search.SCOPE_FIELDS)}


def _rational(value):
    if type(value) not in (str, int):
        raise ValueError('exact rational witness entries required')
    return Q(value)


def _vector(value, length=None):
    if not isinstance(value, list) or length is not None and len(value) != length:
        raise ValueError('rational witness vector has wrong shape')
    return [_rational(entry) for entry in value]


def _clean_witness(witness, source_cell):
    allowed = WITNESS_FIELDS | frozenset(LOCAL_SELECTORS+COUNT_SELECTORS)
    if (not isinstance(witness, dict) or set(witness)-allowed
            or not {'tilt', 'parameters', 'variance_dual'} <= set(witness)
            or _rational(witness['tilt']) != Q(3, 16)):
        raise ValueError('known base-tilt witness parameters required')
    lam, _, mu = _vector(witness['parameters'], 3)
    _vector(witness['variance_dual'], 2)
    if lam <= 0 or mu < 0:
        raise ValueError('positive output tilt and nonnegative occupancy dual required')
    for key in ('regional_exact_zero', 'regional_feedback_uniform_classes',
                'regional_feedback_uniform_replace', 'regional_direct_counts',
                'regional_tilted_atom', 'regional_fine_tilts', 'regional_tilted_variance'):
        if key in witness and type(witness[key]) is not bool:
            raise ValueError('boolean regional witness selector required')
    for key in ('regional_joint_return_through', 'regional_lazy_density_through',
                'regional_feedback_classes_from', 'regional_feedback_classes_through'):
        if key in witness and type(witness[key]) is not int:
            raise ValueError('integer regional witness cutoff required')
    if 'weights_dual' in witness:
        if not isinstance(witness['weights_dual'], list):
            raise ValueError('weight-dual list required')
        for row in witness['weights_dual']:
            _vector(row)
    for name in ('variance_partition', 'regional_count_parts'):
        if name in witness and not isinstance(witness[name], list):
            raise ValueError('witness partitions must be lists')
    for part in witness.get('variance_partition', []):
        if not isinstance(part, dict) or set(part) != {'interval', 'dual'}:
            raise ValueError('variance parts contain only intervals and rational duals')
        _vector(part['interval'], 2)
        _vector(part['dual'], 3)
    for part in witness.get('regional_count_parts', []):
        if not isinstance(part, dict) or set(part) != {'interval', 'dual', 'mgf_witnesses'}:
            raise ValueError('regional parts contain only intervals and rational MGF witnesses')
        if not isinstance(part['mgf_witnesses'], list):
            raise ValueError('regional MGF witness list required')
        for row in part['mgf_witnesses']:
            if (not isinstance(row, dict) or not {'tilt', 'dual'} <= set(row)
                    or set(row)-{'tilt', 'dual', 'tilted_atom', 'tilted_variance_dual'}):
                raise ValueError('known raw MGF witness parameter fields required')
            _rational(row['tilt'])
            _vector(row['dual'], 3)
            if 'tilted_variance_dual' in row:
                _vector(row['tilted_variance_dual'], 3)
            if 'tilted_atom' in row and type(row['tilted_atom']) is not bool:
                raise ValueError('boolean tilted-atom selector required')
    midpoint = sum(source_cell)/2
    if not 0 < midpoint < 1:
        raise ValueError('interior hint midpoint required')
    if 'variance_partition' in witness or 'regional_count_parts' in witness:
        clean = rebase_witness(source_cell, (midpoint, midpoint), witness)
    else:
        clean = copy.deepcopy(witness)
    return dict(mean=str(midpoint), witness=clean)


def proposals(record, scope):
    """Extract only parameter proposals from an exactly matching hill source."""
    expected = canonical_scope(scope)
    if not isinstance(record, dict) or canonical_scope(record.get('scope')) != expected:
        raise ValueError('hint source must match the complete current hill scope')
    if record.get('schema') == hill_cover.POINT_SCHEMA:
        rows = record.get('points')
        if not isinstance(rows, list) or not rows:
            raise ValueError('nonempty checked point list required')
        result = []
        for row in rows:
            if (not isinstance(row, dict) or not isinstance(row.get('checked'), dict)
                    or not isinstance(row['checked'].get('witness'), dict)):
                raise ValueError('each hint point must contain a completed rational witness')
            mean = _rational(row['mean'])
            if (not 0 < mean < 1 or tuple(map(Q, row['cell'])) != (mean, mean)
                    or not Q(scope['root'][0]) <= mean <= Q(scope['root'][1])):
                raise ValueError('matching interior singleton cell within the root required')
            result.append(_clean_witness(row['checked']['witness'], (mean, mean)))
    elif record.get('schema') == hill_cover.SCHEMA:
        cells = hill_cover.validate_record(record)
        result = [_clean_witness(row['witness'], cells[path])
                  for path, row in sorted(record['cover']['leaves'].items())]
        if not result:
            raise ValueError('cover hint source must contain accepted interval witnesses')
    else:
        raise ValueError('hill point or hill interval-cover source required')
    return result


def load_hints(paths, scope):
    paths = list(paths)
    if not paths:
        raise ValueError('at least one hint source required')
    rows, sources = [], []
    seen = set()
    for value in paths:
        path = Path(value).resolve()
        if path in seen:
            raise ValueError('distinct hint sources required')
        seen.add(path)
        raw = path.read_bytes()
        record = json.loads(raw)
        hints = proposals(record, scope)
        rows.extend(hints)
        sources.append(dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest(),
                            schema=record['schema'], hint_count=len(hints)))
    return rows, sources


def check_sources(sources):
    for source in sources:
        if hashlib.sha256(Path(source['path']).read_bytes()).hexdigest() != source['sha256']:
            raise ValueError('hint source changed during search')


def wrap(model, scope, paths, target_bits):
    """Return (search proxy, fresh regional cache, source metadata).

    The model must already be freshly authenticated for scope. The root,
    cutoff, update count, and variance geometry are checked again here.
    Final replay must call hill_cover.fresh_model directly, not this factory.
    """
    canonical_scope(scope)
    if (tuple(model.root) != tuple(map(Q, scope['root']))
            or model.threshold != scope['threshold'] or model.q_min != scope['minimum_groups']
            or model.data.get('updates') != scope['updates']
            or model.variance_bins != scope['variance_bins']):
        raise ValueError('fresh model geometry differs from the hint target scope')
    probes, sources = load_hints(paths, scope)
    cache = RegionalCache()
    proxy = PointReuseModel(model, probes, target_bits, updates=scope['updates'],
        outward=lambda cell, witness: cache.bound(model, cell, witness),
        regional_warm=lambda witness: cache.is_warm(model, witness))
    return proxy, cache, sources
