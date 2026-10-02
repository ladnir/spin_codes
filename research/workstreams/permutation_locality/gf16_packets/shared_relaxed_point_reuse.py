"""Reuse point-search duals on finite cells, with no bound transfer.

This is an isolated search helper. Rescaling a variance partition only
chooses new subintervals: all rational duals remain proposals, and outward()
must rebuild their bounds on the new cell. No search engine imports it.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from time import perf_counter

from flint import arb, ctx
import shared_relaxed_strategy as proof
import variance_partition as variance


def rebase_witness(source_cell, target_cell, witness):
    """Rescale a complete variance partition; preserve its rational duals.

    Unlike a child restriction, the target need not lie in the source.
    The transformed witness carries no numerical upper bound. The local
    operator and each mean/variance/MGF dual must be freshly evaluated.
    """
    source_cell, target_cell = tuple(map(Q, source_cell)), tuple(map(Q, target_cell))
    if len(source_cell) != 2 or len(target_cell) != 2 or not isinstance(witness, dict):
        raise ValueError('two mean cells and a rational witness required')
    source_max = variance.maximum_variance(source_cell)
    target_max = variance.maximum_variance(target_cell)
    result = copy.deepcopy(witness)
    if 'variance_partition' not in witness:
        if 'regional_count_parts' in witness:
            raise ValueError('regional witnesses require a complete variance partition')
        return result
    parts = variance.validate(source_cell, witness['variance_partition'])
    regional = witness.get('regional_count_parts')
    if regional is not None:
        if not isinstance(regional, list) or len(regional) != len(parts):
            raise ValueError('one regional family per variance interval required')
        for saved, (interval, dual) in zip(regional, parts):
            if (not isinstance(saved, dict)
                    or tuple(map(Q, saved.get('interval', ()))) != interval
                    or tuple(map(Q, saved.get('dual', ()))) != dual
                    or not isinstance(saved.get('mgf_witnesses'), list)
                    or len(saved['mgf_witnesses']) > 128
                    or any(not isinstance(v, dict) for v in saved['mgf_witnesses'])):
                raise ValueError('regional MGF families must match the source partition')
    scale = target_max/source_max
    rebased, families = [], []
    for i, (interval, dual) in enumerate(parts):
        interval = [str(v*scale) for v in interval]
        dual = list(map(str, dual))
        rebased.append(dict(interval=interval, dual=dual))
        if regional is not None:
            saved = copy.deepcopy(regional[i])
            saved.update(interval=interval, dual=dual)
            families.append(saved)
    variance.validate(target_cell, rebased)
    result['variance_partition'] = rebased
    if regional is not None:
        result['regional_count_parts'] = families
    return result


def nearest_witness(probes, cell):
    """Choose by mean distance only; saved scores/endpoints are not read."""
    cell = tuple(map(Q, cell))
    variance.maximum_variance(cell)
    if not isinstance(probes, list) or not probes:
        raise ValueError('nonempty list of point proposals required')
    candidates = []
    for index, point in enumerate(probes):
        if not isinstance(point, dict) or not isinstance(point.get('witness'), dict):
            raise ValueError('each point requires its rational witness')
        mean = Q(point['mean'])
        variance.maximum_variance((mean, mean))
        candidates.append((abs(mean-sum(cell)/2), index, mean, point['witness']))
    _, index, mean, witness = min(candidates, key=lambda row: row[:2])
    return index, mean, rebase_witness((mean, mean), cell, witness)


def point_descendant(cell, point, depth):
    """Return one dyadic descendant containing a selected interior point."""
    cell = tuple(map(Q, cell)); point = Q(point)
    variance.maximum_variance(cell)
    if (type(depth) is not int or not 0 <= depth <= 16
            or not cell[0] <= point <= cell[1]):
        raise ValueError('contained point and integer depth in 0..16 required')
    suffix = ''
    for _ in range(depth):
        mid = sum(cell)/2
        if point < mid:
            suffix += '0'; cell = (cell[0], mid)
        else:
            suffix += '1'; cell = (mid, cell[1])
    return suffix, cell


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('points', type=Path)
    parser.add_argument('cover', type=Path)
    parser.add_argument('--paths', nargs='+', required=True)
    parser.add_argument('--descend', type=int, default=0,
        help='Check only a descendant containing the nearest point; never modify the full cover')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--target-bits', type=int, default=40)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if (args.output.exists() or args.output.resolve() in (args.points.resolve(), args.cover.resolve())
            or not 256 <= args.precision <= 1024 or args.target_bits < 20
            or not 0 <= args.descend <= 16):
        parser.error('new output, precision 256..1024, and target margin >=20 required')
    point_raw, cover_raw = args.points.read_bytes(), args.cover.read_bytes()
    points, cover = json.loads(point_raw), json.loads(cover_raw)
    point_row, cover_row = proof.validate_record(points), proof.validate_record(cover)
    for key in ('updates', 'K', 'N', 'minimum_groups', 'base_tilt', 'variance_bins',
                'mixture', 'cap_sha256'):
        if points[key] != cover[key]:
            parser.error('point and cover records require the identical checked comparison')
    if (point_row['threshold'] != cover_row['threshold']
            or tuple(map(Q, point_row['root'])) != tuple(map(Q, cover_row['root']))):
        parser.error('point and cover scopes differ')
    model = proof.build_model(points, args.precision)
    cells = proof.sc.partition(model, cover_row['cover']['leaves'], cover_row['cover']['unresolved'])
    saved = dict(cover_row['cover']['leaves'], **cover_row['cover']['unresolved'])
    if len(set(args.paths)) != len(args.paths) or any(path not in cells for path in args.paths):
        parser.error('distinct existing partition paths required')
    result = dict(schema='shared-gf16-point-reuse-probe-1', precision=args.precision,
        target_bits=args.target_bits, threshold=model.threshold,
        point_sha256=hashlib.sha256(point_raw).hexdigest(),
        cover_sha256=hashlib.sha256(cover_raw).hexdigest(), rows=[],
        note='Fresh checks of selected cells only; no complete certificate or search-state mutation.')
    for path in args.paths:
        cell = cells[path]
        index, mean, candidate = nearest_witness(point_row['probes'], cell)
        original_cell = cell
        suffix = ''
        if args.descend:
            suffix, cell = point_descendant(cell, mean, args.descend)
            candidate = rebase_witness((mean, mean), cell, point_row['probes'][index]['witness'])
        choices = [('nearest-point', candidate)]
        if isinstance(saved[path].get('witness'), dict):
            retained = saved[path]['witness']
            if args.descend:
                from shared_relaxed_parallel import restrict_witness
                retained = restrict_witness(original_cell, cell, retained)
            choices.append(('retained-cell', retained))
        for name, witness in choices:
            ctx.prec = args.precision
            start = perf_counter()
            upper = model.outward(cell, witness)
            elapsed = perf_counter()-start
            if ctx.prec != args.precision or not upper > 0:
                raise ArithmeticError('positive fresh bound at requested precision required')
            row = dict(path=path+suffix, original_path=path, cell=list(map(str, cell)), kind=name,
                point_index=index, point_mean=str(mean), witness=witness,
                upper=proof.endpoint(proof.dyadic([int(v) for v in upper.upper().man_exp()])),
                margin_passed=bool(upper < arb(2)**-args.target_bits), seconds=elapsed)
            result['rows'].append(row)
            print('POINT REUSE', path+suffix, name, 'mean', mean, 'log2 upper', upper.log()/arb(2).log(),
                  'seconds', elapsed, flush=True)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
