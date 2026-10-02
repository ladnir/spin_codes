"""Parallel continuation of a shared GF16 R2 dense partition.

The parent authenticates one exact comparison. Separate workers reconstruct
the actual inner and recheck retained witnesses with those same rational
components. Checkpoints contain a complete partition, including unfinished
cells. A full claim still requires the independent whole-code assembler.
"""
import argparse
import atexit
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager, redirect_stdout
import copy
from fractions import Fraction as Q
import hashlib
import heapq
import json
from math import isfinite
from multiprocessing import get_context
import os
from pathlib import Path
import sys
from types import SimpleNamespace

from flint import arb, ctx
import shared_relaxed_strategy as proof


def component_digest(components):
    rows = [[str(Q(c)), str(Q(p)), a] for c, p, a in components]
    return hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()


def retarget_record(record, distance, index, source_sha256):
    """Change only the selected cutoff claim; retain witnesses as proposals.

    The returned record is a copy. Search must replay every cell at its new
    cutoff, and the whole-code assembler must still regenerate both ranges.
    Old point probes are archived under the provenance, not relabeled.
    """
    old = proof.validate_record(record, index)
    if (not isinstance(source_sha256, str) or len(source_sha256) != 64
            or any(c not in '0123456789abcdef' for c in source_sha256)
            or type(distance) not in (str, int, Q)):
        raise ValueError('exact distance and original source SHA256 required')
    try:
        delta = Q(distance)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('valid exact relative distance required') from error
    if not 0 < delta < Q(1,2):
        raise ValueError('retargeted relative distance must lie in (0,1/2)')
    result = copy.deepcopy(record)
    row = result['results'][index]
    row['distance'] = str(delta)
    row['threshold'] = int(delta*proof.sc.N)
    row['retargeted_from'] = dict(source_sha256=source_sha256, result_index=index,
        distance=old['distance'], threshold=old['threshold'],
        prior_probes=copy.deepcopy(old.get('probes', [])))
    row['probes'] = []
    row.pop('rechecked_precision', None)
    result['proof_status'] = ('Retargeted claim only. Saved cell witnesses are proposals; '
        'every cell must be freshly replayed at the new cutoff before acceptance. '
        'Full sparse/dense assembly is still required.')
    proof.validate_record(result, index)
    return result


def worker_model(components, record, index, precision):
    row = proof.validate_record(record, index)
    if type(precision) is not int or precision < 256:
        raise ValueError('worker precision at least 256 required')
    data = proof.birth_classes.actual(2)
    ctx.prec = precision
    rows = [(str(i), c, tuple(proof.sc.probabilities(p)), active)
            for i, (c, p, active) in enumerate(components)]
    model = proof.sc.Model(rows, data, row['threshold'], record['minimum_groups'],
        Q(record['base_tilt']), inner=proof.birth_classes, variance_shuffle=True,
        variance_bins=record['variance_bins'], regional_count=True)
    if (model.components != components or model.root != tuple(map(Q, row['root']))
            or any(data.get(k) != v for k, v in dict(bits=19, windows=32, updates=2).items())):
        raise ArithmeticError('worker geometry or components differ from checked parent model')
    model.point_probes = row.get('probes', []) if record.get('reuse_point_witnesses', False) else []
    model.split_on_midpoint = record.get('split_on_midpoint', False)
    return model


def initialize(components, record, index, precision, log_directory):
    global _model, _precision, _digest, _log, _region_cache
    _log = open(Path(log_directory)/f'worker-{os.getpid()}.log', 'x')
    with redirect_stdout(_log):
        _model = worker_model(components, record, index, precision)
    _precision = precision
    _digest = component_digest(_model.components)
    _region_cache = None
    entries = record.get('region_cache_entries', 0)
    if entries:
        from shared_relaxed_region_cache import install
        context = install(entries)
        _region_cache = context.__enter__()
        atexit.register(context.__exit__, None, None, None)


def evaluate(model, job, precision, digest):
    path, cell, target_bits, old = job
    upper = None
    witness = old
    retained = None
    split = False
    ctx.prec = precision
    if old is not None:
        value = model.outward(cell, old)
        if ctx.prec != precision or not value > 0 or not value.is_finite():
            raise ArithmeticError('positive fresh retained bound at worker precision required')
        retained = value
        if value < arb(2)**-target_bits:
            upper = value
            score = float(value.log()/arb(2).log())
    if upper is None and getattr(model, 'point_probes', []) and 0 < cell[0] <= cell[1] < 1:
        from shared_relaxed_point_reuse import nearest_witness
        _, _, candidate = nearest_witness(model.point_probes, cell)
        value = model.outward(cell, candidate)
        if ctx.prec != precision or not value > 0 or not value.is_finite():
            raise ArithmeticError('finite positive fresh point-reuse bound required')
        if retained is None or value.upper() < retained.upper():
            retained, witness = value, candidate
        if retained < arb(2)**-target_bits:
            upper = retained
            score = float(retained.log()/arb(2).log())
    if (upper is None and retained is not None and getattr(model, 'split_on_midpoint', False)
            and cell[0] < cell[1]):
        mid = sum(cell)/2
        candidate = restrict_witness(cell, (mid, mid), witness)
        value = model.outward((mid, mid), candidate)
        if ctx.prec != precision or not value > 0 or not value.is_finite():
            raise ArithmeticError('finite positive fresh midpoint bound required')
        if value < arb(2)**-(target_bits+2):
            # A passing midpoint never accepts the interval. Return a failed
            # cell for bisection, retaining only its rational duals.
            split = True
            score = float(retained.log()/arb(2).log())
    if upper is None and not split:
        model.proposal_stop_bits = target_bits+2
        score, witness = model.proposal(cell)
        # Proposals can initialize caches at another precision. The final
        # inequality is always recomputed at the requested worker precision.
        ctx.prec = precision
        if score < -(target_bits+2):
            value = model.outward(cell, witness)
            if ctx.prec != precision or not value > 0:
                raise ArithmeticError('positive fresh proposal bound at worker precision required')
            score = float(value.log()/arb(2).log())
            if value < arb(2)**-target_bits:
                upper = value
    if ctx.prec != precision or not isfinite(score):
        raise ArithmeticError('invalid worker precision or score')
    return dict(path=path, cell=list(map(str, cell)), proposal=float(score), witness=witness,
        upper=None if upper is None else [int(v) for v in upper.upper().man_exp()],
        component_sha256=digest)


def evaluate_worker(job):
    with redirect_stdout(_log):
        result = evaluate(_model, job, _precision, _digest)
        print('SHARED WORKER finished', job[0], result['proposal'], flush=True)
        if _region_cache is not None:
            print('SEARCH REGION CACHE hits', _region_cache.hits, 'misses', _region_cache.misses, flush=True)
        return result


def restrict_witness(parent, child, witness):
    """Restrict saved rational duals; do not transfer a numerical bound.

    The child evaluation rebuilds every bound. Variance intervals cover the
    child's full range, and any regional MGF families keep matching intervals.
    """
    if (not isinstance(witness, dict) or len(parent) != 2 or len(child) != 2):
        raise ValueError('parent, child, and rational witness required')
    parent, child = tuple(map(Q, parent)), tuple(map(Q, child))
    if not 0 <= parent[0] <= child[0] <= child[1] <= parent[1] <= 1:
        raise ValueError('child must lie within the parent mean interval')
    result = copy.deepcopy(witness)
    if 'variance_partition' not in witness:
        if 'regional_count_parts' in witness:
            raise ValueError('regional witnesses require a variance partition')
        return result
    import variance_partition as variance
    parts = variance.validate(parent, witness['variance_partition'])
    upper = variance.maximum_variance(child)
    saved = witness.get('regional_count_parts')
    if saved is not None:
        if not isinstance(saved, list) or len(saved) != len(parts):
            raise ValueError('one regional family per parent variance part required')
        for part, (interval, dual) in zip(saved, parts):
            if (not isinstance(part, dict) or tuple(map(Q, part.get('interval', ()))) != interval
                    or tuple(map(Q, part.get('dual', ()))) != dual
                    or not isinstance(part.get('mgf_witnesses'), list)
                    or len(part['mgf_witnesses']) > 128
                    or any(not isinstance(v, dict) for v in part['mgf_witnesses'])):
                raise ValueError('regional families must match the parent variance partition')
    clipped, regional = [], []
    for index, ((lo, hi), dual) in enumerate(parts):
        if lo >= upper:
            break
        interval = list(map(str, (lo, min(hi, upper))))
        clipped.append(dict(interval=interval, dual=list(map(str, dual))))
        if saved is not None:
            part = copy.deepcopy(saved[index])
            part.update(interval=interval, dual=list(map(str, dual)))
            regional.append(part)
    variance.validate(child, clipped)
    result['variance_partition'] = clipped
    if saved is not None:
        result['regional_count_parts'] = regional
    return result


def search(model, cover, mapping, workers, max_cells, max_depth, target_bits, checkpoint,
           reuse_parent_witness=False):
    if (type(workers) is not int or not 1 <= workers <= 4
            or type(max_cells) is not int or max_cells < 1
            or type(max_depth) is not int or not 1 <= max_depth <= 64
            or type(target_bits) is not int or target_bits < 20
            or type(reuse_parent_witness) is not bool
            or not isinstance(cover, dict) or not isinstance(cover.get('leaves'), dict)
            or not isinstance(cover.get('unresolved'), dict)):
        raise ValueError('valid search limits, margin >=20, and partition dictionaries required')
    cells = proof.sc.partition(model, cover['leaves'], cover['unresolved'])
    previous = cover.get('visited', 0)
    if type(previous) is not int or previous < 0:
        raise ValueError('nonnegative prior work count required')
    saved = dict(cover['leaves'])
    if reuse_parent_witness:
        saved.update(cover['unresolved'])
    pending = [(len(path), path, cell, saved.get(path, {}).get('witness'))
               for path, cell in cells.items()]
    if not pending:
        raise ValueError('nonempty full partition required')
    heapq.heapify(pending)
    leaves, unresolved, visited = {}, {}, 0
    digest = component_digest(model.components)
    while pending and visited < max_cells:
        batch = [heapq.heappop(pending)
                 for _ in range(min(workers, len(pending), max_cells-visited))]
        jobs = [(path, cell, target_bits, old) for _, path, cell, old in batch]
        results = list(mapping(evaluate_worker, jobs))
        if len(results) != len(jobs):
            raise ArithmeticError('missing worker result')
        for (depth, path, cell, _), result in zip(batch, results):
            if (result.get('path') != path or result.get('cell') != list(map(str, cell))
                    or result.get('component_sha256') != digest
                    or not isinstance(result.get('witness'), dict)
                    or not isfinite(result.get('proposal', float('nan')))):
                raise ArithmeticError('worker scope or comparison mismatch')
            visited += 1
            if result['upper'] is not None:
                if not 0 < proof.dyadic(result['upper']) < Q(2)**-target_bits:
                    raise ArithmeticError('invalid worker outward endpoint')
                leaves[path] = dict(witness=result['witness'], proposal=result['proposal'])
            elif depth >= max_depth:
                unresolved[path] = dict(cell=list(map(str, cell)),
                    proposal=result['proposal'], witness=result['witness'])
            else:
                mid = sum(cell)/2
                for suffix, child in (('0', (cell[0], mid)), ('1', (mid, cell[1]))):
                    old = (restrict_witness(cell, child, result['witness'])
                           if reuse_parent_witness else None)
                    heapq.heappush(pending, (depth+1, path+suffix, child, old))
        result = dict(leaves=dict(leaves), unresolved={**unresolved,
            **{p:dict(cell=list(map(str,c)), **({'witness':w} if w is not None else {}))
               for _,p,c,w in pending}},
            visited=previous+visited, component_sha256=digest)
        proof.sc.partition(model, result['leaves'], result['unresolved'])
        checkpoint(result)
        print('SHARED PARALLEL visited', previous+visited, 'accepted', len(leaves),
              'unresolved', len(result['unresolved']), flush=True)
    return result


@contextmanager
def spawn_path():
    original = list(sys.path)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        yield
    finally:
        sys.path[:] = original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('resume', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--index', type=int, default=0)
    parser.add_argument('--workers', type=int, choices=range(1,5), default=4)
    parser.add_argument('--max-cells', type=int, default=500)
    parser.add_argument('--max-depth', type=int, default=24)
    parser.add_argument('--target-bits', type=int, default=40)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--retarget-distance',
        help='Change only the selected distance claim; freshly replay every cell at the new cutoff')
    parser.add_argument('--reuse-parent-witness', action='store_true',
        help='Try restricted parent duals on each child before a new proposal; every bound is rebuilt')
    parser.add_argument('--reuse-point-witnesses', action='store_true',
        help='Try a rebased point witness from this same record; rebuild its finite-cell bound')
    parser.add_argument('--split-on-midpoint', action='store_true',
        help='Bisect instead of reoptimizing when an old witness passes only at the midpoint; never accept from a point')
    parser.add_argument('--region-cache-entries', type=int, choices=range(9), default=0,
        help='Optional bounded worker cache of fresh local polynomials; full assembly remains uncached')
    parser.add_argument('--completion-driven', action='store_true',
        help='Refill idle workers as individual interval jobs finish; preserve in-flight cells in checkpoints')
    parser.add_argument('--retain-search-leaves', action='store_true',
        help='With completion-driven scheduling, retain old leaf labels as unverified search hints; final replay checks all leaves')
    args = parser.parse_args()
    logs = args.output.with_suffix('.workers')
    staging = args.output.with_suffix('.writing.json')
    if (args.output.exists() or logs.exists() or staging.exists()
            or args.resume.resolve() == args.output.resolve()):
        parser.error('use new output and worker-log paths; previous proofs are preserved')
    if (args.max_cells < 1 or not 1 <= args.max_depth <= 64
            or args.target_bits < 20 or args.precision < 256):
        parser.error('valid budget, depth, margin >=20, and precision >=256 required')
    if args.retain_search_leaves and (not args.completion_driven or args.retarget_distance is not None):
        parser.error('retained search hints require completion-driven scheduling without retargeting')
    raw = args.resume.read_bytes()
    record = json.loads(raw)
    if args.retarget_distance is not None:
        try:
            record = retarget_record(record, args.retarget_distance, args.index,
                                     hashlib.sha256(raw).hexdigest())
        except ValueError as error:
            parser.error(str(error))
    row = proof.validate_record(record, args.index)
    record['reuse_point_witnesses'] = args.reuse_point_witnesses
    record['split_on_midpoint'] = args.split_on_midpoint
    record['region_cache_entries'] = args.region_cache_entries
    if args.reuse_point_witnesses and not row.get('probes'):
        parser.error('point reuse requires point proposals in the selected record')
    if not isinstance(row.get('cover'),dict):
        parser.error('a previous full partition, including unresolved cells, is required')
    # Reject malformed paths before rebuilding the expensive actual maps.
    proof.sc.partition(SimpleNamespace(root=tuple(map(Q,row['root']))),
                       row['cover']['leaves'],row['cover']['unresolved'])
    model = proof.build_model(record, args.precision, args.index)
    output_record = copy.deepcopy(record)
    logs.mkdir(parents=True)
    def save(state):
        output_record['results'][args.index]['cover'] = state
        if args.retain_search_leaves:
            output_record['results'][args.index].pop('rechecked_precision', None)
        else:
            output_record['results'][args.index]['rechecked_precision'] = args.precision
        output_record.update(precision=args.precision, target_bits=args.target_bits,
            search_workers=args.workers, parallel_selected_index=args.index,
            search_scheduler='completion' if args.completion_driven else 'batch',
            retain_search_leaves=args.retain_search_leaves,
            reuse_parent_witness=args.reuse_parent_witness,
            resume_sha256=hashlib.sha256(raw).hexdigest(),
            proof_status=('Saved leaf labels are unverified search hints, not fresh bounds. '
                          'Complete fresh dense/sparse assembly must replay EVERY leaf.'
                          if args.retain_search_leaves else
                          'Retained leaves in the selected result were freshly rechecked against one common exact comparison. '
                          'Complete fresh dense/sparse assembly still required.'))
        staging.write_text(json.dumps(output_record,indent=2)+'\n')
        staging.replace(args.output)
        if args.completion_driven:
            print('SHARED COMPLETION visited',state['visited'],'accepted',len(state['leaves']),
                  'unresolved',len(state['unresolved']),flush=True)
    worker_initializer = initialize
    if args.completion_driven:
        import shared_relaxed_completion as completion
        worker_initializer = completion.initialize
    with spawn_path(), ProcessPoolExecutor(max_workers=args.workers, mp_context=get_context('spawn'),
            initializer=worker_initializer, initargs=(model.components,record,args.index,args.precision,logs)) as executor:
        if args.completion_driven:
            state = completion.search(model,row['cover'],executor,args.workers,args.max_cells,
                args.max_depth,args.target_bits,args.precision,save,args.reuse_parent_witness,
                retain_search_leaves=args.retain_search_leaves)
        else:
            state = search(model,row['cover'],executor.map,args.workers,args.max_cells,
                           args.max_depth,args.target_bits,save,args.reuse_parent_witness)
    print('SHARED PARALLEL finished',len(state['leaves']),'accepted',
          len(state['unresolved']),'unresolved',flush=True)


if __name__ == '__main__':
    main()
