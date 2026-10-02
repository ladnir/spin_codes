"""Optional fresh-process dense replay for the selected physical-t64 scope.

Every spawned worker freshly authenticates outer counts, the mixture, and
the selected physical maps. Its complete scope must match the search input.
Workers evaluate only rational cell witnesses; saved bounds are never sent.
The parent validates one result per cell and emits the normal dense replay
schema. This file does not claim a whole-code or implementation certificate.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from fractions import Fraction as Q
from multiprocessing import get_context
import os
from pathlib import Path
import sys
from time import monotonic

from flint import ctx
import dense_t64
import whole_t64 as whole


def prepare(source, proposal, precision, workers):
    whole.require(type(precision) is int and precision >= 128
        and type(workers) is int and 1 <= workers <= 4, 'precision>=128 and workers1..4 required')
    record, source_meta = whole.read_source(source)
    cells = dense_t64.validate_cover_record(record, complete=True)
    scope = record['scope']
    whole.require(all(scope.get(key) == value for key, value in whole.DENSE_SCOPE.items()),
        'explicit selected physical64/macro128 dense scope required')
    _, proposal_meta = whole.read_source(proposal)
    whole.require(record['source'] == proposal_meta, 'immutable mixture proposal source differs')
    whole.claim(scope)
    whole.require(type(scope.get('variance_bins')) is int and 1 <= scope['variance_bins'] <= 64
        and Q(scope['base_tilt']) == Q(3, 16), 'supported dense proof parameters required')
    jobs = [(path, tuple(cells[path]), deepcopy(record['cover']['leaves'][path]['witness']))
            for path in sorted(cells)]
    return record, source_meta, proposal_meta, cells, jobs


def _initialize(proposal, expected_scope, precision):
    global _model, _scope, _precision, _worker
    whole.unchanged(proposal)
    model, scope, source = dense_t64.fresh_model(proposal['path'], precision=precision,
        variance_bins=expected_scope['variance_bins'], distance=Q(expected_scope['distance']),
        birth_density=expected_scope['birth_density'])
    whole.require(scope == expected_scope and source == proposal,
        'spawned worker freshly authenticated a different scope or proposal source')
    whole.require(tuple(model.root) == tuple(map(Q, expected_scope['root']))
        and ctx.prec == precision, 'fresh worker root or precision differs')
    whole.unchanged(proposal)
    _model, _scope, _precision = model, scope, precision
    _worker = dict(pid=os.getpid(), precision=precision,
        scope_sha256=dense_t64.prior.fingerprint(scope), source=source,
        initialization='fresh dense_t64.fresh_model in spawned process')
    print('T64 DENSE worker ready', os.getpid(), flush=True)


def evaluate(model, scope, precision, job):
    """One fresh outward evaluation; only its exact proof row is returned."""
    path, cell, witness = job
    expected = dense_t64.geometry.path_cell(scope['root'], path)
    whole.require(tuple(cell) == expected and tuple(model.root) == tuple(map(Q, scope['root']))
        and isinstance(witness, dict), 'job geometry differs from the fresh model')
    ctx.prec = precision
    upper = model.outward(expected, witness)
    whole.require(ctx.prec == precision, 'worker changed requested Arb precision')
    return dict(path=path, cell=list(map(str, expected)), upper=dense_t64.endpoint(upper))


def _evaluate_worker(job):
    whole.unchanged(_worker['source'])
    row = evaluate(_model, _scope, _precision, job)
    whole.unchanged(_worker['source'])
    return dict(row=row, worker=dict(_worker))


def check_result(message, cells, seen, scope, proposal, precision):
    """Validate result/provenance independently of executor ordering."""
    whole.require(isinstance(message, dict) and set(message) == {'row', 'worker'}, 'worker result envelope required')
    row, worker = message['row'], message['worker']
    whole.require(isinstance(row, dict) and set(row) == {'path', 'cell', 'upper'}
        and isinstance(row.get('path'), str) and row['path'] in cells and row['path'] not in seen,
        'worker result omits, repeats, or adds a cell')
    whole.require(row['cell'] == list(map(str, cells[row['path']])), 'worker cell differs from exact partition')
    whole.dyadic(row['upper'])
    whole.require(isinstance(worker, dict) and type(worker.get('pid')) is int and worker['pid'] > 0
        and type(worker.get('precision')) is int and worker['precision'] == precision
        and worker.get('scope_sha256') == dense_t64.prior.fingerprint(scope)
        and worker.get('source') == proposal
        and worker.get('initialization') == 'fresh dense_t64.fresh_model in spawned process',
        'worker initialization, scope, precision, or source metadata differ')
    return row, worker


def run(source, proposal, output, precision=384, workers=2):
    output = Path(output).resolve()
    whole.require(not output.exists(), 'fresh dense replay output required')
    old, source_meta, proposal_meta, cells, jobs = prepare(source, proposal, precision, workers)
    scope = old['scope']
    inputs = [source_meta, proposal_meta, scope['maps']['source']]
    for item in inputs:
        whole.unchanged(item)
    manifest = whole.proof_source_manifest(inputs)
    start = monotonic()
    result = dict(schema='packed-gl32-t64-s16-dense-replay-1', scope=scope,
        source=proposal_meta, input_source=source_meta, precision=precision,
        workers=workers, worker_start_method='spawn', worker_metadata=[], cells=[],
        complete_dense=False, whole_code_certificate=False,
        proof_status='Fresh dense replay; sparse proof remains separate.')
    output.parent.mkdir(parents=True, exist_ok=True)
    import json
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    def save():
        result['elapsed_seconds'] = monotonic()-start
        output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    total, seen, participants = Q(0), set(), {}
    # Explicit spawn avoids inheriting prepared models or floating state.
    # A temporary import path also supports calls from a test/orchestrator.
    old_path, old_environment = sys.path[:], {key: os.environ.get(key) for key in
        ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'PYTHONDONTWRITEBYTECODE')}
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    try:
        with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn'),
                initializer=_initialize, initargs=(proposal_meta, scope, precision)) as executor:
            for message in executor.map(_evaluate_worker, jobs, chunksize=1):
                row, worker = check_result(message, cells, seen, scope, proposal_meta, precision)
                if worker['pid'] in participants:
                    whole.require(participants[worker['pid']] == worker, 'worker metadata changed between cells')
                participants[worker['pid']] = worker
                whole.require(len(participants) <= workers, 'more worker processes than requested')
                seen.add(row['path'])
                total += whole.dyadic(row['upper'])
                result['cells'].append(row)
                result['worker_metadata'] = [participants[pid] for pid in sorted(participants)]
                result['aggregate_upper'] = whole.compact_endpoint(total, precision+32)
                for item in inputs:
                    whole.unchanged(item)
                save()
                print('T64 PARALLEL REPLAY', len(seen), '/', len(cells), 'path', row['path'], flush=True)
    finally:
        sys.path[:] = old_path
        for key, value in old_environment.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    whole.require(seen == set(cells), 'fresh worker results do not cover the complete dense partition')
    after = whole.proof_source_manifest(inputs)
    whole.require(after == manifest, 'proof sources or inputs changed during dense replay')
    for item in inputs:
        whole.unchanged(item)
    result['cells'].sort(key=lambda row: row['path'])
    result.update(complete_dense=True, aggregate_below_2_minus_40=bool(total < Q(2)**-40),
        aggregation='Exact sum of fresh cell endpoints, rounded upward once to precision+32 mantissa bits.',
        source_manifest_before=manifest, source_manifest_after=after)
    save()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=384)
    parser.add_argument('--workers', type=int, choices=range(1, 5), default=2)
    args = parser.parse_args()
    run(args.input, args.source, args.output, args.precision, args.workers)


if __name__ == '__main__':
    main()
