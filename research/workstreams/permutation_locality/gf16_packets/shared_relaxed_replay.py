"""Fresh dense verification in separate arithmetic contexts.

This is a verifier, not a search: every supplied leaf must form part of a
complete partition and is recomputed from its rational witness. The parent
regenerates the shared-support comparison; each worker reconstructs the
actual inner with those authenticated comparison coefficients.
"""
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from fractions import Fraction as Q
from multiprocessing import get_context
import os
from pathlib import Path

from flint import ctx
import shared_relaxed_parallel as search
import shared_relaxed_strategy as proof


def evaluate(model, job, precision, digest):
    path, cell, witness = job
    ctx.prec = precision
    upper = model.outward(cell, witness)
    if ctx.prec != precision or not upper > 0 or not upper.is_finite():
        raise ArithmeticError('finite positive fresh bound at requested precision required')
    return dict(path=path, cell=list(map(str, cell)), precision=precision,
        component_sha256=digest, output_tilt=str(Q(witness['parameters'][0])),
        upper=[int(v) for v in upper.upper().man_exp()])


def initialize(components, record, index, precision, directory):
    global _model, _precision, _digest, _log
    _log = open(Path(directory)/f'replay-{os.getpid()}.log', 'x', buffering=1)
    with redirect_stdout(_log):
        _model = search.worker_model(components, record, index, precision)
    _precision = precision
    _digest = search.component_digest(_model.components)


def evaluate_worker(job):
    with redirect_stdout(_log):
        row = evaluate(_model, job, _precision, _digest)
        print('DENSE FRESH WORKER checked', row['path'], flush=True)
        return row


def replay(model, leaves, mapping, precision):
    """Consume fresh worker evaluations; reject missing or mismatched cells."""
    if type(precision) is not int or precision < 256:
        raise ValueError('integer precision at least 256 required')
    cells = proof.sc.partition(model, leaves, {})
    digest = search.component_digest(model.components)
    jobs = [(path, cells[path], leaves[path]['witness']) for path in leaves]
    rows = iter(mapping(evaluate_worker, jobs))
    checked, total = [], Q(0)
    for path, cell, witness in jobs:
        row = next(rows, None)
        if (not isinstance(row, dict) or row.get('path') != path
                or row.get('cell') != list(map(str, cell))
                or type(row.get('precision')) is not int or row['precision'] != precision
                or row.get('component_sha256') != digest
                or row.get('output_tilt') != str(Q(witness['parameters'][0]))):
            raise ArithmeticError('missing worker result or mismatched proof scope')
        bound = proof.dyadic(row.get('upper'))
        if bound <= 0 or Q(row['output_tilt']) <= 0:
            raise ArithmeticError('positive outward bound and output tilt required')
        total += bound
        checked.append(dict(path=path, upper=row['upper'], output_tilt=row['output_tilt']))
        print('SHARED DENSE PARALLEL REPLAY', len(checked), 'of', len(jobs),
              'path', path, flush=True)
    if next(rows, None) is not None:
        raise ArithmeticError('unexpected extra worker result')
    return dict(complete=True, leaves=len(checked), upper=proof.endpoint(total),
        checked=checked, threshold=model.threshold, ensemble=proof.ENSEMBLE)


def replay_dense(record, precision=384, index=0, workers=4, logdir=None):
    proof.require_uncached_verification()
    if type(workers) is not int or not 1 <= workers <= 8 or logdir is None:
        raise ValueError('one through eight workers and a new log directory required')
    row = proof.validate_record(record, index, complete=True)
    directory = Path(logdir)
    if directory.exists():
        raise ValueError('use a new worker-log directory; prior evidence is preserved')
    model = proof.build_model(record, precision, index)
    directory.mkdir(parents=True)
    with search.spawn_path(), ProcessPoolExecutor(max_workers=workers,
            mp_context=get_context('spawn'), initializer=initialize,
            initargs=(model.components, record, index, precision, directory)) as executor:
        return replay(model, row['cover']['leaves'], executor.map, precision)
