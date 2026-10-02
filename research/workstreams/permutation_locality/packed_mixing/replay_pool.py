"""Optional, uncached process-parallel replay of authenticated dense witnesses.

The caller must first authenticate the exact positive mixture against fresh
outer bounds and validate the complete partition. This helper verifies only
the supplied cells, not those prerequisites or a whole-code certificate.
Every spawned worker constructs one ordinary scalar_cover.Model. No search
wrapper, regional-operator cache, or saved numerical endpoint is used.
"""
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager, nullcontext, redirect_stdout
from fractions import Fraction as Q
import hashlib
import json
from multiprocessing import get_context
import os
from pathlib import Path
import sys


SCHEMA = 'packed-canonical-gl32-dense-replay-pool-1'


def _exact(value):
    if type(value) not in (int, str, Q):
        raise ValueError('exact integer, rational string, or Fraction required')
    try:
        return Q(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('finite exact rational required') from error


def _plain(value):
    """Copy an exact JSON witness, without introducing floating arithmetic."""
    if type(value) is Q:
        return str(value)
    if value is None or type(value) in (str, int, bool):
        return value
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, dict) and all(type(k) is str for k in value):
        return {k: _plain(v) for k, v in value.items()}
    raise ValueError('witnesses must contain exact JSON data, not floats or objects')


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def prepare(mixture_rows, threshold, q_min, precision, tasks):
    """Validate/copy worker inputs and impose deterministic binary-path order."""
    if (type(threshold) is not int or not 0 <= threshold < 1 << 20
            or type(q_min) is not int or not 1 <= q_min <= 2048
            or type(precision) is not int or precision < 128):
        raise ValueError('valid integer cutoff, minimum groups, and precision >=128 required')
    if not isinstance(mixture_rows, (list, tuple)) or not mixture_rows:
        raise ValueError('nonempty exact positive mixture required')
    mixture = []
    for row in mixture_rows:
        if isinstance(row, dict) and set(row) == {'mass', 'activity'}:
            mass, activity = row['mass'], row['activity']
        elif isinstance(row, (list, tuple)) and len(row) == 2:
            mass, activity = row
        else:
            raise ValueError('mixture rows must be (mass, activity) exact pairs')
        mass, activity = _exact(mass), _exact(activity)
        if mass <= 0 or not 0 < activity <= 1:
            raise ValueError('positive masses and activities in (0,1] required')
        mixture.append(dict(mass=str(mass), activity=str(activity)))
    if len({r['activity'] for r in mixture}) != len(mixture):
        raise ValueError('distinct mixture centers required')
    scope = dict(schema=SCHEMA, ensemble='canonical-gl32-width8-shared4-r2',
        K=1 << 20, N=1 << 21, updates=2, block_width=8, threshold=threshold,
        minimum_groups=q_min, maximum_groups=2048, precision=precision,
        base_tilt='3/16', variance_bins=16, regional_count=True, mixture=mixture)
    scope['scope_sha256'] = _digest(scope)
    jobs = _tasks(tasks)
    return scope, jobs


def _tasks(tasks):
    jobs = []
    for task in tasks:
        if not isinstance(task, (list, tuple)) or len(task) != 3:
            raise ValueError('tasks must be (path, cell, witness) triples')
        path, cell, witness = task
        if (type(path) is not str or len(path) > 64 or set(path)-{'0', '1'}
                or not isinstance(cell, (list, tuple)) or len(cell) != 2
                or not isinstance(witness, dict)):
            raise ValueError('binary path, two exact endpoints, and witness required')
        cell = tuple(map(_exact, cell))
        if not 0 <= cell[0] < cell[1] <= 1:
            raise ValueError('nonempty mean cell in [0,1] required')
        witness = _plain(witness)
        parameters = witness.get('parameters')
        if (not isinstance(parameters, list) or len(parameters) != 3
                or _exact(parameters[0]) <= 0 or _exact(parameters[2]) < 0
                or _exact(witness.get('tilt')) <= 0):
            raise ValueError('positive output/activity tilts and nonnegative group dual required')
        _exact(parameters[1])
        jobs.append((path, cell, witness))
    jobs.sort(key=lambda job: job[0])
    if not jobs or len({job[0] for job in jobs}) != len(jobs):
        raise ValueError('nonempty set of unique task paths required')
    for previous, current in zip(jobs, jobs[1:]):
        if current[0].startswith(previous[0]) or previous[1][1] > current[1][0]:
            raise ValueError('task paths and mean cells must not overlap')
    return jobs


def _endpoint(pair):
    if (not isinstance(pair, list) or len(pair) != 2
            or any(type(v) is not int for v in pair) or pair[0] <= 0):
        raise ValueError('positive exact dyadic endpoint required')
    return list(pair)


def _cell_for_path(root, path):
    lo, hi = map(Q, root)
    for bit in path:
        mid = (lo+hi)/2
        if bit == '0':
            hi = mid
        else:
            lo = mid
    return lo, hi


def evaluate(model, scope, task):
    """Evaluate one rational witness through the ordinary model front door."""
    from flint import ctx
    path, cell, witness = task
    if _cell_for_path(model.root, path) != cell:
        raise ArithmeticError('fresh model domain differs from task geometry')
    ctx.prec = scope['precision']
    upper = model.outward(cell, witness)
    if ctx.prec != scope['precision'] or not upper.is_finite() or not upper > 0:
        raise ArithmeticError('finite positive fresh bound at requested precision required')
    return dict(path=path, cell=list(map(str, cell)), precision=scope['precision'],
        threshold=scope['threshold'], minimum_groups=scope['minimum_groups'],
        scope_sha256=scope['scope_sha256'], witness_sha256=_digest(witness),
        output_tilt=str(_exact(witness['parameters'][0])),
        upper=[int(v) for v in upper.upper().man_exp()])


def validate_results(scope, tasks, rows):
    """Check one-to-one result scopes; the caller aggregates the endpoints."""
    jobs, rows = _tasks(tasks), iter(rows)
    checked = []
    for path, cell, witness in jobs:
        row = next(rows, None)
        expected = dict(path=path, cell=list(map(str, cell)), precision=scope['precision'],
            threshold=scope['threshold'], minimum_groups=scope['minimum_groups'],
            scope_sha256=scope['scope_sha256'], witness_sha256=_digest(witness),
            output_tilt=str(_exact(witness['parameters'][0])))
        if (not isinstance(row, dict) or any(row.get(k) != v for k, v in expected.items())
                or any(type(row.get(k)) is not int for k in ('precision', 'threshold', 'minimum_groups'))):
            raise ArithmeticError('missing, duplicate, reordered, or mismatched worker result')
        checked.append(dict(expected, upper=_endpoint(row.get('upper'))))
    if next(rows, None) is not None:
        raise ArithmeticError('unexpected extra worker result')
    return checked


def _build_model(scope):
    from flint import ctx
    parent = Path(__file__).resolve().parents[1]
    for directory in (parent, parent/'gf16_packets'):
        if str(directory) not in sys.path:
            sys.path.insert(0, str(directory))
    import birth_classes
    import shared_mixture
    import scalar_cover as sc
    data = birth_classes.actual(2)
    # actual() may set its own working precision; restore before model setup.
    ctx.prec = scope['precision']
    mixture = [(Q(r['mass']), Q(r['activity'])) for r in scope['mixture']]
    model = sc.Model(shared_mixture.as_components(mixture), data, scope['threshold'],
        scope['minimum_groups'], Q(3, 16), inner=birth_classes,
        variance_shuffle=True, variance_bins=16, regional_count=True)
    if type(model) is not sc.Model or ctx.prec != scope['precision']:
        raise ArithmeticError('fresh unwrapped model at requested precision required')
    return model


def _initialize(scope, logdir):
    global _model, _scope, _log
    _scope = scope
    _log = (open(Path(logdir)/f'replay-{os.getpid()}.log', 'x', buffering=1)
        if logdir is not None else None)
    with redirect_stdout(_log) if _log is not None else nullcontext():
        _model = _build_model(scope)


def _evaluate_worker(task):
    with redirect_stdout(_log) if _log is not None else nullcontext():
        result = evaluate(_model, _scope, task)
        if _log is not None:
            print('PACKED FRESH REPLAY checked', result['path'], flush=True)
        return result


@contextmanager
def _spawn_path():
    previous = sys.path[:]
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        yield
    finally:
        sys.path[:] = previous


def replay_dense(mixture_rows, threshold, q_min, precision, tasks, workers=3, logdir=None):
    """Spawn 1..4 independent uncached workers; return only checked fresh rows.

    The caller must use the normal ``if __name__ == '__main__'`` guard. An
    optional log directory must be new; existing evidence is never replaced.
    This function does not authenticate outer bounds or assert full coverage.
    """
    if type(workers) is not int or not 1 <= workers <= 4:
        raise ValueError('integer worker count in 1..4 required')
    scope, jobs = prepare(mixture_rows, threshold, q_min, precision, tasks)
    if logdir is not None:
        directory = Path(logdir)
        if directory.exists():
            raise ValueError('a new worker-log directory is required')
        directory.mkdir(parents=True)
        logdir = str(directory.resolve())
    with _spawn_path(), ProcessPoolExecutor(max_workers=workers,
            mp_context=get_context('spawn'), initializer=_initialize,
            initargs=(scope, logdir)) as executor:
        return validate_results(scope, jobs, executor.map(_evaluate_worker, jobs, chunksize=1))
