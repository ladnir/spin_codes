"""Dense coverage and fresh replay for separately scoped GL32 hill contexts.

R2, R3, and R4 are different constructions. The source fixes that choice,
the distance, and the authenticated outer refinement. Point bounds never
become accepted intervals. Search receipts are proposals for a later fresh
replay; even a successful dense replay excludes the sparse prefix.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext, redirect_stdout
import copy
from fractions import Fraction as Q
import hashlib
import json
from multiprocessing import get_context
import os
from pathlib import Path
import sys

from flint import arb, ctx
import cell_search
import search_strategy as geometry

SCHEMA = 'packed-canonical-gl32-hill-dense-cover-1'
CONTEXT_SCHEMA = 'packed-canonical-gl32-hill-dense-context-1'
POINT_SCHEMA = 'packed-canonical-gl32-dense-hill-points-1'
INCIDENCE_SCHEMA = 'packed-canonical-full32-h5-incidence-expected-cdf-1'
EKR_SCHEMA = 'packed-canonical-full32-h5-incidence-ekr-expected-cdf-1'


def validate_scope(scope):
    expected = dict(schema=CONTEXT_SCHEMA, K=1 << 20, N=1 << 21,
        block_width=8, maximum_groups=2048, comparison='direct-expected-shell-majorant',
        last_lp=104, refined=True, regional_count=True)
    if (not isinstance(scope, dict) or any(scope.get(k) != v for k, v in expected.items())
            or type(scope.get('updates')) is not int or scope['updates'] not in (2, 3, 4)
            or scope.get('ensemble') != f"canonical-gl32-width8-shared4-r{scope['updates']}"
            or type(scope.get('minimum_groups')) is not int
            or not 1 <= scope['minimum_groups'] <= 2048
            or type(scope.get('variance_bins')) is not int or not 1 <= scope['variance_bins'] <= 64
            or scope.get('refined') is not True or scope.get('regional_count') is not True
            or Q(scope.get('base_tilt', 0)) != Q(3, 16)
            or type(scope.get('threshold')) is not int
            or scope['threshold'] != cell_search.dense.cutoff(scope.get('distance', 0))
            or not isinstance(scope.get('outer_premises'), dict)
            or scope['outer_premises'].get('schema') not in (INCIDENCE_SCHEMA, EKR_SCHEMA)
            or not isinstance(scope.get('mixture'), list) or not scope['mixture']):
        raise ValueError('exact canonical GL32 R2/R3/R4 refined-count scope required')
    geometry.path_cell(scope['root'], '')
    for key in ('expected_cdf_sha256', 'comparison_caps_sha256'):
        digest = scope.get(key)
        if not isinstance(digest, str) or len(digest) != 64 or set(digest.lower())-set('0123456789abcdef'):
            raise ValueError('exact comparison fingerprints required')
    return scope['outer_premises']['schema'] == EKR_SCHEMA


def validate_cover(scope, cover, *, complete=False):
    validate_scope(scope)
    if (not isinstance(cover, dict) or not isinstance(cover.get('leaves'), dict)
            or not isinstance(cover.get('unresolved'), dict)
            or type(cover.get('visited', 0)) is not int or cover.get('visited', 0) < 0):
        raise ValueError('dense partition and nonnegative work count required')
    cells = geometry.partition(scope['root'], cover['leaves'], cover['unresolved'])
    for path, row in dict(cover['leaves'], **cover['unresolved']).items():
        if 'cell' not in row or tuple(map(Q, row['cell'])) != cells[path]:
            raise ValueError('every interval must retain its exact coordinates')
        if path in cover['leaves'] and not isinstance(row.get('witness'), dict):
            raise ValueError('accepted search intervals need rational witnesses')
    if complete and cover['unresolved']:
        raise ValueError('fresh dense replay requires complete interval coverage')
    return cells


def validate_record(record, *, complete=False):
    if (not isinstance(record, dict) or record.get('schema') != SCHEMA
            or record.get('whole_code_certificate') is not False
            or type(record.get('precision')) is not int or record['precision'] < 128
            or type(record.get('cell_target_bits')) is not int or record['cell_target_bits'] < 20):
        raise ValueError('scoped dense-only hill receipt required')
    return validate_cover(record['scope'], record['cover'], complete=complete)


def from_source(source):
    """Start from a hill context or resume its exact partition, without scores."""
    if source.get('schema') == SCHEMA:
        validate_record(source)
        return copy.deepcopy(source['scope']), copy.deepcopy(source['cover'])
    if source.get('schema') == POINT_SCHEMA:
        scope = source['scope']
    elif source.get('schema') == CONTEXT_SCHEMA:
        scope = source
    else:
        raise ValueError('refined hill context, point receipt, or hill cover required')
    validate_scope(scope)
    scope = {key: copy.deepcopy(scope[key]) for key in ('schema', *cell_search.SCOPE_FIELDS)}
    cover = dict(leaves={}, unresolved={'': dict(cell=copy.deepcopy(scope['root']))}, visited=0)
    return scope, cover


def _authenticate_counts(intersection):
    if intersection:
        from outer_hill_intersection import authenticated_bch_cdf
    else:
        from outer_hill_incidence import authenticated_bch_cdf
    return authenticated_bch_cdf()


def _construct_model(scope, mixture, precision):
    parent = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(parent))
    sys.path.insert(0, str(parent/'gf16_packets'))
    import shared_mixture
    import scalar_cover as sc
    import birth_classes
    data = birth_classes.actual(scope['updates'])
    ctx.prec = precision  # actual() can reset the global working precision.
    return sc.Model(shared_mixture.as_components(mixture), data, scope['threshold'],
        scope['minimum_groups'], Q(3, 16), inner=birth_classes,
        variance_shuffle=True, variance_bins=scope['variance_bins'], regional_count=True)


def fresh_model(scope, precision):
    """Rebuild an unwrapped model; no search cache or old endpoints enter."""
    intersection = validate_scope(scope)
    if type(precision) is not int or precision < 128:
        raise ValueError('working precision >=128 required')
    caps, premises = _authenticate_counts(intersection)
    from monotone import transport_shells
    from local_models import full_block
    comparison = transport_shells(premises['canonical_cdf'], full_block(8))
    if (premises != scope['outer_premises']
            or cell_search.dense.fingerprint(caps) != scope['expected_cdf_sha256']
            or cell_search.dense.fingerprint(comparison) != scope['comparison_caps_sha256']):
        raise ValueError('fresh count premises differ from the hill context')
    # exact_mixture imports the legacy GF16 comparison module before the
    # inner factory runs. A cold CLI must establish both roots here.
    parent = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(parent))
    sys.path.insert(0, str(parent/'gf16_packets'))
    mixture, _ = cell_search.dense.exact_mixture(comparison, scope['mixture'])
    model = _construct_model(scope, mixture, precision)
    if (tuple(map(Q, scope['root'])) != tuple(model.root)
            or model.data.get('updates') != scope['updates']):
        raise ValueError('fresh domain or actual update count differs from scope')
    return model


def select_paths(cover, paths=None, *, shard_index=None, shard_count=None):
    """Select exact pending roots; shards are contiguous in sorted path order."""
    frontier = sorted(cover['unresolved'])
    if paths is not None:
        if (shard_index is not None or shard_count is not None
                or not isinstance(paths, (list, tuple)) or not paths
                or any(type(path) is not str for path in paths)
                or len(set(paths)) != len(paths) or not set(paths) <= set(frontier)):
            raise ValueError('distinct current pending paths, without shard flags, required')
        return sorted(paths)
    if shard_index is None and shard_count is None:
        return frontier
    if (type(shard_count) is not int or type(shard_index) is not int
            or not 1 <= shard_count <= len(frontier) or not 0 <= shard_index < shard_count):
        raise ValueError('valid shard index/count with no empty shards required')
    return frontier[len(frontier)*shard_index//shard_count:
                    len(frontier)*(shard_index+1)//shard_count]


def search(model, scope, cover, *, precision=256, target_bits=52, max_cells=200,
           max_depth=22, checkpoint_every=10, checkpoint=None, selected_paths=None):
    """Search selected pending roots while retaining the entire global partition."""
    validate_cover(scope, cover)
    assigned = select_paths(cover, selected_paths)
    if (type(precision) is not int or precision < 128
            or type(target_bits) is not int or target_bits < 20
            or type(max_cells) is not int or max_cells < 0
            or type(max_depth) is not int or not 1 <= max_depth <= 64
            or type(checkpoint_every) is not int or checkpoint_every < 1):
        raise ValueError('valid integer precision, target, and work limits required')
    if (tuple(model.root) != tuple(map(Q, scope['root']))
            or model.data.get('updates') != scope['updates']):
        raise ValueError('search model differs from the declared construction')
    previous = copy.deepcopy(cover)
    untouched = {path: row for path, row in previous['unresolved'].items() if path not in assigned}
    model.proposal_stop_bits = target_bits+2
    if not previous['unresolved']:
        if checkpoint:
            checkpoint(previous)
        return previous

    def combine(partial):
        cell_search.validate_cover(scope['root'], assigned, partial)
        result = dict(leaves={**copy.deepcopy(previous['leaves']), **partial['leaves']},
            unresolved={**copy.deepcopy(untouched), **partial['unresolved']},
            visited=previous.get('visited', 0)+partial['visited'])
        validate_cover(scope, result)
        return result

    def save(partial):
        if checkpoint:
            checkpoint(combine(partial))

    result = cell_search.search(model, assigned, precision=precision,
        target_bits=target_bits, max_cells=max_cells, max_depth=max_depth,
        checkpoint_every=checkpoint_every, checkpoint=save)
    return combine(result)


def merge_records(snapshot, source, partials):
    """Merge disjoint selected-subtree searches from one frozen source, not proofs.

    ``source`` is the parsed source's absolute path/SHA256 metadata. Callers
    must hash the actual bytes and keep all input bytes unchanged while merging.
    Partial selections need not exhaust the frontier: all other roots remain pending.
    """
    validate_record(snapshot)
    if (not isinstance(source, dict) or set(source) != {'path', 'sha256'}
            or not isinstance(source['path'], str) or not Path(source['path']).is_absolute()
            or not isinstance(source['sha256'], str) or len(source['sha256']) != 64
            or set(source['sha256'])-set('0123456789abcdef')):
        raise ValueError('frozen source absolute path and SHA256 required')
    base = snapshot['cover']
    leaves, pending = copy.deepcopy(base['leaves']), copy.deepcopy(base['unresolved'])
    assigned_all = set()
    visits = base.get('visited', 0)
    selections = []
    search_targets = [snapshot['cell_target_bits']]
    for partial in partials:
        validate_record(partial)
        if (partial.get('source') != source or partial['scope'] != snapshot['scope']
                or partial['precision'] != snapshot['precision']):
            raise ValueError('partial frozen source, exact scope, or precision differs')
        selection = partial.get('search_selection', {})
        roots = select_paths(base, selection.get('assigned_roots', []))
        if selection.get('source_frontier') != sorted(base['unresolved']):
            raise ValueError('partial source frontier differs')
        if assigned_all.intersection(roots):
            raise ValueError('selected source roots overlap across partials')
        fragment = dict(leaves={}, unresolved={}, visited=0)
        for name in ('leaves', 'unresolved'):
            outside = {}
            for path, row in partial['cover'][name].items():
                if any(path.startswith(root) for root in roots):
                    fragment[name][path] = copy.deepcopy(row)
                else:
                    outside[path] = row
            expected = (base['leaves'] if name == 'leaves' else
                {path: row for path, row in base['unresolved'].items() if path not in roots})
            if outside != expected:
                raise ValueError('partial changed rows outside its selected subtrees')
        cell_search.validate_cover(snapshot['scope']['root'], roots, fragment)
        work = partial['cover'].get('visited', 0)-base.get('visited', 0)
        if work < 0:
            raise ValueError('partial work count precedes frozen source')
        visits += work
        assigned_all.update(roots)
        selections.append(roots)
        search_targets.append(partial['cell_target_bits'])
        for root in roots:
            del pending[root]
        leaves.update(fragment['leaves'])
        pending.update(fragment['unresolved'])
    # Old and newly found endpoints/scores are deliberately discarded alike.
    def clean(rows, accepted):
        return {path: dict(cell=list(map(str, geometry.path_cell(snapshot['scope']['root'], path))),
            **({'witness': copy.deepcopy(row['witness'])} if accepted else {}))
            for path, row in rows.items()}
    result = dict(schema=SCHEMA, scope=copy.deepcopy(snapshot['scope']),
        precision=snapshot['precision'], cell_target_bits=min(search_targets),
        whole_code_certificate=False, source=copy.deepcopy(source),
        cover=dict(leaves=clean(leaves, True), unresolved=clean(pending, False), visited=visits),
        merge=dict(assigned_roots=sorted(assigned_all), selections=selections,
            source_and_partial_cell_targets=search_targets, numerical_evaluations=0),
        complete_search_partition=not pending, final_replay_required=True,
        proof_status='Merged search witnesses only; independent fresh replay is mandatory.')
    validate_record(result)
    return result


def _replay_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def _replay_plain(value):
    if value is None or type(value) in (str, int, bool):
        return value
    if isinstance(value, list):
        return [_replay_plain(item) for item in value]
    if isinstance(value, dict) and all(type(key) is str for key in value):
        return {key: _replay_plain(item) for key, item in value.items()}
    raise ValueError('worker witnesses must contain exact JSON parameters only')


def _check_replay_model(model, scope, precision):
    if (ctx.prec != precision or tuple(model.root) != tuple(map(Q, scope['root']))
            or type(model.data.get('updates')) is not int
            or model.data['updates'] != scope['updates'] or model.data.get('windows') != 32
            or model.threshold != scope['threshold'] or model.q_min != scope['minimum_groups']
            or model.tilt != Q(scope['base_tilt']) or model.variance_bins != scope['variance_bins']
            or model.variance_shuffle is not True or model.regional_count is not True):
        raise ArithmeticError('fresh worker model differs from its exact replay scope')


def _initialize_replay_worker(scope, precision, logdir):
    global _worker_model, _worker_scope, _worker_precision, _worker_log
    _worker_scope, _worker_precision = copy.deepcopy(scope), precision
    _worker_log = (open(Path(logdir)/f'worker-{os.getpid()}.log', 'x', buffering=1)
                   if logdir is not None else None)
    with redirect_stdout(_worker_log) if _worker_log is not None else nullcontext():
        # Every spawned process authenticates the outer counts, verifies the
        # rational mixture, and builds its own actual R2/R3/R4 inner operators.
        _worker_model = fresh_model(_worker_scope, precision)
        import scalar_cover
        if type(_worker_model) is not scalar_cover.Model:
            raise ArithmeticError('worker replay requires the ordinary unwrapped model')
        _check_replay_model(_worker_model, _worker_scope, precision)


def _replay_expected(scope, precision, task):
    path, cell, witness = task
    if (type(path) is not str or set(path)-{'0', '1'} or len(path) > 64
            or tuple(map(Q, cell)) != geometry.path_cell(scope['root'], path)
            or not isinstance(witness, dict)):
        raise ValueError('assigned binary path, exact geometry, and witness required')
    return dict(path=path, cell=list(map(str, map(Q, cell))), precision=precision,
        updates=scope['updates'], scope_sha256=_replay_digest(scope),
        witness_sha256=_replay_digest(witness))


def _evaluate_replay(model, scope, precision, task):
    expected = _replay_expected(scope, precision, task)
    cell, witness = tuple(map(Q, task[1])), copy.deepcopy(task[2])
    ctx.prec = precision
    _check_replay_model(model, scope, precision)
    upper = model.outward(cell, witness)
    _check_replay_model(model, scope, precision)
    if (not upper.is_finite() or not upper > 0
            or _replay_digest(witness) != expected['witness_sha256']
            or _replay_digest(scope) != expected['scope_sha256']):
        raise ArithmeticError('finite positive fresh worker bound and unchanged parameters required')
    return dict(expected, upper=[int(value) for value in upper.upper().man_exp()])


def _evaluate_replay_worker(task):
    with redirect_stdout(_worker_log) if _worker_log is not None else nullcontext():
        result = _evaluate_replay(_worker_model, _worker_scope, _worker_precision, task)
        print('FRESH HILL WORKER checked', result['path'], flush=True)
        return result


def _checked_replay_rows(scope, precision, tasks, rows):
    """Validate ordered one-to-one worker results before any endpoint is used."""
    rows = iter(rows)
    for task in tasks:
        expected = _replay_expected(scope, precision, task)
        row = next(rows, None)
        if (not isinstance(row, dict) or any(row.get(key) != value for key, value in expected.items())
                or type(row.get('precision')) is not int or type(row.get('updates')) is not int):
            raise ArithmeticError('missing, duplicate, reordered, or mismatched worker result')
        endpoint = row.get('upper')
        if (not isinstance(endpoint, list) or len(endpoint) != 2
                or any(type(value) is not int for value in endpoint) or endpoint[0] <= 0):
            raise ArithmeticError('positive exact worker dyadic endpoint required')
        yield dict(expected, upper=list(endpoint))
    if next(rows, None) is not None:
        raise ArithmeticError('extra worker replay result')


def _parallel_replay(record, cells, precision, target_bits, checkpoint, workers, logdir):
    scope = _replay_plain(record['scope'])
    tasks = [(path, list(map(str, cells[path])), _replay_plain(row['witness']))
             for path, row in sorted(record['cover']['leaves'].items())]
    if logdir is not None:
        directory = Path(logdir)
        if directory.exists():
            raise ValueError('new worker log directory required')
        directory.mkdir(parents=True)
        logdir = str(directory.resolve())
    result = dict(precision=precision, target_bits=target_bits, checked=[], workers=workers,
        complete_dense=False, passed=False, whole_code_certificate=False, scope=copy.deepcopy(scope),
        note='Fresh independent uncached workers; dense occupancies only, not a whole-code claim.')
    if checkpoint:
        checkpoint(copy.deepcopy(result))
    old_path = sys.path[:]
    thread_keys = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')
    old_environment = {key: os.environ.get(key) for key in thread_keys}
    executor = None
    total = Q(0)
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        os.environ.update({key: '1' for key in thread_keys})
        executor = ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn'),
            initializer=_initialize_replay_worker, initargs=(scope, precision, logdir))
        rows = executor.map(_evaluate_replay_worker, tasks, chunksize=1)
        for row in _checked_replay_rows(scope, precision, tasks, rows):
            total += cell_search.dense.dyadic(row['upper'])
            result['checked'].append(row)
            ctx.prec = precision
            bound = arb(total.numerator)/total.denominator
            result['aggregate_upper'] = [int(value) for value in bound.upper().man_exp()]
            result['log2_upper'] = str(bound.log()/arb(2).log())
            if checkpoint:
                checkpoint(copy.deepcopy(result))
    except BaseException:
        # Python3.14 can terminate running workers as well as cancel queued
        # cells. Older runtimes still cancel pending work and join safely.
        if executor is not None and hasattr(executor, 'terminate_workers'):
            executor.terminate_workers()
        raise
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        sys.path[:] = old_path
        for key, value in old_environment.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    if len(result['checked']) != len(cells):
        raise ArithmeticError('fresh workers did not cover the complete partition')
    result['complete_dense'] = True
    result['passed'] = (total < Q(2)**-target_bits
        and cell_search.dense.dyadic(result['aggregate_upper']) < Q(2)**-target_bits)
    if checkpoint:
        checkpoint(copy.deepcopy(result))
    return result


def replay(record, precision=384, target_bits=41, checkpoint=None, *, workers=1, logdir=None):
    """Authenticate and freshly replay every leaf of a complete dense cover."""
    cells = validate_record(record, complete=True)
    if (type(precision) is not int or precision < 128 or type(target_bits) is not int or target_bits < 1
            or type(workers) is not int or not 1 <= workers <= 4):
        raise ValueError('valid replay precision and aggregate target required')
    if workers > 1:
        return _parallel_replay(record, cells, precision, target_bits, checkpoint, workers, logdir)
    model = fresh_model(record['scope'], precision)
    result = dict(precision=precision, target_bits=target_bits, checked=[],
        complete_dense=False, passed=False, whole_code_certificate=False,
        scope=copy.deepcopy(record['scope']),
        note='Dense occupancies only; a matching sparse-prefix proof is still required.')
    total = Q(0)
    if checkpoint:
        checkpoint(copy.deepcopy(result))
    for path, row in sorted(record['cover']['leaves'].items()):
        ctx.prec = precision
        upper = model.outward(cells[path], copy.deepcopy(row['witness']))
        if ctx.prec != precision or not upper.is_finite() or not upper > 0:
            raise ArithmeticError('finite positive fresh interval endpoint at requested precision required')
        if (tuple(model.root) != tuple(map(Q, record['scope']['root']))
                or model.data.get('updates') != record['scope']['updates']):
            raise ArithmeticError('dense construction changed during replay')
        endpoint = [int(value) for value in upper.upper().man_exp()]
        total += cell_search.dense.dyadic(endpoint)
        result['checked'].append(dict(path=path, cell=list(map(str, cells[path])), upper=endpoint))
        # The exact endpoint sum determines the gate. A compact outward
        # envelope avoids enormous mantissas when cell exponents differ.
        bound = arb(total.numerator)/total.denominator
        result['aggregate_upper'] = [int(value) for value in bound.upper().man_exp()]
        result['log2_upper'] = str(bound.log()/arb(2).log())
        if checkpoint:
            checkpoint(copy.deepcopy(result))
    result['complete_dense'] = True
    result['passed'] = (total < Q(2)**-target_bits
        and cell_search.dense.dyadic(result['aggregate_upper']) < Q(2)**-target_bits)
    if checkpoint:
        checkpoint(copy.deepcopy(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--replay-workers', type=int, default=1,
        help='Fresh independent uncached replay processes, 1 through 4')
    parser.add_argument('--merge', type=Path, nargs='+',
        help='Merge search shards against the positional frozen source; no numerical replay')
    parser.add_argument('--hints', type=Path, nargs='+',
        help='Same-scope witness proposals for search only; saved bounds are ignored')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--cell-target-bits', type=int, default=52)
    parser.add_argument('--target-bits', type=int, default=41, help='Fresh replay aggregate target')
    parser.add_argument('--max-cells', type=int, default=200)
    parser.add_argument('--max-depth', type=int, default=22)
    parser.add_argument('--checkpoint-every', type=int, default=10)
    parser.add_argument('--initial-width', type=Q, help='Geometry-only presplit of unresolved cells')
    parser.add_argument('--paths', nargs='+', help='Exact pending paths after any initial presplit')
    parser.add_argument('--shard-index', type=int, help='Zero-based contiguous frontier shard')
    parser.add_argument('--shard-count', type=int)
    args = parser.parse_args()
    temporary = args.output.with_name(args.output.name+'.tmp')
    if (args.output.exists() or temporary.exists() or args.precision < 128
            or args.cell_target_bits < 20 or args.target_bits < 1
            or args.max_cells < 0 or not 1 <= args.max_depth <= 64 or args.checkpoint_every < 1
            or not 1 <= args.replay_workers <= 4 or not args.replay and args.replay_workers != 1
            or args.initial_width is not None and not 0 < args.initial_width <= 1
            or args.replay and any(value is not None for value in
                (args.initial_width, args.paths, args.shard_index, args.shard_count,
                 args.merge, args.hints))
            or args.merge is not None and any(value is not None for value in
                (args.initial_width, args.paths, args.shard_index, args.shard_count, args.hints))):
        parser.error('fresh output and valid proof/search limits required')
    raw = args.source.read_bytes()
    input_bytes = [(args.source.resolve(), raw)]
    hint_sources = None
    source = json.loads(raw)
    scope, cover = from_source(source)
    if args.replay:
        validate_record(source, complete=True)
    record = dict(schema=SCHEMA, scope=scope, cover=cover, precision=args.precision,
        cell_target_bits=args.cell_target_bits, whole_code_certificate=False,
        source=dict(path=str(args.source.resolve()), sha256=hashlib.sha256(raw).hexdigest()),
        proof_status='Dense search only; accepted witnesses require independent fresh replay.')

    def save():
        def check_inputs():
            for path, expected in input_bytes:
                if path.read_bytes() != expected:
                    raise ValueError('source or shard changed during this invocation')
            if hint_sources is not None:
                import hill_reuse
                hill_reuse.check_sources(hint_sources)
        check_inputs()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(json.dumps(record, indent=2)+'\n')
        check_inputs()
        temporary.replace(args.output)

    if args.merge is not None:
        partials, metadata = [], []
        seen_paths = {args.source.resolve()}
        for path in args.merge:
            path = path.resolve()
            if path in seen_paths:
                parser.error('distinct shard paths, separate from the frozen source, required')
            seen_paths.add(path)
            encoded = path.read_bytes()
            input_bytes.append((path, encoded))
            partials.append(json.loads(encoded))
            metadata.append(dict(path=str(path), sha256=hashlib.sha256(encoded).hexdigest()))
        record = merge_records(source, record['source'], partials)
        record['partial_sources'] = metadata
        save()
        return
    if args.replay:
        def replay_checkpoint(result):
            record['fresh_replay'] = result
            record['proof_status'] = 'Fresh dense replay only; no whole-code certificate.'
            save()
        replay(source, args.precision, args.target_bits, replay_checkpoint,
            workers=args.replay_workers,
            logdir=args.output.with_name(args.output.name+'.workers') if args.replay_workers > 1 else None)
        return
    if args.initial_width is not None:
        cover, record['search_preprocessing'] = geometry.presplit(
            scope['root'], cover, args.initial_width, args.max_depth)
        record['cover'] = cover
    try:
        selected = select_paths(cover, args.paths,
            shard_index=args.shard_index, shard_count=args.shard_count)
    except ValueError as error:
        parser.error(str(error))
    record['search_selection'] = dict(assigned_roots=selected,
        source_frontier=sorted(cover['unresolved']))
    model = fresh_model(scope, args.precision)
    if args.hints:
        import hill_reuse
        model, regional_cache, hint_sources = hill_reuse.wrap(
            model, scope, args.hints, args.cell_target_bits)
        record['hint_sources'] = copy.deepcopy(hint_sources)
        record['hint_policy'] = 'Witness proposals only; every accepted cell requires fresh outward evaluation.'
    def search_checkpoint(result):
        record['cover'] = result
        record['complete_search_partition'] = not result['unresolved']
        save()
    search(model, scope, cover, precision=args.precision, target_bits=args.cell_target_bits,
        max_cells=args.max_cells, max_depth=args.max_depth,
        checkpoint_every=args.checkpoint_every, checkpoint=search_checkpoint,
        selected_paths=selected or None)


if __name__ == '__main__':
    main()
