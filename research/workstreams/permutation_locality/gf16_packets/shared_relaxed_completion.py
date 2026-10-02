"""Optional completion-driven scheduler for shared GF16 interval search.

The caller supplies an executor initialized with this module's initialize
wrapper. It delegates to the unchanged core initializer through its named
module, avoiding duplicate __main__ globals under Windows spawn. Only
scheduling changes: the existing worker evaluates every complete interval,
and the final whole-code verifier is unchanged.
Checkpoints retain pending and in-flight cells as unresolved, so interruption
never removes a cell from the exhaustive partition.
"""
from concurrent.futures import FIRST_COMPLETED, Future, wait
import copy
from fractions import Fraction as Q
import hashlib
import heapq
import json
from math import isfinite

import shared_relaxed_parallel as core


def initialize(components, record, index, precision, log_directory):
    """Initialize the same named module used by evaluate_worker after spawn."""
    core.initialize(components, record, index, precision, log_directory)


def scope(model, precision):
    """Bind an invocation to the checked parent's complete proof setting."""
    if (type(precision) is not int or precision < 256
            or type(model.threshold) is not int or not 0 <= model.threshold < 1 << 21
            or type(model.q_min) is not int or not 1 <= model.q_min <= 2048
            or type(model.variance_bins) is not int or not 1 <= model.variance_bins <= 64
            or model.variance_shuffle is not True or model.regional_count is not True
            or model.inner.__name__ != 'birth_classes'
            or any(type(model.data.get(k)) is not int or model.data[k] != v
                   for k, v in dict(bits=19, windows=32, updates=2).items())):
        raise ValueError('the checked shared GF16 R2 regional model is required')
    root = tuple(map(Q, model.root))
    if len(root) != 2 or not 0 <= root[0] < root[1] <= 1 or Q(model.tilt) <= 0:
        raise ValueError('nondegenerate root interval and positive base tilt required')
    return dict(ensemble=core.proof.ENSEMBLE,
        component_sha256=core.component_digest(model.components),
        precision=precision, threshold=model.threshold, minimum_groups=model.q_min,
        root=list(map(str, root)), base_tilt=str(Q(model.tilt)),
        state_bits=19, packets_per_step=32, updates=2, inner='birth_classes',
        variance_shuffle=True, variance_bins=model.variance_bins, regional_count=True)


def evaluate_worker(request):
    """Tag the unchanged full-interval worker after checking its model."""
    expected = scope(core._model, core._precision)
    if (not isinstance(request, dict) or set(request) != {'id', 'scope', 'job'}
            or type(request['id']) is not int or request['id'] < 1
            or request['scope'] != expected):
        raise ArithmeticError('completion worker request/model mismatch')
    checked_job(request)
    result = core.evaluate_worker(request['job'])
    return dict(result, request_id=request['id'], proof_scope=expected)


def checked_job(request):
    """Reconstruct the full interval from its path, never from point metadata."""
    job = request.get('job')
    if not isinstance(job, (tuple, list)) or len(job) != 4:
        raise ArithmeticError('a complete interval job is required')
    path, cell, target_bits, witness = job
    if (not isinstance(path, str) or len(path) > 64 or set(path)-{'0','1'}
            or type(target_bits) is not int or target_bits < 20
            or (witness is not None and not isinstance(witness, dict))):
        raise ArithmeticError('invalid interval path, margin, or retained witness')
    expected = tuple(map(Q, request['scope']['root']))
    for step in path:
        mid = sum(expected)/2
        expected = (expected[0],mid) if step == '0' else (mid,expected[1])
    if cell != expected:
        raise ArithmeticError('job must cover its complete path interval')


def checked_result(result, request):
    """Validate the job's whole-interval scope and exact positive endpoint."""
    checked_job(request)
    path, cell, target_bits, _ = request['job']
    if (not isinstance(result, dict) or result.get('request_id') != request['id']
            or type(result.get('request_id')) is not int
            or result.get('proof_scope') != request['scope']
            or result.get('path') != path or result.get('cell') != list(map(str, cell))
            or result.get('component_sha256') != request['scope']['component_sha256']
            or not isinstance(result.get('witness'), dict)
            or type(result.get('proposal')) not in (int, float)
            or not isfinite(result['proposal']) or 'upper' not in result):
        raise ArithmeticError('completion worker result/scope mismatch')
    if result['upper'] is not None:
        try:
            value = core.proof.dyadic(result['upper'])
        except (ValueError, TypeError, KeyError, OverflowError) as error:
            raise ArithmeticError('invalid completion worker endpoint') from error
        if not 0 < value < Q(2)**-target_bits:
            raise ArithmeticError('completion worker endpoint misses the requested margin')
    return result


def search(model, cover, executor, workers, max_cells, max_depth, target_bits,
           precision, checkpoint, reuse_parent_witness=False, *,
           retain_search_leaves=False, wait_for=wait):
    """Keep at most workers jobs in flight and consume whichever finishes.

    max_cells bounds submissions in this invocation, including retained-leaf
    replay. visited counts validated completions, as in the existing search.
    A failed or cancelled invocation leaves every unconsumed cell unresolved.
    The executor belongs to the caller and is never shut down by this helper.
    Opt-in retained leaves are saved search hints, not fresh certificates.
    The whole-code verifier must freshly evaluate them before any proof claim.
    """
    if (type(workers) is not int or not 1 <= workers <= 4
            or type(max_cells) is not int or max_cells < 1
            or type(max_depth) is not int or not 1 <= max_depth <= 64
            or type(target_bits) is not int or target_bits < 20
            or type(reuse_parent_witness) is not bool
            or type(retain_search_leaves) is not bool
            or not isinstance(cover, dict) or not isinstance(cover.get('leaves'), dict)
            or not isinstance(cover.get('unresolved'), dict)
            or not callable(checkpoint) or not callable(wait_for)
            or not callable(getattr(executor, 'submit', None))):
        raise ValueError('valid search limits, executor, callbacks, and partition required')
    invocation = scope(model, precision)
    cells = core.proof.sc.partition(model, cover['leaves'], cover['unresolved'])
    previous = cover.get('visited', 0)
    if type(previous) is not int or previous < 0 or not cells:
        raise ValueError('nonnegative previous work and nonempty full partition required')
    retained = {}
    source_sha256 = None
    source_scope = None
    if retain_search_leaves:
        if cover.get('component_sha256') != invocation['component_sha256']:
            raise ValueError('retained search hints require the same comparison digest')
        source_scheduler = cover.get('completion_scheduler', {})
        if not isinstance(source_scheduler, dict):
            raise ValueError('retained completion metadata must be a mapping')
        source_scope = source_scheduler.get('proof_scope')
        if source_scope is not None and source_scope != invocation:
            raise ValueError('retained completion hints require the same recorded proof scope')
        for path, row in cover['leaves'].items():
            if (not isinstance(row, dict) or not isinstance(row.get('witness'), dict)
                    or not row['witness'] or type(row.get('proposal')) not in (int, float)
                    or not isfinite(row['proposal'])):
                raise ValueError('retained source leaves need witnesses and finite search proposals')
            if row['proposal'] < -target_bits:
                # Deliberately discard saved endpoints. This is an untrusted
                # search label, checked afresh only by the final verifier.
                retained[path] = dict(witness=copy.deepcopy(row['witness']),
                    proposal=row['proposal'], retained_search_hint=True)
        source_sha256 = hashlib.sha256(json.dumps(cover,sort_keys=True,
            separators=(',', ':'),allow_nan=False).encode()).hexdigest()
    saved = copy.deepcopy(cover['leaves'])
    if reuse_parent_witness:
        saved.update(copy.deepcopy(cover['unresolved']))
    pending = [(len(path), path, cell, saved.get(path, {}).get('witness'))
               for path, cell in cells.items() if path not in retained]
    heapq.heapify(pending)
    leaves, unresolved, in_flight = dict(retained), {}, {}
    seen_futures = set()
    submitted = completed = 0
    latest = None

    def snapshot():
        waiting = dict(unresolved)
        for _, path, cell, witness in pending:
            if path in waiting:
                raise ArithmeticError('duplicate pending path')
            waiting[path] = dict(cell=list(map(str, cell)),
                                **({'witness':witness} if witness is not None else {}))
        for item, _ in in_flight.values():
            _, path, cell, witness = item
            if path in waiting:
                raise ArithmeticError('in-flight cell overlaps another unresolved cell')
            waiting[path] = dict(cell=list(map(str, cell)),
                                **({'witness':witness} if witness is not None else {}))
        state = dict(leaves=leaves, unresolved=waiting, visited=previous+completed,
            component_sha256=invocation['component_sha256'],
            completion_scheduler=dict(submitted=submitted, completed=completed,
                in_flight=sorted(item[1] for item, _ in in_flight.values()),
                workers=workers, max_cells=max_cells, proof_scope=invocation,
                retain_search_leaves=retain_search_leaves,
                retained_search_leaf_count=len(retained),
                retained_search_leaf_paths=sorted(retained),
                freshly_checked_leaf_count=len(leaves)-len(retained),
                source_cover_sha256=source_sha256, source_visited=previous,
                source_recorded_proof_scope=source_scope,
                status=('Search checkpoint only. Retained leaf labels were NOT freshly '
                    'verified in this invocation; saved endpoints were discarded. The '
                    'whole-code verifier must freshly replay EVERY leaf before a proof claim.'
                    if retained else 'Search checkpoint only. Every accepted leaf was '
                    'freshly checked in this invocation; whole-code replay remains required.')))
        core.proof.sc.partition(model, state['leaves'], state['unresolved'])
        return copy.deepcopy(state)

    def publish():
        nonlocal latest
        latest = snapshot()
        # A user callback must not mutate search state or the returned snapshot.
        checkpoint(copy.deepcopy(latest))

    def refill():
        nonlocal submitted
        while pending and len(in_flight) < workers and submitted < max_cells:
            item = heapq.heappop(pending)
            _, path, cell, old = item
            request = dict(id=submitted+1, scope=copy.deepcopy(invocation),
                           job=(path, cell, target_bits, copy.deepcopy(old)))
            try:
                future = executor.submit(evaluate_worker, request)
                if not isinstance(future, Future) or future in seen_futures:
                    raise ArithmeticError('executor must return one fresh Future per submission')
            except BaseException:
                heapq.heappush(pending, item)
                raise
            in_flight[future] = item, request
            seen_futures.add(future)
            submitted += 1

    try:
        # Default: all prior leaves await replay. Opt-in retained entries are
        # explicitly tagged hints, never advertised as fresh certificates.
        publish()
        refill()
        publish()
        while in_flight:
            done, _ = wait_for(tuple(in_flight), return_when=FIRST_COMPLETED)
            if (not done or not isinstance(done, (set, frozenset))
                    or not done.issubset(in_flight)
                    or any(not future.done() for future in done)):
                raise ArithmeticError('completion wait returned missing, foreign, or unfinished work')
            for future in sorted(done, key=lambda f:in_flight[f][1]['id']):
                item, request = in_flight[future]
                depth, path, cell, _ = item
                result = checked_result(future.result(), request)
                entry = dict(witness=copy.deepcopy(result['witness']), proposal=result['proposal'])
                children = []
                if result['upper'] is None and depth < max_depth:
                    mid = sum(cell)/2
                    for suffix, child in (('0', (cell[0], mid)), ('1', (mid, cell[1]))):
                        old = (core.restrict_witness(cell, child, result['witness'])
                               if reuse_parent_witness else None)
                        children.append((depth+1, path+suffix, child, old))
                # No partition mutation occurs until validation/restriction succeeds.
                del in_flight[future]
                if result['upper'] is not None:
                    leaves[path] = entry
                elif children:
                    for child in children:
                        heapq.heappush(pending, child)
                else:
                    unresolved[path] = dict(entry, cell=list(map(str, cell)))
                completed += 1
                refill()
                publish()
        return latest
    except BaseException:
        # Only this invocation's unconsumed futures are touched. Running jobs
        # may be uncancellable; their cells still remain in the checkpoint.
        for future in in_flight:
            future.cancel()
        try:
            publish()
        except BaseException:
            # Preserve the original error. The last successful checkpoint
            # already contains all work that had not yet been consumed.
            pass
        raise
