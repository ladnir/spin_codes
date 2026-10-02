"""Use old GL32 point witnesses as search hints, never as bound receipts.

The caller authenticates and constructs the current scalar model first.
Only rational witness parameters cross from a point file to that model;
all bounds, including the acceptance test, are evaluated again. The model
must remain unchanged while wrapped. No cache survives the current process.
"""
import copy
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import sys

from flint import arb, ctx

PARENT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PARENT))
sys.path.insert(0, str(PARENT / 'gf16_packets'))
from shared_relaxed_point_reuse import nearest_witness, rebase_witness

POINT_SCHEMA = 'canonical-packed-gl32-dense-points-1'
WITNESS_FIELDS = frozenset((
    'tilt', 'parameters', 'variance_dual', 'weights_dual',
    'variance_partition', 'regional_count_parts', 'regional_tilted_atom',
    'regional_fine_tilts', 'regional_tilted_variance', 'regional_direct_counts',
    'regional_exact_zero', 'regional_lazy_density_through',
    'regional_joint_return_through', 'regional_feedback_classes_from',
    'regional_feedback_classes_through'))


def point_proposals(record):
    """Extract proposals from the known R2 point scope, ignoring old bounds.

    Distance, mixture, old comparison caps, and numerical scores/endpoints
    are intentionally not transferred. They need not match the new model.
    """
    if (not isinstance(record, dict) or record.get('schema') != POINT_SCHEMA
            or type(record.get('updates')) is not int or record['updates'] != 2
            or type(record.get('minimum_groups')) is not int or record['minimum_groups'] != 33
            or record.get('maximum_groups') != 2048
            or record.get('last_lp') != 104 or record.get('refined') is not True
            or Q(record.get('base_tilt', 0)) != Q(3, 16)
            or record.get('variance_shuffle') is not True
            or record.get('variance_bins') != 16 or record.get('regional_count') is not True
            or not 0 < Q(record.get('distance', 0)) < Q(1, 2)):
        raise ValueError('canonical GL32 R2 point scope with occupancies 33..2048 required')
    rows = record.get('probes')
    if not isinstance(rows, list) or not rows:
        raise ValueError('nonempty list of point witnesses required')
    result = []
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get('witness'), dict):
            raise ValueError('each point must provide a witness, not an empty-cell claim')
        witness = row['witness']
        if (set(witness) - WITNESS_FIELDS
                or not {'tilt', 'parameters', 'variance_dual'} <= set(witness)
                or Q(witness['tilt']) != Q(3, 16)):
            raise ValueError('known base-tilt rational witness fields required')
        mean = Q(row['mean'])
        # Even before choosing a target cell, check the source partition.
        clean = rebase_witness((mean, mean), (mean, mean), witness)
        result.append(dict(mean=str(mean), witness=clean))
    return result


def load_points(paths):
    """Return (proposal rows, source hashes) for one or more point files."""
    paths = list(paths)
    if not paths:
        raise ValueError('at least one point source path required')
    probes, sources = [], []
    for path in paths:
        path = Path(path)
        raw = path.read_bytes()
        rows = point_proposals(json.loads(raw))
        probes.extend(rows)
        sources.append(dict(path=str(path.resolve()),
            sha256=hashlib.sha256(raw).hexdigest(), point_count=len(rows)))
    return probes, sources


def _cell(cell):
    cell = tuple(map(Q, cell))
    if len(cell) != 2 or not 0 <= cell[0] <= cell[1] <= 1:
        raise ValueError('ordered mean cell in [0,1] required')
    return cell


def _witness_key(witness):
    return json.dumps(witness, sort_keys=True, separators=(',', ':'), allow_nan=False)


class PointReuseModel:
    """Proxy a fresh scalar model; optimize only proposal selection.

    A successful hint must freshly beat target_bits+2, as required by
    scalar_cover.run. The immediately following outward call can consume
    that same result once, but only for identical cell/witness/precision.
    Any mismatch or subsequent proposal discards this one-entry cache.
    Replay in another process always performs a new outward evaluation.

    By default the last accepted proposal supplies the first hint for the
    next cell, followed by at most one distinct nearest-point hint. Only
    its source cell and copied witness survive, not its numerical bound.
    Standalone replay calls to outward never update this rolling hint.
    An optional regional_warm(witness) predicate permits a rolling regional
    hint to skip floating screening when its exact local operator is cached.
    It never substitutes for fresh exact evaluation on the target cell.
    maximum_optimization_width optionally defers expensive fallback search
    on broader cells: failed hints then request bisection without claiming
    any bound. None preserves the ordinary fallback at every width.
    R2 remains the default. R3/R4 require an explicit updates argument that
    matches the fresh actual-data object; witnesses never select the variant.
    """
    def __init__(self, model, probes, target_bits, *, outward=None, regional_proposal=None,
                 rolling=True, regional_warm=None, maximum_optimization_width=None, updates=2):
        if type(target_bits) is not int or target_bits < 1:
            raise ValueError('positive integer target bits required')
        if type(updates) is not int or updates not in (2, 3, 4):
            raise ValueError('explicit actual inner updates must be 2, 3, or 4')
        if (model.q_min != 33 or model.tilt != Q(3, 16)
                or not model.variance_shuffle or model.variance_bins != 16
                or not model.regional_count or model.inner.__name__ != 'birth_classes'
                or model.data.get('windows') != 32 or type(model.data.get('updates')) is not int
                or model.data['updates'] != updates):
            raise ValueError('fresh canonical GL32 scalar model with matching actual updates required')
        if not isinstance(probes, list) or not probes:
            raise ValueError('nonempty point-proposal list required')
        if outward is not None and not callable(outward):
            raise ValueError('outward override must be a fresh bound evaluator')
        if regional_proposal is not None and not callable(regional_proposal):
            raise ValueError('regional proposal override must be callable')
        if type(rolling) is not bool:
            raise ValueError('boolean rolling-reuse option required')
        if regional_warm is not None and not callable(regional_warm):
            raise ValueError('regional warmth predicate must be callable')
        if maximum_optimization_width is not None:
            if isinstance(maximum_optimization_width, bool):
                raise ValueError('optimization width must be a positive rational at most one')
            maximum_optimization_width = Q(maximum_optimization_width)
            if not 0 < maximum_optimization_width <= 1:
                raise ValueError('optimization width must be a positive rational at most one')
        self._model = model
        self._evaluate = model.outward if outward is None else outward
        self._regional_proposal = regional_proposal
        self._regional_warm = regional_warm
        self._maximum_optimization_width = maximum_optimization_width
        self._probes = copy.deepcopy(probes)
        self._target_bits = target_bits
        self._pending = None
        self._awaiting = None
        self._rolling_enabled = rolling
        self._rolling = None
        self.stats = dict(attempts=0, accepted=0, fallback=0, cache_hits=0,
            regional_screened=0, regional_rejected=0, rolling_attempts=0,
            rolling_accepted=0, rolling_updates=0, duplicate_hints=0,
            warm_exact_attempts=0, warm_exact_accepted=0, warm_exact_rejected=0,
            wide_splits=0)
        # This controls the ordinary fallback search, not proof acceptance.
        model.proposal_stop_bits = target_bits + 2

    def __getattr__(self, name):
        return getattr(self._model, name)

    def _fresh(self, cell, witness):
        precision = ctx.prec
        key = _witness_key(witness)
        upper = self._evaluate(cell, witness)
        if ctx.prec != precision or _witness_key(witness) != key:
            raise ArithmeticError('outward evaluation changed precision or witness')
        if not upper.is_finite() or not upper > 0:
            raise ArithmeticError('positive finite fresh outward bound required')
        return upper

    def _hints(self, cell):
        if self._rolling is not None:
            source, witness = self._rolling
            yield 'rolling', rebase_witness(source, cell, witness)
        # Lazy: a successful rolling hint never needs to prepare this one.
        _, _, witness = nearest_witness(self._probes, cell)
        yield 'point', witness

    def _issued(self, cell, score, witness):
        # Only the matching subsequent outward check can promote a hint.
        # A caller's unrelated checkpoint replay must not seed this slot.
        if math.isfinite(score) and score < -(self._target_bits + 2):
            self._awaiting = (cell, _witness_key(witness), ctx.prec)
        return score, witness

    def _is_warm(self, witness):
        if self._regional_warm is None:
            return False
        precision, key = ctx.prec, _witness_key(witness)
        warm = self._regional_warm(witness)
        if ctx.prec != precision or _witness_key(witness) != key:
            raise ArithmeticError('warmth check changed precision or witness')
        if type(warm) is not bool:
            raise ValueError('regional warmth predicate must return a boolean')
        return warm

    def proposal(self, cell):
        self._pending = None
        self._awaiting = None
        cell = _cell(cell)
        if 0 < cell[0] <= cell[1] < 1:
            seen = set()
            for kind, witness in self._hints(cell):
                key = _witness_key(witness)
                if key in seen:
                    self.stats['duplicate_hints'] += 1
                    continue
                seen.add(key)
                self.stats['attempts'] += 1
                if kind == 'rolling':
                    self.stats['rolling_attempts'] += 1
                warm = False
                if 'regional_count_parts' in witness:
                    warm = kind == 'rolling' and self._is_warm(witness)
                    if warm:
                        # All target-cell count/variance terms are still
                        # checked exactly below; only the float screen goes.
                        self.stats['warm_exact_attempts'] += 1
                    else:
                        # Saved complete MGF families make this a fixed-witness
                        # floating evaluation, not another LP optimization search.
                        from regional_count import propose
                        evaluator = self._regional_proposal or propose
                        precision = ctx.prec
                        score, witness = evaluator(self._model, cell, witness)
                        if ctx.prec != precision or not math.isfinite(score):
                            raise ArithmeticError('finite regional screening score at unchanged precision required')
                        self.stats['regional_screened'] += 1
                        if not score < -(self._target_bits + 2):
                            self.stats['regional_rejected'] += 1
                            continue
                upper = self._fresh(cell, witness)
                if upper < arb(2)**-(self._target_bits + 2):
                    score = float(upper.log()/arb(2).log())
                    # Avoid a threshold-rounding disagreement with scalar_cover.run.
                    if score < -(self._target_bits + 2):
                        self._pending = (cell, _witness_key(witness), ctx.prec, arb(upper))
                        self.stats['accepted'] += 1
                        if kind == 'rolling':
                            self.stats['rolling_accepted'] += 1
                        if warm:
                            self.stats['warm_exact_accepted'] += 1
                        return self._issued(cell, score, witness)
                if warm:
                    self.stats['warm_exact_rejected'] += 1
        if (self._maximum_optimization_width is not None
                and cell[1]-cell[0] > self._maximum_optimization_width):
            self.stats['wide_splits'] += 1
            # This is a scheduler response, not a numerical bound. The
            # scalar driver bisects (or records unresolved at its depth cap).
            return 0., {}
        # Endpoint cells are outside the rebasing helper's strict interior.
        self.stats['fallback'] += 1
        return self._issued(cell, *self._model.proposal(cell))

    def outward(self, cell, witness):
        cell = _cell(cell)
        pending, self._pending = self._pending, None
        awaiting, self._awaiting = self._awaiting, None
        key = (cell, _witness_key(witness), ctx.prec)
        if pending is not None and pending[:3] == key:
            self.stats['cache_hits'] += 1
            upper = arb(pending[3])
        else:
            upper = self._fresh(cell, witness)
        if (self._rolling_enabled and awaiting == key and 0 < cell[0] <= cell[1] < 1
                and upper < arb(2)**-self._target_bits):
            self._rolling = (cell, copy.deepcopy(witness))
            self.stats['rolling_updates'] += 1
        return upper
