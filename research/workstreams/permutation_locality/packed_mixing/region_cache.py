"""One-entry cache of freshly built regional operators, for proof search only.

The caller keeps its authenticated actual-data object immutable. A retained
reference prevents identity reuse. No numerical disk receipt enters this
cache. Cell geometry, outer coefficients, and all MGF inequalities are
reevaluated by outward(); only the cell-independent inner polynomial is reused.
Final independent replay continues to use regional_count.outward directly.
"""
from fractions import Fraction as Q
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'gf16_packets'))

from flint import arb, arb_mat, ctx
import regional_count as rc
import variance_partition as variance

LOCAL_SELECTORS = (
    'regional_exact_zero', 'regional_joint_return_through', 'regional_lazy_density_through',
    'regional_feedback_classes_from', 'regional_feedback_classes_through',
    'regional_feedback_uniform_classes', 'regional_feedback_uniform_replace',
)
COUNT_SELECTORS = (
    'regional_count_parts', 'regional_direct_counts', 'regional_tilted_atom',
    'regional_fine_tilts', 'regional_tilted_variance',
)


def local_key(model, witness):
    """Reject unclassified selectors rather than silently reuse stale operators."""
    if not isinstance(model.data, dict) or not isinstance(witness, dict):
        raise ValueError('actual data dictionary and witness required')
    if any(k.startswith('regional_') and k not in LOCAL_SELECTORS+COUNT_SELECTORS for k in witness):
        raise ValueError('unknown regional selector requires cache-key classification')
    windows = model.data.get('windows')
    if type(windows) is not int or windows < 1:
        raise ValueError('positive local packet geometry required')
    parameters = witness.get('parameters')
    if not isinstance(parameters, list) or len(parameters) != 3 or Q(parameters[0]) <= 0:
        raise ValueError('positive output tilt and two outer duals required')
    for name in ('regional_exact_zero', 'regional_feedback_uniform_classes', 'regional_feedback_uniform_replace'):
        if type(witness.get(name, False)) is not bool:
            raise ValueError('boolean local selectors required')
    for name, limit in (('regional_joint_return_through', min(4, windows)),
                        ('regional_lazy_density_through', windows)):
        if name in witness and (type(witness[name]) is not int or not 0 <= witness[name] <= limit):
            raise ValueError('valid local occupancy cutoff required')
    rc.feedback_class_interval(witness, windows)
    selectors = tuple((name, name in witness, witness.get(name)) for name in LOCAL_SELECTORS)
    # Function identities also invalidate a cache if tests or research code
    # replace the construction functions within the current process.
    return id(model.data), ctx.prec, windows, Q(parameters[0]), selectors, rc.sc.G, rc.local_operators, rc.placement


class RegionalCache:
    def __init__(self):
        self._saved = None
        self.hits = 0
        self.misses = 0

    def clear(self):
        self._saved = None

    def is_warm(self, model, witness):
        """Check the local key without constructing or exposing coefficients.

        This is a search-cost hint, not a bound or an acceptance test. Key
        validation is performed even when cold; malformed selectors fail
        exactly as they would on an attempted polynomial lookup.
        """
        key = local_key(model, witness)
        if self._saved is None or self._saved[0] != key:
            return False
        if self._saved[1] is not model.data:
            raise ArithmeticError('actual-data identity changed')
        return True

    def _polynomial(self, model, witness):
        key = local_key(model, witness)
        if self._saved is not None and self._saved[0] == key:
            if self._saved[1] is not model.data:
                raise ArithmeticError('actual-data identity changed')
            self.hits += 1
            return self._saved[2]
        precision = ctx.prec
        polynomial = rc.placement(rc.local_operators(model, witness))
        if ctx.prec != precision or local_key(model, witness) != key:
            raise ArithmeticError('operator construction changed precision or selectors')
        if len(polynomial) != rc.sc.G+1 or not polynomial:
            raise ArithmeticError('complete regional polynomial required')
        size = polynomial[0].nrows()
        if size < 1 or any(m.nrows() != size or m.ncols() != size for m in polynomial):
            raise ArithmeticError('matching square regional matrices required')
        if any(not m[i, j].is_finite() or m[i, j] < 0
                for m in polynomial for i in range(size) for j in range(size)):
            raise ArithmeticError('finite nonnegative outward coefficients required')
        # The cache owns its matrices. The private evaluation path below only
        # uses nonmutating arithmetic and never exposes these references.
        self._saved = key, model.data, tuple(arb_mat(m) for m in polynomial)
        self.misses += 1
        return self._saved[2]

    def polynomial(self, model, witness):
        """Return independent copies when an external caller requests coefficients."""
        return [arb_mat(m) for m in self._polynomial(model, witness)]

    def bound(self, model, cell, witness):
        """Scalar-cover adapter preserving its front-door witness validation."""
        if not isinstance(witness, dict):
            raise ValueError('witness dictionary required')
        if any(k.startswith('regional_') and k not in LOCAL_SELECTORS+COUNT_SELECTORS for k in witness):
            raise ValueError('unknown regional selector requires cache-key classification')
        if 'regional_count_parts' not in witness:
            return model.outward(cell, witness)
        if len(witness['parameters']) != 3:
            raise ValueError('output tilt and two outer duals required')
        lam, eta, mu = map(Q, witness['parameters'])
        tilt = Q(witness['tilt'])
        if lam <= 0 or mu < 0:
            raise ValueError('positive tilt and nonnegative occupancy dual required')
        model.family(tilt)
        if tilt != model.tilt and 'weights_dual' not in witness:
            raise ValueError('alternate tilt missing weight duals')
        if model.variance_shuffle and tilt == model.tilt and 'variance_dual' not in witness:
            raise ValueError('base tilt missing variance dual')
        if not model.regional_count:
            raise ValueError('regional count witness requires enabled model')
        return self.outward(model, cell, witness)[0]

    def outward(self, model, cell, witness):
        """Same directed arithmetic as regional_count.outward, with cached placement."""
        parts, checked = rc.prepare_witness(model, cell, witness)
        lam = Q(witness['parameters'][0])
        if lam <= 0:
            raise ValueError('positive output tilt required')
        precision = ctx.prec
        region = self._polynomial(model, witness)
        direct = witness.get('regional_direct_counts', False)
        weights, _ = model.weights(cell, model.tilt)
        scale = Q(1) if direct else sum(weights)
        masses = ([rc.aq(model.tilt)**(-j) for j in range(rc.sc.G+1)] if direct else
                  rc.binomial_masses(rc.sc.G, weights[1]/scale))
        cs, _, _ = model.family(model.tilt)
        total = arb(0)
        size = region[0].nrows()
        for (interval, dual), part in zip(parts, checked['regional_count_parts']):
            eta, mu, gamma = dual
            cap = variance.factor(model, cell, interval[0], witness['variance_dual'])
            ratios = (rc.count_mass_caps if direct else rc.count_ratios)(
                model.features, model.active, cell, interval, model.q_min, rc.sc.G,
                rc.aq(cap), part['mgf_witnesses'])
            matrix = sum((mass*ratio*value for mass, ratio, value in zip(masses, ratios, region)), arb_mat(size, size))
            power = matrix**rc.sc.REGIONS
            moment = sum((power[0, j] for j in range(size)), arb(0))
            count = sum((rc.aq(c)*rc.aq(eta*f+mu*a+gamma*f*(1-f)).exp()
                         for c, f, a in zip(cs, model.features, model.active)), arb(0))
            exponent = rc.sc.G*(min(eta*x for x in cell)+min(gamma*v for v in interval))+mu*model.q_min
            term = (rc.sc.G*count.log()-rc.aq(exponent)+moment.log()).exp()
            total = rc.up(total+term)
        upper = rc.up((total.log()+rc.sc.PACKETS*rc.aq(scale).log()+rc.aq(lam)*model.threshold).exp())
        if ctx.prec != precision:
            raise ArithmeticError('precision changed during cached regional evaluation')
        return upper, checked
