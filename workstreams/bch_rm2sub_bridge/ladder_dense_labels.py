"""Label-sum dense bound with the all-one label retained as a unit term.

Discovery and outward box evaluation share witnesses, not arithmetic.
The iid comparison costs (L+1)^256 and preserves the actual state process.
"""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
from flint import arb
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

import ladder_instance as identity
import activation_density_arb as transfer
import frontier_dense as intervals
import k30_dense_unrestricted_screen as row_search
import constant_split_grid as split

core, base = identity.core, identity.core.base
aa = intervals.aa


def log_label_sum(ps, costs, eta):
    """Exact input rationals; outward Arb normalizer with one constant term."""
    core.require(len(ps) == len(costs) and all(0 < p < 1 for p in ps), 'Ordinary probabilities required')
    e = aa(eta)
    return ((-e).exp()+sum((aa(g)*(-e*aa(p)).exp() for p, g in zip(ps, costs)), arb(0))).log()


def positive_terminal(matrix, power):
    """Discovery only, scaled before each multiply to avoid initial underflow."""
    matrix = matrix.copy()
    maximum = float(matrix.max())
    if not maximum > 0:
        return -math.inf
    matrix /= maximum
    scale = math.log(maximum)
    for _ in range(power):
        matrix = matrix@matrix
        maximum = float(matrix.max())
        if not maximum > 0:
            return -math.inf
        matrix /= maximum
        scale = 2*scale+math.log(maximum)
    total = float(matrix[0].sum())
    return math.log(total)+scale if total > 0 else -math.inf


class Checker:
    def __init__(self, spec, ps, bands=None):
        self.spec = spec
        self.rows, self.cutoff = spec['rows'], spec['cutoff']
        self.t, self.s, self.ac, self.kernel = core.inputs.load(spec['configuration'])
        self.power = spec['message_exponent']-6
        core.require(self.t == 128 and 1 << self.power == 256*self.rows//self.t, 'Wrong epoch power')
        self.ps = tuple(F(p) for p in ps)
        self.bands = tuple(tuple(b) for b in (split.ordinary_bands() if bands is None else bands))
        caps = core.inputs.caps_module.caps()
        core.require(len(self.ps) == len(self.bands) and all(0 < p < 1 for p in self.ps) and
            sorted(w for b in self.bands for w in b) == sorted(w for w in caps if w != 256), 'Invalid ordinary partition')
        self.pmin = min(self.ps)
        self.costs = tuple(split.sparse.costs_for(self.bands, self.ps, caps))
        self.float_p = np.array([float(p) for p in self.ps]+[1.])
        self.float_g = np.array([math.log(g) for g in self.costs]+[0.])
        self.epoch = lru_cache(maxsize=256)(self._epoch)
        self.label_sum = lru_cache(maxsize=4096)(lambda eta: log_label_sum(self.ps, self.costs, eta))

    def _epoch(self, tilt):
        # Integer tilt units are 1/40 in log(lambda).
        core.require(type(tilt) is int and -400 <= tilt <= 160, 'Invalid tilt index')
        lam = (arb(tilt)/40).exp()
        coefficients = transfer.epoch(self.t, self.s, self.ac, self.kernel, (-lam).exp())
        parts = [intervals.monotone_parts([intervals.rational(row[k]) for row in coefficients]) for k in range(16)]
        ascending = [[aa(v) for v in a] for a, d in parts]
        descending = [[aa(v) for v in d] for a, d in parts]
        floating = np.array([[float(v) for v in row] for row in coefficients]).reshape(-1, 4, 4)
        return lam, ascending, descending, floating

    def density(self, vlo, vhi):
        return self.pmin+(1-self.pmin)*vlo, self.pmin+(1-self.pmin)*vhi

    def eta(self, nu):
        # Only a proposal: rounded to an exact rational before certification.
        def objective(e):
            return float(logsumexp(self.float_g-e*self.float_p)+e*nu)
        fit = minimize_scalar(objective, bounds=(-8192, 8192), method='bounded', options={'xatol': 1e-5})
        core.require(fit.success, 'Auxiliary slope search failed')
        return F(round(4*float(fit.x)), 4)

    def witness(self, lo, hi, vlo, vhi):
        nlo, nhi = self.density(vlo, vhi)
        nu = float((nlo+nhi)/2)
        r = float(F(lo+hi, 2*self.rows))*nu
        weights = np.array([math.comb(self.t,j)*r**j*(1-r)**(self.t-j) for j in range(self.t+1)])
        prediction = round(40*math.log((lo+hi)/(2*self.rows)))
        candidates = sorted(set(max(-400, min(100, v)) for v in range(prediction-60, prediction+113, 12)))
        best = None
        for tilt in candidates:
            lam, _, _, co = self.epoch(tilt)
            value = positive_terminal(np.einsum('i,ijk->jk', weights, co), self.power)+self.cutoff*float(lam)
            if best is None or value < best[0]:
                best = value, tilt
        return dict(tilt=best[1], eta=base.encode(self.eta(nu)))

    def bound(self, lo, hi, vlo, vhi, witness):
        core.require(1 <= lo <= hi <= self.rows and 0 <= vlo <= vhi <= 1, 'Invalid box')
        eta = base.decode(witness['eta'])
        core.require(abs(eta) <= 10000, 'Invalid auxiliary slope')
        nlo, nhi = self.density(vlo, vhi)
        rlo, rhi = F(lo, self.rows)*nlo, F(hi, self.rows)*nhi
        lam, ascending, descending, _ = self.epoch(witness['tilt'])
        low = intervals.bernstein_weights(self.t, aa(rlo))
        high = intervals.bernstein_weights(self.t, aa(rhi))
        matrix = tuple(max(arb(0), (sum((w*c for w,c in zip(high,a)), arb(0))-
                    sum((w*c for w,c in zip(low,d)), arb(0))).upper()) for a,d in zip(ascending, descending))
        result = transfer.power(matrix, 1 << self.power)
        # Q*(log S + eta*nu), bounded jointly at the four rectangle corners.
        log_s = self.label_sum(eta)
        label = max((q*(log_s+aa(eta*nu))).upper() for q in (lo,hi) for nu in (nlo,nhi))
        location = F(1,2) if 2*lo <= self.rows <= 2*hi else F(hi if 2*hi < self.rows else lo, self.rows)
        count = self.rows*intervals.entropy(location)+256*arb(self.rows+1).log()
        # This bound sums every Q in the box, not merely one occupancy.
        return (count+label+self.cutoff*lam+sum(result[:4], arb(0)).log()+arb(hi-lo+1).log()).upper()

    def point_screen(self, q, v):
        witness = self.witness(q, q, v, v)
        return float(-self.bound(q, q, v, v, witness)/arb(2).log()), witness


def row_probabilities(slope):
    ps, _ = row_search.row_witnesses(core.inputs.caps_module.caps(), slope)
    core.require(ps[-1] == 1, 'Wrong constant-row witness')
    return ps[:-1]


def singleton_probabilities(slope):
    bands = [(w,) for w in sorted(core.inputs.caps_module.caps()) if w != 256]
    ps = []
    for (w,) in bands:
        r = w/256
        p = 2*r/(math.sqrt((1-slope)**2+4*slope*r)+1-slope)
        ps.append(F.from_float(p))
    return ps, bands
