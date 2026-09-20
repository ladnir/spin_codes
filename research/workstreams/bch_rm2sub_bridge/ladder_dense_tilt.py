"""Two-tilt dense witnesses; fixed label coordinates, exact all-one factor.

The input-count tilt changes the reference law, not the encoder. Discovery
uses binary64; bound() independently reconstructs all constants with Arb.
"""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
from flint import arb
from scipy.optimize import brentq
from scipy.special import logsumexp

import ladder_dense_labels as labels

core, base, aa = labels.core, labels.base, labels.aa


def tilted_costs(bands, ps, caps, log_x):
    """Gamma'_g=max_w A_w/binom(n,w)/q^w/(1-q)^(n-w) * x^w."""
    ordinary = []
    for band, p in zip(bands, ps):
        values = [aa(F(caps[w], math.comb(256, w))/(p**w*(1-p)**(256-w))).log()
                  +w*log_x for w in band]
        ordinary.append(max(v.upper() for v in values))
    # The all-one density is one before tilting; its normalizer is x^256.
    return ordinary+[256*log_x]


def batched_terminal(matrices, power):
    """Positive fixed-width powers, discovery only."""
    matrices = matrices.copy()
    scales = matrices.max(axis=(1, 2))
    matrices /= scales[:, None, None]
    logs = np.log(scales)
    for _ in range(power):
        matrices = matrices@matrices
        scales = matrices.max(axis=(1, 2))
        matrices /= scales[:, None, None]
        logs = 2*logs+np.log(scales)
    totals = matrices[:, 0].sum(axis=1)
    # An underflowed zero is not evidence of a tiny moment. Reject that
    # discovery candidate; every chosen witness is still recomputed in Arb.
    result = np.full(len(totals), np.inf)
    valid = totals > 0
    result[valid] = np.log(totals[valid])+logs[valid]
    return result


def scalar_bound(rows, lo, hi, nlo, nhi, log_s, eta, x):
    """Concave-log tangent keeps input and label normalizers coupled."""
    theta0 = F(lo+hi,2*rows)*(nlo+nhi)/2
    a = 1/x-1
    b0 = 1+a*aa(theta0)
    gradient = a/b0
    intercept = 256*rows*(b0.log()-gradient*aa(theta0))
    slope = aa(eta)+256*gradient
    return max((intercept+q*(log_s+slope*aa(nu))).upper()
               for q in (lo,hi) for nu in (nlo,nhi))


class Checker(labels.Checker):
    def __init__(self, spec, ps, bands=None):
        super().__init__(spec, ps, bands)
        self.caps = core.inputs.caps_module.caps()
        self.shells = [np.array([math.log(self.caps[w])-math.log(math.comb(256,w))
                      -w*math.log(float(p))-(256-w)*math.log1p(-float(p))
                      for w in band]) for band,p in zip(self.bands,self.ps)]
        self.input_costs = lru_cache(maxsize=2048)(self._input_costs)
        self.normalizer = lru_cache(maxsize=4096)(self._normalizer)

    def _input_costs(self, tilt):
        core.require(type(tilt) is int and -400 <= tilt <= 400, 'Invalid input tilt')
        return tilted_costs(self.bands, self.ps, self.caps, arb(tilt)/40)

    def _normalizer(self, tilt, eta):
        costs = self.input_costs(tilt)
        return sum(((g-aa(eta*p)).exp() for g,p in zip(costs,self.ps+(F(1),))), arb(0)).log()

    def witness(self, lo, hi, vlo, vhi):
        nlo, nhi = self.density(vlo, vhi)
        nu = float((nlo+nhi)/2)
        q = (lo+hi)/2
        theta = q/self.rows*nu
        input_tilts = np.arange(-120, 121, 8)
        logs_x = input_tilts/40
        costs = np.array([[max(a+np.array(b)*x) for a,b in zip(self.shells,self.bands)]+[256*x]
                          for x in logs_x])
        etas = []
        label_cost = []
        for g in costs:
            def derivative(e):
                values = g-e*self.float_p
                return nu-float(np.exp(values-logsumexp(values))@self.float_p)
            low, high = -10000., 10000.
            eta = low if derivative(low) >= 0 else high if derivative(high) <= 0 else brentq(derivative,low,high)
            eta = round(eta*4)/4
            etas.append(eta)
            label_cost.append(q*(float(logsumexp(g-eta*self.float_p))+eta*nu))
        factors = 1-theta+theta*np.exp(-logs_x)
        r = theta*np.exp(-logs_x)/factors
        weights = np.array([[math.comb(self.t,j)*v**j*(1-v)**(self.t-j)
                             for j in range(self.t+1)] for v in r])
        common = np.array(label_cost)+256*self.rows*np.log(factors)
        prediction = round(40*math.log(q/self.rows))
        candidates = sorted(set(max(-400,min(120,v)) for v in range(prediction-80,prediction+161,12)))
        best = None
        for tilt in candidates:
            lam, _, _, epoch = self.epoch(tilt)
            values = common+batched_terminal(np.einsum('bi,ijk->bjk',weights,epoch),self.power)+self.cutoff*float(lam)
            i = int(np.argmin(values))
            if best is None or values[i] < best[0]:
                best = (values[i], dict(tilt=tilt,input_tilt=int(input_tilts[i]),eta=base.encode(F(etas[i]))))
        return best[1]

    def bound(self, lo, hi, vlo, vhi, witness):
        core.require(1 <= lo <= hi <= self.rows and 0 <= vlo <= vhi <= 1, 'Invalid box')
        eta = base.decode(witness['eta'])
        core.require(abs(eta) <= 10000, 'Invalid auxiliary slope')
        input_tilt = witness['input_tilt']
        # Validate before evaluating exp or the normalizer.
        self.input_costs(input_tilt)
        x = (arb(input_tilt)/40).exp()
        nlo, nhi = self.density(vlo,vhi)
        theta_lo, theta_hi = F(lo,self.rows)*nlo, F(hi,self.rows)*nhi
        rlo = (aa(theta_lo)/(x+(1-x)*aa(theta_lo))).lower()
        rhi = (aa(theta_hi)/(x+(1-x)*aa(theta_hi))).upper()
        # Exact endpoint theta=1 avoids division-induced intervals above one.
        rlo, rhi = max(arb(0),rlo), min(arb(1),rhi)
        lam, ascending, descending, _ = self.epoch(witness['tilt'])
        low = labels.intervals.bernstein_weights(self.t,rlo)
        high = labels.intervals.bernstein_weights(self.t,rhi)
        matrix = tuple(max(arb(0),(sum((w*c for w,c in zip(high,a)),arb(0))-
                       sum((w*c for w,c in zip(low,d)),arb(0))).upper()) for a,d in zip(ascending,descending))
        result = labels.transfer.power(matrix,1 << self.power)
        log_s = self.normalizer(input_tilt,eta)
        scalar = scalar_bound(self.rows,lo,hi,nlo,nhi,log_s,eta,x)
        location = F(1,2) if 2*lo <= self.rows <= 2*hi else F(hi if 2*hi < self.rows else lo,self.rows)
        count = self.rows*labels.intervals.entropy(location)+256*arb(self.rows+1).log()+arb(hi-lo+1).log()
        return (count+scalar+self.cutoff*lam+sum(result[:4],arb(0)).log()).upper()
