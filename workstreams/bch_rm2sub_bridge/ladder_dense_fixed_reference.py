"""Fixed iid probability per box, paying for distinct label-coordinate sums.

The input tilt may depend on theta. This requires the explicit composition
count below; it cannot reuse the fixed-x label sum without that extra cost.
"""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
from flint import arb
from scipy.optimize import brentq
from scipy.special import logsumexp

import ladder_dense_tilt as tilted

labels,core,base,aa = tilted.labels,tilted.core,tilted.base,tilted.aa


def scalar_secant(rows,lo,hi,nlo,nhi,r,eta,alpha,beta):
    """Secant of log S followed by concave-log tangents; exact corner cover."""
    theta0 = (F(lo,rows)*nlo+F(hi,rows)*nhi)/2
    core.require(0 < theta0 < 1 and 0 <= beta <= 256,'Invalid tangent or secant slope')
    t0 = aa(theta0)
    log_r,log_not_r = aa(r).log(),(1-aa(r)).log()
    values = []
    for q in (lo,hi):
        for nu in (nlo,nhi):
            theta = aa(F(q,rows)*nu)
            log_theta = t0.log()+(theta-t0)/t0
            log_not_theta = (1-t0).log()-(theta-t0)/(1-t0)
            value = q*alpha+aa(eta*q*nu)-q*aa(beta)*(log_r-log_not_r)-256*rows*log_not_r
            value += q*aa(beta)*log_theta+aa(256*rows-q*beta)*log_not_theta
            values.append(value.upper())
    return max(values)


class Checker(tilted.Checker):
    def __init__(self,spec,ps,bands=None):
        super().__init__(spec,ps,bands)
        self.shell_logs = [[aa(F(self.caps[w],math.comb(256,w))/(p**w*(1-p)**(256-w))).log()
                            for w in band] for band,p in zip(self.bands,self.ps)]
        self.fixed_epoch = lru_cache(maxsize=256)(self._fixed_epoch)
        self.fixed_moment = lru_cache(maxsize=1024)(self._fixed_moment)

    def _fixed_epoch(self,index):
        core.require(type(index) is int and -400 <= index <= 160,'Invalid output tilt')
        lam = (arb(index)/40).exp()
        return lam,labels.transfer.epoch(self.t,self.s,self.ac,self.kernel,(-lam).exp())

    def _fixed_moment(self,index,r):
        lam,coefficients = self.fixed_epoch(index)
        weights = labels.intervals.bernstein_weights(self.t,aa(r))
        matrix = tuple(sum((w*row[k] for w,row in zip(weights,coefficients)),arb(0)).upper() for k in range(16))
        result = labels.transfer.power(matrix,1 << self.power)
        return (self.cutoff*lam+sum(result[:4],arb(0)).log()).upper()

    def log_sum(self,log_x,eta):
        terms = []
        for band,p,logs in zip(self.bands,self.ps,self.shell_logs):
            terms.append((max((g+w*log_x).upper() for g,w in zip(logs,band))-aa(eta*p)).exp())
        terms.append((256*log_x-aa(eta)).exp())
        return sum(terms,arb(0)).log().upper()

    def witness(self,lo,hi,a,b):
        result = super().witness(lo,hi,a,b)
        nlo,nhi = self.density(a,b)
        nu,q = float((nlo+nhi)/2),(lo+hi)/2
        theta = q/self.rows*nu
        # Fine local discovery avoids declaring a point unresolved merely
        # because the inherited grid steps by .2/.3 in the log tilts.
        xs = list(range(max(-400,result['input_tilt']-8),min(400,result['input_tilt']+8)+1))
        logs_x = np.array(xs)/40
        costs = np.array([[max(g+np.array(band)*x) for g,band in zip(self.shells,self.bands)]+[256*x]
                          for x in logs_x])
        etas,common = [],[]
        for g in costs:
            def derivative(e):
                values = g-e*self.float_p
                return nu-float(np.exp(values-logsumexp(values))@self.float_p)
            eta = -10000. if derivative(-10000.) >= 0 else 10000. if derivative(10000.) <= 0 else brentq(derivative,-10000.,10000.)
            eta = round(4*eta)/4
            etas.append(eta)
            common.append(q*(float(logsumexp(g-eta*self.float_p))+eta*nu))
        factors = 1-theta+theta*np.exp(-logs_x)
        rs = theta*np.exp(-logs_x)/factors
        weights = np.array([[math.comb(self.t,j)*v**j*(1-v)**(self.t-j) for j in range(self.t+1)] for v in rs])
        common = np.array(common)+256*self.rows*np.log(factors)
        best = None
        for tilt in range(max(-400,result['tilt']-10),min(160,result['tilt']+10)+1,2):
            lam,_,_,epoch = self.epoch(tilt)
            values = common+tilted.batched_terminal(np.einsum('bi,ijk->bjk',weights,epoch),self.power)+self.cutoff*float(lam)
            i = int(np.argmin(values))
            if math.isfinite(values[i]) and (best is None or values[i] < best[0]):
                best = values[i],dict(tilt=tilt,input_tilt=xs[i],eta=base.encode(F(etas[i])))
        if best is not None:
            result = best[1]
        x = math.exp(result['input_tilt']/40)
        r = theta/(x+(1-x)*theta)
        # Choosing a nearby rational r is harmless; it only proposes a law.
        result['fixed_r'] = base.encode(F.from_float(min(1-1e-12,max(1e-12,r))))
        return result

    def bound(self,lo,hi,a,b,witness):
        core.require(1 <= lo <= hi <= self.rows and 0 <= a <= b <= 1,'Invalid box')
        r,eta = base.decode(witness['fixed_r']),base.decode(witness['eta'])
        core.require(0 < r < 1 and abs(eta) <= 10000,'Invalid fixed reference witness')
        nlo,nhi = self.density(a,b)
        theta_lo,theta_hi = F(lo,self.rows)*nlo,F(hi,self.rows)*nhi
        # Production boxes have positive density width; this also excludes
        # the isolated theta=1 point where a separate limiting formula is needed.
        core.require(0 < theta_lo <= theta_hi <= 1 and theta_lo < 1,'Degenerate reference box')
        log_odds_r = aa(r).log()-(1-aa(r)).log()
        u0 = aa(theta_lo).log()-(1-aa(theta_lo)).log()-log_odds_r
        y0 = self.log_sum(u0,eta)
        if theta_hi == 1:
            # All exponents are at most 256, so log S(u)-256u decreases.
            beta = F(256)
            alpha = (y0-256*u0).upper()
        else:
            u1 = aa(theta_hi).log()-(1-aa(theta_hi)).log()-log_odds_r
            y1 = self.log_sum(u1,eta)
            slope = 0. if theta_lo == theta_hi else float((y1-y0)/(u1-u0))
            beta = F(min(1024,max(0,round(4*slope))),4) if math.isfinite(slope) else F(128)
            # An arbitrary slope in [0,256] is valid when its line dominates
            # both endpoints of the convex log S function.
            alpha = max((y0-aa(beta)*u0).upper(),(y1-aa(beta)*u1).upper())
        scalar = scalar_secant(self.rows,lo,hi,nlo,nhi,r,eta,alpha,beta)
        location = F(1,2) if 2*lo <= self.rows <= 2*hi else F(hi if 2*hi < self.rows else lo,self.rows)
        count = self.rows*labels.intervals.entropy(location)+256*arb(self.rows+1).log()+arb(hi-lo+1).log()
        # At most this many distinct nu values at each Q. Paying for all
        # compositions allows x(theta) to vary inside a retained box.
        count += arb(math.comb(hi+len(self.bands),len(self.bands))).log()
        return (count+scalar+self.fixed_moment(witness['tilt'],r)).upper()
