"""Uniform Bernoulli moments from endpoint bounds and logit convexity."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
from types import SimpleNamespace

from flint import arb,ctx
import fast_combined_dense as prior

fixed = prior.high_band_dense.parent.parent.fixed_input
model = fixed.model


def interpolate(log_lower,log_upper,lower,upper,bits):
    assert 0 < lower <= upper < 1 and type(bits) is int and bits > 0
    width = model.number(upper*(1-lower)/(lower*(1-upper))).log()
    # For z=logit(p), log E[f(X)] = log(sum_x f(x)e^(z|x|))
    # - N log(1+e^z). The first term is convex; the second term's
    # secant defect is at most N*(z1-z0)^2/32, since h'' <= 1/4.
    return model.up(max(log_lower.upper(),log_upper.upper())+bits*width**2/32)


class Checker(prior.Checker):
    def __init__(self,exponent):
        super().__init__(exponent)
        self.zero_checker.fixed_input_bound = self.fixed_input_bound

    def fixed_input_bound(self,lo,hi,a,b,witness):
        old = fixed.Checker.fixed_input_bound(self,lo,hi,a,b,witness)
        if old is None:
            return None
        nlo,nhi = self.density(a,b)
        theta_lo,theta_hi = F(lo,self.rows)*nlo,F(hi,self.rows)*nhi
        theta = (theta_lo+theta_hi)/2
        r,eta = (model.base.decode(witness[k]) for k in ('fixed_r','eta'))
        x = theta*(1-r)/(r*(1-theta))
        lower = theta_lo/(x+(1-x)*theta_lo)
        upper = theta_hi/(x+(1-x)*theta_hi)
        lam = (arb(witness['tilt'])/40).exp()
        moments = []
        for p in sorted({lower,upper}):
            probability = model.number(p)
            matrix = self.engine.bernoulli(probability,lam)
            activated = fixed.bernoulli_activation.activated(matrix,self.engine,probability,lam)
            values = [model.independent.terminal(m,self.engine.n,1 << self.power)
                      for m in (matrix,activated)]
            moments.append(min(model.up(value.log()) for value in values))
        moment = interpolate(moments[0],moments[-1],lower,upper,self.engine.output_bits)
        costs = fixed.scalar.tilted_costs(self.bands,self.ps,self.caps,model.number(x).log())
        log_s = sum(((g-model.number(eta*p)).exp()
                     for g,p in zip(costs,self.ps+(F(1),))),arb(0)).log().upper()
        scalar = fixed.scalar.scalar_bound(self.rows,lo,hi,nlo,nhi,log_s,eta,model.number(x))
        locations = sum(math.comb(self.rows,q) for q in range(lo,hi+1))
        comparison = self.comparison(lo,hi,a,b,ctx.prec)
        improved = model.up(arb(locations).log()+256*comparison.log()+scalar+self.cutoff*lam+moment)
        return min(old,improved)


def run(*args):
    driver = prior.high_band_dense.parent.parent.budget_dense
    original = driver.short_dense
    driver.short_dense = SimpleNamespace(Checker=Checker,
        mixed_dense=prior.high_band_dense.parent.parent.short_dense.mixed_dense)
    try:
        return driver.run(*args)
    finally:
        driver.short_dense = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--nodes',type=int,default=300)
    p.add_argument('--seconds',type=float,default=300)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.output.resolve(),16,64,a.seed.resolve(),a.nodes,a.seconds,42,a.verify)
