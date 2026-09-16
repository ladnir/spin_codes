"""Keep the input normalizer coupled to the Bernoulli moment over a box."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
from types import SimpleNamespace

from flint import arb,ctx
import convex_input_dense as prior

fixed,model = prior.fixed,prior.model


def coupled_majorant(values,lower,upper,bits,slope):
    assert 0 < lower <= upper < 1 and bits > 0
    width = model.number(upper*(1-lower)/(lower*(1-upper))).log()
    # log(D(theta)^N M(r(theta))) is a Bernoulli log moment in theta.
    # Its logit second derivative is >= -N/4. The additional linear
    # term slope*theta contributes >= -abs(slope)/4.
    return model.up(max(v.upper() for v in values)
                    +model.number(bits+abs(slope))*width**2/32)


class Checker(prior.prior.Checker):
    def coupled_bound(self,lo,hi,a,b,witness):
        if hi-lo > 7:
            return None
        nlo,nhi = self.density(a,b)
        if nhi == 1:
            return None
        center = (F(lo,self.rows)*nlo+F(hi,self.rows)*nhi)/2
        r,eta = (model.base.decode(witness[k]) for k in ('fixed_r','eta'))
        x = center*(1-r)/(r*(1-center))
        lam = (arb(witness['tilt'])/40).exp()
        costs = fixed.scalar.tilted_costs(self.bands,self.ps,self.caps,model.number(x).log())
        terms = [(cost-model.number(eta*p)).exp()
                 for cost,p in zip(costs,self.ps+(F(1),))]
        low = sum((v for v,p in zip(terms,self.ps+(F(1),)) if p <= F(3,4)),arb(0))
        high = sum((v for v,p in zip(terms,self.ps+(F(1),)) if p > F(3,4)),arb(0))
        xi = prior.prior.high_band_dense.XI
        total = low+high
        log_totals = [total.log().upper(),low.log().upper(),
                      (low+high*(-model.number(xi)).exp()).log().upper()]
        branch_sums = [arb(0),arb(0),arb(0)]
        for q in range(lo,hi+1):
            lower,upper = F(q,self.rows)*nlo,F(q,self.rows)*nhi
            values = []
            for theta in sorted({lower,upper}):
                probability = model.number(theta/(x+(1-x)*theta))
                matrix = self.engine.bernoulli(probability,lam)
                activated = fixed.bernoulli_activation.activated(matrix,self.engine,probability,lam)
                moment = min(model.up(model.independent.terminal(m,self.engine.n,1 << self.power).log())
                             for m in (matrix,activated))
                normalizer = model.number(1-theta+theta/x)
                values.append(self.engine.output_bits*normalizer.log()+model.number(eta*self.rows*theta)+moment)
            coupled = coupled_majorant(values,lower,upper,self.engine.output_bits,eta*self.rows)
            common = arb(math.comb(self.rows,q)).log()+self.cutoff*lam+coupled
            for index,log_total in enumerate(log_totals):
                branch_sums[index] += (common+q*log_total).exp()
        unrestricted = self.comparison(lo,hi,a,b,ctx.prec)**256
        unsplit = (branch_sums[0]*unrestricted).log()
        exceptional = branch_sums[2]*unrestricted*model.number(xi).exp()
        if nlo > max(p for p in self.ps if p <= F(3,4)):
            split = exceptional.log()
        else:
            restricted = self.constant_comparison(lo,hi,a,b,ctx.prec)**256
            split = (branch_sums[1]*restricted+exceptional).log()
        return min(model.up(unsplit),model.up(split))

    def bound(self,lo,hi,a,b,witness):
        old = super().bound(lo,hi,a,b,witness)
        if old < -80*arb(2).log():
            return old
        improved = self.coupled_bound(lo,hi,a,b,witness)
        return old if improved is None else min(old,improved)


def run(*args):
    driver = prior.prior.high_band_dense.parent.parent.budget_dense
    original = driver.short_dense
    driver.short_dense = SimpleNamespace(Checker=Checker,
        mixed_dense=prior.prior.high_band_dense.parent.parent.short_dense.mixed_dense)
    try:
        return driver.run(*args)
    finally:
        driver.short_dense = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--nodes',type=int,default=400)
    p.add_argument('--seconds',type=float,default=600)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.output.resolve(),16,64,a.seed.resolve(),a.nodes,a.seconds,42,a.verify)
