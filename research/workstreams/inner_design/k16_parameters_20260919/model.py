"""Explicit-step adapter for the existing BCH-256 outward bound engine.

Historical producers stay unchanged. The exact same fixed-j and independent-map
Fourier formulas now take t from the instance. Only the outer stays fixed at 256.
"""
from fractions import Fraction as F
from functools import lru_cache
import math
from pathlib import Path
import sys
from types import FunctionType, SimpleNamespace

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'finite_migration'))
import adaptive_length_study as study
import fixed_input
import dense_bch
import density_probe
from flint import arb, arb_poly

base = study.core.mixing.model
up, number = base.up, base.number
single = base.independent.single


@lru_cache(maxsize=16)
def maps(t, s):
    record = study.record_for(256, 128) if (t, s) == (128, 19) else study.core.inner(t, s)
    a, b = record['expansion_columns'], record['feedback_columns']
    exact = base.independent.g.tv.fixed.maps
    sa = exact.spectrum(exact.generators(a, s))
    sb = exact.spectrum(exact.generators(b, s))
    kernel = exact.dual_spectrum(sb, t, s)
    assert sa[0] == sb[0] == 1 and sa.get(t, 0) == 0
    assert sum(sa.values()) == sum(sb.values()) == 1 << s
    assert sum(kernel.values()) == 1 << (t-s)
    assert len(a) == len(b) == len(set(b)) == t and all(a) and all(b)
    k = [kernel.get(j, 0) for j in range(t+1)]
    assert k[:3] == [1, 0, 0]
    spectrum = {w: n for w, n in sb.items() if w}
    caps = base.independent.g.fiber_caps(t, s, spectrum, k)
    low = base.independent.search.low_cancellation(a, b, s)
    return record, spectrum, k, caps, low


class Engine(base.Engine):
    def __init__(self, t, s, exponent=16):
        record, self.b_spectrum, self.kernel, self.caps, self.low = maps(t, s)
        self.record = record
        self.t, self.s = t, s
        self.a_columns, self.columns = record['expansion_columns'], record['feedback_columns']
        self.spectrum = record['spectrum']
        self.levels = sorted(self.spectrum)
        self.n, self.m = len(self.levels)+2, (1 << s)-1
        self.exponent, self.length = exponent, 1 << (exponent-7)
        self.outer_length, self.output_bits = 256, 256*self.length
        self.cutoff = self.output_bits//10
        assert self.length % t == 0

    def identity(self):
        result = super().identity()
        result['inner'] = dict(t=self.t, s=self.s, transvection_rounds=1,
            expansion_columns=self.a_columns, feedback_columns=self.columns)
        return result

    @lru_cache(maxsize=24)
    def epoch(self, log_lam):
        t, n, m = self.t, self.n, self.m
        z = (-number(log_lam).exp()).exp()
        powers = [z**j for j in range(2*t+1)]
        def hist(items):
            return sum((count*powers[w] for w, count in items), arb(0))
        output = []
        for j in range(t+1):
            total = math.comb(t, j)
            nonzero = total-self.kernel[j]
            nk, r = arb(nonzero)/total, arb(int(self.caps[j]['cap']))/total
            moments = [up(sum((math.comb(w,h)*math.comb(t-w,j-h)*powers[w+j-2*h]
                for h in range(max(0,j-t+w),min(w,j)+1)), arb(0))/total) for w in self.levels]
            d = max(moments)
            cd = min(up(d),up(nk),up(r*powers[min(abs(w-j) for w in self.levels)]))
            if j in self.low:
                cd = min(cd,max((up(hist(p)/total) for p in self.low[j]['patterns']),default=arb(0)))
            matrix = [arb(0) for _ in range(n*n)]
            matrix[0], matrix[1] = arb(self.kernel[j])/total*powers[j], nk*powers[j]
            matrix[n], matrix[n+1] = cd/2+min(up(d),up(nk))/(2*m), d/2
            for i,v in enumerate(self.levels):
                matrix[n+i+2] = d*self.spectrum[v]/(2*m)
                mass = min(nonzero,self.spectrum[v]*int(self.caps[j]['cap']))
                cancel = min(up(moments[i]),up(arb(mass)*powers[abs(v-j)]/(self.spectrum[v]*total)))
                if j in self.low:
                    cancel = min(cancel,up(hist(self.low[j]['by_weight'].get(v,{}).items())/(self.spectrum[v]*total)))
                matrix[(i+2)*n] = cancel/2+min(up(moments[i]),up(nk))/(2*m)
                for k,w in enumerate(self.levels):
                    matrix[(i+2)*n+k+2] = moments[i]*self.spectrum[w]/(2*m)
                if nonzero:
                    matrix[(i+2)*n+1] = moments[i]/2
                else:
                    matrix[(i+2)*n+i+2] += moments[i]/2
            output.append(tuple(map(up,matrix)))
        return output

    @lru_cache(maxsize=8)
    def region(self, log_lam, maximum):
        assert 0 <= maximum <= self.length
        rows, n, t = self.epoch(log_lam), self.n, self.t
        a = tuple(arb_poly([row[k]*math.comb(t,j)
            for j,row in enumerate(rows[:min(t,maximum)+1])]) for k in range(n*n))
        current = tuple(arb_poly([int(i==j)]) for i in range(n) for j in range(n))
        remaining = self.length//t
        while remaining:
            if remaining & 1:
                current = self.poly_mul(current,a,maximum)
            remaining >>= 1
            if remaining:
                a = self.poly_mul(a,a,maximum)
        return [tuple(max(arb(0),up(p[j]/math.comb(self.length,j))) for p in current)
                for j in range(maximum+1)]

    def q1(self, log_lam, linear=False):
        zero, one = self.epoch(log_lam)[:2]
        h = self.length//self.t
        if linear:
            rz, ra = single.identity(self.n), (arb(0),)*(self.n*self.n)
            for _ in range(h):
                ra = single.add(single.mul(ra,zero,self.n),single.mul(rz,one,self.n))
                rz = single.mul(rz,zero,self.n)
            region = rz, tuple(v/h for v in ra)
        else:
            region = single.regions(zero,one,self.n,h)
        return single.moments(*region,self.n,256)

    def bernoulli(self, theta, lam):
        assert theta > 0 and theta < 1 and lam > 0
        t,n,m = self.t,self.n,self.m
        z = (-lam).exp()
        g0,g1 = 1-theta+theta*z, theta+(1-theta)*z
        rho0 = min(arb(1),up(abs(1-2*theta*z/g0)))
        rho1 = min(arb(1),up(abs(1-2*(1-theta)*z/g1)))
        entries = []
        for v in self.levels:
            cap = arb(1)
            for w,count in self.b_spectrum.items():
                ends = (max(0,v+w-t),min(v,w))
                cap += count*max(up(rho0**(w-h)*rho1**h) for h in ends)
            entries.append(up(g0**(t-v)*g1**v*(cap/(2*(m+1))+arb(1)/(2*m))))
        result = [arb(0) for _ in range(n*n)]
        result[0] = sum((self.kernel[j]*(theta*z)**j*(1-theta)**(t-j) for j in range(t+1)),arb(0))
        result[1] = sum(((math.comb(t,j)-self.kernel[j])*(theta*z)**j*(1-theta)**(t-j) for j in range(t+1)),arb(0))
        for i,entry in [(1,max(entries))]+[(i+2,e) for i,e in enumerate(entries)]:
            result[i*n] = entry
            for k,w in enumerate(self.levels):
                result[i*n+k+2] = entry*self.spectrum[w]
        return tuple(map(up,result))


def activated(matrix, engine, theta, lam):
    z = (-lam).exp()
    g0 = 1-theta+theta*z
    rho = min(arb(1),up(abs(1-2*theta*z/g0)))
    cap = up(g0**engine.t*(1+sum((count*rho**w for w,count in engine.b_spectrum.items()),arb(0)))/(engine.m+1))
    result = list(matrix)
    result[1] = arb(0)
    for k,w in enumerate(engine.levels):
        result[k+2] = up(cap*engine.spectrum[w])
    return tuple(result)


class Checker(fixed_input.Checker):
    def __init__(self, t, s, exponent=16):
        # Retain its outer-only counting banks and empty caches, replacing
        # the default engine before evaluating any transfer or moment.
        dense_bch.Checker.__init__(self, exponent)
        self.engine = Engine(t,s,exponent)
        self.t,self.s = t,s
        self.power = exponent+1-(t.bit_length()-1)
        self.comparison = lru_cache(maxsize=2048)(self._comparison)
        assert (1 << self.power) == self.engine.output_bits//t

    def _fixed_moment(self, index, r):
        lam,theta = (arb(index)/40).exp(),number(r)
        coefficients = self.engine.epoch(F(index,40))
        weights = [math.comb(self.t,j)*theta**j*(1-theta)**(self.t-j) for j in range(self.t+1)]
        matrix = tuple(up(sum((w*row[k] for w,row in zip(weights,coefficients)),arb(0)))
            for k in range(self.engine.n*self.engine.n))
        direct = self.engine.bernoulli(theta,lam)
        matrices = (matrix,direct,activated(direct,self.engine,theta,lam))
        value = min(up(base.independent.terminal(a,self.engine.n,1 << self.power)) for a in matrices)
        return up(self.cutoff*lam+value.log())

    # Reuse the scalar reduction verbatim, changing only its call to the
    # activation transfer. Its private globals do not mutate any old module.
    fixed_input_bound = FunctionType(fixed_input.Checker.fixed_input_bound.__code__,
        {**fixed_input.Checker.fixed_input_bound.__globals__,
         'bernoulli_activation': SimpleNamespace(activated=activated)})

    def witness(self, *args):
        # The direct search already infers t from the kernel spectrum.
        return study.core.grid.ladder.candidate.dense_search.Checker.witness(self,*args)
