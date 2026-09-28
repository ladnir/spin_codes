"""Scalar domination retaining a band of the total outer input weight.

For one shuffled row of length B, let d[w] = Abar[w]/C(B,w), with w in
the permitted interval I. A particular array x of n rows has capped mass
product_i d[wt(x_i)]. Define the max-product convolution

    D_n(W) = max_{w_1+...+w_n=W, w_i in I} product_i d[w_i].

For W in [L,U], its ratio to the iid Bernoulli(p) array measure is at most

    max_{L <= W <= U} D_n(W)/(p**W * (1-p)**(Bn-W)).

There is no binomial or row-label factor in D_n: it bounds the mass of
each individual array, not the sum of all arrays of a given weight.

The exact convolution below is deliberately limited to small checks.
Production bounds use an affine majorant supplied by a rational witness s:

    D_n(W) <= gamma_I(s)**n * s**W * (1-s)**(Bn-W).

The ratio of this majorant to Bernoulli(p) is monotone in W. Its maximum
on the band is attained at U when s >= p, and at L otherwise. Minimizing
over any finite set of witnesses remains an upper bound. The choice s=p
recovers the usual gamma_I(p)**n bound. Different disjoint bands can use
different reference probabilities p and independent witness choices s.

Multiply each band factor by its reference expectation of the same
nonnegative encoder function, then sum the bands. Each reference moment
may include inputs outside its assigned band; that only enlarges its term.
Group-location and active-row-label multiplicities belong outside this
scalar bound, once each. This module supplies no inner distance bound.
"""
from fractions import Fraction as Q
from math import comb, log

from flint import arb

from row_counts import row_gamma_exact, row_gamma_function


def _rational(value):
    value = Q(value)
    return arb(value.numerator)/value.denominator


def _up(value):
    if not value.is_finite():
        raise ArithmeticError('finite outward aggregate factor required')
    return arb(value.upper())


class AggregateWeights:
    """One row-weight interval repeated for a positive number of rows.

    Binary64 methods are proposals only. majorant_arb replays exact
    rational witnesses with outward arithmetic at the caller's precision.
    """
    def __init__(self, caps, lo, hi, rows, *, exclude_zero=True):
        self.caps = tuple(Q(x) for x in caps)
        self.length = len(self.caps)-1
        if not self.caps or any(x < 0 for x in self.caps):
            raise ValueError('nonempty nonnegative spectrum caps required')
        if (not isinstance(lo,int) or not isinstance(hi,int)
                or not 0 <= lo <= hi <= self.length):
            raise ValueError('invalid row-weight interval')
        if not isinstance(rows,int) or rows < 1:
            raise ValueError('a positive integer row count is required')
        self.lo, self.hi, self.rows = lo, hi, rows
        self.exclude_zero = exclude_zero
        self.total_bits = self.length*rows
        self.weights = tuple(w for w in range(max(lo,1 if exclude_zero else 0),hi+1)
                             if self.caps[w])
        self.density = {w:self.caps[w]/comb(self.length,w) for w in self.weights}
        self.minimum = rows*min(self.weights) if self.weights else 0
        self.maximum = rows*max(self.weights) if self.weights else -1
        self._gammas = {}
        self._log_gamma = row_gamma_function(self.caps,lo,hi,exclude_zero=exclude_zero)

    @staticmethod
    def _probability(value):
        value = Q(value)
        if not 0 < value < 1:
            raise ValueError('an interior reference probability is required')
        return value

    def _band(self, lower, upper):
        if (not isinstance(lower,int) or not isinstance(upper,int)
                or not 0 <= lower <= upper <= self.total_bits):
            raise ValueError('invalid total input-weight band')
        lower, upper = max(lower,self.minimum),min(upper,self.maximum)
        return (lower,upper) if lower <= upper else None

    def gamma(self, witness):
        witness = self._probability(witness)
        if witness not in self._gammas:
            self._gammas[witness] = row_gamma_exact(self.caps,self.lo,self.hi,witness,
                                                   exclude_zero=self.exclude_zero)
        return self._gammas[witness]

    def _parameters(self, p, lower, upper, witness):
        p, witness = self._probability(p),self._probability(witness)
        band = self._band(lower,upper)
        if band is None:
            return p,witness,None
        return p,witness,band[int(witness >= p)]

    def majorant_exact(self, p, lower, upper, witness):
        """Exact rational factor; intended for bounded checks, not large n."""
        p,s,endpoint = self._parameters(p,lower,upper,witness)
        if endpoint is None:
            return Q(0)
        return self.gamma(s)**self.rows*(s/p)**endpoint*((1-s)/(1-p))**(self.total_bits-endpoint)

    def majorant_arb(self, p, lower, upper, witness):
        """Outward factor for exact rational p and affine-majorant witness."""
        p,s,endpoint = self._parameters(p,lower,upper,witness)
        if endpoint is None:
            return arb(0)
        return _up(_rational(self.gamma(s))**self.rows
                   *_rational(s/p)**endpoint
                   *_rational((1-s)/(1-p))**(self.total_bits-endpoint))

    def best_arb(self, p, lower, upper, witnesses):
        """Minimum of independently valid outward factors, with its witness."""
        witnesses = tuple(self._probability(s) for s in witnesses)
        if not witnesses:
            raise ValueError('at least one affine-majorant witness is required')
        choices = [(self.majorant_arb(p,lower,upper,s),s) for s in witnesses]
        return min(choices,key=lambda item:item[0])

    def majorant_log(self, p, lower, upper, witness):
        """Natural-log binary64 proposal; never a certificate factor."""
        p,s = float(p),float(witness)
        if not 0 < p < 1 or not 0 < s < 1:
            raise ValueError('interior proposal probabilities required')
        band = self._band(lower,upper)
        if band is None:
            return float('-inf')
        endpoint = band[int(s >= p)]
        return (self.rows*self._log_gamma(s)+endpoint*log(s/p)
                +(self.total_bits-endpoint)*log((1-s)/(1-p)))

    def exact_products(self, *, max_operations=200_000):
        """Small exact max-product convolution, with an explicit work cap.

        Reaching this cap fails before the next update; it does not return
        a truncated product table. Production users should use majorant_arb.
        """
        if not isinstance(max_operations,int) or max_operations < 1:
            raise ValueError('a positive operation cap is required')
        current, operations = {0:Q(1)},0
        for _ in range(self.rows):
            operations += len(current)*len(self.weights)
            if operations > max_operations:
                raise ValueError('exact convolution work cap exceeded; use the affine majorant')
            updated = {}
            for total,value in current.items():
                for w,density in self.density.items():
                    weight, candidate = total+w,value*density
                    if candidate > updated.get(weight,Q(0)):
                        updated[weight] = candidate
            current = updated
        return current

    def exact_band(self, p, lower, upper, *, max_operations=200_000):
        """Optimal scalar for the capped measure, by bounded exact DP."""
        p = self._probability(p)
        band = self._band(lower,upper)
        if band is None:
            return Q(0)
        products = self.exact_products(max_operations=max_operations)
        return max((value/(p**w*(1-p)**(self.total_bits-w))
                    for w,value in products.items() if band[0] <= w <= band[1]),default=Q(0))
