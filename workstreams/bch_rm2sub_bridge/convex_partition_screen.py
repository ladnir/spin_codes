"""Stable scalar partition screen for an affine-extended log envelope.

Infinite extension gives a positive-coefficient upper bound. This helper
uses floating arithmetic and is not an outward evaluator.
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.special import gammaln, logsumexp
from scipy.stats import poisson


def poisson_tail_series(start, log_lambda):
    lam = math.exp(log_lambda)
    if lam < .9*start:
        term = total = 1.
        weighted = 0.
        for j in range(1, 10000):
            term *= lam/(start+j)
            total += term
            weighted += j*term
            ratio = lam/(start+j+1)
            if term*ratio/(1-ratio) <= 1e-16*total:
                break
        else:
            raise ArithmeticError('Tail series did not converge')
        return start*log_lambda-gammaln(start+1)+math.log(total), start+weighted/total
    logsf = float(poisson.logsf(start-1, lam))
    assert math.isfinite(logsf)
    mean = lam*math.exp(float(poisson.logsf(start-2, lam))-logsf)
    return lam+logsf, mean


class Partition:
    def __init__(self, length):
        self.length = length
        self.k = np.arange(length)
        self.factorial = gammaln(self.k+1)

    def evaluate(self, psi, slope, theta):
        tail, tail_mean = poisson_tail_series(self.length, theta+slope)
        tail += psi[-1]-slope*(self.length-1)
        logs = np.r_[psi-self.factorial+self.k*theta, tail]
        logtotal = float(logsumexp(logs))
        weights = np.exp(logs-logtotal)
        mean = float(weights[:-1]@self.k+weights[-1]*tail_mean)
        gradient = weights[:-1].copy()
        gradient[-1] += weights[-1]
        slope_gradient = weights[-1]*(tail_mean-(self.length-1))
        return logtotal, mean, gradient, float(slope_gradient)

    def fit(self, psi, slope, total):
        theta = brentq(lambda t:self.evaluate(psi, slope, t)[1]-total/256, -20, 12, xtol=1e-12)
        value, mean, gradient, slope_gradient = self.evaluate(psi, slope, theta)
        objective = gammaln(total+1)-total*(math.log(256)+theta)+256*value
        return dict(theta=theta, value=float(objective), gradient=gradient,
                    slope_gradient=slope_gradient, mean=mean)
