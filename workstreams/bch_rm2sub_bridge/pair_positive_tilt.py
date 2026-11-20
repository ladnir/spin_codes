"""Fast discovery-only positive pair-kernel tilts with analytic gradients."""
import numpy as np
from scipy.optimize import minimize
from screen_pair_type_bound import load_pairs, probabilities, SIGNS


class PairTilt:
    def __init__(self):
        powers, self.counts = load_pairs()
        self.exponents = 4*powers

    def value_gradient(self, theta, m):
        p = probabilities(theta)
        linear = SIGNS@p
        terms = self.counts*np.prod(linear[None, :]**self.exponents, axis=1)
        value = np.sum(terms)
        assert value > 0
        derivative = np.divide(terms[:, None]*self.exponents, linear[None, :],
                               out=np.zeros_like(self.exponents, dtype=float),
                               where=linear[None, :] != 0).sum(axis=0)
        gradient_p = SIGNS.T@derivative/value
        objective = 64*np.log(value)-m@np.log(p)
        gradient = (64*p*gradient_p-m)[1:]
        return float(objective), gradient

    def fit(self, m, initial=None):
        m = np.asarray(m, dtype=float)
        assert m.shape == (4,) and min(m) >= 0 and sum(m) == 8192
        if initial is None:
            initial = np.log(np.maximum(m[1:], .25)/m[0])
        fit = minimize(self.value_gradient, initial, args=(m,), jac=True,
                       method='L-BFGS-B', bounds=[(-32, 8)]*3,
                       options={'ftol':1e-14, 'gtol':1e-7, 'maxiter':200})
        return probabilities(fit.x), float(fit.fun), fit.x
