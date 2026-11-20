"""LP proposals and exact rational duals for functions of T80 overlap."""
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
import bridge as base
from johnson_overlap import hahn


class MomentProblem:
    def __init__(self, maximum_degree=24):
        receipt = base.read(base.HERE/'generated/weight80_fourth.json')
        for name, digest in receipt['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
        self.caps = {80-r['difference_weight']//2:base.decode(r['probability_upper']) for r in receipt['rows']}
        self.atoms = sorted(self.caps)
        self.norms, self.exact_matrix = [], []
        for degree in range(2, maximum_degree+1):
            row = [hahn(256, 80, degree, 80-a) for a in self.atoms]
            norm = max(abs(v)*self.caps[a] for a, v in zip(self.atoms, row))
            self.norms.append(norm)
            self.exact_matrix.append([v/norm for v in row])
        self.matrix = np.array([[float(v*self.caps[a]) for a, v in zip(self.atoms, row)] for row in self.exact_matrix])
        self.scales = np.array([float(self.caps[a]) for a in self.atoms])
        self.equations = np.array([self.scales, np.array(self.atoms)*self.scales])

    def exact_bound(self, values, degree, intercept, slope, multipliers):
        assert len(values) == len(self.atoms) and len(multipliers) == degree-1
        assert min(multipliers, default=F(0)) >= 0
        upper = intercept+25*slope
        for i, a in enumerate(self.atoms):
            affine = intercept+slope*a-sum((mult*self.exact_matrix[j][i] for j, mult in enumerate(multipliers)), F(0))
            upper += self.caps[a]*max(F(0), values[i]-affine)
        return upper

    def solve(self, values, degree):
        assert degree <= len(self.norms)+1
        objective = np.array([float(v*self.caps[a]) for a, v in zip(self.atoms, values)])
        scale = max(abs(objective))
        assert scale > 0
        fit = linprog(-objective/scale, A_ub=-self.matrix[:degree-1], b_ub=np.zeros(degree-1),
                      A_eq=self.equations, b_eq=[1, 25], bounds=[(0, 1)]*len(self.atoms), method='highs')
        assert fit.success
        intercept, slope = [F.from_float(float(-v*scale)) for v in fit.eqlin.marginals]
        multipliers = [max(F(0), F.from_float(float(-v*scale))) for v in fit.ineqlin.marginals]
        upper = self.exact_bound(values, degree, intercept, slope, multipliers)
        return dict(degree=degree, upper=base.encode(upper), intercept=base.encode(intercept),
                    slope=base.encode(slope), multipliers=[base.encode(v) for v in multipliers],
                    diagnostic_primal_objective=float(-fit.fun*scale))

    def replay(self, values, row):
        upper = self.exact_bound(values, row['degree'], base.decode(row['intercept']),
                                 base.decode(row['slope']), list(map(base.decode, row['multipliers'])))
        assert upper == base.decode(row['upper'])
        return upper
