"""Two-dimensional moment tables guide cold dense witness optimization.

The table is shared by both BCH block sizes at a fixed output length and
inner map. It proposes witnesses only; all retained bounds use the exact
positive transfer, with log-domain evaluation where required.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp

import refine_bch_dense_v1 as base
import close_bch_dense_v2 as cover

HERE = Path(__file__).resolve().parent


class FastRefiner(base.Refiner):
    def __init__(self, counts, block, t, s, ac, kernel, length, bands):
        super().__init__(counts, block, t, s, ac, kernel, length, bands)
        self.xgrid = np.arange(-9., 1.001, .04)
        self.ygrid = np.arange(-14., 8.001, .04)
        self.size = block*length
        identity = hashlib.sha256(json.dumps([self.size, t, s, ac, kernel,
                   hashlib.sha256(Path(__file__).read_bytes()).hexdigest()], sort_keys=True).encode()).hexdigest()
        directory = HERE/'dense_moment_tables_v1'; directory.mkdir(exist_ok=True)
        path = directory/f'{identity}.json'
        if path.exists():
            data = json.loads(path.read_text())
            if data['identity'] != identity:
                raise ValueError('moment table identity changed')
            self.table = np.array(data['values'])
        else:
            theta = 1/(1+np.exp(-self.ygrid))
            values = []
            for x in self.xgrid:
                epoch = self.epochs.at(math.exp(x))
                values.append(base.transfer.terminal_logs(base.transfer.epoch_mixture(epoch, t, theta), self.size//t)/self.size)
            self.table = np.array(values)
            path.write_text(json.dumps(dict(identity=identity, values=self.table.tolist()))+'\n')

    def lookup(self, x, y):
        i = int(np.clip(np.searchsorted(self.xgrid, x)-1, 0, len(self.xgrid)-2))
        j = int(np.clip(np.searchsorted(self.ygrid, y)-1, 0, len(self.ygrid)-2))
        u, v = (x-self.xgrid[i])/.04, (y-self.ygrid[j])/.04
        a, b = self.table[i, j], self.table[i+1, j]
        c, d = self.table[i, j+1], self.table[i+1, j+1]
        value = (1-u)*(1-v)*a+u*(1-v)*b+(1-u)*v*c+u*v*d
        dx = ((1-v)*(b-a)+v*(d-c))/.04
        dy = ((1-u)*(c-a)+u*(d-b))/.04
        return value, dx, dy

    def refine(self, box):
        corners = base.transfer.typed.vertices(box['lower'], box['upper'], self.length)
        penalty = base.transfer.typed.lattice_log_count(box['lower'], box['upper'])
        witness = box['witness']; probabilities = np.array(witness['probabilities'])
        costs = self.costs(probabilities); initial = np.array(witness['proposal'])
        tilt = float(witness['log_surprisal'])
        baseline = float(max(self.value(corners, tilt, initial, probabilities, costs)))+penalty
        if abs(baseline-box['own_log_bound']) > 2e-6:
            raise ArithmeticError('fast-refinement baseline failed replay')
        anchor = int(np.argmax(initial)); free = np.array([j for j in range(len(initial)) if j != anchor])
        center = float(np.clip(tilt, -9., 1.))
        start = np.r_[center, np.log(initial/initial[anchor])[free]]
        log_mult = gammaln(self.length+1)-gammaln(corners+1).sum(axis=1)
        def objective(coordinates):
            logits = np.zeros(len(initial)); logits[free] = coordinates[1:]
            lp = logits-logsumexp(logits); proposal = np.exp(lp)
            theta = float(proposal@probabilities); y = math.log(theta)-math.log1p(-theta)
            value, dx, dy = self.lookup(float(coordinates[0]), y)
            terms = corners@(costs-self.block*lp)-(self.block-1)*log_mult
            aggregate = float(logsumexp(terms)); mean = np.exp(terms-aggregate)@corners
            gradient = proposal*(1+dy*(probabilities-theta)/(theta*(1-theta)))-mean/self.length
            lam = math.exp(coordinates[0])
            return (value+(aggregate+penalty+self.cutoff*lam)/self.size,
                    np.r_[dx+self.cutoff*lam/self.size, gradient[free]])
        result = minimize(objective, start, jac=True, method='L-BFGS-B',
                          bounds=[(max(-9., center-1.), min(1., center+1.))]+[(-28., 28.)]*len(free),
                          options={'maxiter': 100, 'ftol': 1e-13, 'gtol': 1e-9, 'maxls': 30})
        logits = np.zeros(len(initial)); logits[free] = result.x[1:]
        proposal = np.exp(logits-logsumexp(logits)); z = float(result.x[0])
        values = self.value(corners, z, proposal, probabilities, costs)
        bound = float(max(values))+penalty
        selected = dict(box)
        if bound < baseline:
            selected = dict(box, own_log_bound=bound,
                            witness=dict(log_surprisal=z, proposal=proposal.tolist(), probabilities=probabilities.tolist(),
                                         log_density_costs=costs.tolist(), maximum_vertex=corners[int(np.argmax(values))].tolist()))
        predicted = objective(result.x)[0]*self.size
        # Exact polishing is useful near closure, or when interpolation
        # crosses a sharp state-transition regime. Coarse geometric boxes
        # with a large deficit need subdivision, not expensive polishing.
        if abs(selected['own_log_bound']) < 2000*math.log(2) or bound-predicted > 200*math.log(2):
            polished = super().refine(selected)
            if polished['own_log_bound'] < selected['own_log_bound']:
                selected = dict(box, own_log_bound=polished['own_log_bound'], witness=polished['witness'])
        return selected


class FastCoverRefiner(cover.CoverRefiner, FastRefiner):
    """CoverRefiner's cold optimizer dispatches to FastRefiner via this MRO."""
