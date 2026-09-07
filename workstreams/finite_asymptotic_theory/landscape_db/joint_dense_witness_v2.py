"""Optimize tilt, type proposal, and band counting probabilities jointly."""
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp

import close_bch_dense_v2 as cover
import fast_dense_witness_v1 as fast


class JointRefiner(fast.FastRefiner):
    def objective(self, corners, penalty, anchor, probabilities, temperature=.001):
        free = np.array([j for j in range(len(probabilities)) if j != anchor])
        variable = np.array([j for j in range(1, len(probabilities)) if probabilities[j] < 1.])
        weights = [np.array(band) for band in self.bands]
        shell_logs = [np.array([math.log(self.counts[w])-math.log(math.comb(self.block, w)) for w in band])
                      for band in self.bands]
        log_mult = gammaln(self.length+1)-gammaln(corners+1).sum(axis=1)

        def decode(coordinates):
            logits = np.zeros(len(probabilities)); logits[free] = coordinates[1:1+len(free)]
            proposal = np.exp(logits-logsumexp(logits))
            p = np.array(probabilities)
            p[variable] = 1/(1+np.exp(-coordinates[1+len(free):]))
            return proposal, p

        def evaluate(coordinates):
            proposal, p = decode(coordinates)
            costs = np.zeros(len(p)); derivatives = np.zeros(len(p))
            for j in range(1, len(p)):
                if p[j] == 1.:
                    costs[j] = math.log(self.counts[self.block])
                else:
                    candidates = shell_logs[j-1]-weights[j-1]*math.log(p[j])-(self.block-weights[j-1])*math.log1p(-p[j])
                    normalizer = logsumexp(candidates/temperature)
                    costs[j] = temperature*normalizer
                    derivatives[j] = self.block*p[j]-float(np.exp(candidates/temperature-normalizer)@weights[j-1])
            theta = float(proposal@p)
            y = math.log(theta)-math.log1p(-theta)
            moment, dx, dy = self.lookup(float(coordinates[0]), y)
            terms = corners@(costs-self.block*np.log(proposal))-(self.block-1)*log_mult
            total = float(logsumexp(terms)); mean = np.exp(terms-total)@corners
            lam = math.exp(coordinates[0]); theta_derivative = dy/(theta*(1-theta))
            g_proposal = proposal*(1+theta_derivative*(p-theta))-mean/self.length
            g_probability = (theta_derivative*proposal*p*(1-p)+mean*derivatives/self.size)[variable]
            value = moment+(total+penalty+self.cutoff*lam)/self.size
            return value, np.r_[dx+self.cutoff*lam/self.size, g_proposal[free], g_probability]
        return evaluate, decode, free, variable

    def refine(self, box):
        corners = fast.base.transfer.typed.vertices(box['lower'], box['upper'], self.length)
        penalty = fast.base.transfer.typed.lattice_log_count(box['lower'], box['upper'])
        witness = box['witness']; p = np.array(witness['probabilities']); pi = np.array(witness['proposal'])
        baseline = float(max(self.value(corners, witness['log_surprisal'], pi, p, self.costs(p))))+penalty
        if abs(baseline-box['own_log_bound']) > 2e-6:
            raise ArithmeticError('joint baseline failed direct replay')
        anchor = int(np.argmax(pi))
        objective, decode, free, variable = self.objective(corners, penalty, anchor, p)
        initial = np.r_[np.clip(witness['log_surprisal'], -9., 1.), np.log(pi/pi[anchor])[free],
                        np.clip(np.log(p[variable]/(1-p[variable])), -12., 12.)]
        for temperature in (.1, .01, .001):
            objective, decode, free, variable = self.objective(corners, penalty, anchor, p, temperature)
            result = minimize(objective, initial, jac=True, method='L-BFGS-B',
                              bounds=[(-9., 1.)]+[(-28., 28.)]*len(free)+[(-12., 12.)]*len(variable),
                              options={'maxiter': 160, 'ftol': 1e-14, 'gtol': 1e-10, 'maxls': 40})
            initial = result.x
        proposal, probabilities = decode(result.x); z = float(result.x[0]); costs = self.costs(probabilities)
        values = self.value(corners, z, proposal, probabilities, costs)
        value = float(max(values))+penalty
        selected = dict(box)
        if value < baseline:
            selected = dict(box, own_log_bound=value,
                            witness=dict(log_surprisal=z, proposal=proposal.tolist(), probabilities=probabilities.tolist(),
                                         log_density_costs=costs.tolist(), maximum_vertex=corners[int(np.argmax(values))].tolist()))
        # The interpolated objective only proposes a fixed counting measure.
        # Near closure, polish its tilt and type proposal with direct moments.
        if abs(selected['own_log_bound']) < 2000*math.log(2):
            polished = fast.base.Refiner.refine(self, selected)
            if polished['own_log_bound'] < selected['own_log_bound']:
                selected = dict(box, own_log_bound=polished['own_log_bound'], witness=polished['witness'])
        return selected


class JointCoverRefiner(cover.CoverRefiner, JointRefiner):
    """Preserve warm witnesses and disjoint subdivision while replacing the search."""
