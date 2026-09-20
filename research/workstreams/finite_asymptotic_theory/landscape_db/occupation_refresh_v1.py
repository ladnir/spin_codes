"""Four-state occupation bounds, including multi-bit zero syndromes.

Separate producer version: no receipt-bound historical evaluator is changed.
All arithmetic is nearest binary64; witnesses are not outward certificates.
"""
import math

import numpy as np
from scipy.special import logsumexp
from scipy.optimize import minimize

import activation_q1_refresh as q1
import activation_occupation as old
import balanced_occupation as balanced
import composition_boxes as boxes
import composition_occupation as composition
import typed_dense_boxes as typed


class Epochs:
    """Cache the hypergeometric laws shared by every tilt of one inner map."""
    def __init__(self, t, s, spectrum, kernel):
        # Reuse the historical input validator, including the kernel mass.
        old.epoch_logs(t, s, spectrum, kernel, 1., 0)
        self.t, self.s = t, s
        self.spectrum, self.kernel = spectrum, kernel
        self.weights = sorted(spectrum)
        self.mass = np.array([math.log(spectrum[w]) - math.log((1 << s)-1)
                              for w in self.weights])
        self.hyper = []
        for j in range(t+1):
            law = np.full((len(self.weights), t+1), -np.inf)
            total = math.log(math.comb(t, j))
            for i, w in enumerate(self.weights):
                for v in range(max(0, j-t+w), min(w, j)+1):
                    law[i, w+j-2*v] = math.log(math.comb(w, v)*math.comb(t-w, j-v))-total
            self.hyper.append(law)

    def at(self, lam, maximum=None):
        if not math.isfinite(lam) or lam <= 0:
            raise ValueError('positive finite tilt required')
        maximum = self.t if maximum is None else min(maximum, self.t)
        m = (1 << self.s)-1
        log_m, log_kappa = math.log(m), math.log1p(1/(m-1))
        survive = math.log1p(-1/m)
        out = np.full((maximum+1, 4, 4), -np.inf)
        zero, _ = q1.epoch_logs(self.t, self.s, self.spectrum, [lam])
        out[0] = zero[0]
        for j in range(1, maximum+1):
            moments = logsumexp(self.hyper[j]-lam*np.arange(self.t+1), axis=1)
            arbitrary = float(max(moments))
            uniform = float(logsumexp(moments+self.mass))
            total = math.comb(self.t, j)
            k = self.kernel[j]
            beta = math.log(k)-math.log(total) if k else -math.inf
            nonkernel = math.log(total-k)-math.log(total) if k < total else -math.inf
            # All zero-syndrome inputs are retained, including nonzero ones.
            out[j, 0, 0] = beta-lam*j
            out[j, 0, 1] = nonkernel-lam*j
            for state, moment in ((1, arbitrary), (2, uniform),
                                  (3, min(arbitrary, log_kappa+uniform))):
                # For any live prestate, fresh multiplication cancels a
                # nonzero syndrome with probability 1/M. Kernel inputs
                # never cancel and refresh exactly to U. To avoid needing
                # a joint kernel/output enumerator, dominate the combined
                # surviving mass in L by (1-1/M)*moment + kernel_upper/M.
                kernel_upper = min(moment, beta-lam*max(0, min(self.weights)-j))
                nonkernel_upper = min(moment, nonkernel-lam*max(0, min(self.weights)-j))
                out[j, state, 0] = nonkernel_upper-log_m
                out[j, state, 3] = np.logaddexp(moment+survive, kernel_upper-log_m)
        return out


def polynomial_product(left, right, maximum):
    limit = min(maximum, len(left)+len(right)-2)
    result = np.full((limit+1, 4, 4), -np.inf)
    for degree in range(limit+1):
        lo, hi = max(0, degree-len(right)+1), min(degree, len(left)-1)
        products = q1.matrix_product(left[lo:hi+1], right[degree-hi:degree-lo+1][::-1])
        result[degree] = np.logaddexp.reduce(products, axis=0)
    return result


def region_logs(epoch, t, length, maximum):
    if length < t or length % t or maximum > length or len(epoch) != min(t, maximum)+1:
        raise ValueError('incomplete or non-native region')
    power = epoch+np.array([math.log(math.comb(t, j)) for j in range(len(epoch))])[:, None, None]
    current = np.full((1, 4, 4), -np.inf)
    current[0, np.arange(4), np.arange(4)] = 0.
    exponent = length//t
    while exponent:
        if exponent & 1:
            current = polynomial_product(current, power, maximum)
        exponent >>= 1
        if exponent:
            power = polynomial_product(power, power, maximum)
    return current-np.array([math.log(math.comb(length, j)) for j in range(maximum+1)])[:, None, None]


def terminal_logs(matrices, count):
    if count < 1:
        raise ValueError('positive product length required')
    power = matrices
    current = np.full_like(power, -np.inf)
    current[:, np.arange(4), np.arange(4)] = 0.
    while count:
        if count & 1:
            current = q1.matrix_product(current, power)
        count >>= 1
        if count:
            power = q1.matrix_product(power, power)
    return np.logaddexp.reduce(current[:, 0], axis=1)


def epoch_mixture(epoch, t, thetas):
    theta = np.asarray(thetas)
    if len(epoch) != t+1 or np.any((theta <= 0) | (theta >= 1)):
        raise ValueError('complete epochs and interior probabilities required')
    matrix = np.full((len(theta), 4, 4), -np.inf)
    for j in range(t+1):
        mass = math.log(math.comb(t, j))+j*np.log(theta)+(t-j)*np.log1p(-theta)
        np.logaddexp(matrix, mass[:, None, None]+epoch[j], out=matrix)
    return matrix


class SparseModel(composition.SparseComposition):
    def components(self, regions, cutoff, lam, shift):
        if shift not in self.auxiliary_cache:
            roots, active, inactive = old.density_roots(self.counts, self.block, self.bands, shift)
            self.auxiliary_cache[shift] = (
                composition.probability_logs(self.indices, active, inactive),
                self.block*np.sum(roots[self.indices], axis=1))
        probabilities, cost = self.auxiliary_cache[shift]
        matrix = np.full((len(probabilities), 4, 4), -np.inf)
        for j in range(probabilities.shape[1]):
            np.logaddexp(matrix, probabilities[:, j, None, None]+regions[j], out=matrix)
        return np.minimum(self.trivial, terminal_logs(matrix, self.block)+cost+cutoff*lam)


class CompositionBoxes(boxes.CompositionBoxes):
    def bound(self, regions, lower, upper, total, length, cutoff, lam):
        matrix = self.matrix(regions, lower, total)
        return (math.log(math.comb(length, total))
                +boxes.log_assignment_count_upper(lower, upper, total)
                +cutoff*lam+float(terminal_logs(matrix[None], self.block)[0]))


class DenseModel(typed.TypedDense):
    """Adaptive complete type cover, with four-state moments and fixed witnesses."""
    def __init__(self, counts, block, t, s, a_counts, kernel, length, tilts,
                 probability_scales=(.5, .75, 1.), bands=None):
        # Keep the inherited partition/search, replacing its transfer and
        # evaluator. Initialize the small counting bank without old epochs.
        super().__init__(counts, block, t, s, a_counts, kernel, length, [],
                         probability_scales, False, bands)
        for shift in (1., 2., 3.):
            p = self.ps.copy(); costs = self.log_gammas.copy()
            for g, band in enumerate(self.bands, 1):
                if band == [block]:
                    continue
                eta = math.log(p[g])-math.log1p(-p[g])+shift
                lp, ln = -np.logaddexp(0., -eta), -np.logaddexp(0., eta)
                p[g] = math.exp(lp)
                costs[g] = max(math.log(counts[w])-math.log(math.comb(block, w))-w*lp-(block-w)*ln for w in band)
            self.probability_banks.append((10.+shift, p, costs))
        self.evaluated_boxes = 0
        self.tilts = list(tilts)
        prepared = Epochs(t, s, a_counts, kernel)
        self.epochs = [prepared.at(math.exp(z)) for z in tilts]
        # A table guides cold proposal searches. Final reported witnesses
        # always use the direct positive transfer, not interpolation.
        self.eta_grid = np.arange(-16., 16.001, .05)
        theta = 1/(1+np.exp(-self.eta_grid))
        size = self.block*self.length
        self.proposal_tables = [terminal_logs(epoch_mixture(epoch, t, theta), size//t)/size
                                +np.logaddexp(0., self.eta_grid) for epoch in self.epochs]

    def optimize_proposal(self, corners, probabilities, costs, index, initial):
        table = self.proposal_tables[index]
        anchor = int(np.argmax(initial))
        free = np.array([j for j in range(len(initial)) if j != anchor])
        mult = typed.gammaln(self.length+1)-typed.gammaln(corners+1).sum(axis=1)
        def objective(coordinates):
            logits = np.zeros(len(initial)); logits[free] = coordinates
            lp = logits-logsumexp(logits); proposal = np.exp(lp)
            theta = float(proposal@probabilities)
            eta = math.log(theta)-math.log1p(-theta)
            cell = int(np.clip(np.searchsorted(self.eta_grid, eta)-1, 0, len(table)-2))
            density = (table[cell+1]-table[cell])/.05
            value = table[cell]+(eta-self.eta_grid[cell])*density-np.logaddexp(0., eta)
            terms = corners@(costs-self.block*lp)-(self.block-1)*mult
            total = float(logsumexp(terms))
            mean = np.exp(terms-total)@corners
            derivative = (density-theta)/(theta*(1-theta))
            gradient = proposal*(1+derivative*(probabilities-theta))-mean/self.length
            return value+total/(self.block*self.length), gradient[free]
        solution = minimize(objective, np.log(initial/initial[anchor])[free], jac=True,
                            method='L-BFGS-B', bounds=[(-28., 28.)]*len(free),
                            options={'maxiter': 60, 'ftol': 1e-12, 'gtol': 1e-8})
        logits = np.zeros(len(initial)); logits[free] = solution.x
        return np.exp(logits-logsumexp(logits))

    def evaluate(self, lower, upper):
        corners = typed.vertices(lower, upper, self.length)
        if not len(corners):
            return None
        initial = typed.proposal_for(corners, self.length)
        candidates = []
        for shift in (-1., 0., 1., 2.):
            proposal = initial.copy()
            proposal[1:] *= math.exp(shift)
            proposal /= proposal.sum()
            for scale, probabilities, costs in self.probability_banks:
                candidates.append((proposal, scale, probabilities, costs))
        thetas = np.array([proposal@p for proposal, scale, p, costs in candidates])
        eta = np.log(thetas)-np.log1p(-thetas)
        penalty = typed.lattice_log_count(lower, upper)
        best, witness = math.inf, None
        ranked = []
        for tilt_index, tilt in enumerate(self.tilts):
            moments = (np.interp(eta, self.eta_grid, self.proposal_tables[tilt_index])
                       -np.logaddexp(0., eta))*(self.block*self.length)
            for (proposal, scale, p, costs), moment in zip(candidates, moments):
                values = typed.point_logs(corners, self.length, self.block, costs,
                                          proposal, float(moment), self.cutoff, math.exp(tilt))
                bound = float(max(values))+penalty
                ranked.append((bound, tilt_index, proposal, scale, p, costs))
        # Interpolation only ranks candidate proposals. Re-evaluate both
        # the initial and optimized versions of three distinct candidates.
        seen = set()
        for _, index, initial, scale, p, costs in sorted(ranked, key=lambda x: x[0]):
            key = (index, scale)
            if key in seen:
                continue
            seen.add(key)
            optimized = self.optimize_proposal(corners, p, costs, index, initial)
            for proposal in (initial, optimized):
                moment = float(terminal_logs(epoch_mixture(self.epochs[index], self.t, [proposal@p]),
                                              self.block*self.length//self.t)[0])
                values = typed.point_logs(corners, self.length, self.block, costs, proposal,
                                          moment, self.cutoff, math.exp(self.tilts[index]))
                bound = float(max(values))+penalty
                if bound < best:
                    best = bound
                    witness = dict(log_surprisal=self.tilts[index], proposal=proposal.tolist(),
                                   probabilities=p.tolist(), log_density_costs=costs.tolist(),
                                   probability_scale=scale, proposal_optimized=proposal is optimized,
                                   maximum_vertex=corners[int(np.argmax(values))].tolist())
            if len(seen) == 3:
                break
        self.evaluated_boxes += 1
        if self.evaluated_boxes % 256 == 0:
            print(f'dense boxes evaluated: {self.evaluated_boxes}', flush=True)
        return dict(lower=lower.tolist(), upper=upper.tolist(), corners=corners,
                    own_log_bound=best, best_log_bound=best, witness=witness,
                    parent=None, children=[])
