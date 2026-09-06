"""Joint optimization of the fixed-total multiplier and a convex envelope.

The LP separation oracle checks every retained (overlap,marginal-sum)
cell. Only numerical proposals are produced; no proof gate is implied.
"""
import argparse
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.sparse import coo_matrix, diags, hstack, vstack, csr_matrix
import bridge as base


class EnvelopeOracle:
    def __init__(self, table, sums, right):
        self.table, self.delta, self.right = table, sums-1631.25, right
        self.length = len(table)
        n = self.length
        curvature = diags([-np.ones(n-2), 2*np.ones(n-2), -np.ones(n-2)],
                          [0, 1, 2], shape=(n-2, n), format='csr')
        slopes = csr_matrix(([1., -1., -1., 1.], ([0, 0, 1, 1], [0, 1, n-2, n-1])), shape=(2, n))
        self.fixed = hstack([vstack([curvature, slopes]), csr_matrix((n, 1))], format='csr')
        self.fixed_rhs = np.r_[np.zeros(n-2), .7, right]
        self.cuts = set(zip(range(n), np.argmax(table, axis=1).tolist()))

    def solve(self, weights):
        n = self.length
        objective = np.r_[weights+1e-14, 0.]
        for _ in range(40):
            pairs = sorted(self.cuts)
            kk, ss = np.array(pairs).T
            rows = np.repeat(np.arange(len(pairs)), 2)
            columns = np.column_stack([kk, np.full(len(pairs), n)]).ravel()
            values = np.column_stack([-np.ones(len(pairs)), -self.delta[ss]]).ravel()
            matrix = coo_matrix((values, (rows, columns)), shape=(len(pairs), n+1)).tocsr()
            fit = linprog(objective, A_ub=vstack([self.fixed, matrix]),
                          b_ub=np.r_[self.fixed_rhs, -self.table[kk, ss]],
                          bounds=[(None, None)]*n+[(-4, 0)], method='highs')
            assert fit.success
            x = fit.x
            residual = self.table-x[-1]*self.delta[None, :]-x[:-1, None]
            maxima = np.max(residual, axis=1)
            indices = np.argmax(residual, axis=1)
            if max(maxima) < 1e-6:
                assert max(self.fixed@x-self.fixed_rhs) < 1e-6
                return x
            added = {(int(k), int(indices[k])) for k in np.flatnonzero(maxima > 1e-7)}-self.cuts
            assert added, ('Numerical separation residual', float(max(maxima)))
            self.cuts.update(added)
        raise AssertionError('Separation iteration limit')


def run(tag, totals, slopes, iterations):
    output = base.HERE/'generated'/f'joint_weighted_{tag}_screen.json'
    assert not output.exists()
    source = base.HERE/'generated/marginal_sum_overlap_screen.json'
    saved = base.read(source)
    for name, digest in saved['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    table = np.array([[v if v is not None else -np.inf for v in row] for row in saved['log_bounds']])
    sums = np.array(saved['marginal_sums'])
    n = len(table)
    k = np.arange(208801)
    factorial = gammaln(k+1)
    law_path = base.HERE/'generated/johnson_convex_law.json'
    law = base.read(law_path)
    atoms = np.array([r['overlap'] for r in law['rows']])
    logmass = np.log([float(base.decode(r['probability'])) for r in law['rows']])
    r = np.arange(2611)
    def logcomb(n, j):
        return gammaln(n+1)-gammaln(j+1)-gammaln(n-j+1)
    support_logmass = logcomb(2610, r)+logcomb(5579, 2610-r)-logcomb(8189, 2610)
    rows = []
    for a in totals:
        assert n <= a <= 208800
        def outer(theta):
            return logsumexp(support_logmass+r*logsumexp(logmass+atoms*theta))-a*theta
        mass = minimize_scalar(outer, bounds=(-6, 6), method='bounded')
        trials, best = [], None
        for right in slopes:
            oracle = EnvelopeOracle(table, sums, right)
            mean = a/256
            x = oracle.solve(np.exp(-mean+k[:n]*math.log(mean)-factorial[:n]))
            def extended(x):
                return np.r_[x[:-1], x[-2]+right*np.arange(1, a+2-n)]
            def value(x, theta):
                return gammaln(a+1)-a*(math.log(256)+theta)+256*logsumexp(extended(x)-factorial[:a+1]+k[:a+1]*theta)
            def best_tilt(x):
                return minimize_scalar(lambda theta:value(x, theta), bounds=(-12, 12), method='bounded')
            for iteration in range(iterations):
                tilt = best_tilt(x)
                logs = extended(x)-factorial[:a+1]+k[:a+1]*tilt.x
                weights = np.exp(logs-logsumexp(logs))
                gradient = weights[:n].copy()
                gradient[-1] += sum(weights[n:])
                vertex = oracle.solve(gradient)
                gap = 256*float(gradient@(x[:-1]-vertex[:-1]))
                if gap < 1e-6:
                    break
                direction = vertex-x
                step = minimize_scalar(lambda eta:value(x+eta*direction, tilt.x), bounds=(0, 1), method='bounded')
                x += step.x*direction
            tilt = best_tilt(x)
            residual = float(np.max(table-x[-1]*(sums-1631.25)[None, :]-x[:-1, None]))
            assert residual < 1e-5
            trial = dict(multiplier=float(x[-1]), extension_slope=right, inner_log_bound=float(tilt.fun),
                         inner_log_tilt=float(tilt.x), maximum_constraint_residual=residual,
                         final_iteration=iteration, conditional_gradient_gap=float(gap), cuts=len(oracle.cuts))
            trials.append(trial)
            if best is None or trial['inner_log_bound'] < best['inner_log_bound']:
                best = dict(trial, proposed_log_envelope=x[:-1].tolist())
            print('Joint weighted A', a, 'slope', right, 'b', float(x[-1]),
                  'combined bits', (tilt.fun+min(0, mass.fun))/math.log(2), flush=True)
        rows.append(dict(total_core_overlap=a, best=best, trials=trials,
                         actual_mass_log_bound=float(min(0, mass.fun)), mass_log_tilt=float(mass.x),
                         combined_summand_bound_bits_screen=float((best['inner_log_bound']+min(0, mass.fun))/math.log(2))))
    base.write_new(output, dict(status='JOINT_FIXED_TOTAL_MULTIPLIER_ENVELOPE_SCREEN_ONLY', rows=rows,
        multiplier_range=[-4, 0], outer_actual_mgf_comparison='nu_J and exact hypergeometric support intersection',
        outward_certified=False, full_second_moment_certified=False, actual_second_moment_estimated=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__), source, law_path]}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag', required=True)
    parser.add_argument('--totals', nargs='+', type=int, required=True)
    parser.add_argument('--slopes', nargs='+', type=float, default=[0, 1, 2, 3, 4, 6])
    parser.add_argument('--iterations', type=int, default=12)
    args = parser.parse_args()
    assert args.iterations > 0
    run(args.tag, args.totals, args.slopes, args.iterations)
