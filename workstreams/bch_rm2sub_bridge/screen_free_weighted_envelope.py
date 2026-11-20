"""Joint convex screen over envelope values, total-weight multiplier and slope."""
import argparse
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.sparse import coo_matrix, diags, hstack, vstack, csr_matrix
import bridge as base
import convex_partition_screen as partition


class Oracle:
    def __init__(self, table, sums):
        self.table, self.delta = table, sums-1631.25
        n = self.n = len(table)
        curvature = diags([-np.ones(n-2), 2*np.ones(n-2), -np.ones(n-2)],
                          [0, 1, 2], shape=(n-2, n), format='csr')
        curvature = hstack([curvature, csr_matrix((n-2, 2))])
        slope_rows = csr_matrix(([1., -1., -1., 1., -1.],
                                ([0, 0, 1, 1, 1], [0, 1, n-2, n-1, n+1])), shape=(2, n+2))
        self.fixed = vstack([curvature, slope_rows])
        self.fixed_rhs = np.r_[np.zeros(n-2), .7, 0.]
        self.cuts = set(zip(range(n), np.argmax(table, axis=1).tolist()))

    def solve(self, gradient):
        n = self.n
        objective = gradient.copy()
        objective[:n] += 1e-14
        objective[-1] += 1e-14
        for _ in range(40):
            kk, ss = np.array(sorted(self.cuts)).T
            rows = np.repeat(np.arange(len(kk)), 2)
            columns = np.column_stack([kk, np.full(len(kk), n)]).ravel()
            values = np.column_stack([-np.ones(len(kk)), -self.delta[ss]]).ravel()
            matrix = coo_matrix((values, (rows, columns)), shape=(len(kk), n+2)).tocsr()
            fit = linprog(objective, A_ub=vstack([self.fixed, matrix]),
                          b_ub=np.r_[self.fixed_rhs, -self.table[kk, ss]],
                          bounds=[(None, None)]*n+[(-4, 0), (0, 8)], method='highs')
            assert fit.success
            x = fit.x
            residual = self.table-x[n]*self.delta[None, :]-x[:n, None]
            maxima, indices = np.max(residual, axis=1), np.argmax(residual, axis=1)
            if max(maxima) < 1e-6:
                assert max(self.fixed@x-self.fixed_rhs) < 1e-6
                return x
            added = {(int(k), int(indices[k])) for k in np.flatnonzero(maxima > 1e-7)}-self.cuts
            assert added, ('Numerical separation residual', float(max(maxima)))
            self.cuts.update(added)
        raise AssertionError('Separation iteration limit')


def run(tag, totals, iterations):
    output = base.HERE/'generated'/f'free_weighted_{tag}_screen.json'
    assert not output.exists()
    source = base.HERE/'generated/marginal_sum_overlap_screen.json'
    saved = base.read(source)
    for name, digest in saved['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    table = np.array([[v if v is not None else -np.inf for v in row] for row in saved['log_bounds']])
    sums = np.array(saved['marginal_sums'])
    n = len(table)
    part = partition.Partition(n)
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
        oracle = Oracle(table, sums)
        mean = a/256
        weights = np.exp(-mean+part.k*math.log(mean)-part.factorial)
        x = oracle.solve(np.r_[weights, 0., 0.])
        def evaluate(x):
            return part.fit(x[:n], x[-1], a)
        constant = np.r_[np.full(n, np.max(table)), 0., 0.]
        if evaluate(constant)['value'] < evaluate(x)['value']:
            x = constant
        def outer(theta):
            return logsumexp(support_logmass+r*logsumexp(logmass+atoms*theta))-a*theta
        mass = minimize_scalar(outer, bounds=(-6, 6), method='bounded')
        history = []
        for iteration in range(iterations+1):
            result = evaluate(x)
            gradient = np.r_[result['gradient'], 0., result['slope_gradient']]
            vertex = oracle.solve(gradient)
            gap = 256*float(gradient@(x-vertex))
            if iteration % 20 == 0 or iteration == iterations or gap < 1e-5:
                history.append(dict(iteration=iteration, gap=gap, multiplier=float(x[n]), extension_slope=float(x[-1]),
                                    combined_summand_bound_bits_screen=float((result['value']+min(0, mass.fun))/math.log(2))))
                print('Free weighted A', a, history[-1], flush=True)
            if gap < 1e-5 or iteration == iterations:
                break
            direction = vertex-x
            step = minimize_scalar(lambda eta:evaluate(x+eta*direction)['value'], bounds=(0, 1),
                                   method='bounded', options={'xatol':1e-6, 'maxiter':40})
            x += step.x*direction
        residual = float(np.max(table-x[n]*(sums-1631.25)[None, :]-x[:n, None]))
        assert residual < 1e-5 and max(oracle.fixed@x-oracle.fixed_rhs) < 1e-5
        rows.append(dict(total_core_overlap=a, multiplier=float(x[n]), extension_slope=float(x[-1]),
                         proposed_log_envelope=x[:n].tolist(), inner_log_bound=result['value'], inner_log_tilt=result['theta'],
                         actual_mass_log_bound=float(min(0, mass.fun)), mass_log_tilt=float(mass.x),
                         combined_summand_bound_bits_screen=float((result['value']+min(0, mass.fun))/math.log(2)),
                         conditional_gradient_gap=gap, maximum_constraint_residual=residual, history=history))
    base.write_new(output, dict(status='JOINT_FREE_SLOPE_FIXED_TOTAL_ENVELOPE_SCREEN_ONLY', rows=rows,
        infinite_affine_extension=True, multiplier_range=[-4, 0], slope_range=[0, 8],
        outward_certified=False, full_second_moment_certified=False, actual_second_moment_estimated=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__), Path(partition.__file__), source, law_path]}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag', required=True)
    parser.add_argument('--totals', nargs='+', type=int, required=True)
    parser.add_argument('--iterations', type=int, default=200)
    args = parser.parse_args()
    run(args.tag, args.totals, args.iterations)
