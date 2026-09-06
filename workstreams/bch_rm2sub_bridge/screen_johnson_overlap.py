"""Use all Johnson positivity constraints to screen sharper row moments."""
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import bridge as base
from johnson_overlap import hahn


def problem():
    source = base.HERE/'generated/weight80_fourth.json'
    receipt = base.read(source)
    caps = {80-r['difference_weight']//2:base.decode(r['probability_upper']) for r in receipt['rows']}
    atoms = sorted(caps)
    scales = np.ones(len(atoms))
    scales[-1] = float(caps[80])
    matrix, normalizers = [], []
    for degree in range(2, 81):
        values = [hahn(256, 80, degree, 80-a) for a in atoms]
        norm = max(abs(v) for v in values[:-1])
        normalizers.append(norm)
        matrix.append([float(v/norm)*s for v, s in zip(values, scales)])
    bounds = [(0, float(caps[a])) for a in atoms[:-1]]+[(0, 1)]
    equations = np.array([scales, np.array(atoms)*scales])
    return atoms, caps, scales, np.array(matrix), normalizers, bounds, equations


def run():
    output = base.HERE/'generated/johnson_overlap_screen.json'
    assert not output.exists()
    atoms, caps, scales, matrix, normalizers, bounds, equations = problem()
    objective = np.array([(a-25)**2 for a in atoms])*scales
    rows = []
    for degree in [1, 3, 5, 10, 20, 40, 80]:
        fit = linprog(-objective, A_ub=-matrix[:degree-1] if degree>1 else None,
                      b_ub=np.zeros(degree-1) if degree>1 else None,
                      A_eq=equations, b_eq=[1, 25], bounds=bounds, method='highs')
        row = dict(maximum_degree=degree, success=bool(fit.success), message=fit.message)
        if fit.success:
            row.update(variance_upper_screen=float(-fit.fun),
                       equality_dual=(-fit.eqlin.marginals).tolist(),
                       positivity_multipliers=(-fit.ineqlin.marginals).tolist() if degree>1 else [],
                       primal_mass=(fit.x*scales).tolist())
        rows.append(row)
        print('Johnson overlap LP', degree, row.get('variance_upper_screen'), row['success'], flush=True)
    base.write_new(output, dict(status='JOHNSON_ROW_OVERLAP_LP_SCREEN_ONLY', rows=rows,
        atoms=atoms, normalizers=[base.encode(x) for x in normalizers],
        mathematical_source='https://arxiv.org/html/2405.07666v2#S4.SS2',
        outward_certified=False, full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), base.HERE/'johnson_overlap.py', base.HERE/'generated/weight80_fourth.json']}))


if __name__ == '__main__':
    run()
