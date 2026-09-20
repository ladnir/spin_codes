"""Atom-cap scaled Johnson LP; every proposed dual is checked exactly."""
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import bridge as base
from johnson_overlap import hahn


def run():
    output = base.HERE/'generated/johnson_scaled_variance.json'
    assert not output.exists()
    source = base.HERE/'generated/weight80_fourth.json'
    saved = base.read(source)
    for name, digest in saved['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    caps = {80-r['difference_weight']//2:base.decode(r['probability_upper']) for r in saved['rows']}
    atoms = sorted(caps)
    exact_matrix, norms = [], []
    for degree in range(2, 81):
        values = [hahn(256, 80, degree, 80-a) for a in atoms]
        norm = max(abs(v)*caps[a] for a, v in zip(atoms, values))
        exact_matrix.append([v/norm for v in values])
        norms.append(norm)
    matrix = np.array([[float(v*caps[a]) for a, v in zip(atoms, row)] for row in exact_matrix])
    scales = np.array([float(caps[a]) for a in atoms])
    equations = np.array([scales, np.array(atoms)*scales])
    objective = np.array([(a-25)**2 for a in atoms])*scales
    rows = []
    for degree in [1, 3, 5, 10, 20, 40, 80]:
        fit = linprog(-objective, A_ub=-matrix[:degree-1] if degree>1 else None,
                      b_ub=np.zeros(degree-1) if degree>1 else None,
                      A_eq=equations, b_eq=[1, 25], bounds=[(0, 1)]*len(atoms), method='highs')
        row = dict(maximum_degree=degree, optimizer_success=bool(fit.success))
        if fit.success:
            intercept, slope = [F.from_float(float(-v)) for v in fit.eqlin.marginals]
            multipliers = [max(F(0), F.from_float(float(-v))) for v in fit.ineqlin.marginals] if degree>1 else []
            upper = intercept+25*slope
            for i, a in enumerate(atoms):
                affine = intercept+slope*a-sum((mult*exact_matrix[j][i] for j, mult in enumerate(multipliers)), F(0))
                upper += caps[a]*max(F(0), F((a-25)**2)-affine)
            row.update(variance_upper=base.encode(upper), diagnostic_lp_variance=float(-fit.fun),
                       equality_dual=[base.encode(intercept), base.encode(slope)],
                       positivity_multipliers=[base.encode(v) for v in multipliers])
        rows.append(row)
        print('Scaled Johnson degree', degree, 'LP', row.get('diagnostic_lp_variance'),
              'EXACT upper', float(base.decode(row['variance_upper'])) if fit.success else None, flush=True)
    base.write_new(output, dict(status='EXACT_DUAL_BOUNDS_FROM_CAP_SCALED_JOHNSON_LP',
        atoms=atoms, rows=rows, normalizers=[base.encode(v) for v in norms],
        exact_dual_hinge_bounds=True, full_second_moment_certified=False,
        mathematical_source='https://arxiv.org/html/2405.07666v2#S4.SS2',
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), base.HERE/'johnson_overlap.py', source]}))


if __name__ == '__main__':
    run()
