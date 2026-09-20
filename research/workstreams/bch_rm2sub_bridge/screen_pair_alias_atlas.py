"""A three-tilt screen for a wider shared-region pair-type atlas.

The Gaussian-looking overlap envelope is only tested pointwise here. This
screen makes no outward claim and no distributional or second-moment claim.
"""
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from flint import arb, arb_poly, ctx
import bridge as base
import tightened_occupancy as tight
from screen_pair_type_bound import load_pairs, probabilities
from certify_pair_type_fourier_tracked import evaluate


def run():
    output = base.HERE/'generated/pair_alias_wide_atlas_screen.json'
    assert not output.exists()
    powers, counts = load_pairs()
    grid = (256, 256, 128)
    total = math.prod(grid)
    roots = [np.exp(2j*np.pi*np.arange(n)/n) for n in grid]
    ctx.prec = 192
    kernel = tight.kernel_spectrum()
    poly = arb_poly([kernel.get(j, 0) for j in range(129)])**64
    singles = {j:float((poly[j]/math.comb(8192, j)).log()) for j in range(740, 901, 2)}
    tiles = [(list(range(740, 819, 2)), list(range(740, 819, 2)), 780, 780),
             (list(range(740, 819, 2)), list(range(820, 901, 2)), 780, 860),
             (list(range(820, 901, 2)), list(range(820, 901, 2)), 860, 860)]
    rows = []
    for first, second, cx, cy in tiles:
        center_overlap = round(cx*cy/8192)
        center = np.array([8192-cx-cy+center_overlap, cy-center_overlap, cx-center_overlap, center_overlap])
        def objective(theta):
            p = probabilities(theta)
            return 64*math.log(float(evaluate(p[:, None], powers, counts)[0]))-float(center@np.log(p))
        opt = minimize(objective, np.log(center[1:]/center[0]), method='BFGS')
        p = probabilities(opt.x)
        table = np.empty(total, dtype=np.complex128)
        for start in range(0, total, 2048):
            indices = np.arange(start, min(start+2048, total))
            i = indices//(grid[1]*grid[2])
            j = indices//grid[2] % grid[1]
            k = indices % grid[2]
            inputs = np.vstack([np.full(len(indices), p[0], dtype=np.complex128),
                                p[1]*roots[0][i], p[2]*roots[1][j], p[3]*roots[2][k]])
            value = evaluate(inputs, powers, counts)*2**30
            for _ in range(6):
                value = value*value
            table[start:start+len(indices)] = value
        aliases = np.fft.fftn(table.reshape(grid))/total
        del table
        x, y = np.meshgrid(first, second, indexing='ij')
        single = np.array([[singles[a]+singles[b] for b in second] for a in first])
        maximum, envelope, worst, ew = -math.inf, -math.inf, None, None
        maximag, minimum_atom = 0., math.inf
        for overlap in range(40, 127):
            m = [8192-x-y+overlap, y-overlap, x-overlap, np.full(x.shape, overlap)]
            atom = aliases[m[1] % grid[0], m[2] % grid[1], m[3] % grid[2]]
            minimum_atom = min(minimum_atom, float(atom.real.min()))
            if not np.all(atom.real > 0):
                bad = np.unravel_index(int(np.argmin(atom.real)), x.shape)
                base.write_new(output, dict(status='SCREEN_INDETERMINATE_NONPOSITIVE_ATOM', grid=list(grid),
                    completed_tiles=rows, failed_center_type=center.tolist(),
                    failed_type=[int(v[bad]) for v in m], minimum_real_atom=minimum_atom,
                    local_sha256={str(Path(__file__).relative_to(base.HERE)):base.sha(Path(__file__))}))
                print('Wide atlas screen indeterminate: nonpositive numerical atom at',
                      [int(v[bad]) for v in m], flush=True)
                return
            logden = gammaln(8193)-sum(gammaln(v+1) for v in m)+sum(v*math.log(prob) for v, prob in zip(m, p))
            logratio = np.log(atom.real)-1920*math.log(2)-logden-single
            adjusted = logratio-(overlap-x*y/8192)**2/20000
            for field, label in [(logratio, 'raw'), (adjusted, 'envelope')]:
                index = np.unravel_index(int(np.argmax(field)), x.shape)
                value = float(field[index])
                if label == 'raw' and value > maximum:
                    maximum, worst = value, [int(v[index]) for v in m]
                if label == 'envelope' and value > envelope:
                    envelope, ew = value, [int(v[index]) for v in m]
            maximag = max(maximag, float(np.max(np.abs(atom.imag)/atom.real)))
        row = dict(first_weights=first, second_weights=second, center_type=center.tolist(),
                   tilt_probabilities=p.tolist(), maximum_ratio=math.exp(maximum), worst_type=worst,
                   overlap_quadratic_coefficient=1/20000, envelope_intercept=envelope, envelope_worst_type=ew,
                   max_relative_imaginary=maximag, minimum_atom=minimum_atom)
        rows.append(row)
        print('Wide tile', cx, cy, 'raw ratio', row['maximum_ratio'], 'overlap envelope intercept', envelope,
              'at', ew, 'max relative imaginary', maximag, flush=True)
    base.write_new(output, dict(status='FLOATING_WIDE_ALIAS_ATLAS_SCREEN_ONLY', grid=list(grid), rows=rows,
        covered_with_message_swap_symmetry=True, first_weights=list(range(740, 901, 2)),
        second_weights=list(range(740, 901, 2)), overlaps=list(range(40, 127)), type_count=81*81*87,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),
            base.HERE/'screen_pair_type_bound.py', base.HERE/'certify_pair_type_fourier_tracked.py',
            base.HERE/'generated/t128_s15_pair_spectrum.json', Path(tight.__file__)]}))


if __name__ == '__main__':
    run()
