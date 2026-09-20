"""Screen a marginal-coordinate Fourier atlas with an overlap envelope.

Variables mark x=wt(U), y=wt(V), k=overlap, instead of three nonzero cell
counts. This reduces alias losses at mixed marginal/overlap corners.
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


def tiles():
    ranges = [list(range(740, 819, 2)), list(range(820, 901, 2))]
    for i, j in [(0, 0), (0, 1), (1, 1)]:
        for lower, upper, ck in [(40, 79, 60), (80, 119, 100), (120, 160, 140)]:
            yield ranges[i], ranges[j], list(range(lower, upper+1)), 780+80*i, 780+80*j, ck


def run():
    path = base.HERE/'generated/marginal_alias_atlas_screen.json'
    assert not path.exists()
    powers, counts = load_pairs()
    grid = (128, 128, 128)
    total = math.prod(grid)
    roots = np.exp(2j*np.pi*np.arange(128)/128)
    ctx.prec = 192
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    singles = {j:float((kernel[j]/math.comb(8192, j)).log()) for j in range(740, 901, 2)}
    rows = []
    for xs, ys, ks, cx, cy, ck in tiles():
        center = np.array([8192-cx-cy+ck, cy-ck, cx-ck, ck])
        def objective(theta):
            p = probabilities(theta)
            return 64*math.log(float(evaluate(p[:, None], powers, counts)[0].real))-float(center@np.log(p))
        opt = minimize(objective, np.log(center[1:]/center[0]), method='BFGS')
        p = probabilities(opt.x)
        values = np.empty(total, dtype=np.complex128)
        for start in range(0, total, 2048):
            indices = np.arange(start, min(start+2048, total))
            a = indices//16384
            b = indices//128 % 128
            c = indices % 128
            inputs = np.vstack([np.full(len(indices), p[0], dtype=np.complex128),
                                p[1]*roots[b], p[2]*roots[a], p[3]*roots[(a+b+c) % 128]])
            z = evaluate(inputs, powers, counts)*2**30
            for _ in range(6):
                z = z*z
            values[start:start+len(indices)] = z
        aliases = np.fft.fftn(values.reshape(grid))/total
        del values
        x, y = np.meshgrid(xs, ys, indexing='ij')
        single = np.array([[singles[a]+singles[b] for b in ys] for a in xs])
        maximum, envelope, worst, ew = -math.inf, -math.inf, None, None
        maximag = 0.
        for k in ks:
            m = [8192-x-y+k, y-k, x-k, np.full(x.shape, k)]
            atom = aliases[x % 128, y % 128, k % 128]
            assert np.all(atom.real > 0)
            denominator = gammaln(8193)-sum(gammaln(v+1) for v in m)+sum(v*math.log(prob) for v, prob in zip(m, p))
            logratio = np.log(atom.real)-1920*math.log(2)-denominator-single
            adjusted = logratio-0.00012*(k-x*y/8192)**2
            loc = np.unravel_index(int(np.argmax(logratio)), x.shape)
            if logratio[loc] > maximum:
                maximum, worst = float(logratio[loc]), [int(v[loc]) for v in m]
            loc = np.unravel_index(int(np.argmax(adjusted)), x.shape)
            if adjusted[loc] > envelope:
                envelope, ew = float(adjusted[loc]), [int(v[loc]) for v in m]
            maximag = max(maximag, float(np.max(np.abs(atom.imag)/atom.real)))
        row = dict(first_weights=xs, second_weights=ys, overlaps=ks, center_type=center.tolist(),
                   tilt_probabilities=p.tolist(), maximum_ratio=math.exp(maximum), maximum_type=worst,
                   overlap_quadratic_coefficient=0.00012, envelope_intercept=envelope,
                   envelope_maximum_type=ew, maximum_relative_imaginary=maximag)
        rows.append(row)
        print('Marginal atlas', cx, cy, ck, 'max ratio', math.exp(maximum), 'envelope intercept', envelope,
              'imaginary relative', maximag, flush=True)
    base.write_new(path, dict(status='MARGINAL_COORDINATE_ALIAS_ATLAS_SCREEN_ONLY', grid=list(grid),
        rows=rows, first_weights=list(range(740, 901, 2)), second_weights=list(range(740, 901, 2)),
        overlaps=list(range(40, 161)), type_count=81*81*121, pair_swap_symmetry=True,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),
            base.HERE/'screen_pair_type_bound.py', base.HERE/'certify_pair_type_fourier_tracked.py',
            base.HERE/'generated/t128_s15_pair_spectrum.json', Path(tight.__file__)]}))


if __name__ == '__main__':
    run()
