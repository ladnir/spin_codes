"""Parameterized proof-first screens. No floating result is a certificate."""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb, ctx
import bridge as base
import larger_state_maps as maps
import activation_bridge as q1
import polynomial_regions as poly
import tightened_occupancy as tight
import occupation_three as groups
import screen_exponential_modes as cap_loader
import screen_refined_bands as hull
import scaled_adaptive as scaled
from screen_constant_row_split import moment


def q1_screen(name):
    t, s, spectrum, _ = maps.load(name)
    best = np.full(257, np.inf)
    witness = np.zeros(257, dtype=int)
    for tenth in range(-120, 1):
        lam = math.exp(tenth/10)
        value = q1.log_coefficients(*q1.log_regions(t, s, spectrum, lam)) + base.CUTOFF*lam
        value = math.log(base.ROWS) + np.minimum(0., value)
        witness[value < best] = tenth
        best = np.minimum(best, value)
    coefficients = {w: F.from_float(math.exp(max(-700., float(best[w])))) for w in base.WEIGHTS}
    upper, _, _ = base.bch_bound(coefficients)
    margin = math.log2(upper.denominator)-math.log2(upper.numerator)
    base.write_new(base.HERE/'generated'/f'larger_{name}_q1_screen.json', dict(
        status='LARGER_STATE_Q1_SCREEN_ONLY', configuration=name, margin_bits_diagnostic=margin,
        coefficient_rows={str(w): dict(log_coefficient=float(best[w]), witness_tenth=int(witness[w])) for w in base.WEIGHTS},
        local_sha256={**maps.sources(), Path(__file__).name: base.sha(Path(__file__)),
                      Path(q1.__file__).name: base.sha(Path(q1.__file__))}))
    print(name, 'Q1 screen margin', margin, flush=True)


def range_screen(name, qs, tilt, tag):
    path = base.HERE/'generated'/f'larger_{name}_{tag}_screen.json'
    assert not path.exists()
    t, s, spectrum, kernel = maps.load(name)
    caps, cap_sources = cap_loader.latest_caps()
    bands = groups.BANDS
    weights = [np.array(b) for b in bands]
    logs = [np.array([math.log(caps[w])-math.log(math.comb(256, w)) for w in b]) for b in bands]
    ctx.prec = 256
    lam = (arb(tilt)/10).exp()
    coefficients = poly.regions(t, s, spectrum, kernel, (-lam).exp(), max(qs))
    region = np.array([[float(v.log()) if v > 0 else -math.inf for v in row] for row in coefficients]).reshape(-1, 3, 3)
    print(name, 'region coefficients ready', tilt, flush=True)
    rows = []
    for q in qs:
        selected = region[:q+1]
        ps = []
        pure = []
        offset = base.CUTOFF*float(lam)+math.log(math.comb(base.ROWS, q))+q*math.log(len(bands))
        for w, v in zip(weights, logs):
            def objective(theta):
                p = 1/(1+math.exp(-theta))
                gamma = float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
                return moment(selected, q, p)+q*gamma
            opt = minimize_scalar(objective, bounds=(-6, 14), method='bounded', options={'xatol':1e-6})
            ps.append(1/(1+math.exp(-float(opt.x))))
            pure.append(-(float(opt.fun)+offset)/math.log(2))
        ps = np.array(ps)
        gamma = np.array([float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))) for v,w,p in zip(logs,weights,ps)])
        roots = np.exp(gamma/256)
        keep = hull.active_lines(ps, roots)
        mantissas, exponents = scaled.initial(coefficients[:q+1], outward=False)
        value = scaled.terminal_log(mantissas, exponents, ps[keep], roots[keep])+offset
        row = dict(occupation=q, witness_tenth=tilt, p=[base.encode(F.from_float(float(p))) for p in ps],
                   margin_bits_diagnostic=-value/math.log(2), pure_margins_with_assignment_cost=pure)
        rows.append(row)
        print(name, 'Q', q, 'margin', row['margin_bits_diagnostic'], 'worst pure', min(pure), flush=True)
    dependencies = [Path(__file__), Path(poly.__file__), Path(tight.__file__), Path(groups.__file__),
                    Path(cap_loader.__file__), Path(hull.__file__), Path(scaled.__file__),
                    base.HERE/'general_occupancy.py', base.HERE/'activation_bridge.py',
                    base.HERE/'screen_constant_row_split.py'] + cap_sources
    base.write_new(path, dict(status='LARGER_STATE_RANGE_SCREEN_ONLY', configuration=name,
        bands=bands, rows=rows, used_caps={str(w):caps[w] for w in base.WEIGHTS},
        local_sha256={**maps.sources(), **{str(p.relative_to(base.HERE)):base.sha(p) for p in dependencies}},
        all_occupations_certified=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=('q1', 'range'))
    parser.add_argument('--configuration', choices=maps.NAMES, required=True)
    parser.add_argument('--occupancies', nargs='+', type=int, default=[128, 1024, 2048, 4096, 8192])
    parser.add_argument('--tilt', type=int, default=-35)
    parser.add_argument('--tag', default='first')
    args = parser.parse_args()
    assert args.tag.isidentifier() and all(1 <= q <= 8192 for q in args.occupancies)
    if args.mode == 'q1':
        q1_screen(args.configuration)
    else:
        range_screen(args.configuration, args.occupancies, args.tilt, args.tag)
