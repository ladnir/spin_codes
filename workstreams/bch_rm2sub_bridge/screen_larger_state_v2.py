"""Larger-state discovery with the all-one band kept exact (p=1, cost=1)."""
import argparse
import math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb, ctx
import bridge as base
import larger_state_maps as maps
import polynomial_regions as poly
import occupation_three as groups
import screen_exponential_modes as cap_loader
import scaled_adaptive as scaled
import positive_line_hull as hull
from screen_constant_row_split import moment


def region_logs(name, tilt):
    path = base.HERE/'generated'/f'larger_v2_{name}_region_{tilt}.npy'
    receipt = path.with_suffix('.json')
    if receipt.exists():
        saved = base.read(receipt)
        assert base.sha(path) == saved['array_sha256']
        for source, digest in saved['local_sha256'].items():
            assert base.sha(base.HERE/source) == digest
        return np.load(path, allow_pickle=False)
    assert not path.exists()
    t, s, spectrum, kernel = maps.load(name)
    ctx.prec = 256
    lam = (arb(tilt)/10).exp()
    coefficients = poly.regions(t, s, spectrum, kernel, (-lam).exp(), 8192)
    values = np.array([[float(v.log()) if v > 0 else -math.inf for v in row] for row in coefficients]).reshape(-1, 3, 3)
    with path.open('xb') as stream:
        np.save(stream, values, allow_pickle=False)
    files = [Path(__file__), Path(poly.__file__), base.HERE/'tightened_occupancy.py',
             base.HERE/'general_occupancy.py', base.HERE/'activation_bridge.py']
    base.write_new(receipt, dict(status='LARGER_STATE_DISCOVERY_LOG_CACHE_ONLY', configuration=name,
        tilt=tilt, array_sha256=base.sha(path), local_sha256={**maps.sources(),
        **{str(p.relative_to(base.HERE)):base.sha(p) for p in files}}))
    return values


def screen(name, qs, tilt, tag):
    output = base.HERE/'generated'/f'larger_v2_{name}_{tag}_screen.json'
    assert not output.exists()
    caps, cap_sources = cap_loader.latest_caps()
    bands = groups.BANDS
    assert bands[-1] == (256,) and caps[256] == 1
    weights = [np.array(b) for b in bands[:-1]]
    logs = [np.array([math.log(caps[w])-math.log(math.comb(256, w)) for w in b]) for b in bands[:-1]]
    region = region_logs(name, tilt)
    print(name, 'v2 region ready', tilt, flush=True)
    rows = []
    for q in qs:
        selected = region[:q+1]
        ps, gamma, pure = [], [], []
        offset = base.CUTOFF*math.exp(tilt/10)+math.log(math.comb(8192, q))+q*math.log(len(bands))
        for w, v in zip(weights, logs):
            def objective(theta):
                p = 1/(1+math.exp(-theta))
                g = float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
                return moment(selected, q, p)+q*g
            opt = minimize_scalar(objective, bounds=(-6, 14), method='bounded', options={'xatol':1e-6})
            p = 1/(1+math.exp(-float(opt.x)))
            ps.append(p)
            gamma.append(float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))))
            pure.append(-(float(opt.fun)+offset)/math.log(2))
        ps.append(1.)
        gamma.append(0.)
        pure.append(-(moment(selected, q, 1.)+offset)/math.log(2))
        ps = np.array(ps)
        roots = np.exp(np.array(gamma)/256)
        left, right = roots*(1-ps), roots*ps
        keep = hull.indices(left, right)
        exponents = np.floor(selected.max(axis=(1,2))/math.log(2)).astype(np.int64)
        mantissas = np.exp(selected-exponents[:,None,None]*math.log(2))
        value = scaled.terminal_log(mantissas, exponents, ps[keep], roots[keep])+offset
        row = dict(occupation=q, witness_tenth=tilt, p=[base.encode(F.from_float(float(p))) for p in ps],
                   margin_bits_diagnostic=-value/math.log(2), pure_margins_with_assignment_cost=pure)
        rows.append(row)
        print(name, 'Q', q, 'tilt', tilt, 'margin', row['margin_bits_diagnostic'], 'worst pure', min(pure), flush=True)
    files = [Path(__file__), Path(poly.__file__), Path(groups.__file__), Path(cap_loader.__file__),
             Path(scaled.__file__), Path(hull.__file__), base.HERE/'screen_constant_row_split.py',
             base.HERE/'general_occupancy.py', base.HERE/'tightened_occupancy.py',
             base.HERE/f'generated/larger_v2_{name}_region_{tilt}.json'] + cap_sources
    base.write_new(output, dict(status='LARGER_STATE_EXACT_CONSTANT_BAND_SCREEN_ONLY', configuration=name,
        bands=bands, rows=rows, used_caps={str(w):caps[w] for w in base.WEIGHTS},
        cap_receipts=[str(p.relative_to(base.HERE)) for p in cap_sources],
        local_sha256={**maps.sources(), **{str(p.relative_to(base.HERE)):base.sha(p) for p in files}},
        all_occupations_certified=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--configuration', choices=maps.NAMES, required=True)
    parser.add_argument('--occupancies', nargs='+', type=int, required=True)
    parser.add_argument('--tilt', type=int, required=True)
    parser.add_argument('--tag', required=True)
    args = parser.parse_args()
    assert args.tag.isidentifier() and all(1 <= q <= 8192 for q in args.occupancies)
    screen(args.configuration, args.occupancies, args.tilt, args.tag)
