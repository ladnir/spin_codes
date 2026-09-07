"""Bounded, serial dense-occupation screens for fixed BCH-256 inner maps.

Uses the activation-aware four-state transfer and the existing square-root
region-density domination. Samples occupancies and type proportions only:
neither a passing nor a weak screen is a full distance conclusion.
"""
import argparse
import math
from pathlib import Path
import time

import numpy as np
from scipy.special import logsumexp, xlogy, xlog1py

import bridge as base
import dual_track_q1 as q1
import larger_state_maps as maps
import k30_sparse as caps_module
import k30_dense_unrestricted_screen as rows_module
import poisson_density_factor as density
import occupation_refresh_v1 as refresh


def load(name):
    if name in maps.NAMES:
        return maps.load(name)
    t, s, spectrum = base.load_map(name)
    data = base.read(base.HERE / 'inputs' / f'{name}_b_kernel_spectrum.json')
    kernel = {r['total_weight']: r['kernel_words'] for r in data['by_total_weight']}
    return t, s, spectrum, kernel


def binomial_logs(t, rates):
    rates = np.asarray(rates, dtype=float)
    if rates.ndim != 1 or np.any(~np.isfinite(rates)) or np.any((rates < 0) | (rates > 1)):
        raise ValueError('A finite vector of probabilities is required')
    j = np.arange(t + 1)
    return (np.array([math.log(math.comb(t, int(v))) for v in j])[None, :]
            + xlogy(j[None, :], rates[:, None])
            + xlog1py((t-j)[None, :], -rates[:, None]))


def mixture(epoch, weights):
    return logsumexp(weights[:, :, None, None] + epoch[None, :, :, :], axis=1)


def source_hashes():
    files = {Path(__file__), Path(base.__file__), Path(maps.__file__),
             Path(caps_module.__file__), Path(rows_module.__file__),
             Path(density.__file__), Path(refresh.__file__),
             Path(q1.refresh.__file__), Path(refresh.old.__file__),
             base.HERE/'k30_dense_screen.py', base.HERE/'general_batch_certificate.py',
             base.HERE/'occupation_three.py', base.HERE/'POISSON_DENSITY_REFINEMENT.md',
             base.HERE/'inputs/manifest.json', maps.DIRECTORY/'manifest.json',
             base.HERE/'generated/christoffel_oa29_caps.json'}
    files.update((base.HERE/'inputs').glob('*.json'))
    files.update(maps.DIRECTORY.glob('*.json'))
    files.update((base.HERE/'generated').glob('joint_shell_*/cap.json'))
    return {p.relative_to(base.ROOT).as_posix(): base.sha(p) for p in sorted(files)}


def run(name, exponent, occupations, output, seconds, grid_size, tilt_count):
    if output.exists():
        raise FileExistsError('Use a new output path; screens are immutable')
    t, s, spectrum, kernel = load(name)
    if exponent < 7:
        raise ValueError('Message dimension must contain complete outer rows')
    length = 1 << (exponent - 7)
    if length < t or length % t or not occupations or any(q < 2 or q > length for q in occupations):
        raise ValueError('Invalid epochs or occupancies')
    if seconds <= 0 or grid_size < 3 or tilt_count < 3:
        raise ValueError('Positive time and nontrivial witness grids required')
    cutoff = 256 * length // 10
    counts = caps_module.caps()
    prepared = refresh.Epochs(t, s, spectrum, [kernel.get(j, 0) for j in range(t+1)])
    banks = []
    for slope in (0., .75, 1.25):
        ps, hull = rows_module.row_witnesses(counts, slope)
        grid = np.unique(np.r_[np.linspace(min(map(float, ps)), 1., grid_size),
                               [x for x, _ in hull]])
        costs = np.interp(grid, [x for x, _ in hull], [y for _, y in hull])
        banks.append((slope, ps, hull, grid, costs))
    results = []
    started = time.monotonic()
    for q in occupations:
        if time.monotonic()-started >= seconds:
            break
        fraction = q / length
        tilts = np.geomspace(max(1e-7, fraction/100), min(6., 100*fraction), tilt_count)
        weights = [binomial_logs(t, fraction * grid) for _, _, _, grid, _ in banks]
        best = [np.full(len(grid), np.inf) for _, _, _, grid, _ in banks]
        witnesses = [np.zeros(len(grid)) for _, _, _, grid, _ in banks]
        for lam in tilts:
            epoch = prepared.at(float(lam))
            for i, law in enumerate(weights):
                value = refresh.terminal_logs(mixture(epoch, law), 256*length//t) + cutoff*lam
                better = value < best[i]
                best[i][better] = value[better]
                witnesses[i][better] = lam
        count_cost = (math.log(math.comb(length, q)) + q*math.log(len(caps_module.BANDS))
                      + 256*math.log(density.density_factor(length)))
        comparisons = []
        for i, (slope, _, _, grid, costs) in enumerate(banks):
            margins = -(count_cost + q*costs + best[i])/math.log(2)
            if not np.isfinite(margins).all():
                raise ArithmeticError('Nonfinite diagnostic')
            index = int(np.argmin(margins))
            comparisons.append(dict(slope=slope, margin_bits=float(margins[index]),
                worst_sample_active_density=float(grid[index]),
                tilt=float(witnesses[i][index]), sampled_densities=len(grid)))
        selected = max(comparisons, key=lambda r: r['margin_bits'])
        results.append(dict(occupation=q, **selected, slope_comparisons=comparisons))
        print(name, 'Q', q, 'dense sample margin', round(selected['margin_bits'], 4),
              'seconds', round(time.monotonic()-started, 2), flush=True)
    base.write_new(output, dict(status='BCH256_INNER_DENSE_SAMPLED_DIAGNOSTIC',
        configuration=name, message_exponent=exponent, rows=length, cutoff=cutoff,
        minimum_A_weight=min(spectrum), minimum_kernel_weight=min(j for j,n in kernel.items() if j and n),
        requested_occupancies=occupations, results=results, full_distance_proved=False,
        completed_requested_samples=len(results)==len(occupations),
        density_factor=density.density_factor(length), tilt_count=tilt_count,
        seconds=time.monotonic()-started, source_sha256=source_hashes(),
        shell_caps={str(w):n for w,n in counts.items()},
        row_banks=[dict(slope=slope, probabilities=[base.encode(p) for p in ps], hull=hull)
                   for slope,ps,hull,_,_ in banks],
        limitations=['Nearest binary64, not an outward certificate.',
                     'No coverage between sampled occupancies or active densities.',
                     'Weak upper bounds do not prove a first-moment obstruction or encoder failure.',
                     'The budget is checked between occupations; one occupation may exceed it.',
                     'Fixed selected maps, fresh independent multipliers, no performance claim.']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--configuration', choices=(*maps.NAMES, *base.CONFIGS), required=True)
    parser.add_argument('--m', type=int, default=20)
    parser.add_argument('--occupancies', nargs='+', type=int, default=[512, 1024, 2620, 4096, 6144, 8192])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=120)
    parser.add_argument('--grid-size', type=int, default=65)
    parser.add_argument('--tilt-count', type=int, default=41)
    args = parser.parse_args()
    run(args.configuration, args.m, args.occupancies, args.output, args.seconds, args.grid_size, args.tilt_count)
