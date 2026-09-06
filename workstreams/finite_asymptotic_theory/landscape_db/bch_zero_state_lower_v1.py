"""Lower bounds on BCH bad-word first moments from all-zero state paths.

This can diagnose a first-moment obstruction, not the probability of code
failure. Exact formulas are evaluated in nearest binary64 arithmetic.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def positive_product(left, right, maximum):
    out = np.full(min(maximum+1, len(left)+len(right)-1), -np.inf)
    for k in range(len(out)):
        lo, hi = max(0, k-len(right)+1), min(k, len(left)-1)
        out[k] = np.logaddexp.reduce(left[lo:hi+1]+right[k-hi:k-lo+1][::-1])
    return out


def kernel_region_logs(kernel, epochs, maximum):
    # All inputs used here have even regional weight. The studied inner
    # maps include a parity check, so odd kernel coefficients vanish.
    if any(kernel[1::2]):
        raise ValueError('even-kernel shortcut requires the parity check')
    power = np.array([math.log(n) if n else -math.inf for n in kernel[::2]])
    current = np.array([0.]); remaining = epochs
    while remaining:
        if remaining & 1:
            current = positive_product(current, power, maximum//2)
        remaining >>= 1
        if remaining:
            power = positive_product(power, power, maximum//2)
    length = (len(kernel)-1)*epochs
    return current-np.array([math.log(math.comb(length, 2*j)) for j in range(len(current))])


def kl(a, p):
    if a == 0:
        return -math.log1p(-p)
    if a == 1:
        return -math.log(p)
    return a*math.log(a/p)+(1-a)*math.log((1-a)/(1-p))


def concentrated_range(q, p, block):
    """Chernoff interval whose outside mass is negligible against 2^(1-B)."""
    target = (block+5)*LN2+math.log(block)
    lo, hi = 0, int(q*p)
    while lo < hi:
        mid = (lo+hi+1)//2
        if q*kl((mid-1)/q, p) >= target:
            lo = mid
        else:
            hi = mid-1
    lower = lo
    lo, hi = math.ceil(q*p), q
    while lo < hi:
        mid = (lo+hi)//2
        if q*kl((mid+1)/q, p) >= target:
            hi = mid
        else:
            lo = mid+1
    return lower, lo, (1-block)*LN2+math.log1p(-1/32)


def convex_lower_value(xs, ys, mean):
    """Evaluate the lower convex hull, retaining its supporting segment."""
    hull = []
    for x, y in zip(xs, ys):
        if not math.isfinite(y):
            return -math.inf, None
        while len(hull) >= 2:
            a, b = hull[-2:]
            if (b[1]-a[1])*(x-b[0]) < (y-b[1])*(b[0]-a[0]):
                break
            hull.pop()
        hull.append((int(x), float(y)))
    if not hull or not hull[0][0] <= mean <= hull[-1][0]:
        return -math.inf, None
    for a, b in zip(hull, hull[1:]):
        if a[0] <= mean <= b[0]:
            value = a[1]+(mean-a[0])*(b[1]-a[1])/(b[0]-a[0])
            return value, [list(a), list(b)]
    return hull[0][1], [list(hull[0])]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--block', type=int, choices=(64, 128), default=128)
    parser.add_argument('--step', type=int, choices=(64, 128, 256), default=256)
    parser.add_argument('--state', type=int, default=20)
    parser.add_argument('--exponent', type=int, default=20)
    args = parser.parse_args()
    _, spectra, maps, dependencies = study.load_inputs()
    block, t, s = args.block, args.step, args.state
    length = (1 << args.exponent)//(block//2); cutoff = block*length//10
    counts, config = spectra[block], maps[t, s]
    if length % t:
        raise ValueError('whole epochs required')
    queries = []
    for weight, count in counts.items():
        maximum_q = min(length, cutoff//weight)
        for fraction in (.35, .5, .65, .8, .9, 1.):
            q = 2*int(maximum_q*fraction/2)
            if q < 2 or weight == block:
                continue
            p = weight/block; lo, hi, log_good = concentrated_range(q, p, block)
            first = lo+(lo % 2); last = hi-(hi % 2)
            if first < 4 or first > last:
                continue
            queries.append((q, weight, count, first, last, log_good))
    maximum = max(hi for q, w, n, lo, hi, g in queries)
    logs = kernel_region_logs(config['kernel_counts'], length//t, maximum)
    results = []
    for q, weight, count, lo, hi, log_good in queries:
        xs = np.arange(lo, hi+1, 2); ys = logs[lo//2:hi//2+1]
        value, segment = convex_lower_value(xs, ys, q*weight/block)
        if not math.isfinite(value):
            continue
        bound = math.log(math.comb(length, q))+q*math.log(count)+log_good+block*value
        results.append(dict(occupation=q, outer_weight=weight, output_weight=q*weight,
                            regional_weight_range=[lo, hi], lower_convex_segment=segment,
                            good_parity_log_lower=log_good, log_first_moment_lower=bound,
                            first_moment_lower_bits=bound/LN2))
    results.sort(key=lambda r: r['log_first_moment_lower'], reverse=True)
    for path in (Path(__file__),):
        dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
    payload = dict(status='BINARY64_FIRST_MOMENT_LOWER_DIAGNOSTIC', arguments=vars(args),
                   source_sha256=dependencies, kernel_coefficient_maximum=maximum,
                   best=results[0], queries=results,
                   limitations=['A large expected bad-word count does not imply a large failure probability.',
                                'Nearest binary64; no outward lower certificate.',
                                'Each lower bound restricts to even Q and a single exact outer weight shell.'])
    path = HERE/f'bch_zero_state_lower_b{block}_t{t}_s{s}_e{args.exponent}.json'
    path.write_text(json.dumps(payload, indent=2)+'\n')
    print(json.dumps(payload['best'], indent=2), flush=True)


if __name__ == '__main__':
    main()
