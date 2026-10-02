"""Support reduction through cheap hyperplanes of binary subcodes.

An h-space has m=2^h-1 hyperplanes. Each active coordinate disappears in
exactly one hyperplane, so their average support is (m-1)/m times its
support. Count the necessarily small hyperplanes, then bound extensions
inside shortened codes. All arithmetic is integral; no LP solver is used.
"""
import argparse
import json
from math import comb, log2
from pathlib import Path
import shared_support
from bch_joint_support import rank_total
from count_refinements import refine_pair


def bound(caps, dimensions, g, h, u, last=None):
    n = len(dimensions)-1
    last = n if last is None else last
    if not 2 <= h <= g == len(caps) or not 0 <= u <= last <= n:
        raise ValueError('valid rank and support range required')
    lower = caps[h-2]
    if not any(lower):
        return 0, None
    minimum = next(v for v, count in enumerate(lower) if count)
    if minimum > u:
        return 0, None
    spaces = [x//rank_total(h-1, g, h-1) for x in lower]
    if any(a > b for a, b in zip(spaces, spaces[1:])):
        raise ValueError('nondecreasing lower-rank CDF caps required')
    multiplicity = rank_total(h, g, h)
    m = (1 << h)-1
    best, witness = caps[h-1][u], None
    for stop in range(max(minimum, (m-1)*u//m), u+1):
        numerator = m*(stop+1)-(m-1)*u
        denominator = stop+1-minimum
        good = m if stop == u else min(m, (numerator+denominator-1)//denominator)
        if good <= 0:
            continue
        for t in range(u, last+1):
            extensions = (1 << (dimensions[t]-h+1))-1 if dimensions[t] >= h else 0
            # A decreasing weight, truncated after stop, permits CDF
            # summation by parts. Differences are not treated as shell caps.
            moment = sum((spaces[v]-(spaces[v-1] if v else 0))*comb(n-v, t-v)
                         for v in range(stop+1))
            upper = multiplicity*extensions*moment//(good*comb(n-u, t-u))
            if upper < best:
                best, witness = upper, dict(stop=stop, containing_set=t, guaranteed_hyperplanes=good)
    return best, witness


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supports', type=int, nargs='+', default=[104, 108, 112, 116, 119, 120, 128])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    _, raw = shared_support.shared_counts(True)
    spectrum = shared_support.authenticated_caps()
    caps, _, dims, _ = refine_pair(raw, spectrum)
    rows = []
    for u in args.supports:
        upper, witness = bound(caps, dims, 4, 4, u, last=192)
        gain = log2(caps[3][u])-log2(upper) if upper else None
        print('HYPERPLANE support', u, 'gain bits', gain, 'witness', witness, flush=True)
        rows.append(dict(support=u, prior=str(caps[3][u]), upper=str(upper), gain_bits=gain, witness=witness))
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(dict(schema='shared-gf16-hyperplane-counts-1', rows=rows), indent=2)+'\n')


if __name__ == '__main__':
    main()
