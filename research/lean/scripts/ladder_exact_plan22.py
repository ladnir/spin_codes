"""Seeded exact-piece planner for the s22 ladder tail (untrusted producer).

For each q-range [lo, hi] it searches a small tilt/alpha window and splits the range
(geometric midpoint) until max(exact_term(lo), exact_term(hi)) meets the target.
The seed is the empirical optimum lambda*d ~= 40*lo seen at G20/G22 (independent of
L), searched over 10*2^(k/4) within +-1 octave, alpha in {0.55..1}. Lean checks every
q of every piece exactly; this only proposes pieces.

Usage: python scripts/ladder_exact_plan22.py e out.json lo:hi [lo:hi ...]
"""
import json
import sys
from fractions import Fraction as Q
from math import log2

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import ladder_tail22 as T  # noqa: E402

ALS = [Q(11, 20), Q(3, 5), Q(13, 20), Q(7, 10), Q(4, 5), Q(1)]


def candidates(lo):
    k0 = round(4 * log2(max(40 * lo, 10) / 10))
    return sorted({max(0, k) for k in range(k0 - 4, k0 + 5)})


def best_for(e, lo, hi, D):
    best = (1e18, None)
    for k in candidates(lo):
        lam = T.lam_of(round(10 * 2 ** (k / 4)), D)
        f = T.fam(lam)
        for al in ALS:
            v = T.exact_term(e, lo, f, float(al))
            if v >= best[0]:
                continue          # the max over endpoints cannot beat the current best
            v = max(v, T.exact_term(e, hi, f, float(al)))
            if v < best[0]:
                best = (v, (lam, al))
    return best


def main():
    e, out = int(sys.argv[1]), sys.argv[2]
    work = [tuple(map(int, r.split(':'))) for r in sys.argv[3:]][::-1]
    L, NP, D = T.geom(e)
    cfg = json.loads((T.ROOT / f'.lake/turbo/ladder/g{e}/config.json').read_text())
    target = -(cfg['q1_bits'] + 20)
    pieces = []
    while work:
        lo, hi = work.pop()
        v, arg = best_for(e, lo, hi, D)
        if v > target and hi > lo:
            mid = min(max(int(round((lo * (hi + 1)) ** 0.5)), lo + 1), hi)
            print(f'  split exact {lo}..{hi} (worst {v:.1f}) at {mid}', flush=True)
            work += [(mid, hi), (lo, mid - 1)]
            continue
        assert v <= target, (lo, hi, v)
        lam, al = arg
        pieces.append(dict(kind='exact', tilt=f'{lam.numerator}/{lam.denominator}',
                           alpha=f'{al.numerator}/{al.denominator}', p=None, qlo=lo, qhi=hi))
        print(f'  exact {lo}..{hi}: lambda*d {float(lam) * D:.1f} alpha {al} worst {v:.1f}', flush=True)
        with open(out, 'w') as fh:
            json.dump(sorted(pieces, key=lambda p: p['qlo']), fh, indent=1)
    print('EXACTDONE', len(pieces), flush=True)


if __name__ == '__main__':
    main()
