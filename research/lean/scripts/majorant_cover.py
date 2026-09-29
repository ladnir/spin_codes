"""Generate replayable covers for the paper's refined BA majorant.

All accepted bounds use exact fixed-point operations mirrored from MajorantDefs.
Floating point is used only to select the positive generating-function witness.
Every rational endpoint is rounded OUTWARDS, including segment intersections.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import time

import fixmirror as F
import majorant as M

S = F.SCALE
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / 'workstreams/inner_design/imt_asymptotic/d11/OUTER_REFINED.json'


def centered(lg, u, A, B, W, slope, intercept):
    if not (F.feasOK(A, B) and F.feasOK(B, W)):
        return None
    am, bm, wm = F.midPt(A), F.midPt(B), F.midPt(W)
    vals = [F.fgObjL(lg, u, am), F.fpiEvalL(lg, am, bm),
            F.fpiEvalL(lg, bm, wm), F.flogIWL(lg, u),
            F.fDpaL(lg, A, B), F.fDpcL(lg, A, B),
            F.fDpaL(lg, B, W), F.fDpcL(lg, B, W)]
    if any(v is None for v in vals):
        return None
    g, p, q, lu, pa, pb, qb, qw = vals
    c = F.sub(F.add(F.add(g, p), q), F.add(F.mul(slope, wm), intercept))
    return F.add(F.add(F.add(c, F.mul(F.sub(pa, lu), F.sub(A, am))),
                       F.mul(F.add(pb, qb), F.sub(B, bm))),
                 F.mul(F.sub(qw, slope), F.sub(W, wm)))


def plain(lg, u, A, B, W, slope, intercept):
    vals = [F.fgObjL(lg, u, A), F.fpiBestClampL(lg, A, B),
            F.fpiBestClampL(lg, B, W)]
    if any(v is None for v in vals):
        return None
    g, p, q = vals
    return F.sub(F.add(F.add(g, p), q), F.add(F.mul(slope, W), intercept))


@lru_cache(maxsize=65536)
def witness(A):
    return M.pick_u(A)


def build_segment(seg, terms=10, cap=2000000):
    lg = lru_cache(maxsize=500000)(F.directLog(terms))
    slope, intercept = M.frac_fix(seg['slope']), M.frac_fix(seg['intercept'])
    A0, B0 = (0, S), (0, S)
    W0 = (M.frac_fix(seg['w0'])[0], M.frac_fix(seg['w1'])[1])
    state = dict(nodes=0, leaves=0, deepest=0, centered=0, plain=0, infeasible=0)
    start = last = time.monotonic()

    def visit(A, B, W, depth):
        nonlocal last
        state['nodes'] += 1
        state['deepest'] = max(state['deepest'], depth)
        if state['nodes'] > cap or depth > 100:
            raise RuntimeError((state, A, B, W))
        now = time.monotonic()
        if now-last > 15:
            print(dict(state, seconds=round(now-start, 1)), flush=True)
            last = now
        if F.infeasible(A, B, W):
            state['infeasible'] += 1
            state['leaves'] += 1
            return ['leaf', 1]
        u = witness(A)
        v = centered(lg, F.sc(u), A, B, W, slope, intercept)
        kind = 'centered'
        if v is None or v[1] >= 0:
            v = plain(lg, F.sc(u), A, B, W, slope, intercept)
            kind = 'plain'
        if v is not None and v[1] < 0:
            state[kind] += 1
            state['leaves'] += 1
            return ['leaf', u if kind == 'centered' else -u]
        # The original generator uses absolute widths and double weight width.
        d = max(range(3), key=lambda j: (A[1]-A[0], B[1]-B[0], 2*(W[1]-W[0]))[j])
        box = [A, B, W]
        lo, hi = box[d]
        m = (lo+hi)//2
        if not lo < m < hi:
            raise RuntimeError(('indivisible', box, v))
        left, right = box.copy(), box.copy()
        left[d], right[d] = (lo, m), (m, hi)
        return ['split', d, m, visit(*left, depth+1), visit(*right, depth+1)]

    tree = visit(A0, B0, W0, 0)
    return dict(tree=tree, A=A0, B=B0, W=W0, slope=slope, intercept=intercept,
                terms=terms, stats=dict(state, seconds=time.monotonic()-start))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('segment', help='index or all')
    parser.add_argument('--cap', type=int, default=2000000)
    parser.add_argument('--output', type=Path, default=HERE/'majorant_data')
    args = parser.parse_args()
    data = json.loads(SOURCE.read_text())
    segs = [dict(slope=M.Fr(s['slope']), intercept=M.Fr(s['intercept']),
                 w0=M.Fr(s['weight'][0]), w1=M.Fr(s['weight'][1]))
            for s in data['left_segments']]
    args.output.mkdir(exist_ok=True)
    for i in range(len(segs)) if args.segment == 'all' else [int(args.segment)]:
        print('segment', i, flush=True)
        result = build_segment(segs[i], cap=args.cap)
        result['segment'] = i
        result['source'] = str(SOURCE)
        result['rational_segment'] = data['left_segments'][i]
        (args.output/f'segment{i}.json').write_text(json.dumps(result)+'\n')
        print(result['stats'], flush=True)


if __name__ == '__main__':
    main()
