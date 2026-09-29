"""Constraint propagation at every node.  (T7c)

Diagnostic finding (`majorant_why.py`): 81% of deep failures are boxes that
*straddle* the feasibility boundary.  Their bound must cover infeasible points,
where the inequality is not even claimed, so they cannot certify however small
they get -- only shrinking until they fall wholly inside or wholly outside
helps.  That means resolving a 2-D surface at extreme resolution, which is the
astronomical box count.

Feasibility is `a/2 < c < 1-a/2` and `c/2 < x < 1-c/2`.  Each is a pair of
linear constraints, so a box can be *narrowed* against them before bounding:

    c <= 2x, c <= 2(1-x), a <= 2c, a <= 2(1-c), x >= c/2, x <= 1-c/2

Narrowing is sound: it discards only points that violate feasibility, and the
claim says nothing about those.  Iterated to a fixed point, a straddling box
becomes a smaller interior box rather than a subdivision cascade.

Usage:  python -B -u scripts/majorant_narrow.py <segment> [cap]
"""
import sys, os, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixmirror as F
import majorant as M

S = F.SCALE


def narrow(A, C, X):
    """Shrink (A,C,X) to its feasible part.  None if provably empty."""
    a0, a1 = A
    c0, c1 = C
    x0, x1 = X
    for _ in range(8):
        pa0, pa1, pc0, pc1, px0, px1 = a0, a1, c0, c1, x0, x1
        # a/2 < c < 1 - a/2   =>   a <= 2c and a <= 2(1-c); c >= a/2, c <= 1-a/2
        a1 = min(a1, 2 * c1, 2 * (S - c0))
        c0 = max(c0, (a0 + 1) // 2)
        c1 = min(c1, S - (a0 + 1) // 2)
        # c/2 < x < 1 - c/2   =>   c <= 2x and c <= 2(1-x); x >= c/2, x <= 1-c/2
        c1 = min(c1, 2 * x1, 2 * (S - x0))
        x0 = max(x0, (c0 + 1) // 2)
        x1 = min(x1, S - (c0 + 1) // 2)
        if a1 < a0 or c1 < c0 or x1 < x0:
            return None
        if (a0, a1, c0, c1, x0, x1) == (pa0, pa1, pc0, pc1, px0, px1):
            break
    return (a0, a1), (c0, c1), (x0, x1)


def run(seg, maxn):
    lg = F.directLog(M.TERMS)
    X0 = (int(seg["w0"] * S), int(seg["w1"] * S))
    root = narrow((0, S), (0, S), X0)
    assert root is not None
    A0, C0, Xr = root
    wA = max(1, A0[1] - A0[0]); wC = max(1, C0[1] - C0[0]); wX = Xr[1] - Xr[0]
    V0 = wA * wC * max(1, wX)
    st = {"nodes": 0, "vol": 0, "deep": 0, "leaves": 0}
    t0 = time.time()

    def report():
        el = time.time() - t0
        fr = st["vol"] / V0
        eta = (el / fr - el) if fr > 0 else float("inf")
        print(f"    nodes={st['nodes']:>9d}  resolved={fr:.9f}  "
              f"depth<={st['deep']}  elapsed={el/60:.1f}m  eta={eta/3600:.2f}h",
              flush=True)

    def build(box, depth):
        st["nodes"] += 1
        st["deep"] = max(st["deep"], depth)
        if st["nodes"] > maxn:
            raise RuntimeError("node cap")
        if st["nodes"] % 5000 == 0:
            report()
        A, C, X = box
        vol = (A[1]-A[0]) * (C[1]-C[0]) * (X[1]-X[0])
        nb = narrow(A, C, X)
        if nb is None:                      # wholly infeasible
            st["vol"] += vol; st["leaves"] += 1
            return
        A, C, X = nb
        # the sliver narrowing discarded is infeasible, hence resolved: count it
        # now, or the fraction can never reach 1 even on success
        vol = (A[1]-A[0]) * (C[1]-C[0]) * (X[1]-X[0])
        st["vol"] += (box[0][1]-box[0][0]) * (box[1][1]-box[1][0])             * (box[2][1]-box[2][0]) - vol
        v = M.leaf_value(lg, M.pick_u(A), A, C, X, seg["slope"], seg["intercept"])
        if v is not None and v[1] < 0:
            st["vol"] += vol; st["leaves"] += 1
            return
        if depth >= 200:
            raise RuntimeError(f"depth cap at {A} {C} {X}")
        wa = (A[1]-A[0]) / wA; wc = (C[1]-C[0]) / wC
        wx = (X[1]-X[0]) / wX if wX else 0.0
        if wa >= wc and wa >= wx:
            m = F.fdiv(A[0]+A[1], 2); lo, hi = ((A[0], m), C, X), ((m, A[1]), C, X)
        elif wc >= wx:
            m = F.fdiv(C[0]+C[1], 2); lo, hi = (A, (C[0], m), X), (A, (m, C[1]), X)
        else:
            m = F.fdiv(X[0]+X[1], 2); lo, hi = (A, C, (X[0], m)), (A, C, (m, X[1]))
        build(lo, depth+1); build(hi, depth+1)

    sys.setrecursionlimit(200000)
    build((A0, C0, Xr), 0)
    report()
    return st


if __name__ == "__main__":
    segs, _ = M.load_segments()
    i = int(sys.argv[1])
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 5000000
    s = segs[i]
    print(f"segment {i}: slope {float(s['slope']):.6g} "
          f"x in [{float(s['w0']):.6f}, {float(s['w1']):.6f}]", flush=True)
    try:
        st = run(s, cap)
        print(f"  closed: {st['nodes']} nodes, {st['leaves']} leaves, "
              f"deepest {st['deep']}", flush=True)
    except RuntimeError as e:
        print(f"  FAILED: {e}", flush=True)
