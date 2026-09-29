"""Will the segment-0 cover actually finish?  (T7c)

`majorant_split.py` reports node counts, which say nothing about how much of
the search is left.  The build is a DFS over a 3-D box, so the *resolved
volume* -- the total volume of boxes that have been closed or pruned -- is an
exact progress measure, and elapsed/fraction is a real ETA.

Also records the deepest box reached, because the run aborts at MAXDEPTH
regardless of how much time is available.

Usage:  python -B -u scripts/majorant_progress.py <segment> [cap]
"""
import sys, os, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixmirror as F
import majorant as M

S = F.SCALE
TERMS = M.TERMS
MAXDEPTH = int(os.environ.get("MAJ_MAXDEPTH", M.MAXDEPTH))


def run(seg, thr, maxn):
    lg = F.directLog(TERMS)
    X0 = (int(seg["w0"] * S), int(seg["w1"] * S))
    # same feasibility-tightened start as majorant_split.py
    cmax = min(S, 2 * X0[1])
    amax = min(S, 2 * cmax)
    A0 = (0, amax)
    C0 = (0, cmax)
    wA, wC, wX = A0[1] - A0[0], C0[1] - C0[0], X0[1] - X0[0]
    V0 = wA * wC * wX
    st = {"nodes": 0, "vol": 0, "deep": 0}
    t0 = time.time()

    def vol(A, C, X):
        return (A[1] - A[0]) * (C[1] - C[0]) * (X[1] - X[0])

    def report():
        el = time.time() - t0
        fr = st["vol"] / V0
        eta = (el / fr - el) if fr > 0 else float("inf")
        print(f"    nodes={st['nodes']:>9d}  resolved={fr:.9f}  "
              f"depth<={st['deep']}  elapsed={el/60:.1f}m  "
              f"eta={eta/3600:.2f}h", flush=True)

    def build(A, C, X, depth):
        st["nodes"] += 1
        st["deep"] = max(st["deep"], depth)
        if st["nodes"] > maxn:
            raise RuntimeError("node cap")
        if st["nodes"] % 20000 == 0:
            report()
        if M.infeasible(A, C, X):
            st["vol"] += vol(A, C, X)
            return 1
        v = M.leaf_value(lg, M.pick_u(A), A, C, X, seg["slope"], seg["intercept"])
        if v is not None and v[1] < thr:
            st["vol"] += vol(A, C, X)
            return 1
        if depth >= MAXDEPTH:
            raise RuntimeError(f"depth cap at {A} {C} {X}")
        # widest relative to initial width
        wa = (A[1] - A[0]) / wA
        wc = (C[1] - C[0]) / wC
        wx = (X[1] - X[0]) / wX if wX else 0.0
        if wa >= wc and wa >= wx:
            d = 0
        elif wc >= wx:
            d = 1
        else:
            d = 2
        if d == 0:
            m = F.fdiv(A[0] + A[1], 2)
            lo, hi = ((A[0], m), C, X), ((m, A[1]), C, X)
        elif d == 1:
            m = F.fdiv(C[0] + C[1], 2)
            lo, hi = (A, (C[0], m), X), (A, (m, C[1]), X)
        else:
            m = F.fdiv(X[0] + X[1], 2)
            lo, hi = (A, C, (X[0], m)), (A, C, (m, X[1]))
        return build(*lo, depth + 1) + build(*hi, depth + 1)

    sys.setrecursionlimit(200000)
    leaves = build(A0, C0, X0, 0)
    report()
    return st["nodes"], leaves, st["deep"]


if __name__ == "__main__":
    segs, _ = M.load_segments()
    i = int(sys.argv[1])
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 20000000
    s = segs[i]
    print(f"segment {i}: slope {float(s['slope']):.6g} "
          f"x in [{float(s['w0']):.6f}, {float(s['w1']):.6f}]", flush=True)
    try:
        n, lv, dp = run(s, 0, cap)
        print(f"  closed: {n} nodes, {lv} leaves, deepest {dp}", flush=True)
    except RuntimeError as e:
        print(f"  FAILED: {e}", flush=True)
