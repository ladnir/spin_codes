"""Why do boxes fail to certify?  (T7c diagnostic)

A box is split again whenever `leaf_value` does not return a strictly negative
upper bound.  That happens for two very different reasons:

  unevaluable : `leaf_value` returned None -- some log argument straddles zero,
                so the bound does not exist on this box at all.  Subdividing
                such a box produces children that also straddle, forever.
  positive    : the bound exists but is not yet negative.  Subdividing this
                *does* make progress.

`infeasible` only prunes boxes lying wholly outside the feasible region, so a
box straddling the boundary is neither pruned nor evaluable.  If the deep
failures are overwhelmingly `unevaluable`, the cascade is a boundary artifact
and not a looseness in the bound.

Usage:  python -B -u scripts/majorant_why.py <segment> [cap] [mindepth]
"""
import sys, os, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixmirror as F
import majorant as M

S = F.SCALE


def run(seg, maxn, mindepth):
    lg = F.directLog(M.TERMS)
    X0 = (int(seg["w0"] * S), int(seg["w1"] * S))
    A0, C0 = (0, S), (0, S)
    wA, wC, wX = A0[1] - A0[0], C0[1] - C0[0], X0[1] - X0[0]
    st = {"nodes": 0}
    tally = collections.Counter()
    samples = []

    def classify(A, C, X):
        """Which piece is missing, if any."""
        u = M.pick_u(A)
        if F.fgObjL(lg, F.sc(u), A) is None:
            return "unevaluable:g"
        if F.fpiBestClampL(lg, A, C) is None:
            return "unevaluable:pi(a,c)"
        if F.fpiBestClampL(lg, C, X) is None:
            return "unevaluable:pi(c,x)"
        return "positive"

    def straddles(A, C, X):
        """Does the box cross the feasibility boundary rather than lie inside?"""
        inside = (2 * C[0] >= A[1] and 2 * S - A[1] >= 2 * C[1]
                  and 2 * X[0] >= C[1] and 2 * S - C[1] >= 2 * X[1])
        return not inside

    def build(A, C, X, depth):
        st["nodes"] += 1
        if st["nodes"] > maxn:
            raise StopIteration
        if M.infeasible(A, C, X):
            return
        v = M.leaf_value(lg, M.pick_u(A), A, C, X, seg["slope"], seg["intercept"])
        if v is not None and v[1] < 0:
            return
        if depth >= mindepth:
            kind = classify(A, C, X)
            tally[(kind, straddles(A, C, X))] += 1
            if len(samples) < 6 and kind.startswith("unevaluable"):
                samples.append((depth, kind,
                                (A[0] / S, A[1] / S), (C[0] / S, C[1] / S),
                                (X[0] / S, X[1] / S)))
        if depth >= M.MAXDEPTH:
            return
        wa = (A[1] - A[0]) / wA
        wc = (C[1] - C[0]) / wC
        wx = (X[1] - X[0]) / wX if wX else 0.0
        if wa >= wc and wa >= wx:
            m = F.fdiv(A[0] + A[1], 2)
            build((A[0], m), C, X, depth + 1); build((m, A[1]), C, X, depth + 1)
        elif wc >= wx:
            m = F.fdiv(C[0] + C[1], 2)
            build(A, (C[0], m), X, depth + 1); build(A, (m, C[1]), X, depth + 1)
        else:
            m = F.fdiv(X[0] + X[1], 2)
            build(A, C, (X[0], m), depth + 1); build(A, C, (m, X[1]), depth + 1)

    sys.setrecursionlimit(200000)
    try:
        build(A0, C0, X0, 0)
    except StopIteration:
        pass
    return st["nodes"], tally, samples


if __name__ == "__main__":
    segs, _ = M.load_segments()
    i = int(sys.argv[1])
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 150000
    mind = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    n, tally, samples = run(segs[i], cap, mind)
    print(f"segment {i}: {n} nodes, failures at depth >= {mind}\n")
    tot = sum(tally.values())
    print(f"{'reason':<24} {'straddles':<10} {'count':>9}  {'share':>7}")
    for (kind, strad), c in tally.most_common():
        print(f"{kind:<24} {str(strad):<10} {c:>9}  {100*c/tot:>6.2f}%")
    print(f"{'TOTAL':<24} {'':<10} {tot:>9}")
    if samples:
        print("\nsample unevaluable boxes:")
        for d, k, A, C, X in samples:
            print(f"  depth {d:>3} {k:<22} a=[{A[0]:.6f},{A[1]:.6f}] "
                  f"c=[{C[0]:.6f},{C[1]:.6f}] x=[{X[0]:.6f},{X[1]:.6f}]")
