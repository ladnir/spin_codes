"""Split-heuristic experiment for the concave-majorant cover (T7c).

`majorant.py` splits the dimension that is widest relative to its own initial
width.  That is the same mis-scaling that inflated the dense-tail cover, and
the certificate's own count (101,940 boxes for the whole claim, 15 left
segments) says the verifier is far more efficient than the branch-and-bound
here.

This compares two rules on the same segment:

  width : the current one, widest relative to initial width
  trial : split each dimension, score by the resulting max upper bound,
          and take the dimension that lowers it most

Usage:  python -B scripts/majorant_split.py <segment> <width|trial> [cap]
"""
import sys, os, math, json
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixmirror as F
import majorant as M

S = F.SCALE
TERMS = M.TERMS
MAXDEPTH = M.MAXDEPTH


def run(seg, thr, maxn, rule):
    lg = F.directLog(TERMS)
    X0 = (int(seg["w0"] * S), int(seg["w1"] * S))
    # Feasibility (c/2 < x < 1-c/2 and a/2 < c < 1-a/2) caps `c` at 2*max(x)
    # and `a` at 2*max(c).  Starting from [0,1]^2 instead makes the search
    # subdivide a large infeasible region before `infeasible` prunes it --
    # which is what stalls the narrow-x segments.
    cmax = min(S, 2 * X0[1])
    amax = min(S, 2 * cmax)
    A0 = (0, amax)
    C0 = (0, cmax)
    wA, wC, wX = A0[1] - A0[0], C0[1] - C0[0], X0[1] - X0[0]
    st = {"nodes": 0, "evals": 0}

    def upper(A, C, X):
        st["evals"] += 1
        v = M.leaf_value(lg, M.pick_u(A), A, C, X, seg["slope"], seg["intercept"])
        return None if v is None else v[1]

    def halves(A, C, X, d):
        if d == 0:
            m = F.fdiv(A[0] + A[1], 2)
            return ((A[0], m), C, X), ((m, A[1]), C, X)
        if d == 1:
            m = F.fdiv(C[0] + C[1], 2)
            return (A, (C[0], m), X), (A, (m, C[1]), X)
        m = F.fdiv(X[0] + X[1], 2)
        return (A, C, (X[0], m)), (A, C, (m, X[1]))

    BIG = 10 ** 40

    def choose(A, C, X):
        if rule == "width":
            wa = (A[1] - A[0]) / wA
            wc = (C[1] - C[0]) / wC
            wx = (X[1] - X[0]) / wX if wX else 0.0
            if wa >= wc and wa >= wx:
                return 0
            return 1 if wc >= wx else 2
        if rule == "trial":
            best, bestscore = 0, None
            for d in (0, 1, 2):
                if d == 2 and not wX:
                    continue
                lo, hi = halves(A, C, X, d)
                ul, uh = upper(*lo), upper(*hi)
                sl = -BIG if ul is None else ul
                sh = -BIG if uh is None else uh
                sc = max(sl, sh)
                if bestscore is None or sc < bestscore:
                    best, bestscore = d, sc
            return best
        # "collapse": how much of the overestimate is attributable to each
        # dimension, measured by pinning it to its midpoint.  Unlike the
        # greedy trial rule this does not starve a dimension whose immediate
        # gain is small but whose total contribution is large.
        base = upper(A, C, X)
        if base is None:
            return 0
        best, bestdrop = 0, None
        for d in (0, 1, 2):
            if d == 2 and not wX:
                continue
            if d == 0:
                m = F.fdiv(A[0] + A[1], 2); box = ((m, m), C, X)
            elif d == 1:
                m = F.fdiv(C[0] + C[1], 2); box = (A, (m, m), X)
            else:
                m = F.fdiv(X[0] + X[1], 2); box = (A, C, (m, m))
            u = upper(*box)
            drop = BIG if u is None else base - u
            if bestdrop is None or drop > bestdrop:
                best, bestdrop = d, drop
        return best

    def build(A, C, X, depth):
        st["nodes"] += 1
        if st["nodes"] > maxn:
            raise RuntimeError("node cap")
        if st["nodes"] % 20000 == 0:
            print(f"    nodes={st['nodes']} evals={st['evals']}", flush=True)
        if M.infeasible(A, C, X):
            return 1
        v = M.leaf_value(lg, M.pick_u(A), A, C, X, seg["slope"], seg["intercept"])
        if v is not None and v[1] < thr:
            return 1
        if depth >= MAXDEPTH:
            raise RuntimeError(f"depth cap at {A} {C} {X}")
        d = choose(A, C, X)
        lo, hi = halves(A, C, X, d)
        return build(*lo, depth + 1) + build(*hi, depth + 1)

    sys.setrecursionlimit(200000)
    leaves = build(A0, C0, X0, 0)
    return st["nodes"], leaves, st["evals"]


if __name__ == "__main__":
    segs, _ = M.load_segments()
    i = int(sys.argv[1])
    rule = sys.argv[2]
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 400000
    s = segs[i]
    print(f"segment {i} rule={rule}: slope {float(s['slope']):.6g} "
          f"x in [{float(s['w0']):.6f}, {float(s['w1']):.6f}]", flush=True)
    try:
        n, lv, ev = run(s, 0, cap, rule)
        print(f"  closed: {n} nodes, {lv} leaves, {ev} extra evals", flush=True)
    except RuntimeError as e:
        print(f"  FAILED: {e}", flush=True)
