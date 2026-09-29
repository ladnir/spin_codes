"""Cost probe for the concave-majorant cover (`eq:ba-spectrum-majorant`).

The claim is `a_BA(x) ≤ â_BA(x)` on `W = [13/125, 112/125]`, where `â_BA` is the
lower envelope of affine supports.  Since the envelope is a minimum, dominating
the *active* (smallest) line on an x-box dominates all of them, so each segment
contributes one 3-D cover of

    g(a) + π(a,c) + π(c,x) - (s·x + intercept) < 0

over `a, c ∈ [0,1]` feasible and `x` in that segment's active interval.

Driven by `fixmirror`, i.e. by exactly the bound Lean checks.

Usage:  python -B scripts/majorant.py <segment-index|all> [limit-nodes]
"""
import sys, os, math, json
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixmirror as F

S = F.SCALE
TERMS = 10
MAXDEPTH = 60

CERT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "workstreams", "paper_architecture", "certificates",
    "single_sampled_ba_rm2sub", "golay_ba3_concave_majorant.json")


def load_segments():
    d = json.load(open(CERT))
    segs = []
    for s in d["left_segments"]:
        segs.append({
            "slope": Fr(s["slope"]),
            "intercept": Fr(s["intercept"]),
            "w0": Fr(s["weight"][0]),
            "w1": Fr(s["weight"][1]),
        })
    return segs, Fr(d["central_constant_upper"])


def frac_fix(q):
    """A Fix point enclosing the rational q."""
    return F.ofFrac(q.numerator, q.denominator)


def gobj_real(u, a):
    G = 1 + 759 * u ** 8 + 2576 * u ** 12 + 759 * u ** 16 + u ** 24
    return math.log(G) / 24 - a * math.log(u)


def pick_u(A):
    a = (A[0] + A[1]) / 2 / S
    lo, hi = 1e-4, 3.0
    for _ in range(200):
        m1 = lo + (hi - lo) / 3
        m2 = hi - (hi - lo) / 3
        if gobj_real(m1, a) < gobj_real(m2, a):
            hi = m2
        else:
            lo = m1
    return max(1, int(round((lo + hi) / 2 * S)))


def infeasible(A, C, X):
    return (2 * C[1] < A[0] or 2 * S - A[0] < 2 * C[0]
            or 2 * X[1] < C[0] or 2 * S - C[0] < 2 * X[0])


def leaf_value(lg, u, A, C, X, slope, intercept):
    g = F.fgObjL(lg, F.sc(u), A)
    if g is None:
        return None
    p1 = F.fpiBestClampL(lg, A, C)
    if p1 is None:
        return None
    p2 = F.fpiBestClampL(lg, C, X)
    if p2 is None:
        return None
    line = F.add(F.mul(frac_fix(slope), X), frac_fix(intercept))
    return F.sub(F.add(F.add(g, p1), p2), line)


def run_segment(seg, thr, maxn):
    lg = F.directLog(TERMS)
    A0 = (0, S)
    C0 = (0, S)
    X0 = (int(seg["w0"] * S), int(seg["w1"] * S))
    wA, wC, wX = A0[1] - A0[0], C0[1] - C0[0], X0[1] - X0[0]
    state = {"nodes": 0, "maxacc": None}

    def widest(A, C, X):
        wa = (A[1] - A[0]) / wA
        wc = (C[1] - C[0]) / wC
        wx = (X[1] - X[0]) / wX if wX else 0.0
        if wa >= wc and wa >= wx:
            return 0
        return 1 if wc >= wx else 2

    def build(A, C, X, depth):
        state["nodes"] += 1
        if state["nodes"] > maxn:
            raise RuntimeError("node cap")
        if state["nodes"] % 4000 == 0:
            print(f"    nodes={state['nodes']}", flush=True)
        if infeasible(A, C, X):
            return 1
        u = pick_u(A)
        v = leaf_value(lg, u, A, C, X, seg["slope"], seg["intercept"])
        if v is not None and v[1] < thr:
            if state["maxacc"] is None or v[1] > state["maxacc"]:
                state["maxacc"] = v[1]
            return 1
        if depth >= MAXDEPTH:
            raise RuntimeError(f"depth cap at {A} {C} {X}")
        d = widest(A, C, X)
        if d == 0:
            m = F.fdiv(A[0] + A[1], 2)
            return build((A[0], m), C, X, depth + 1) + build((m, A[1]), C, X, depth + 1)
        if d == 1:
            m = F.fdiv(C[0] + C[1], 2)
            return build(A, (C[0], m), X, depth + 1) + build(A, (m, C[1]), X, depth + 1)
        m = F.fdiv(X[0] + X[1], 2)
        return build(A, C, (X[0], m), depth + 1) + build(A, C, (m, X[1]), depth + 1)

    sys.setrecursionlimit(200000)
    leaves = build(A0, C0, X0, 0)
    return state["nodes"], leaves, state["maxacc"]


if __name__ == "__main__":
    segs, central = load_segments()
    which = sys.argv[1]
    maxn = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    thr = 0          # the certificate's own margin is -2.2e-9; start at 0
    idxs = range(len(segs)) if which == "all" else [int(which)]
    for i in idxs:
        s = segs[i]
        print(f"segment {i}: slope {float(s['slope']):.6g} "
              f"intercept {float(s['intercept']):.9g} "
              f"x in [{float(s['w0']):.6f}, {float(s['w1']):.6f}]", flush=True)
        try:
            n, lv, mx = run_segment(s, thr, maxn)
            print(f"  closed: {n} nodes, {lv} leaves, "
                  f"max accepted {mx / S if mx is not None else None}", flush=True)
        except RuntimeError as e:
            print(f"  FAILED: {e} after {maxn} nodes", flush=True)
