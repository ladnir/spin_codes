"""Generate the dense-tail certificate: a split tree and a logarithm table.

The search is driven by `fixmirror`, which is the same bound Lean checks, so
every leaf this emits is a leaf Lean accepts.  The original mpmath verifier is
not used here — its bound is different (both are sound, neither dominates), and
a partition chosen by one is not guaranteed to satisfy the other.

Usage:  python -B scripts/gen_dense_tail.py [max_nodes] [terms]
"""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixmirror as F
if os.environ.get('NO_MEET'):
    F.USE_MEET = False

S = F.SCALE
TERMS = int(sys.argv[2]) if len(sys.argv) > 2 else 10
THR = -76800000000000000000000          # -7.68e-8, the paper's constant
MAXN = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
MAXDEPTH = 60

# the root box: alpha in [1/125,1], beta in [0,1], omega in [1/500,13/125]
A0 = (S // 125, S)
B0 = (0, S)
W0 = (S // 500, 13 * S // 125)
W_A, W_B, W_W = A0[1] - A0[0], B0[1] - B0[0], W0[1] - W0[0]

requested = set()


def lg(p):
    requested.add(p)
    return F.flogQ(p, S, TERMS)


def gobj_real(u, a):
    G = 1 + 759 * u ** 8 + 2576 * u ** 12 + 759 * u ** 16 + u ** 24
    return math.log(G) / 24 - a * math.log(u)


def pick_u(A):
    """Witness minimising the objective at the box's midpoint alpha."""
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


def widest(A, B, W):
    # the verifier's normalisation: alpha and beta raw, omega scaled by 125/99.
    # Normalising omega by its own initial width instead subdivides it about
    # ten times harder than it needs, which multiplies the whole subtree.
    wa = A[1] - A[0]
    wb = B[1] - B[0]
    ww = (W[1] - W[0]) * 125 // 99
    if wa >= wb and wa >= ww:
        return 0
    return 1 if wb >= ww else 2


nodes = 0
leaves = 0
deepest = 0
maxacc = None
stalled = []


def build(A, B, W, depth):
    """Return a tree node, or None if this subtree could not be closed."""
    global nodes, leaves, deepest, maxacc
    nodes += 1
    deepest = max(deepest, depth)
    if nodes % 2000 == 0:
        print(f"  nodes={nodes} leaves={leaves} deepest={deepest} "
              f"logs={len(requested)}", flush=True)
    if nodes > MAXN:
        raise RuntimeError(f"node cap {MAXN} exceeded")
    if F.infeasible(A, B, W):
        leaves += 1
        return ("leaf", 1)
    u = pick_u(A)
    if 0 < u:
        v = F.leafValue(lg, u, A, B, W)
        if v is not None and v[1] < THR:
            leaves += 1
            if maxacc is None or v[1] > maxacc:
                maxacc = v[1]
            return ("leaf", u)
    if depth >= MAXDEPTH:
        stalled.append((A, B, W, depth))
        return None
    d = widest(A, B, W)
    if d == 0:
        m = F.fdiv(A[0] + A[1], 2)
        l = build((A[0], m), B, W, depth + 1)
        r = build((m, A[1]), B, W, depth + 1)
    elif d == 1:
        m = F.fdiv(B[0] + B[1], 2)
        l = build(A, (B[0], m), W, depth + 1)
        r = build(A, (m, B[1]), W, depth + 1)
    else:
        m = F.fdiv(W[0] + W[1], 2)
        l = build(A, B, (W[0], m), depth + 1)
        r = build(A, B, (m, W[1]), depth + 1)
    if l is None or r is None:
        return None
    return ("split", d, m, l, r)


if __name__ == "__main__":
    sys.setrecursionlimit(100000)
    try:
        tree = build(A0, B0, W0, 0)
    except RuntimeError as e:
        print(json.dumps({"status": "node-cap", "error": str(e), "nodes": nodes,
                          "leaves": leaves, "deepest": deepest}))
        raise SystemExit(1)
    print(json.dumps({
        "status": "closed" if tree is not None else "stalled",
        "nodes": nodes, "leaves": leaves, "deepest": deepest,
        "distinct_logs": len(requested),
        "max_accepted_upper": (maxacc / S) if maxacc is not None else None,
        "threshold": THR / S,
        "stalled_boxes": len(stalled),
    }))
    if tree is None:
        raise SystemExit(1)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dense_tail_tree.json")
    with open(out, "w") as f:
        json.dump({"terms": TERMS, "thr": str(THR),
                   "root": [[str(x) for x in A0], [str(x) for x in B0],
                            [str(x) for x in W0]],
                   "logs": sorted(str(p) for p in requested)}, f)
    print("wrote", out)
