"""Emit an interval cover as a set of self-contained Lean modules.

The monolithic form does not scale: one 23,541-entry table elaborates to a
246 MB `.olean`, and one `decide` over 13,513 leaves reached 16 GB before being
killed.  So the tree is cut into parts of bounded leaf count, and each part
gets its own module holding *its own* table, its own subtree, and its own
finished claim.

Splitting at the level of the claim rather than the check is what makes the
tables independent: `claim_split_fst/snd/thd` combine conclusions, so two parts
never have to share an oracle.

Usage:  python -B scripts/emit_cover.py [--small] [limit] [terms]
"""
import sys, os, json, math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixmirror as F

S = F.SCALE
SMALL = "--small" in sys.argv
argv = [a for a in sys.argv[1:] if not a.startswith("--")]
LIMIT = int(argv[0]) if argv else 900          # max leaves per part
TERMS = int(argv[1]) if len(argv) > 1 else 10
THR = -76800000000000000000000
MAXDEPTH = 60

if SMALL:
    A0 = (10 * S // 1000, 200 * S // 1000)
    B0 = (100 * S // 1000, 300 * S // 1000)
    W0 = (50 * S // 1000, 104 * S // 1000)
    TAG = "Small"
else:
    A0 = (S // 125, S)
    B0 = (0, S)
    W0 = (S // 500, 13 * S // 125)
    TAG = ""

plain_lg = F.directLog(TERMS)


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


def widest(A, B, W):
    wa = A[1] - A[0]
    wb = B[1] - B[0]
    ww = (W[1] - W[0]) * 125 // 99
    if wa >= wb and wa >= ww:
        return 0
    return 1 if wb >= ww else 2


nodes = 0


def build(A, B, W, depth):
    global nodes
    nodes += 1
    if nodes % 4000 == 0:
        print(f"  search nodes={nodes}", flush=True)
    if F.infeasible(A, B, W):
        return ("leaf", 1)
    u = pick_u(A)
    v = F.leafValue(plain_lg, u, A, B, W)
    if v is not None and v[1] < THR:
        return ("leaf", u)
    if depth >= MAXDEPTH:
        raise RuntimeError(f"depth cap at {A} {B} {W}")
    d = widest(A, B, W)
    if d == 0:
        m = F.fdiv(A[0] + A[1], 2)
        return ("split", 0, m, build((A[0], m), B, W, depth + 1),
                build((m, A[1]), B, W, depth + 1))
    if d == 1:
        m = F.fdiv(B[0] + B[1], 2)
        return ("split", 1, m, build(A, (B[0], m), W, depth + 1),
                build(A, (m, B[1]), W, depth + 1))
    m = F.fdiv(W[0] + W[1], 2)
    return ("split", 2, m, build(A, B, (W[0], m), depth + 1),
            build(A, B, (m, W[1]), depth + 1))


def leaves_of(t):
    if t[0] == "leaf":
        return 1
    return leaves_of(t[3]) + leaves_of(t[4])


parts = []          # (tree, A, B, W)


def cut(t, A, B, W):
    """Replace every subtree of at most LIMIT leaves by a part reference."""
    if leaves_of(t) <= LIMIT:
        parts.append((t, A, B, W))
        return ("part", len(parts) - 1)
    _, d, m, l, r = t
    if d == 0:
        return ("split", 0, m, cut(l, (A[0], m), B, W), cut(r, (m, A[1]), B, W))
    if d == 1:
        return ("split", 1, m, cut(l, A, (B[0], m), W), cut(r, A, (m, B[1]), W))
    return ("split", 2, m, cut(l, A, B, (W[0], m)), cut(r, A, B, (m, W[1])))


def collect(t, A, B, W, out):
    if t[0] == "leaf":
        def rec(p):
            out.add(p)
            return F.flogQ(p, S, TERMS)
        F.leafOK(rec, THR, t[1], A, B, W)
        return
    _, d, m, l, r = t
    if d == 0:
        collect(l, (A[0], m), B, W, out); collect(r, (m, A[1]), B, W, out)
    elif d == 1:
        collect(l, A, (B[0], m), W, out); collect(r, A, (m, B[1]), W, out)
    else:
        collect(l, A, B, (W[0], m), out); collect(r, A, B, (m, W[1]), out)


def emit_tree(t, buf):
    if t[0] == "leaf":
        buf.append(f".leaf {t[1]}")
        return
    buf.append(f".split {t[1]} {t[2]} (")
    emit_tree(t[3], buf)
    buf.append(") (")
    emit_tree(t[4], buf)
    buf.append(")")


def emit_logs(keys, vals, lo, hi, buf):
    if lo > hi:
        buf.append(".leaf")
        return
    mid = (lo + hi) // 2
    k = keys[mid]
    v = vals[k]
    buf.append(f".node {k} ⟨{v[0]}, {v[1]}⟩ (")
    emit_logs(keys, vals, lo, mid - 1, buf)
    buf.append(") (")
    emit_logs(keys, vals, mid + 1, hi, buf)
    buf.append(")")


def fix(x):
    return f"⟨{x[0]}, {x[1]}⟩"


def emit_claim(node, A, B, W):
    """The boxes must be passed explicitly: `setHi ?A m` does not unify with a
    part's own box literal, so leaving them implicit leaves metavariables."""
    if node[0] == "part":
        return f"part{TAG}{node[1]}_claim"
    _, d, m, l, r = node
    fn = ["claim_split_fst", "claim_split_snd", "claim_split_thd"][d]
    if d == 0:
        la, lb, lw = (A[0], m), B, W
        ra, rb, rw = (m, A[1]), B, W
    elif d == 1:
        la, lb, lw = A, (B[0], m), W
        ra, rb, rw = A, (m, B[1]), W
    else:
        la, lb, lw = A, B, (W[0], m)
        ra, rb, rw = A, B, (m, W[1])
    return (f"({fn} (A := {fix(A)}) (B := {fix(B)}) (W := {fix(W)}) (m := {m})"
            f" ({emit_claim(l, la, lb, lw)}) ({emit_claim(r, ra, rb, rw)}))")


HEADER = ("set_option maxRecDepth 8000000" + chr(10)
          + "set_option maxHeartbeats 0" + chr(10) + chr(10))

if __name__ == "__main__":
    sys.setrecursionlimit(200000)
    tree = build(A0, B0, W0, 0)
    total_leaves = leaves_of(tree)
    print(f"search closed: {nodes} nodes, {total_leaves} leaves")
    skeleton = cut(tree, A0, B0, W0)
    print(f"parts: {len(parts)}  (limit {LIMIT} leaves each)")

    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    cov = os.path.join(root, "SpinCodes", "Cover")
    os.makedirs(cov, exist_ok=True)

    total_logs = 0
    for i, (t, A, B, W) in enumerate(parts):
        args = set()
        collect(t, A, B, W, args)
        keys = sorted(args)
        vals = {k: F.flogQ(k, S, TERMS) for k in keys}
        assert all(v is not None for v in vals.values())
        total_logs += len(keys)
        lb, tb = [], []
        emit_logs(keys, vals, 0, len(keys) - 1, lb)
        emit_tree(t, tb)
        name = f"Part{TAG}{i}"
        with open(os.path.join(cov, f"{name}.lean"), "w", encoding="utf-8") as f:
            f.write(
                f"/- Generated by `scripts/emit_cover.py`.  Do not edit.\n\n"
                f"Part {i} of the cover: {leaves_of(t)} leaves, {len(keys)} logarithms.\n"
                f"Self-contained — its own table, its own subtree, its own claim. -/\n"
                "import SpinCodes.Cover.Params\n"
                "import SpinCodes.Numeric.BAEval\n\n"
                + HEADER +
                "namespace Spin.Cover\n\nopen Spin.Numeric Spin.Numeric.Fix\n\n"
                f"def table{TAG}{i} : LogTree :=\n  " + "".join(lb) + "\n\n"
                f"def tree{TAG}{i} : BoxTree :=\n  " + "".join(tb) + "\n\n"
                f"def a{TAG}{i} : Fix := {fix(A)}\n"
                f"def b{TAG}{i} : Fix := {fix(B)}\n"
                f"def w{TAG}{i} : Fix := {fix(W)}\n\n"
                f"theorem logs{TAG}{i}_ok : (table{TAG}{i}).check terms = true := by decide\n\n"
                f"theorem tree{TAG}{i}_ok :\n"
                f"    checkTree (leafOK (treeLog table{TAG}{i}) thr) tree{TAG}{i}"
                f" a{TAG}{i} b{TAG}{i} w{TAG}{i} = true := by decide\n\n"
                f"theorem part{TAG}{i}_claim (x y z : ℝ) (hx : Mem a{TAG}{i} x)\n"
                f"    (hy : Mem b{TAG}{i} y) (hz : Mem w{TAG}{i} z) :\n"
                f"    DenseClaim thr x y z :=\n"
                f"  denseTail_of_checkTree (treeLog_oracle logs{TAG}{i}_ok)"
                f" tree{TAG}{i} tree{TAG}{i}_ok hx hy hz\n\n"
                "end Spin.Cover\n")

    imports = "\n".join(f"import SpinCodes.Cover.Part{TAG}{i}" for i in range(len(parts)))
    with open(os.path.join(cov, f"DenseTail{TAG}.lean"), "w", encoding="utf-8") as f:
        f.write(
            "/- Generated by `scripts/emit_cover.py`.  Do not edit.\n\n"
            f"Assembles {len(parts)} independently checked parts into the cover of\n"
            f"`eq:structured-ba-dense-tail` over `α ∈ [{A0[0]}/scale, {A0[1]}/scale]` etc.\n"
            "Each `claim_split` is one application of `le_or_gt` at a cut. -/\n"
            + imports + "\n\n" + HEADER +
            "namespace Spin.Cover\n\nopen Spin.Numeric Spin.Numeric.Fix\n\n"
            f"def rootA{TAG} : Fix := {fix(A0)}\n"
            f"def rootB{TAG} : Fix := {fix(B0)}\n"
            f"def rootW{TAG} : Fix := {fix(W0)}\n\n"
            f"/-- **The dense-tail gap**, assembled from {len(parts)} parts. -/\n"
            f"theorem denseTail{TAG} (x y z : ℝ) (hx : Mem rootA{TAG} x)\n"
            f"    (hy : Mem rootB{TAG} y) (hz : Mem rootW{TAG} z) :\n"
            f"    DenseClaim thr x y z :=\n"
            f"  {emit_claim(skeleton, A0, B0, W0)} x y z hx hy hz\n\n"
            "end Spin.Cover\n")

    with open(os.path.join(here, f"cover_meta{TAG}.json"), "w") as f:
        json.dump({"nodes": nodes, "leaves": total_leaves, "parts": len(parts),
                   "logs_total_with_duplication": total_logs, "terms": TERMS,
                   "limit": LIMIT}, f)
    print(f"wrote {len(parts)} part modules + DenseTail{TAG}.lean; "
          f"{total_logs} logarithms across parts")
