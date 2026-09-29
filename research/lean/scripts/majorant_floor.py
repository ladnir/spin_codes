"""Does the `Fix` bound resolve the majorant margin at all?  (T7c)

The certificate's largest accepted residual is `-2.2e-9`, and it works at
`interval_dps = 70`.  The Lean `Fix` layer carries 30 decimal digits.  If the
overestimate of the residual bound plateaus above `-2.2e-9` as boxes shrink,
then no split heuristic can close the cover and the blocker is *precision*,
not search.

This finds the float argmax of the residual on a segment, then shrinks a box
around it and prints the `Fix` upper bound at each scale, against the true
value.

Usage:  python -B scripts/majorant_floor.py <segment>
"""
import sys, os, math
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixmirror as F
import majorant as M

S = F.SCALE


def xlogx(t):
    return 0.0 if t <= 0 else t * math.log(t)


def pi_real(a, c):
    """The paper's π(a,c), via the x log x form."""
    return (xlogx(c - a / 2) + xlogx(1 - c - a / 2) + a * math.log(2.0)
            - xlogx(1 - a) - 2 * xlogx(a / 2) - 0.0) if False else None


def residual_real(a, c, x, slope, intercept):
    """Float residual g(a) + π(a,c) + π(c,x) - (s x + b), via the Fix layer at
    a degenerate (point) box, which is the same formula Lean uses."""
    A = (int(a * S), int(a * S))
    C = (int(c * S), int(c * S))
    X = (int(x * S), int(x * S))
    lg = F.directLog(M.TERMS)
    v = M.leaf_value(lg, M.pick_u(A), A, C, X, slope, intercept)
    if v is None:
        return None
    return (v[0] + v[1]) / 2 / S


def main():
    segs, _ = M.load_segments()
    i = int(sys.argv[1])
    seg = segs[i]
    slope, intercept = seg["slope"], seg["intercept"]
    x0, x1 = float(seg["w0"]), float(seg["w1"])
    print(f"segment {i}: slope {float(slope):.6g}  x in [{x0:.6f}, {x1:.6f}]")

    # coarse float search for the argmax of the residual
    best = None
    N = 26
    for ia in range(1, N):
        a = ia / N
        for ic in range(1, N):
            c = a / 2 + ic / N * (1 - a - 1e-9)
            if not (a / 2 < c < 1 - a / 2):
                continue
            for ix in range(N + 1):
                x = x0 + (x1 - x0) * ix / N
                if not (c / 2 < x < 1 - c / 2):
                    continue
                r = residual_real(a, c, x, slope, intercept)
                if r is None:
                    continue
                if best is None or r > best[0]:
                    best = (r, a, c, x)
    if best is None:
        print("  no feasible point found")
        return
    r, a, c, x = best
    print(f"  coarse argmax ~ a={a:.4f} c={c:.4f} x={x:.6f}   residual {r:.6e}")

    # shrink a box around it and watch the Fix upper bound
    lg = F.directLog(M.TERMS)
    print("   half-width        Fix upper bound        width of enclosure")
    for k in range(2, 40, 3):
        h = 2.0 ** (-k)
        A = (max(0, int((a - h) * S)), min(S, int((a + h) * S)))
        C = (max(0, int((c - h) * S)), min(S, int((c + h) * S)))
        X = (max(0, int((x - h) * S)), min(S, int((x + h) * S)))
        v = M.leaf_value(lg, M.pick_u(A), A, C, X, slope, intercept)
        if v is None:
            print(f"   2^-{k:<3d}  {h:.3e}   (unevaluable)")
            continue
        print(f"   2^-{k:<3d}  {h:.3e}   upper = {v[1]/S:+.6e}"
              f"   enclosure width = {(v[1]-v[0])/S:.3e}")


if __name__ == "__main__":
    main()
