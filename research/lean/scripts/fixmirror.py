"""A Python mirror of the Lean `Fix` layer.

The branch-and-bound that generates the dense-tail certificate has to be
driven by *exactly* the bound Lean will check, not by the mpmath bound of the
original verifier.  Otherwise a box accepted here could fail there, and the
failure would only surface at the end of a twenty-minute kernel run.

So this is a line-by-line port of `SpinCodes/Numeric/FixedDefs.lean` and
`SpinCodes/Numeric/BAEvalDefs.lean`.  It is checked against Lean by
`scripts/check_mirror.py`, which compares integer numerators exactly — the
port is only useful if it agrees to the last digit.

Nothing here is trusted: it chooses the partition and the witnesses.  Lean
re-derives every bound from the emitted data.
"""

PREC = 30
SCALE = 10 ** PREC

# ---------------------------------------------------------------- arithmetic

def fdiv(n, d):
    # Int.ediv for positive d is floor division, which is Python's //
    return n // d


def cdiv(n, d):
    return -((-n) // d)


def sc(n):
    return (n, n)


def ofInt(n):
    return (n * SCALE, n * SCALE)


def ofFrac(p, q):
    return (fdiv(p * SCALE, q), cdiv(p * SCALE, q))


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def neg(a):
    return (-a[1], -a[0])


def sub(a, b):
    return (a[0] - b[1], a[1] - b[0])


def mul(a, b):
    p1 = a[0] * b[0]
    p2 = a[0] * b[1]
    p3 = a[1] * b[0]
    p4 = a[1] * b[1]
    lo = min(min(p1, p2), min(p3, p4))
    hi = max(max(p1, p2), max(p3, p4))
    return (fdiv(lo, SCALE), cdiv(hi, SCALE))


def divInt(a, k):
    return (fdiv(a[0], k), cdiv(a[1], k))


def powf(a, n):
    r = ofInt(1)
    for _ in range(n):
        r = mul(r, a)
    return r


def pm(r):
    return (-r[1], r[1])


ZERO = (0, 0)

# ------------------------------------------------------------------ the log

def oddGo(z2, m, i):
    if m == 0:
        return ZERO
    return add(ofFrac(1, i), mul(z2, oddGo(z2, m - 1, i + 2)))


def oddSeries(m, z):
    return mul(z, oddGo(mul(z, z), m, 1))


def remOdd(c, m):
    return mul(ofInt(3), powf(ofFrac(1, c), 2 * m + 1))


LOG2 = add(mul(ofInt(2), oddSeries(21, ofFrac(1, 3))), pm(remOdd(3, 21)))


def flogA(P, Q, m):
    return add(mul(ofInt(2), oddSeries(m, ofFrac(P - Q, P + Q))), pm(remOdd(5, m)))


def ipow2(k):
    return 2 ** k


def shiftGo(num, den, fuel, k):
    while fuel > 0:
        if 2 * den <= 3 * (num * ipow2(k)):
            return k
        k += 1
        fuel -= 1
    return k


def shiftDownGo(p, q, fuel, j):
    while fuel > 0:
        if 2 * p <= 3 * (q * ipow2(j)):
            return j
        j += 1
        fuel -= 1
    return j


def shiftPair(p, q):
    if 2 * q <= 3 * p:
        return (0, shiftDownGo(p, q, 128, 0))
    return (shiftGo(p, q, 128, 0), 0)


def flogKJ(p, q, k, j, n):
    P = p * ipow2(k)
    Q = q * ipow2(j)
    if not (0 < P and 0 < Q and 2 * Q <= 3 * P and 2 * P <= 3 * Q):
        return None
    return add(sub(flogA(P, Q, n), mul(ofInt(k), LOG2)), mul(ofInt(j), LOG2))


def flogQ(p, q, n):
    k, j = shiftPair(p, q)
    return flogKJ(p, q, k, j, n)


def directLog(n):
    return lambda p: flogQ(p, SCALE, n)


def flogIWL(lg, a):
    L = lg(a[0])
    if L is None:
        return None
    H = lg(a[1])
    if H is None:
        return None
    return (L[0], H[1])

# ---------------------------------------------------------------- x log x

ONE = ofInt(1)


def xlPointL(lg, p):
    if p == 0:
        return ZERO
    l = lg(p)
    if l is None:
        return None
    return mul(sc(p), l)


EINV = 367879441171442321595523770161


def tangentPt(a):
    return max(max(a[0], min(EINV, a[1])), 1)


def fxlogxL(lg, a):
    L = xlPointL(lg, a[0])
    if L is None:
        return None
    H = xlPointL(lg, a[1])
    if H is None:
        return None
    if a[1] == 0:
        return ZERO
    m = tangentPt(a)
    lm = lg(m)
    if lm is None:
        return None
    lo = sub(mul(a, add(lm, ONE)), sc(m))[0]
    return (lo, max(L[1], H[1]))


def fpiEvalL(lg, a, c):
    ts = []
    for x in (c, sub(ONE, c), sub(c, divInt(a, 2)),
              sub(sub(ONE, c), divInt(a, 2)), sub(ONE, a)):
        t = fxlogxL(lg, x)
        if t is None:
            return None
        ts.append(t)
    t1, t2, t3, t4, t5 = ts
    return add(sub(sub(add(add(mul(a, LOG2), t1), t2), t3), t4), t5)

# ------------------------------------------------------------- the centered bound

def midPt(A):
    return sc(fdiv(A[0] + A[1], 2))


def fDpaL(lg, A, C):
    l1 = flogIWL(lg, sub(C, divInt(A, 2)))
    l2 = flogIWL(lg, sub(sub(ONE, C), divInt(A, 2)))
    l3 = flogIWL(lg, sub(ONE, A))
    if l1 is None or l2 is None or l3 is None:
        return None
    return sub(add(LOG2, divInt(add(l1, l2), 2)), l3)


def fDpcL(lg, A, C):
    l1 = flogIWL(lg, C)
    l2 = flogIWL(lg, sub(ONE, C))
    l3 = flogIWL(lg, sub(C, divInt(A, 2)))
    l4 = flogIWL(lg, sub(sub(ONE, C), divInt(A, 2)))
    if l1 is None or l2 is None or l3 is None or l4 is None:
        return None
    return add(sub(sub(l1, l2), l3), l4)


def fpiCenteredL(lg, A, C):
    p0 = fpiEvalL(lg, midPt(A), midPt(C))
    if p0 is None:
        return None
    da = fDpaL(lg, A, C)
    if da is None:
        return None
    dc = fDpcL(lg, A, C)
    if dc is None:
        return None
    return add(add(p0, mul(da, sub(A, midPt(A)))), mul(dc, sub(C, midPt(C))))


def feasOK(A, C):
    return (0 < 2 * C[0] - A[1] and 0 < 2 * SCALE - 2 * C[1] - A[1]
            and A[1] < SCALE and 0 < C[0] and C[1] < SCALE)


def fpiCenteredG(lg, A, C):
    return fpiCenteredL(lg, A, C) if feasOK(A, C) else None


def meet(a, b):
    return (max(a[0], b[0]), min(a[1], b[1]))


def fpiBestL(lg, A, C):
    p = fpiEvalL(lg, A, C)
    q = fpiCenteredG(lg, A, C)
    if p is not None and q is not None:
        return meet(p, q)
    return p if p is not None else q

def clampLo(X):
    return (max(0, X[0]), X[1])


def fpiEvalClampL(lg, A, C):
    xs = [C, sub(ONE, C), clampLo(sub(C, divInt(A, 2))),
          clampLo(sub(sub(ONE, C), divInt(A, 2))), sub(ONE, A)]
    ts = []
    for x in xs:
        t = fxlogxL(lg, x)
        if t is None:
            return None
        ts.append(t)
    t1, t2, t3, t4, t5 = ts
    return add(sub(sub(add(add(mul(A, LOG2), t1), t2), t3), t4), t5)


USE_MEET = True   # measurement toggle only; Lean always meets


def fpiBestClampL(lg, A, C):
    q = fpiCenteredG(lg, A, C)
    if q is not None and not USE_MEET:
        return q
    p = fpiEvalClampL(lg, A, C)
    if p is not None and q is not None:
        return meet(p, q)
    return p if p is not None else q


# ------------------------------------------------------------- the Golay part

def fGolayG(u):
    u2 = mul(u, u)
    v = mul(u2, u2)
    v2 = mul(v, v)
    return add(ONE, mul(v2, add(ofInt(759),
        mul(v, add(ofInt(2576), mul(v, add(ofInt(759), v2)))))))


def fgObjL(lg, u, a):
    lgG = flogIWL(lg, fGolayG(u))
    if lgG is None:
        return None
    lu = flogIWL(lg, u)
    if lu is None:
        return None
    return sub(divInt(lgG, 24), mul(a, lu))

# ------------------------------------------------------------------ the leaf

def infeasible(A, B, W):
    return (2 * B[1] < A[0] or 2 * SCALE - A[0] < 2 * B[0]
            or 2 * W[1] < B[0] or 2 * SCALE - B[0] < 2 * W[0])


def leafValue(lg, u, A, B, W):
    """The objective enclosure a leaf is judged on, or None."""
    g = fgObjL(lg, sc(u), A)
    if g is None:
        return None
    p1 = fpiBestClampL(lg, A, B)
    if p1 is None:
        return None
    p2 = fpiBestClampL(lg, B, W)
    if p2 is None:
        return None
    return add(add(g, p1), p2)


def leafOK(lg, thr, u, A, B, W):
    if infeasible(A, B, W):
        return True
    if not (0 < u):
        return False
    v = leafValue(lg, u, A, B, W)
    return v is not None and v[1] < thr
