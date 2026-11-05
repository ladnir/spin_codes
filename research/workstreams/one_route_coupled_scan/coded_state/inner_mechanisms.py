"""Bounded exact mechanism checks; no native timing or whole-code claim.

All numerics are serial, file-free, and guarded. The field quartet is an
actual distribution, not an assumption that two pair layers are uniform.
"""
import argparse
import itertools
import math
import os
from pathlib import Path
import random
import shutil
import threading
import time
from fractions import Fraction


POLY = 0x10000008D


def mul(a, b, polynomial=POLY):
    degree = polynomial.bit_length()-1
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a >> degree:
            a ^= polynomial
    return out


def remainder(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def trace(a, bits=32, polynomial=POLY):
    out, x = 0, a
    for _ in range(bits):
        out ^= x
        x = mul(x, x, polynomial)
    assert out in (0, 1)
    return out


def power(a,n):
    out = 1
    while n:
        if n & 1:
            out = mul(out,a)
        a = mul(a,a)
        n >>= 1
    return out


def embedded_delta():
    # Embed the existing AES-field element 0x02 in this proof representation.
    generator = next(g for base in range(2,256)
                     if (g := power(base,((1<<32)-1)//255)) != 0
                     and all(power(g,255//p) != 1 for p in (3,5,17)))
    a = 1
    for _ in range(255):
        a = mul(a,generator)
        if power(a,8)^power(a,4)^power(a,3)^a^1 == 0:
            assert a not in (0,1) and power(a,256) == a
            return a
    raise AssertionError('AES-field embedding not found')


def verify_base():
    x = 2
    for i in range(1, 33):
        x = mul(x, x)
        if i == 16:
            assert gcd(x ^ 2, POLY) == 1
    assert x == 2
    alpha = next(1 << i for i in range(32) if trace(1 << i) == 1)
    return alpha


def add(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def quadratic(a, b, mu, base_mul):
    # Direct reference: u^2=u+mu. Do not use the optimized circuit here.
    p = base_mul(a[1], b[1])
    return (base_mul(a[0], b[0]) ^ base_mul(mu, p),
            base_mul(a[0], b[1]) ^ base_mul(a[1], b[0]) ^ p)


def field_ops(alpha, base_mul=mul):
    def m2(a, b):
        return quadratic(a, b, alpha, base_mul)

    def m4(a, b):
        a0, a1, b0, b1 = a[:2], a[2:], b[:2], b[2:]
        p = m2(a1, b1)
        return add(m2(a0, b0), m2((0, alpha), p)) + add(add(m2(a0, b1), m2(a1, b0)), p)

    def prep2(c):
        return c[0], base_mul(alpha, c[1]), c[0] ^ c[1]

    def prep4(c):
        return prep2(c[:2]), prep2(m2((0, alpha), c[2:])), prep2(add(c[:2], c[2:]))

    def run2(c, x):
        p0 = base_mul(c[0], x[0])
        p1 = base_mul(c[1], x[1])
        p2 = base_mul(c[2], x[0] ^ x[1])
        return p0 ^ p1, p0 ^ p2

    def run4(c, x):
        p0 = run2(c[0], x[:2])
        p1 = run2(c[1], x[2:])
        p2 = run2(c[2], add(x[:2], x[2:]))
        return add(p0, p1) + add(p0, p2)

    return m2, m4, prep4, run4


def determinant(matrix, product=mul):
    # Characteristic two: signs disappear. Independent of inversion/rank code.
    n = len(matrix)
    out = 0
    for p in itertools.permutations(range(n)):
        term = 1
        for i in range(n):
            term = product(term, matrix[i][p[i]])
        out ^= term
    return out


def matrix_for(c, multiply):
    columns = [multiply(c, tuple(int(i == j) for i in range(4))) for j in range(4)]
    return [[columns[j][i] for j in range(4)] for i in range(4)]


def minors(matrix):
    return [(rows, cols, determinant([[matrix[i][j] for j in cols] for i in rows]))
            for size in range(1, 5)
            for rows in itertools.combinations(range(4), size)
            for cols in itertools.combinations(range(4), size)]


def shifted(matrix, delta=2):
    return [[v ^ (delta if i == j else 0) for j, v in enumerate(row)]
            for i, row in enumerate(matrix)]


def apply_matrix(matrix, vector):
    out = []
    for row in matrix:
        v = 0
        for c, x in zip(row, vector):
            v ^= mul(c, x)
        out.append(v)
    return tuple(out)


def sweep(q, forward, backward=()):
    r = list(q)
    for i, a in enumerate(forward, 1):
        r[i] ^= mul(a, r[i-1])
    for i in reversed(range(len(backward))):
        r[i] ^= mul(backward[i], r[i+1])
    return tuple(r)


def inverse_sweep(v, forward, backward=()):
    r = list(v)
    # Undo backward sweep first; use its unmodified input on each edge.
    for i in range(len(backward)):
        r[i] ^= mul(backward[i], r[i+1])
    for i in reversed(range(len(forward))):
        r[i+1] ^= mul(forward[i], r[i])
    return tuple(r)


def structural_checks():
    alpha = verify_base()
    delta = embedded_delta()
    m2, m4, prepare, fast = field_ops(alpha)
    rng = random.Random(20261026)
    checks = 0
    for _ in range(128):
        c = tuple(rng.getrandbits(32) for _ in range(4))
        x = tuple(rng.getrandbits(32) for _ in range(4))
        assert fast(prepare(c), x) == m4(c, x)
        assert m4(c, x) == m4(x, c)
        checks += 2
    # One actual-field witness proves every minor polynomial is nonzero.
    c = (0xC423F3E3,0xB01BA9E8,0x5AA87C13,0x915276C9)
    a = matrix_for(c, m4)
    values = minors(a) + minors(shifted(a,delta))
    assert all(value for _, _, value in values)
    assert len(values) == 138
    print(dict(field32_polynomial=hex(POLY), alpha=hex(alpha), delta=hex(delta), quartet_key=list(map(hex,c)),
               all_138_minors_nonzero=True, frozen_witness=True), flush=True)
    checks += len(values)
    calls = [0]
    def counted(a,b):
        calls[0] += 1
        return mul(a,b)
    _, _, prep_counted, run_counted = field_ops(alpha,counted)
    key = prep_counted(c)
    calls[0] = 0
    assert run_counted(key,(1,2,3,4)) == m4(c,(1,2,3,4))
    assert calls[0] == 9
    # Negative controls: accepting arbitrary structured maps would be unsound.
    assert not all(v for _,_,v in minors(matrix_for((1,0,0,0),m4)))
    assert determinant(shifted(matrix_for((delta,0,0,0),m4),delta)) == 0
    checks += 3
    # Exhaustive fixed-input law in GF4^4, independent small field fixture.
    bm = lambda x,y: mul(x,y,7)
    assert trace(2,2,7) == 1
    _, small, _, _ = field_ops(2,bm)
    points = list(itertools.product(range(4), repeat=4))
    for v in points[1:]:
        assert len({small(c,v) for c in points}) == 256
        checks += 1
    # Boundary and sparse inverse obstructions are identities, not codewords.
    for _ in range(32):
        f = [rng.getrandbits(32) or 1 for _ in range(3)]
        b = [rng.getrandbits(32) or 1 for _ in range(3)]
        for j in range(4):
            e = tuple(int(i==j) for i in range(4))
            v = inverse_sweep(e,f,b)
            assert sweep(v,f,b) == e
            assert sum(x != 0 for x in v) <= (2 if j in (0,3) else 3)
            checks += 2
        assert sweep((0,0,0,1),f,b)[3] == 1
        checks += 1
    # XOR butterfly leaves the all-equal vector on a single coordinate.
    q = [1,1,1,1]
    for bit in (1,2):
        for i in range(4):
            if i & bit:
                q[i] ^= q[i ^ bit]
    assert sum(v != 0 for v in q) == 1
    # Rank-one refresh I+u v^T fixes ker(v^T), also when v^T u=0.
    u, v, x = (1,1,0,0), (1,1,1,1), (0,0,1,1)
    assert sum(v[i]*u[i] for i in range(4)) % 2 == 0
    assert sum(v[i]*x[i] for i in range(4)) % 2 == 0
    assert add(x, tuple(ui * (sum(v[i]*x[i] for i in range(4)) % 2) for ui in u)) == x
    checks += 3
    # Raw scalar is uniform on GF(q^4), not restricted to the base field.
    rejection_degree = 2*sum(j*math.comb(4,j)**2 for j in range(1,5))
    assert rejection_degree == 280
    inflation = 1/(1-Fraction((1<<21)//4*rejection_degree,1<<32))
    assert inflation == Fraction(1024,989) and inflation < 2
    print(dict(exact_checks=checks, density_bound=str(inflation),
               variable_base_products_per_quartet=9, runtime_base_xors_per_quartet=15,
               scope='local algebra and distribution transfer; not whole-code'),flush=True)


def concentrated_refresh_checks():
    # Exact raw rank-one law over GF4 in dimension three.
    bm = lambda x,y: mul(x,y,7)
    points = list(itertools.product(range(4),repeat=3))
    a = (1,0,0)
    histogram = {x:0 for x in points}
    for u,v in itertools.product(points,repeat=2):
        dot = 0
        for vi,ai in zip(v,a):
            dot ^= bm(vi,ai)
        result = add(a,tuple(bm(ui,dot) for ui in u))
        histogram[result] += 1
    for x,count in histogram.items():
        assert Fraction(count,4**6) == Fraction(1,4)*int(x==a)+Fraction(3,4**4)
    # Fixed state-position support in every batch. All-lazy constraints are
    # homogeneous in the message and number at most rank*(L-1), for ANY keys.
    # Per-row inclusion across regions is independent Bernoulli(w/512).
    for w in (224,232,235,236,238,240,245,250,251,256):
        p = Fraction(w,512)
        mean = 8192*sum(Fraction(math.comb(256,j))*p**j*(1-p)**(256-j)*(j-128)
                        for j in range(129,257))
        ranks = []
        for rank in (1,2,3,4,5):
            cap = rank*4095
            lower = max(Fraction(0),(mean-cap)/(w*4096//2-cap))
            ranks.append(float(lower))
            if w == 240 and rank == 1:
                assert lower > Fraction(27,10000)
            if w == 251 and rank == 4:
                assert lower > Fraction(15,10000)
        print(dict(state_support=w, output_fraction=w/512,
                   expected_shortened_dimension=float(mean),
                   rank_one_exception_constraints=4095,
                   rigorous_failure_lower_by_rank_1_to_5=ranks),flush=True)
    print(dict(rank_one_exact_law_checks=64, fixed_support_probe=True,
               caveat='zero lower bound is inconclusive; positive lower bound is an obstruction'),flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay',action='store_true')
    args = parser.parse_args()
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free >= 20*(1<<30)
    timer = threading.Timer(60,lambda:os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        structural_checks()
        concentrated_refresh_checks()
        if args.replay:
            from diagnose_pairing_sufficiency import rational_quartet_replay
            rational_quartet_replay()
        print(dict(seconds=time.monotonic()-started, whole_code=False),flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    main()
