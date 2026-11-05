r"""Refine residual terms of the frozen coset return moment.

Fix the actual GF256 row v(a)=(1,a,b(a),a*b(a)) and an inside-zero
set T. The symmetric matrix M_T=Id+A_T^t A_T has rank at least two.
Indeed M00 is0 or1, Mii=1+M0i^2 for i>0, and M03=M12. A nonzero
symmetric rank-one matrix over GF256 is xx^t: write it as lambda*yy^t
and absorb the unique square root of lambda into y. If M00=1, x0=1
forces Mii=M0i^2, a contradiction. If M00=0, its first row vanishes,
while the three other diagonal entries force x1=x2=x3=1. Then M12=1
but M03=0, again a contradiction. Rank zero is excluded by the same
diagonal identities. This argument uses the implemented field and row
structure; it does not change the distribution of state refreshes.

For |I|=ell>=4, put n=16-ell. Let K=I\T have size k and let V be
v requested zeros outside I. The equations are
C_K x_K=M_T u, A_V u=0, where u=C_I x. If k>=4, C_K is onto and
the normalized nonzero-syndrome mass is exactly
q^-(ell-k)*(q^-v-q^-4). If k<4, C_K is injective. The syndrome space
has dimension at most4-v by outside MDS constraints and at most k+2
by rank M_T>=2. The map u to x to Sx is injective, and Sx is supported
on n+k-v remaining output coordinates. Thus its dimension is at most
d=min(4-v,k+2,n+k-v), and its normalized mass is at most q^-ell*(q^d-1).

Sum these positive coefficients with weights
binom(ell,k)binom(n,v)z^(n+k-v)(1-z)^(ell-k+v), 0<=v<=min(3,n).
The old coset expression instead uses dimension4-v for k<4. The rank-M
bound alone improves (k,v)=(0,0),(0,1),(1,0); the output-support bound
can further improve ell>=13. Intersect with the full frozen coset D
cap. Preserve Q, the genuine D1/D2/D3, and every non-E local channel.
Replace raw E, then its one local alpha power; rebuild pair corrections.

This module exports augment_caps and a composition-compatible build_caps.
It has no route, source, count, or global probability calculation. Exact
GF8 RS[8,4] tests check the rank lemma and zero-subset counts, including
one independent enumeration of4096 input vectors. All numerical runs
belong to the parent.180s/256MiB/20GiB; stdout only, no generated files.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations, product
import json
from math import comb
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
COSET_PIN = "8fd5db13b0aca7ebb4149772827985dfba322a5ca02e17ea6fd2931a8199bd06"
Q, LENGTH, DIMENSION = 256, 16, 4
SECONDS, MEMORY = 180, 256*(1 << 20)


def load_coset():
    path = Path(__file__).with_name("cheap_lane_coset_return_refinement.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != COSET_PIN:
        raise RuntimeError("frozen coset-return helper hash mismatch")
    spec = importlib.util.spec_from_file_location("residual_coset_return",path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def positive_moment(ell,z,q=Q,length=LENGTH,rank_refinement=True):
    """Positive rational sum; rank_refinement=False replays q^4 Qell Bn."""
    if not isinstance(z,F) or not 0<z<=1 or not 4<=ell<=length:
        raise ValueError("require rational0<z<=1 and4<=ell<=length")
    n = length-ell
    total = F(0)
    for k in range(ell+1):
        t = ell-k
        for v in range(min(3,n)+1):
            weight = comb(ell,k)*comb(n,v)*z**(n+k-v)*(1-z)**(t+v)
            if k>=4:
                mass = F(1,q**t)*(F(1,q**v)-F(1,q**4))
            else:
                dimension = min(4-v,k+2,n+k-v) if rank_refinement else 4-v
                mass = F(q**dimension-1,q**ell)
            assert mass>=0
            total += weight*mass
    return total


def augment_caps(caps,z,alpha,guard):
    """Refine a baseline or coset map, without mutating either input map."""
    coset = load_coset()
    if "coset_receipt" not in caps:
        caps = coset.augment_caps(caps,z,alpha,guard)
    if not isinstance(alpha,F) or not 0<alpha<=1:
        raise ValueError("require rational0<alpha<=1")
    v,core,small = caps["verifier"],caps["core"],caps["small"]
    if core.Q!=Q or core.Z!=z or core.ALPHA!=alpha:
        raise ValueError("caps do not match the actual field and witnesses")
    quiet = caps["quiet_exact"]
    previous = caps["nonquiet_exact"]
    refined = list(previous)
    rows = []
    for ell in range(4,17):
        guard()
        candidate = positive_moment(ell,z)
        unrefined = positive_moment(ell,z,rank_refinement=False)
        expected = Q**4*quiet[ell]*coset.outside_nonzero(16-ell,z)
        assert unrefined==expected and 0<=candidate<=unrefined
        refined[ell] = min(previous[ell],candidate)
        rows.append({"ell":ell,"old_coset_D":str(previous[ell]),
                     "positive_rank_candidate":str(candidate),"new_D":str(refined[ell])})
    assert refined[:4]==list(previous[:4])
    raw,powered = {},{}
    exponent = v.exact(alpha)
    changed = 0
    largest_ratio = F(1)
    for occ,row in caps["raw"].items():
        guard()
        before = small.return_cap(occ,quiet,previous,core.KAPPA)
        after = small.return_cap(occ,quiet,refined,core.KAPPA)
        assert 0<=after<=before
        point = v.fmpq(before.numerator,before.denominator)
        assert row[2].lower().fmpq()<=point<=row[2].upper().fmpq()
        new_row,new_powered = list(row),list(caps["powered"][occ])
        new_row[2] = v.exact(after)
        new_powered[2] = v.positive(new_row[2])**exponent
        raw[occ],powered[occ] = tuple(new_row),tuple(new_powered)
        for channel in (0,1,3,4,5):
            assert raw[occ][channel] is row[channel]
            assert powered[occ][channel] is caps["powered"][occ][channel]
        if after<before:
            changed += 1
            assert after>0
            largest_ratio = max(largest_ratio,before/after)
    assert raw.keys()==powered.keys()==caps["raw"].keys()==caps["powered"].keys()
    receipt = {"q":Q,"rank_M_T_lower":2,"rank_M_only_terms":[[0,0],[0,1],[1,0]],
               "residual_dimension_upper":"min(4-v,k+2,n+k-v)",
               "ell1_2_3_unchanged":True,"quiet_unchanged":True,
               "changed_channel":"raw E and its single local alpha power only",
               "requires_rebuilt_pair_corrections":True,"changed_shapes":changed,
               "shape_count":len(raw),"maximum_raw_E_ratio":str(largest_ratio),"blocks":rows}
    result = dict(caps)
    report = dict(caps["block_report"])
    report["residual_rank_return"] = receipt
    result.update({"raw":raw,"powered":powered,"nonquiet_exact":tuple(refined),
                   "block_report":report,"tests":dict(caps["tests"]),
                   "residual_rank_receipt":receipt})
    return result


def build_caps(z,alpha,triple,anchor,guard,deadline):
    """Drop-in composition builder; actual ell1/2/3 and coset caps retained."""
    coset = load_coset()
    caps = coset.build_caps(z,alpha,triple,anchor,guard,deadline)
    return augment_caps(caps,z,alpha,guard)


def gf8_mul(a,b):
    value = 0
    while b:
        if b&1:
            value ^= a
        b >>= 1
        a <<= 1
        if a&8:
            a ^= 11
    return value


def gf8_dot(a,b):
    result = 0
    for x,y in zip(a,b):
        result ^= gf8_mul(x,y)
    return result


def gf8_rank(rows,columns):
    matrix = [list(row) for row in rows]
    rank = 0
    for j in range(columns):
        pivot = next((i for i in range(rank,len(matrix)) if matrix[i][j]),None)
        if pivot is None:
            continue
        matrix[rank],matrix[pivot] = matrix[pivot],matrix[rank]
        inverse = next(x for x in range(1,8) if gf8_mul(x,matrix[rank][j])==1)
        matrix[rank] = [gf8_mul(x,inverse) for x in matrix[rank]]
        for i in range(len(matrix)):
            if i!=rank and matrix[i][j]:
                scale = matrix[i][j]
                matrix[i] = [x^gf8_mul(scale,y) for x,y in zip(matrix[i],matrix[rank])]
        rank += 1
        if rank==len(matrix):
            break
    return rank


def exact_tests(guard):
    """Independent GF8 matrix counts and direct vector enumeration."""
    coset = load_coset()
    q,length = 8,8
    a = [(1,x,gf8_mul(x,x),gf8_mul(gf8_mul(x,x),x)) for x in range(q)]
    c = tuple(zip(*a))
    assert [[gf8_dot(row,col) for col in zip(*a)] for row in c]==[[0]*4 for _ in range(4)]
    s = [[gf8_dot(a[i],a[j])^int(i==j) for j in range(length)] for i in range(length)]
    identity = [[int(i==j) for j in range(length)] for i in range(length)]
    assert [[gf8_dot(row,col) for col in zip(*s)] for row in s]==identity
    rank_checks = 0
    for mask in range(1<<length):
        guard()
        chosen = [a[i] for i in range(length) if mask&(1<<i)]
        matrix = [[int(i==j) for j in range(4)] for i in range(4)]
        for row in chosen:
            for i in range(4):
                for j in range(4):
                    matrix[i][j] ^= gf8_mul(row[i],row[j])
        assert matrix[0][0] in (0,1) and matrix[0][3]==matrix[1][2]
        assert all(matrix[i][i]==(1^gf8_mul(matrix[0][i],matrix[0][i])) for i in range(1,4))
        assert gf8_rank(matrix,4)>=2
        rank_checks += 1
    # Rank-one tensors exhaust the symmetric rank-one forms over GF8.
    for x in product(range(q),repeat=4):
        matrix = [[gf8_mul(x[i],x[j]) for j in range(4)] for i in range(4)]
        fits = matrix[0][0] in (0,1) and matrix[0][3]==matrix[1][2]
        fits = fits and all(matrix[i][i]==(1^gf8_mul(matrix[0][i],matrix[0][i])) for i in range(1,4))
        assert not fits
    selected = ((0,1,2,3),(0,2,5,7),(0,1,2,3,4),(0,1,3,5,6,7),tuple(range(8)))
    subset_checks = moment_checks = 0
    direct_histogram = [0]*(length+1)
    for x in product(range(q),repeat=4):
        u = tuple(gf8_dot(tuple(row[i] for i in selected[0]),x) for row in c)
        if any(u):
            y = tuple(gf8_dot(tuple(row[i] for i in selected[0]),x) for row in s)
            direct_histogram[sum(value!=0 for value in y)] += 1
    assert sum(direct_histogram)==q**4-1
    for points in selected:
        guard()
        ell,n = len(points),length-len(points)
        ci = [tuple(row[i] for i in points) for row in c]
        assert gf8_rank(ci,ell)==4
        terms = []
        for mask in range(1<<length):
            zeros = tuple(i for i in range(length) if mask&(1<<i))
            rows = [tuple(s[j][i] for i in points) for j in zeros]
            r0,r1 = gf8_rank(rows,ell),gf8_rank(rows+ci,ell)
            actual = F(1,q**r0)-F(1,q**r1)
            t = sum(i in points for i in zeros)
            v,k = len(zeros)-t,ell-t
            if v>=4:
                upper = F(0)
            elif k>=4:
                upper = F(1,q**t)*(F(1,q**v)-F(1,q**4))
                assert actual==upper
            else:
                upper = F(q**min(4-v,k+2,n+k-v)-1,q**ell)
            assert 0<=actual<=upper
            terms.append((len(zeros),actual))
            subset_checks += 1
        for z in (F(1,4),F(1,2),F(3,4),F(1)):
            actual = sum((z**(length-size)*(1-z)**size*mass for size,mass in terms),F(0))
            upper = positive_moment(ell,z,q,length)
            old = positive_moment(ell,z,q,length,False)
            expected = q**4*coset.quiet_moment(ell,z,q,4)*coset.outside_nonzero(n,z,q,4)
            assert old==expected and 0<=actual<=upper<=old
            if ell==length:
                assert upper==(z+(1-z)/q)**ell-coset.quiet_moment(ell,z,q,4)
            if points==selected[0]:
                direct = sum((count*z**weight for weight,count in enumerate(direct_histogram)),F(0))/q**ell
                assert actual==direct
            moment_checks += 1
    rational_checks = 0
    for z in (F(1,4),F(1,2),F(3,4)):
        for ell in range(4,17):
            guard()
            quiet = coset.quiet_moment(ell,z)
            old = Q**4*quiet*coset.outside_nonzero(16-ell,z)
            assert positive_moment(ell,z,rank_refinement=False)==old
            assert 0<=positive_moment(ell,z)<=old
            if ell==16:
                assert positive_moment(ell,z)==(z+(1-z)/Q)**16-quiet
            rational_checks += 1
    return {"GF8_modulus":"X^3+X+1","CA_zero":True,"S_squared_identity":True,
            "all_subset_rank_checks":rank_checks,"rank_one_tensors":q**4,
            "zero_subset_count_checks":subset_checks,"positive_moment_checks":moment_checks,
            "direct_input_vectors":q**4,"GF256_rational_checks":rational_checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(1,4))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--tests-only",action="store_true")
    args = parser.parse_args()
    if not 0<args.z<1 or not 0<args.alpha<=1:
        parser.error("require0<z<1,0<alpha<=1")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    deadline = started+SECONDS
    coset = load_coset()
    anchor = coset.load("cheap_lane_cone_cauchy_anchors.py",coset.ANCHOR_PIN,"residual_return_anchor")
    peak = anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline:
            raise TimeoutError("180-second residual-rank deadline; no retry")
        if peak()>MEMORY:
            raise MemoryError("256MiB residual-rank budget exceeded")
    tests = exact_tests(guard)
    result = {"scope":"local positive return refinement only; no global bound",
              "z":str(args.z),"alpha":str(args.alpha),"tests":tests,
              "coset_pin":COSET_PIN,"generated_bytes":0,"whole_code_certificate":False}
    if not args.tests_only:
        composition = coset.load_composition()
        triple = composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,
                                  "residual_return_triple")
        from flint import ctx
        ctx.prec = 256
        caps = build_caps(args.z,args.alpha,triple,anchor,guard,deadline)
        result.update({"actual_field_tests":caps["tests"],
                       "coset_return":caps["coset_receipt"],
                       "residual_rank_return":caps["residual_rank_receipt"],"arb_precision":256})
    guard()
    result.update({"seconds":time.monotonic()-started,"peak_process_bytes":peak(),
                   "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()
