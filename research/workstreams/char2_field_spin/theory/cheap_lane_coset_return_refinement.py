"""Retain inside-coordinate weights in the nonzero-syndrome return cap.

Fix an actual GF256 input-point set I, |I|=ell>=4, in one emission
block. Write C=A^t, CA=0, S=Id+AC, and
Q_ell=q^-ell sum_{Cx=0} z^wt(x). The frozen Q_ell is exact.
Put n=16-ell and
B_n=sum_{v=0}^{min(3,n)} binom(n,v) z^(n-v)(1-z)^v(q^-v-q^-4).
This positive sum is the outside-coordinate moment for uniform nonzero
syndromes, with normalization q^-4 (not (q^4-1)^-1).

Since C_I is onto, conditional on u=C_I x the inside output x+A_I u
is uniform on an affine coset of ker C_I. Positive zero-subset expansion
bounds each such coset moment by the homogeneous moment q^4 Q_ell.
This bound holds for every u; it does not assume independence of inside
and outside output weights. Consequently
q^-ell sum_{C_I x!=0} z^wt(S_I x) <= q^4 Q_ell B_(16-ell).
Also q^4 Q_ell<=1 and B_n<=(1-q^-4)z^max(n-3,0). Thus this branch
never exceeds the old floor branch, but we still intersect the entire old
D_ell cap. The ell1/2/3 refinements and every quiet Q_ell are preserved.

augment_caps copies the frozen composition's maps and changes only raw E
and its single local fractional power. Subsequent pair correction/Jensen
tables must be rebuilt from the matching raw and powered maps. This file
computes no route, source, count, regional matrix, or global probability.
Standalone execution checks exact small-field cosets and rational bounds,
then builds local caps only. Parent owns all runs.180s/256MiB/20GiB;
stdout only, no generated files, no frozen-source edits.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations, combinations_with_replacement, product
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
COMPOSITION_PIN = "57c3abe184ce6ee643779c48106b79f1d39c85bfb4175e28404ec86b88bc7169"
ANCHOR_PIN = "319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
Q, DIMENSION, LENGTH = 256, 4, 16
SECONDS, MEMORY = 180, 256*(1 << 20)


def load(name, digest, module_name):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError("frozen helper hash mismatch: "+name)
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_composition():
    return load("cheap_lane_pair_concavity_triple.py", COMPOSITION_PIN,
                "coset_return_composition")


def outside_nonzero(n, z, q=Q, dimension=DIMENSION):
    """Exact normalized outside moment of a punctured MDS evaluation map."""
    if not isinstance(z, F) or not 0 < z <= 1 or not 0 <= n or not 1 <= dimension < q:
        raise ValueError("require rational0<z<=1,n>=0,1<=dimension<q")
    return sum((comb(n, v)*z**(n-v)*(1-z)**v
                *(F(1, q**v)-F(1, q**dimension))
                for v in range(min(dimension-1, n)+1)), F(0))


def quiet_moment(ell, z, q=Q, dimension=DIMENSION):
    """Positive MDS zero-subset formula, used independently in exact tests."""
    return sum((comb(ell, t)*z**t*((1-z)/q)**(ell-t)/q**min(dimension, t)
                for t in range(ell+1)), F(0))


def refine_block_caps(z, quiet, previous):
    """Return exact rational D0..D16 and a compact, reproducible receipt."""
    if not isinstance(z, F) or not 0 < z < 1:
        raise ValueError("require rational0<z<1")
    if len(quiet) != LENGTH+1 or len(previous) != LENGTH+1:
        raise ValueError("expected17 block caps")
    refined = list(previous)
    rows = []
    for ell in range(DIMENSION, LENGTH+1):
        b = outside_nonzero(LENGTH-ell, z)
        factor = Q**DIMENSION*quiet[ell]
        candidate = factor*b
        floor = (Q**DIMENSION-1)*quiet[ell]*z**max(0, 13-ell)
        assert 0 <= factor <= 1 and 0 <= candidate <= b and candidate <= floor
        refined[ell] = min(previous[ell], candidate)
        rows.append({"ell":ell, "inside_factor":str(factor),
                     "outside_moment":str(b), "candidate":str(candidate),
                     "old_D":str(previous[ell]), "new_D":str(refined[ell])})
    assert refined[:DIMENSION] == list(previous[:DIMENSION])
    assert all(0 <= new <= old for new, old in zip(refined, previous))
    return refined, {"q":Q, "dimension":DIMENSION, "length":LENGTH,
                     "normalization":"q^-ell on input; q^-4 on outside syndrome",
                     "ell1_2_3_unchanged":True, "quiet_unchanged":True,
                     "refined_ell":[row["ell"] for row in rows
                                    if F(row["new_D"]) < F(row["old_D"])],
                     "blocks":rows}


def augment_caps(caps, z, alpha, guard):
    """Copy composition caps, replace E only, and retain the same interface.

    The input must be the pinned ell1/2/3 composition at exactly z,alpha.
    No mutation of its maps, rows, tests, or reports is performed. All
    non-E Arb objects are reused unchanged. No pair correction is reused.
    """
    if not isinstance(z, F) or not isinstance(alpha, F) or not 0 < alpha <= 1:
        raise ValueError("require rational witnesses and0<alpha<=1")
    v, core, small = caps["verifier"], caps["core"], caps["small"]
    if core.Q != Q or core.Z != z or core.ALPHA != alpha:
        raise ValueError("caps do not match actual GF256 and requested witnesses")
    _, quiet, old_d, _, _ = core.rational_caps()
    assert all(value == quiet_moment(ell, z) for ell, value in enumerate(quiet))
    previous = list(old_d)
    block_report = dict(caps["block_report"])
    for ell, key in ((1,"new_D1"),(2,"new_D2"),(3,"refined_D3")):
        previous[ell] = F(block_report[key])
    refined, receipt = refine_block_caps(z, quiet, previous)
    raw, powered = {}, {}
    exponent = v.exact(alpha)
    changed = 0
    largest_ratio = F(1)
    for occ, row in caps["raw"].items():
        guard()
        before = small.return_cap(occ, quiet, previous, core.KAPPA)
        after = small.return_cap(occ, quiet, refined, core.KAPPA)
        assert 0 <= after <= before
        # Exact rational endpoints, not midpoint or Arb-ball containment.
        exact_before = v.fmpq(before.numerator, before.denominator)
        assert row[2].lower().fmpq() <= exact_before <= row[2].upper().fmpq()
        new_row, new_powered = list(row), list(caps["powered"][occ])
        new_row[2] = v.exact(after)
        new_powered[2] = v.positive(new_row[2])**exponent
        raw[occ], powered[occ] = tuple(new_row), tuple(new_powered)
        for channel in (0,1,3,4,5):
            assert raw[occ][channel] is row[channel]
            assert powered[occ][channel] is caps["powered"][occ][channel]
        if after < before:
            changed += 1
            assert after > 0
            largest_ratio = max(largest_ratio, before/after)
    assert raw.keys() == powered.keys() == caps["raw"].keys() == caps["powered"].keys()
    receipt.update({"changed_shapes":changed, "shape_count":len(raw),
                    "maximum_raw_E_ratio":str(largest_ratio),
                    "changed_channel":"raw E and its single local alpha power only",
                    "requires_rebuilt_pair_corrections":True})
    result = dict(caps)
    block_report["coset_return"] = receipt
    result.update({"raw":raw, "powered":powered, "block_report":block_report,
                   "tests":dict(caps["tests"]), "coset_receipt":receipt,
                   "quiet_exact":tuple(quiet), "nonquiet_exact":tuple(refined)})
    return result


def build_caps(z, alpha, triple, anchor, guard, deadline):
    """Drop-in replacement for the frozen composition's build_caps."""
    composition = load_composition()
    baseline = composition.build_caps(z, alpha, triple, anchor, guard, deadline)
    return augment_caps(baseline, z, alpha, guard)


def gf4_mul(a, b):
    """Independent GF4 arithmetic, modulus X^2+X+1, for exhaustive toys."""
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return value


def gf4_dot(a, b):
    result = 0
    for x, y in zip(a, b):
        result ^= gf4_mul(x, y)
    return result


def exact_tests(guard):
    """Exhaust all GF4 cosets and validate independent rational inequalities."""
    q, k, length = 4, 2, 4
    a = [(1, point) for point in range(q)]
    c = tuple(zip(*a))
    assert [[gf4_dot(row, col) for col in zip(*a)] for row in c] == [[0,0],[0,0]]
    s = [[gf4_dot(a[i], a[j]) ^ int(i == j) for j in range(length)]
         for i in range(length)]
    identity = [[int(i == j) for j in range(length)] for i in range(length)]
    assert [[gf4_dot(row, col) for col in zip(*s)] for row in s] == identity
    syndromes = tuple(product(range(q), repeat=k))
    coset_checks = moment_checks = 0
    for ell in range(k, length+1):
        for points in combinations(range(length), ell):
            guard()
            outside = tuple(i for i in range(length) if i not in points)
            groups = {u:[] for u in syndromes}
            for x in product(range(q), repeat=ell):
                u = tuple(gf4_dot(tuple(row[i] for i in points), x) for row in c)
                inside = tuple(x[j] ^ gf4_dot(a[i], u) for j, i in enumerate(points))
                full_y = tuple(gf4_dot(tuple(row[i] for i in points), x) for row in s)
                assert sum(value != 0 for value in full_y) == (
                    sum(value != 0 for value in inside)
                    + sum(gf4_dot(a[i], u) != 0 for i in outside))
                groups[u].append((x, inside, full_y))
            assert all(len(vectors) == q**(ell-k) for vectors in groups.values())
            for z in (F(1,4),F(1,2),F(3,4),F(1)):
                quiet = sum((z**sum(t != 0 for t in x)
                             for x, _, _ in groups[(0,0)]), F(0))/q**ell
                assert quiet == quiet_moment(ell,z,q,k)
                homogeneous = q**k*quiet
                for u, vectors in groups.items():
                    inside_moment = sum((z**sum(t != 0 for t in inside)
                                         for _, inside, _ in vectors), F(0))/len(vectors)
                    assert inside_moment <= homogeneous
                    coset_checks += 1
                b = outside_nonzero(length-ell,z,q,k)
                direct_b = sum((z**sum(gf4_dot(a[i],u) != 0 for i in outside)
                                for u in syndromes if any(u)), F(0))/q**k
                assert b == direct_b
                direct_d = sum((z**sum(t != 0 for t in y)
                                for u, vectors in groups.items() if any(u)
                                for _, _, y in vectors), F(0))/q**ell
                candidate = homogeneous*b
                floor = (q**k-1)*quiet*z**max(0,length-ell-k+1)
                g = z+(1-z)/q
                assert 0 <= direct_d <= candidate <= b
                assert candidate <= floor and direct_d <= g**ell-quiet
                moment_checks += 1
    rational_checks = 0
    for z in (F(1,4),F(1,2),F(3,4)):
        guard()
        quiet = [quiet_moment(ell,z) for ell in range(LENGTH+1)]
        g = z+(1-z)/Q
        old = [min(g**ell-quiet[ell],
                   (Q**min(DIMENSION,ell)-1)*quiet[ell]*z**max(0,13-ell))
               for ell in range(LENGTH+1)]
        refined, _ = refine_block_caps(z,quiet,old)
        assert refined[:4] == old[:4]
        assert refined[16] == old[16]
        # Independent subset expansion checks the positive OR recurrence.
        for occ in ((0,0,0,0),(1,2,3,4),(4,4,4,4),(4,8,12,16),(16,16,16,16)):
            def active(ds):
                qp, dp = F(1), F(0)
                for ell in occ:
                    qp, dp = qp*quiet[ell], dp*(quiet[ell]+ds[ell])+qp*ds[ell]
                direct = F(0)
                for mask in range(1,16):
                    term = F(1)
                    for i, ell in enumerate(occ):
                        term *= ds[ell] if mask & (1 << i) else quiet[ell]
                    direct += term
                assert direct == dp
                return dp
            assert 0 <= active(refined) <= active(old)
            rational_checks += 1
        rational_checks += 13
    return {"GF4_modulus":"X^2+X+1", "GF4_CA_zero":True,
            "GF4_S_squared_identity":True, "GF4_input_sets":11,
            "exact_coset_checks":coset_checks, "exact_moment_checks":moment_checks,
            "GF256_rational_monotonicity_checks":rational_checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", type=F, default=F(1,4))
    parser.add_argument("--alpha", type=F, default=F(119,1280))
    parser.add_argument("--tests-only", action="store_true")
    args = parser.parse_args()
    if not 0 < args.z < 1 or not 0 < args.alpha <= 1:
        parser.error("require0<z<1,0<alpha<=1")
    if max(n.bit_length() for x in (args.z,args.alpha)
           for n in (x.numerator,x.denominator)) > 64:
        parser.error("rational witnesses limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    deadline = started+SECONDS
    anchor = load("cheap_lane_cone_cauchy_anchors.py", ANCHOR_PIN, "coset_return_anchor")
    peak = anchor.peak_memory_reader()
    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("180-second coset-return deadline; no retry")
        if peak() > MEMORY:
            raise MemoryError("256MiB coset-return budget exceeded")
    tests = exact_tests(guard)
    result = {"scope":"local return caps only; no route/source/count or global bound",
              "z":str(args.z), "alpha":str(args.alpha), "tests":tests,
              "composition_pin":COMPOSITION_PIN, "anchor_pin":ANCHOR_PIN,
              "generated_bytes":0, "whole_code_certificate":False}
    if not args.tests_only:
        composition = load_composition()
        triple = composition.load("cheap_lane_exact_triple_return.py",
                                  composition.TRIPLE_PIN,"coset_return_triple")
        from flint import ctx
        ctx.prec = 256
        caps = build_caps(args.z,args.alpha,triple,anchor,guard,deadline)
        result.update({"actual_field_tests":caps["tests"],
                       "coset_return":caps["coset_receipt"], "arb_precision":256})
    guard()
    result.update({"seconds":time.monotonic()-started, "peak_process_bytes":peak(),
                   "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(json.dumps(result,indent=2),flush=True)


if __name__ == "__main__":
    main()
