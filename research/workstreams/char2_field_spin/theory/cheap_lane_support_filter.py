"""Exact finite-field checks for a local support-filter inequality.

Let E have q elements. Fix a matrix A with n columns and U=ker_E(A).
Fix p<q E-linear forms l_1,...,l_p on E^n. Repeated forms are allowed.
Let F mean that every l_i restricts to a nonzero form on U. Then

    1_F q^(-rank A)
      <= q/(q-p) Pr_{v uniform in E^n}[Av=0 and all l_i(v)!=0].

If F fails, the right side and the filtered left side are both zero.
Otherwise each forbidden hyperplane has conditional probability 1/q in U.
The union bound leaves at least (q-p)/q of U. No independence is asserted.
For p=0 the empty conjunction is true, including when U={0}.

The equations and forms are fixed before applying this comparison. For
an extinction transition, incoming-state forms may remain required while
the next-state equations force zero. A terminal has no next-state equations.
These are local facts; this helper makes no claim about a chain or a code.

Tests enumerate kernels and multisets of projective forms over small fields.
Named examples cover extinction, terminal equations, duplicate/scaled forms,
zero-dimensional kernels, and a pencil of p=q-1 hyperplanes in E^2.
That pencil leaves exactly 2(q-1) vectors: it shows an order-q necessary
penalty, not equality in the coarse q/(q-p) bound. The p=1 bound is exact.

Standard library only. No generated files, caches, random sampling, or search.
Run with python -B and --self-test. The test budget is five seconds and
requires at least 20 GiB free on the source volume before enumeration.
"""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
import json
import math
from pathlib import Path
import shutil
import time


class Field:
    def __init__(self, q):
        if q not in (2, 3, 4, 8):
            raise ValueError("only the explicitly tested tiny fields are supported")
        self.q = q
        self.add = tuple(tuple((a+b) % 3 if q == 3 else a ^ b
                               for b in range(q)) for a in range(q))
        self.mul = tuple(tuple(self._multiply(a, b) for b in range(q))
                         for a in range(q))
        self.inv = (0,) + tuple(next(b for b in range(1, q)
                                    if self.mul[a][b] == 1)
                                for a in range(1, q))

    def _multiply(self, a, b):
        if self.q == 3:
            return a*b % 3
        modulus = {2: 3, 4: 7, 8: 11}[self.q]
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= modulus
        return out

    def dot(self, a, b):
        out = 0
        for x, y in zip(a, b):
            out = self.add[out][self.mul[x][y]]
        return out

    def rank(self, rows, n):
        a = [list(row) for row in rows]
        pivot = 0
        for col in range(n):
            found = next((i for i in range(pivot, len(a)) if a[i][col]), None)
            if found is None:
                continue
            a[pivot], a[found] = a[found], a[pivot]
            inv = self.inv[a[pivot][col]]
            a[pivot] = [self.mul[inv][x] for x in a[pivot]]
            for i in range(pivot+1, len(a)):
                factor = a[i][col]
                if factor:
                    negative = (3-factor) % 3 if self.q == 3 else factor
                    a[i] = [self.add[x][self.mul[negative][y]]
                            for x, y in zip(a[i], a[pivot])]
            pivot += 1
            if pivot == len(a):
                break
        return pivot


def budget(deadline):
    if time.monotonic() > deadline:
        raise TimeoutError("five-second enumeration budget exceeded")


def vectors(q, n):
    return tuple(product(range(q), repeat=n))


def projective_forms(field, n):
    normalized = set()
    for form in vectors(field.q, n):
        first = next((x for x in form if x), None)
        if first is not None:
            inv = field.inv[first]
            normalized.add(tuple(field.mul[inv][x] for x in form))
    return tuple(sorted(normalized))


def kernels(field, n, deadline):
    """Enumerate every kernel via a basis chosen among projective rows."""
    ambient = vectors(field.q, n)
    rows = projective_forms(field, n)
    seen = {}
    for count in range(n+1):
        for matrix in combinations(rows, count):
            budget(deadline)
            kernel = tuple(v for v in ambient
                           if all(field.dot(row, v) == 0 for row in matrix))
            if kernel not in seen:
                rank = field.rank(matrix, n)
                assert len(kernel) == field.q**(n-rank)
                seen[kernel] = (matrix, rank)
    return tuple((matrix, rank, kernel)
                 for kernel, (matrix, rank) in seen.items())


def check_cached(q, n, rank, kernel, form_indices, zero_masks):
    p = len(form_indices)
    if p >= q:
        raise ValueError("the support-filter comparison requires p<q")
    full = (1 << len(kernel))-1
    good = full
    possible = True
    for index in form_indices:
        mask = zero_masks[index]
        possible &= mask != full
        good &= full ^ mask
    good_count = good.bit_count()
    assert possible == (good_count > 0)
    unfiltered = Fraction(1, q**rank)
    assert unfiltered == Fraction(len(kernel), q**n)
    lhs = unfiltered if possible else Fraction(0)
    rhs = Fraction(q, q-p)*Fraction(good_count, q**n)
    assert lhs <= rhs, (q, n, rank, form_indices, lhs, rhs)
    return possible, lhs == rhs, good_count


def check_case(field, n, matrix, forms):
    if any(len(row) != n for row in tuple(matrix)+tuple(forms)):
        raise ValueError("matrix rows and forms must have n coordinates")
    kernel = tuple(v for v in vectors(field.q, n)
                   if all(field.dot(row, v) == 0 for row in matrix))
    rank = field.rank(matrix, n)
    masks = tuple(sum(1 << i for i, v in enumerate(kernel)
                      if field.dot(form, v) == 0) for form in forms)
    possible, equality, good = check_cached(
        field.q, n, rank, kernel, tuple(range(len(forms))), masks)
    return {"rank": rank, "kernel_size": len(kernel), "p": len(forms),
            "possible": possible, "equality": equality, "good_vectors": good}


def self_test():
    started = time.monotonic()
    deadline = started + 5.0
    counts = Counter()
    exhaustive = []
    # All subspaces and all form multisets of size p<q, including zero forms.
    # Projective representatives suffice: scaling a form does not change zeros.
    for q, n, expected_kernels in ((2, 2, 5), (2, 3, 16),
                                   (3, 2, 6), (4, 2, 7)):
        field = Field(q)
        forms = ((0,)*n,) + projective_forms(field, n)
        spaces = kernels(field, n, deadline)
        assert len(spaces) == expected_kernels
        start_count = counts["comparisons"]
        for matrix, rank, kernel in spaces:
            masks = tuple(sum(1 << i for i, v in enumerate(kernel)
                              if field.dot(form, v) == 0) for form in forms)
            for p in range(q):
                for indices in combinations_with_replacement(range(len(forms)), p):
                    budget(deadline)
                    possible, equality, _ = check_cached(q, n, rank, kernel,
                                                         indices, masks)
                    counts["comparisons"] += 1
                    counts["possible"] += possible
                    counts["impossible"] += not possible
                    counts["equalities"] += equality
        exhaustive.append({"q": q, "ambient_dimension": n,
                           "kernels": len(spaces),
                           "comparisons": counts["comparisons"]-start_count})

    e4 = Field(4)
    # Variables are (x0,x1,u0,u1), with output and pre-refresh state x+u.
    u0, u1 = (0, 0, 1, 0), (0, 0, 0, 1)
    v0, v1 = (1, 0, 1, 0), (0, 1, 0, 1)
    named = {
        "extinction_with_live_incoming": check_case(e4, 4, (v0, v1), (u0, u1)),
        "impossible_incoming_support": check_case(e4, 4, (u0,), (u0,)),
        "terminal_output_only": check_case(e4, 4, (v0,), (u0, u1)),
        "terminal_without_equations": check_case(e4, 4, (), (u0, u1)),
        "duplicate_and_scaled_forms": check_case(e4, 4, (v0,),
                                                  (u0, u0, (0, 0, 2, 0))),
        "zero_kernel_no_forms": check_case(e4, 2, ((1, 0), (0, 1)), ()),
        "zero_kernel_required_form": check_case(e4, 2,
                                                 ((1, 0), (0, 1)), ((1, 0),)),
    }
    assert named["extinction_with_live_incoming"]["possible"]
    assert not named["impossible_incoming_support"]["possible"]
    assert named["terminal_output_only"]["possible"]
    assert named["zero_kernel_no_forms"]["equality"]
    assert not named["zero_kernel_required_form"]["possible"]

    pencils = []
    for q in (2, 3, 4, 8):
        field = Field(q)
        forms = tuple((1, a) for a in range(q-1))
        result = check_case(field, 2, (), forms)
        assert result["good_vectors"] == 2*(q-1)
        necessary = Fraction(q*q, 2*(q-1))
        assert necessary <= q
        pencils.append({"q": q, "p": q-1, "good_vectors": 2*(q-1),
                        "necessary_penalty": str(necessary),
                        "lemma_penalty": q})
    # p=1 attains q/(q-1), independently of the pencil's large-p behavior.
    assert check_case(Field(8), 2, (), ((1, 0),))["equality"]
    # A union of all q+1 projective hyperplanes covers E^2. This is outside
    # the lemma: a nonzero restriction alone is insufficient without p<q.
    boundary = Field(4)
    all_forms = projective_forms(boundary, 2)
    assert all(any(boundary.dot(f, v) == 0 for f in all_forms)
               for v in vectors(4, 2))
    try:
        check_case(boundary, 2, (), all_forms[:4])
    except ValueError:
        pass
    else:
        raise AssertionError("p=q must be rejected")
    budget(deadline)
    return {"status": "passed", "arithmetic": "exact Fraction and field operations",
            "exhaustive": exhaustive, "totals": dict(counts), "named": named,
            "pencils": pencils, "p_at_least_q_rejected": True,
            "q256_p32_penalty": "8/7",
            "q256_p32_penalty_log2_float": math.log2(8/7),
            "whole_code_certificate": False,
            "elapsed_seconds": time.monotonic()-started,
            "generated_data_bytes": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", required=True)
    parser.parse_args()
    volume = Path(__file__).resolve().anchor
    free = shutil.disk_usage(volume).free
    if free < 20*(1 << 30):
        raise SystemExit("requires at least 20 GiB free on source volume")
    result = self_test()
    result["free_space_gib"] = free/(1 << 30)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
