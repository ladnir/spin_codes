"""Exact finite-field checks of the prospective weight-40/42 endpoint test.

No spectrum cap is accepted by this script. Exhaustive small-field checks and
known BCH witnesses validate the reconstruction, not the unknown shell sizes.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import random
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True
from analyze_endpoint_sampling import (ceil_fraction, log_integer,
                                       neg_log_one_minus, ratio_record)
from bch_quotient import verify_quotient_algebra

ROOT = Path(__file__).resolve().parents[1]


class Field:
    def __init__(self, q, modulus):
        self.q = q
        self.modulus = modulus
        def multiply(a, b):
            result = 0
            while b:
                if b & 1:
                    result ^= a
                a <<= 1
                if a & q:
                    a ^= modulus
                b >>= 1
            return result
        self.mul = [bytes(multiply(a, b) for b in range(q)) for a in range(q)]
        self.inv = [0] + [self.power(x, q - 2) for x in range(1, q)]
        assert all(self.mul[x][self.inv[x]] == 1 for x in range(1, q))

    def power(self, a, e):
        result = 1
        while e:
            if e & 1:
                result = self.mul[result][a]
            a = self.mul[a][a]
            e >>= 1
        return result

    def evaluate(self, coefficients, x):
        result = 0
        row = self.mul[x]
        for c in reversed(coefficients):
            result = row[result] ^ c
        return result

    def product(self, points, locator=False):
        result = [1]
        for x in points:
            updated = [0] * (len(result) + 1)
            for j, c in enumerate(result):
                updated[j] ^= c if locator else self.mul[x][c]
                updated[j + 1] ^= self.mul[x][c] if locator else c
            result = updated
        return result

    def solve(self, matrix, rhs):
        n = len(rhs)
        rows = [row[:] + [v] for row, v in zip(matrix, rhs)]
        for j in range(n):
            pivot = next((i for i in range(j, n) if rows[i][j]), None)
            if pivot is None:
                return None
            rows[j], rows[pivot] = rows[pivot], rows[j]
            scale = self.mul[self.inv[rows[j][j]]]
            rows[j] = [scale[v] for v in rows[j]]
            for i in range(n):
                if i == j or rows[i][j] == 0:
                    continue
                scale = self.mul[rows[i][j]]
                rows[i] = [a ^ scale[b] for a, b in zip(rows[i], rows[j])]
        return [row[-1] for row in rows]

    def remainders(self, nodes, exponent, r):
        # Newton interpolation of u^exponent, also producing P_B.
        coefficients = [0] * len(nodes)
        product = [1]
        for x in nodes:
            delta = self.power(x, exponent) ^ self.evaluate(coefficients, x)
            scale = self.mul[delta][self.inv[self.evaluate(product, x)]]
            for j, c in enumerate(product):
                coefficients[j] ^= self.mul[scale][c]
            updated = [0] * (len(product) + 1)
            for j, c in enumerate(product):
                updated[j] ^= self.mul[x][c]
                updated[j + 1] ^= c
            product = updated
        result = [coefficients]
        for _ in range(r):
            previous = result[-1]
            leading = previous[-1]
            shifted = [0] + previous[:-1]
            result.append([c ^ self.mul[leading][p]
                           for c, p in zip(shifted, product)])
        return result

    def reconstruct(self, remainders, sigma, a, r):
        n, m = a + 2*r, a + r
        constraint_indices = [0] + list(range(m + 1, n))
        assert len(constraint_indices) == r
        matrix = [[remainders[j][i] for j in range(1, r + 1)]
                  for i in constraint_indices]
        rhs = [(int(i == 0) ^ self.mul[sigma][remainders[0][i]])
               for i in constraint_indices]
        solution = self.solve(matrix, rhs)
        if solution is None:
            return None
        h = [sigma] + solution
        g = [0] * n
        for coefficient, remainder in zip(h, remainders):
            for i, v in enumerate(remainder):
                g[i] ^= self.mul[coefficient][v]
        assert g[0] == 1 and all(v == 0 for v in g[m + 1:])
        return g[:m + 1], h


def endpoint_roots(field, g, h, exponent):
    return tuple(x for x in range(1, field.q)
                 if field.evaluate(g, x) == field.mul[field.power(x, exponent)][
                     field.evaluate(h, x)])


def exhaustive_small_field():
    field = Field(16, 0x13)
    a, delta, exponent = 1, 3, 9
    reports = []
    for r in (1, 2):
        d, n, w = delta + 2*r, a + 2*r, delta + 2*r + 1
        expected = {}
        # Independent enumeration of supports, using the parity-check p_1=0.
        for support in itertools.combinations(range(1, 16), d):
            syndrome1 = 0
            for x in support:
                syndrome1 ^= x
            if syndrome1:
                continue
            sigma = 0
            for x in support:
                sigma ^= field.power(x, delta)
            roots = tuple(sorted(field.power(x, 13) for x in support))
            expected[(sigma, roots)] = 0
        found = Counter()
        singular_queries = 0
        for nodes in itertools.combinations(range(1, 16), n):
            remainders = field.remainders(nodes, exponent, r)
            for sigma in range(16):
                candidate = field.reconstruct(remainders, sigma, a, r)
                if candidate is None:
                    singular_queries += 1
                    continue
                g, h = candidate
                roots = endpoint_roots(field, g, h, exponent)
                if len(roots) == d:
                    assert h[-1] != 0
                    found[(sigma, roots)] += 1
        assert set(found) == set(expected)
        assert set(found.values()) == {math.comb(d, n)}
        # Independent enumeration of even extended supports verifies puncturing.
        extended_count = 0
        for support in itertools.combinations(range(16), w):
            syndrome1 = 0
            for x in support:
                syndrome1 ^= x
            extended_count += syndrome1 == 0
        assert len(expected) * 16 == w * extended_count
        reports.append(dict(field_size=16, locator_degree=d, base_size=n,
                            queries=16*math.comb(15, n),
                            endpoint_supports=len(expected),
                            endpoint_queries=sum(found.values()),
                            expected_queries_per_endpoint=math.comb(d, n),
                            singular_queries=singular_queries,
                            extended_shell_count=extended_count,
                            exhaustive_identity_passed=True))
        print(json.dumps(reports[-1]), flush=True)
    return reports


def witness_words(weight, family):
    if weight == 40:
        path = ROOT / 'generated' / f'wambach_{family}_orbits_r5_w40.json'
        records = json.loads(path.read_text())['shells']['40']['orbits']
        return path, [int(row['canonical_hex'], 16) for row in records[:8]]
    path = ROOT / 'generated' / f'wambach_{family}_orbits_r5_w42.tsv'
    with path.open(newline='') as stream:
        rows = list(itertools.islice(csv.DictReader(stream, delimiter='\t'), 8))
    return path, [sum(int(row[f'limb{i}'], 16) << (64*i) for i in range(4))
                  for row in rows]


def bch_positive_controls():
    field = Field(256, 0x14d)
    powers = [field.power(2, i) for i in range(255)]
    assert len(set(powers)) == 255
    rng = random.Random(20260904)  # Structural test only, never a certificate.
    reports = []
    for weight in (40, 42):
        r, a, exponent = (weight - 38)//2, 18, 146
        d, n = weight - 1, weight - 20
        for family in ('p', 'q'):
            path, words = witness_words(weight, family)
            sigma_counts = Counter()
            checks = 0
            for word in words:
                support = [powers[i] for i in range(255) if (word >> i) & 1]
                if (word >> 255) & 1:
                    support.append(0)
                assert len(support) == weight
                anchor = support[0]
                punctured = [x ^ anchor for x in support if x != anchor]
                locator = field.product(punctured, locator=True)
                assert all(locator[j] == 0 for j in range(1, 37, 2))
                # Place nonzero syndrome in the representative subspace S.
                sigma = locator[37]
                if sigma:
                    scale = field.power(field.inv[sigma], pow(37, -1, 255))
                    punctured = [field.mul[scale][x] for x in punctured]
                    locator = field.product(punctured, locator=True)
                sigma = locator[37]
                assert sigma in (0, 1)
                for j in range(1, 38):
                    syndrome = 0
                    for x in punctured:
                        syndrome ^= field.power(x, j)
                    assert syndrome == (sigma if j == 37 else 0)
                if family == 'q':
                    assert sigma == 0
                sigma_counts[sigma] += 1
                g, h = locator[::2], locator[37::2]
                roots = tuple(sorted(field.power(x, 253) for x in punctured))
                assert endpoint_roots(field, g, h, exponent) == roots
                for _ in range(64):
                    nodes = rng.sample(roots, n)
                    remainders = field.remainders(nodes, exponent, r)
                    candidate = field.reconstruct(remainders, sigma, a, r)
                    assert candidate == (g, h), (weight, family, sigma)
                    # Independent full n-by-n solve checks the reduced system.
                    matrix = [[field.power(x, j) for j in range(1, a+r+1)] +
                              [field.power(x, exponent+j) for j in range(1, r+1)]
                              for x in nodes]
                    rhs = [1 ^ field.mul[sigma][field.power(x, exponent)]
                           for x in nodes]
                    assert field.solve(matrix, rhs) == g[1:] + h[1:]
                    checks += 1
            reports.append(dict(weight=weight, family=family, witnesses=len(words),
                                sigma_counts=dict(sigma_counts),
                                full_and_reduced_reconstructions=checks,
                                source=str(path.relative_to(ROOT)),
                                source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                passed=True))
            print(json.dumps(reports[-1]), flush=True)
    return reports


def random_query_cross_checks():
    field = Field(256, 0x14d)
    rng = random.Random(20260905)
    reports = []
    for r in (1, 2):
        a, exponent = 18, 146
        n = a + 2*r
        singular = 0
        for _ in range(256):
            nodes = rng.sample(range(1, 256), n)
            sigma = rng.randrange(32)
            candidate = field.reconstruct(field.remainders(nodes, exponent, r),
                                          sigma, a, r)
            matrix = [[field.power(x, j) for j in range(1, a+r+1)] +
                      [field.power(x, exponent+j) for j in range(1, r+1)]
                      for x in nodes]
            rhs = [1 ^ field.mul[sigma][field.power(x, exponent)] for x in nodes]
            full = field.solve(matrix, rhs)
            assert (candidate is None) == (full is None)
            if candidate is None:
                singular += 1
            else:
                g, h = candidate
                assert full == g[1:] + h[1:]
        reports.append(dict(weight=38+2*r, queries=256,
                            singular_queries=singular, passed=True))
    return reports


def sample_plans():
    plans = []
    for weight, cap in ((40, 50_000_000_000_000), (42, 5_000_000_000_000_000)):
        per_word = Fraction(weight * math.comb(weight-1, weight-20),
                            8192 * math.comb(255, weight-20))
        p_bad = (cap + 1) * per_word
        lo, hi = neg_log_one_minus(p_bad)
        log_lo, log_hi = log_integer(1 << 41)
        samples = ceil_fraction(log_hi / lo)
        assert samples * lo >= log_hi
        assert (samples - 1) * hi < log_lo
        plans.append(dict(weight=weight, target_cap=cap,
                          first_bad_integer_conservative=cap+1,
                          per_word_probability=ratio_record(per_word),
                          bad_shell_probability_lower=ratio_record(p_bad),
                          fixed_zero_hit_sample_budget=samples,
                          ideal_IID_false_accept_bound='2^-41',
                          integer_budget_certified_by_rational_log_enclosures=True,
                          execution_status='prospective; not run'))
    return plans


def main():
    payload = dict(classification='exact structural checks and prospective sample budgets; no shell cap established',
                   quotient_algebra=verify_quotient_algebra(),
                   small_field=exhaustive_small_field(),
                   bch_positive_controls=bch_positive_controls(),
                   random_query_cross_checks=random_query_cross_checks(),
                   prospective_plans=sample_plans(),
                   familywise_error='Existing weight-38 test at 2^-40 plus two new tests at 2^-41 each: union bound 2^-39 under ideal IID bytes; not SPIN setup failure probability.',
                   script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    path = ROOT / 'generated' / 'higher_shell_endpoint_reduction_checks.json'
    path.write_text(json.dumps(payload, indent=2) + '\n')
    print(path, flush=True)
    print(json.dumps(payload['prospective_plans'], indent=2), flush=True)


if __name__ == '__main__':
    main()
