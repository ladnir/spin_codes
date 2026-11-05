"""Independent directed verifier for one cheap-lane source-support class.

Fixed witness: q=256, z=21/25, alpha=119/1280, distance threshold734003.
Cover all state paths but ONLY the8192 cyclic consecutive256-group supports
(4096 in each alternating matching orientation). This is not a whole-code
certificate. No frozen helper is imported and no construction is changed.

Exact rational quiet/activation/return caps and paired-route counts feed a
two-state positive cone. Original output-only emission caps use directed
Arb min/max over gamma=0,1/2,1 and compatible block cardinalities. Raise
each local cone cap to alpha BEFORE averaging routes. Positive bivariate
polynomial products extract the exact regional profiles; the final region
has127 updates and one output-only terminal. Fixed variable scales improve
conditioning and are exactly removed in the coefficient normalizers.

Arb intervals enclose the selected majorant, not just floating proposals.
Known-nonnegative polynomial coefficients may be intersected with [0,inf).
No midpoint decisions or FFT coefficients are silently discarded. stdout
only,48-second cooperative guard,20GiB free, estimated peak128MiB.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
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
Q = 256
Z, ALPHA = F(21, 25), F(119, 1280)
USCALE, VSCALE = F(1, 32), F(1, 2048)
KAPPA = F(Q**3, 255**2*241)
ETA, BETA = F(256, 255)**2, F(256**160, 255**32)
LIMIT, MAX_FULL = 734003, 128


def route_records(single_degree):
    """Exact oriented half-slot enumeration, with no omitted relevant tail."""
    half = []
    for f in range(17):
        for s in range(min(single_degree, 16-f)+1):
            for left in range(s+1):
                weight = F(comb(16, f)*comb(16-f, s)*comb(s, left), 2**s)
                half.append((f, s, (f+left, f+s-left), weight))
    records = defaultdict(F)
    for f, s, occ, weight in half:
        for ff, ss, other, ww in half:
            if s+ss <= single_degree:
                records[f+ff, s+ss, tuple(sorted(occ+other))] += weight*ww
    totals = defaultdict(F)
    for (f, s, _), weight in records.items():
        totals[f, s] += weight
    checks = 0
    for s in range(single_degree+1):
        for f in range(33-s):
            assert totals[f, s] == comb(32, f)*comb(32-f, s)
            checks += 1
    return records, checks


def rational_caps():
    g = Z+(1-Z)/Q
    quiet, nonquiet = [], []
    for ell in range(17):
        value = sum(F(comb(ell, k))*Z**k*((1-Z)/Q)**(ell-k)
                    *F(1, Q**min(4, k)) for k in range(ell+1))
        other = min(g**ell-value,
                    (Q**min(4, ell)-1)*value*Z**max(0, 13-ell))
        assert 0 <= value <= g**ell and other >= 0
        quiet.append(value); nonquiet.append(other)
    p = [F(comb(16, h))*min(F(1), F(Q, 255**(17-h))) for h in range(1, 17)]
    activation = [g**ell-((1-Z)/Q)**ell if ell else F(0) for ell in range(65)]
    assert quiet[0] == 1 and nonquiet[0] == 0 and p[-1] == 1
    assert activation[0] == 0 and activation[1] == Z
    return g, quiet, nonquiet, p, activation


class Verifier:
    def __init__(self, single_degree, deadline):
        from flint import arb, arb_poly, fmpq
        self.arb, self.poly, self.fmpq = arb, arb_poly, fmpq
        self.single, self.deadline = single_degree, deadline
        self.intersections = 0

    def guard(self):
        if time.monotonic() > self.deadline:
            raise TimeoutError("48-second directed-verifier budget exceeded")

    def exact(self, value):
        value = F(value)
        return self.arb(value.numerator)/value.denominator

    def positive(self, value):
        if not value.is_finite() or value.upper() < 0:
            raise ArithmeticError("nonfinite or negative positive-expression enclosure")
        if not value >= 0:
            self.intersections += 1
        return value.nonnegative_part()

    def blocks(self):
        """Original output-only block interface, not a later terminal refinement."""
        a = self.arb
        z = self.exact(Z)
        root = z.sqrt()
        tables = []
        for floor_mark, b in ((a(1), z), (root, root), (z, a(1))):
            g = b+(1-b)/Q
            table = []
            for ell in range(17):
                outside = 16-ell
                rv = [comb(outside, v)*b**(outside-v)*(1-b)**v
                      for v in range(min(3, outside)+1)]
                retained = sum(rv, a(0))
                row = [g**ell]
                for h in range(1, 5):
                    rank = sum((w/Q**max(0, h-4+v) for v, w in enumerate(rv)), a(0))
                    split = retained/Q**h+(1-self.exact(F(1, Q**h)))*b**max(0, outside-3)
                    value = g**ell*rank.min(split)
                    if h == 1:
                        value = value.min(g**ell*b**max(0, outside-2))
                    row.append(self.positive(value*floor_mark**max(0, 13-ell)))
                table.append(row)
            tables.append(table)
        return tables

    def local_caps(self, shapes):
        a = self.arb
        g, qb, db, p, activation = rational_caps()
        tables = self.blocks()
        # Cache the compatible two-block max-product tables, never mixing
        # placements inside one candidate product.
        pairs = []
        for table in tables:
            cache = {}
            for x in range(17):
                for y in range(x, 17):
                    row = []
                    for h in range(9):
                        candidates = [table[x][i]*table[y][h-i]
                                      for i in range(max(0, h-4), min(4, h)+1)]
                        value = candidates[0]
                        for candidate in candidates[1:]:
                            value = value.max(candidate)
                        row.append(self.positive(value))
                    cache[x, y] = row
            pairs.append(cache)
        pa = [self.exact(x) for x in p]
        result = {}
        for occ in shapes:
            self.guard()
            th = [None]*17
            for cache in pairs:
                left, right = cache[occ[:2]], cache[occ[2:]]
                for h in range(1, 17):
                    candidates = [left[i]*right[h-i]
                                  for i in range(max(0, h-8), min(8, h)+1)]
                    value = candidates[0]
                    for candidate in candidates[1:]:
                        value = value.max(candidate)
                    th[h] = value if th[h] is None else th[h].min(value)
            u = self.positive(sum((pa[h-1]*th[h] for h in range(1, 17)), a(0)))
            quiet, live = F(1), F(0)
            for ell in occ:
                live, quiet = live*(qb[ell]+db[ell])+quiet*db[ell], quiet*qb[ell]
            # Independent positive-OR identity test for this rational table.
            product = F(1)
            for ell in occ:
                product *= qb[ell]+db[ell]
            assert quiet+live == product
            ell = sum(occ)
            result[occ] = (self.exact(quiet), self.exact(activation[ell]),
                           self.exact(KAPPA*live/Q**16), u, self.exact(g**ell), u)
        return result

    def polynomials(self, records, caps):
        a, poly = self.arb, self.poly
        degree = MAX_FULL+1
        values = [[[a(0) for _ in range(degree)] for _ in range(self.single+1)]
                  for _ in range(6)]
        power = self.exact(ALPHA)
        powered = {occ: tuple(self.positive(x)**power for x in cap)
                   for occ, cap in caps.items()}
        scales = {(f, s): self.exact(USCALE**f*VSCALE**s)
                  for f in range(33) for s in range(self.single+1)}
        for (f, s, occ), weight in records.items():
            factor = self.exact(weight)*scales[f, s]
            for j, value in enumerate(powered[occ]):
                values[j][s][f] += factor*value
        entries = [[poly(row) for row in item] for item in values]
        return [[entries[0], entries[1]], [entries[2], entries[3]]], [[entries[4]], [entries[5]]]

    def multiply(self, left, right, maximum=MAX_FULL):
        """Positive matrix product with a separate singleton-degree axis."""
        assert len(left[0]) == len(right)
        out = []
        for row in left:
            outrow = []
            for j in range(len(right[0])):
                entries = []
                for s in range(self.single+1):
                    self.guard()
                    value = self.poly([])
                    for k, entry in enumerate(row):
                        for t in range(s+1):
                            if len(entry[t]) and len(right[k][j][s-t]):
                                value += (entry[t]*right[k][j][s-t]).truncate(maximum+1)
                    entries.append(self.poly([self.positive(value[f])
                                              for f in range(min(maximum+1, len(value)))]))
                outrow.append(entries)
            out.append(outrow)
        return out

    def regions(self, regular, terminal):
        n, updates = 1, 0
        for _ in range(7):
            terminal, regular = self.multiply(regular, terminal), self.multiply(regular, regular)
            updates, n = n+updates, 2*n
        assert n == 128 and updates == 127
        return regular, terminal

    def extract(self, table, full, single):
        if single > self.single or full > MAX_FULL:
            raise ValueError("requested profile exceeds retained degrees")
        denominator = self.exact(F(comb(4096, full)*comb(4096-full, single))
                                 *USCALE**full*VSCALE**single)
        return [[self.positive(entry[single][full]/denominator) for entry in row] for row in table]

    def matmul(self, left, right):
        return [[sum((left[i][k]*right[k][j] for k in range(len(right))), self.arb(0))
                 for j in range(len(right[0]))] for i in range(len(left))]

    def matpow(self, value, exponent):
        out = [[self.arb(int(i == j)) for j in range(2)] for i in range(2)]
        while exponent:
            if exponent & 1:
                out = self.matmul(out, value)
            exponent //= 2
            if exponent:
                value = self.matmul(value, value)
        return out

    def chain(self, regular, terminal):
        results = []
        for orientation in (0, 1):
            x, y = regular[orientation], regular[1-orientation]
            value = self.matmul(self.matpow(self.matmul(x, y), 127), x)
            value = self.matmul(value, terminal[1-orientation])[0][0]
            results.append(self.positive(value))
        assert 127*2+1 == 255 and 255*128+127 == 32767
        source = (self.exact(ETA)*self.exact(BETA)**256*self.exact(Z)**(-LIMIT))**self.exact(ALPHA)
        return [4096*source*x for x in results]

    def toy_test(self):
        """Independent exact bivariate one-state convolution containment."""
        raw = {(0, 0): F(1, 3), (1, 0): F(2, 5), (0, 1): F(3, 7)}
        entry = [self.poly([self.exact(raw[0, 0]), self.exact(raw[1, 0])]),
                 self.poly([self.exact(raw[0, 1])])]
        entry += [self.poly([]) for _ in range(self.single-1)]
        actual = self.multiply([[entry]], [[entry]], 2)[0][0]
        expected = defaultdict(F)
        for (f, s), x in raw.items():
            for (ff, ss), y in raw.items():
                if s+ss <= self.single:
                    expected[f+ff, s+ss] += x*y
        for (f, s), x in expected.items():
            assert actual[s][f].contains(self.fmpq(x.numerator, x.denominator))
        # Interval min/max enclose the mathematical choices; no midpoint branch.
        x, y = self.exact(F(1, 3)), self.exact(F(2, 3))
        assert x.min(y).contains(self.fmpq(1, 3))
        assert x.max(y).contains(self.fmpq(2, 3))
        return len(expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(192, 256, 384, 512), default=192)
    parser.add_argument("--single-degree", type=int, choices=(2, 4), default=2)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    from flint import arb, ctx
    ctx.prec = args.precision
    v = Verifier(args.single_degree, started+48)
    toy = v.toy_test()
    records, checks = route_records(args.single_degree)
    v.guard()
    shapes = sorted({key[2] for key in records})
    caps = v.local_caps(shapes)
    regular, terminal = v.regions(*v.polynomials(records, caps))
    profiles = ((128, 0), (127, 2))
    rr = [v.extract(regular, *profile) for profile in profiles]
    tt = [v.extract(terminal, *profile) for profile in profiles]
    contributions = v.chain(rr, tt)
    bound = sum(contributions, arb(0))
    v.guard()
    passed = bool(bound < arb(2)**(-200))
    print(json.dumps({"scope": "all state paths; consecutive256-group supports only",
        "directed": True, "whole_code_certificate": False,
        "parameters": {"q": Q, "z": str(Z), "alpha": str(ALPHA), "distance_threshold": LIMIT,
            "source_count_per_orientation": 4096, "profiles": profiles,
            "full_degree": MAX_FULL, "singleton_degree": args.single_degree,
            "full_scale": str(USCALE), "singleton_scale": str(VSCALE)},
        "precision_bits": args.precision, "proved_below_2_pow_minus_200": passed,
        "log2_bound_interval": str(bound.log()/arb(2).log()),
        "orientation_log2_intervals": [str(x.log()/arb(2).log()) for x in contributions],
        "bound_interval": str(bound), "route_records": len(records), "occupancy_shapes": len(shapes),
        "tests": {"exact_route_normalizers": checks, "rational_quiet_return": "passed",
            "positive_OR_recurrences": len(shapes), "toy_convolution_coefficients": toy,
            "terminal_updates": 32767, "directed_minmax": "passed"},
        "nonnegative_intersections": v.intersections, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))
    if not passed:
        raise SystemExit("directed interval does not prove the fixed200-bit target")


if __name__ == "__main__":
    main()
